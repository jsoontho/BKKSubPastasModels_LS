# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:29:28 2024

@author: jtsoonthornran

REAL CASE:

Code to run ESMDA with state vector including pastas parameters,
pumping time series, and subsidence multiplier parameters

Best info is basin-wide pumping

Code for plotting is separate

For well nest with one well (ie one time series)
"""

# ##############################################################################

###############################################################################
# import statements
###############################################################################

import os
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import pickle
import pastas as ps
import sys
import warnings

# Bangkok Subsidence Model Package
import bkk_sub_gw

# Ignoring Pastas warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

# Changing current directory to locaiton of python script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

############################################################################
# Generate ensemble of pumping time series settings
###############################################################################


def generate_pumping_ens(ann_pump, n, option):
    """Generates ensemble of time series of groundwater pumping.

    Input:
    ann_pump - annual pumping rates (initial, mean) + std
    n - number of ensemble members
    option - AR for error or completely random within normal dist

    Returns:
    interp_pump - list_dates with interpolated pumping m3/day
    """

    # Number of values for pumping time series to randomly choose
    n_pump = len(ann_pump.Pump2)

    # ann_dates has Pump2 for each t (mean) and std which will be
    # used to generate normal dist at each t
    # New list for all values about to be randomly chosen and
    # interpolated
    ann_pump.index = pd.to_datetime(ann_pump.index)
    df = pd.DataFrame(index=ann_pump.index)

    # For each ensemble member (each consisting of a time series)
    for i in range(n):

        temp_list = []

        # If taking from normal distribution
        if option[0] == "normal":

            # For each t
            for t in range(n_pump):

                # Mean, std
                temp = np.random.normal(ann_pump.iloc[t, 1],
                                        ann_pump.iloc[t, 2])

                # Make sure temp is > 0
                if temp < 0:

                    temp_list.append(0)

                else:

                    temp_list.append(temp)

        # Formating
        mat = pd.DataFrame(temp_list, index=ann_pump.index,
                           columns=[i])
        df = pd.concat([df, mat], join="outer",
                       keys=["Date", "Date"], axis=1)
        df.columns = df.columns.droplevel()

    return df


# %%#########################################################################
# Pastas settings
###############################################################################

# Importing model

# Folder to save/import graph and model
modelpath = os.path.abspath("models//")

# Total path
tot_path = os.path.abspath("inputs")

# Importing model
# Model files
modelfiles = os.listdir(modelpath)

wellnestlist = ["LCBKK005"]

Wellnest_name = wellnestlist[0]

###############################################################################
# Creating/Importing and Plotting Pastas Model
###############################################################################

# Saving obs and indices
# Creating observations for ESMDA
dobs = pd.Series(np.empty(1, dtype=object))
# Saving gw obs index
gw_obs_indices = []

# Only open figures for each well nest
plt.close("all")

# Saving results for each well nest
models = []
time_mins = []
time_maxs = []
well_names = []

# Reading in groundwater data
full_path = os.path.join(tot_path, Wellnest_name + ".xlsx")
data = pd.read_excel(full_path, skiprows=3)

# BAR PLOT preparation
daterange = pd.date_range(dt.datetime(1978, 12, 31), periods=43,
                          freq="Y").tolist()
df = pd.DataFrame(daterange, columns=["date"])

loc = os.path.join(os.path.abspath("inputs"), "SurveyingLevels.xlsx")

subdata = pd.read_excel(loc, sheet_name=Wellnest_name + "_Leveling",
                        index_col=3)
subdata = pd.DataFrame(subdata)
subdata.index = pd.to_datetime(subdata.index)

# Getting rid of benchmarks outside time period
subdata = subdata[(subdata.Year <= 2020)]

# Benchmarks should start at 0 at the first year.
bench = subdata.loc[:, subdata.columns.str.contains("Land")]

# IMPORTANT INFO
# For benchmark measurements, the first year is 0, the second year
# is the compaction rate over that first year.
# For implicit Calc, the first year has a compaction rate over that
# year, so to shift benchmarks value to the previouse year to match
# Index has the right years
bench.index = bench.index.shift(-1, freq="D")
bench["date"] = bench.index

# Gets the last date of each year
lastdate = bench.groupby(pd.DatetimeIndex(bench["date"]).year,
                         as_index=False).agg(
                             {"date": max}).reset_index(drop=True)
bench = bench.loc[lastdate.date]
bench[bench == 0] = np.nan  # 0's to nan

# Subsidence observations
dobs = pd.concat([dobs, bench.iloc[:, 0].dropna()])

# Saving obs indices
gw_obs_indices.append(bench.iloc[:, 0].dropna().index)

# Groundwater observations
# For all wells in well nest
for wells in data.columns[-(len(data.columns)-2):]:

    # Name of well as a string
    well_name = wells

    #######################################################################
    # Creating Pastas Model
    #######################################################################

    # If file exists:
    try:

        # Load existing model
        wellmodel = [s for s in modelfiles
                     if np.logical_and(Wellnest_name in s,
                                       well_name in s)][0]
        model = ps.io.load(modelpath + "/" + wellmodel)

    except:

        sys.exit("Model doesn't exist")

    # Saving groundwater observations and indices
    gw_obs_indices.append(model.observations().index)
    dobs = pd.concat([dobs, model.observations()])

# Taking away first values
dobs = dobs[1:]

# Pumping
# Annual pumping data (mean), std
pumppath = os.path.join(os.path.abspath("inputs"),
                        "BasinPumping_Annual_ESMDA.xlsx")
pumpsheet = "EstTotalPump_54-60"
annual_pump = pd.read_excel(pumppath, sheet_name=pumpsheet,
                            index_col=0, parse_dates=["Date"])

# List of daily dates to be interpolated
pumpsheet = "EstTotalPump_54-60_Int50"
listdaily_pump = pd.read_excel(pumppath, sheet_name=pumpsheet,
                               index_col=0, parse_dates=["Date"])

annual_pump2 = annual_pump.copy()

# Only until 2024
annual_pump = annual_pump[annual_pump.index <= "2023"]
listdaily_pump = listdaily_pump[listdaily_pump.index <= "2023"]

# %% ESMDA Implementation

tmin = "1978"
tmax = "2021"

# Reading in thickness and storage data
path = os.path.join(os.path.abspath("inputs"),
                    "SUBParametersPriortoManual.xlsx")
# path = os.path.join(os.path.abspath("inputs"), "SUBParameters.xlsx")
Thick_data = pd.read_excel(path, sheet_name="Thickness",
                           index_col=0)  # Thickness
Sskv_data = pd.read_excel(path,
                          sheet_name="Sskv",
                          index_col=0)  # Sskv
Sske_data = pd.read_excel(path,
                          sheet_name="Sske",
                          index_col=0)  # Ssk
K_data = pd.read_excel(path,
                       sheet_name="K",
                       index_col=0)  # K

# Mode can be "raw" as in raw groundwater data vs "Pastas" for importing Pastas
# simulated groundwater in the aquifers
mode = "Pastas"

# If mode is Pastas, need model path
if mode == "Pastas":

    mpath = os.path.abspath("models//")

# Pumping flag, for PASTAS, if changing pumping scenario
pumpflag = 1
# If changing pumping scenario, need pumping sheet/path
if pumpflag == 1:

    ppath = os.path.join(os.path.abspath("inputs"), "BasinPumping.xlsx")

    psheet = "EstTotalPump_54-60_Int50"

# Convergence criteria
CC = 1 * 10**-5

# Number of nodes in clay
node_num = 10

# Using available heads as proxy for missing
proxyflag = 1

ic_run = True

# CALIBRATION SETTINGS
return_sub = False
p_multop = [True, "SsK"]

par_error = [4, 2, 4]

# Ensemble size (Geir used for the figures 1e7)
# (Reduce to speed up; costly is the pd-estimation for plotting,
# not ES-MDA)
ne = int(1000)
na = 8  # Number of assimilation steps
na_win = 1
obs_error = 3.5
par_error.extend([30, 100, 100, 30])
dist = ["norm", "norm", "norm", "norm"]

# ESMDA mode
esmdaflag = "my_gwparam_subSELECT_pump"

# Index of interested parameters
param_index = np.array([0, 1, 2, 3])

# Preallocation of each ensemble result
ens_covmmresults = {}
ens_mresults = {}
ens_CIresults = {}
ens_covmmresults["cov_mm"] = []
ens_mresults["mprior"] = []
ens_CIresults["CI"] = []

# Generate pumping ensemble
option = ["normal", .99]

# Pumping error in prior
pump_err = .5
annual_pump["Std"] = annual_pump['Pump'] * pump_err
pumping_ens = generate_pumping_ens(annual_pump, ne, option)

# Path to import old subsidence results
path = os.path.abspath("models")

# Reload object from file
file2 = open(path + "\\" + Wellnest_name + "_500.pkl", "rb")
model_sub = pickle.load(file2)
file2.close()

# Another option for clay heads
h_ic = [x[7][:, -1] for x in model_sub["all_results"]]
# Gets hidden clay groundwater state and time
t_ic = [x[6] for x in model_sub["all_results"]]

# If number of assimilation windows > 1
# Preallocation and separating obs into different windows
if na_win > 1:

    # Dataframe of obs to split
    gwobs_df = pd.DataFrame(dobs[len(gw_obs_indices[0]):],
                            index=gw_obs_indices[1])
    subobs_df = pd.DataFrame(dobs[:len(gw_obs_indices[0])],
                             index=gw_obs_indices[0])

    # Subsets
    gwobs_df_sub = np.array_split(gwobs_df, na_win)
    subobs_df_sub = np.array_split(subobs_df, na_win)

# Number of ensembles
n = 1
for n_ens in range(n):

    # Number of data assimilation windows
    if na_win == 1:
        sub_m, \
            models, \
                well_names = bkk_sub_gw.bkk_sub.bkk_subsidence(wellnestlist,
                                                               mode, tmin,
                                                               tmax,
                                                               Thick_data,
                                                               K_data,
                                                               Sskv_data,
                                                               Sske_data,
                                                               CC=CC,
                                                               Nz=node_num,
                                                               ic_run=True,
                                                               califlag=True,
                                                               p_multop=p_multop,
                                                               na=na, na_win=na_win,
                                                               ne=ne,
                                                               obs_error=obs_error,
                                                               par_error=par_error,
                                                               return_sub=return_sub,
                                                               proxyflag=proxyflag,
                                                               pumpflag=pumpflag,
                                                               pump_path=ppath,
                                                               pump_sheet=psheet,
                                                               model_path=mpath,
                                                               esmdaflag=esmdaflag,
                                                               esmdaindex=param_index,
                                                               dist=dist,
                                                               user_obs=dobs,
                                                               user_obs_indices=gw_obs_indices,
                                                               pump_ens=pumping_ens,
                                                               annual_pump=annual_pump,
                                                               listdaily_pump=listdaily_pump,
                                                               init_state=[t_ic, h_ic])

        # Get the approximated parameters
        averages = np.average(sub_m["mprior"][-1, :, :], axis=1)
        # For each parameter ensemble, save CI
        well_CI = []
        n_param = len(param_index)
        for i in range(n_param):

            CI = stats.t.interval(0.90, df=ne-1,
                                  loc=averages[i],
                                  scale=stats.sem(
                                      sub_m["mprior"][-1, i, :]))
            well_CI.append(np.insert(CI, 1, averages[i]))

        ens_CIresults["CI"].append(well_CI)
        ens_mresults["mprior"].append(averages)

    else:

        # Where subsets of obs kept
        obs_df = []

        # Where results kept
        sub_m = []

        for win in range(na_win):

            # Temp data
            temp_data = pd.concat([subobs_df_sub[win],
                                   gwobs_df_sub[win]], ignore_index=False)

            # Appending subsets of obs
            obs_df.append(temp_data.iloc[:, 0])

            # New gw_obs_indices for window
            gw_obs_indices = [subobs_df_sub[win].index,
                              gwobs_df_sub[win].index]

            temp_sub_m, \
                models, \
                    well_names = bkk_sub_gw.bkk_sub.bkk_subsidence(wellnestlist,
                                                                   mode, tmin,
                                                                   tmax,
                                                                   Thick_data,
                                                                   K_data,
                                                                   Sskv_data,
                                                                   Sske_data,
                                                                   CC=CC,
                                                                   Nz=node_num,
                                                                   ic_run=True,
                                                                   califlag=True,
                                                                   p_multop=p_multop,
                                                                   na=na, na_win=na_win,
                                                                   ne=ne,
                                                                   obs_error=obs_error,
                                                                   par_error=par_error,
                                                                   return_sub=return_sub,
                                                                   proxyflag=proxyflag,
                                                                   pumpflag=pumpflag,
                                                                   pump_path=ppath,
                                                                   pump_sheet=psheet,
                                                                   model_path=mpath,
                                                                   esmdaflag=esmdaflag,
                                                                   esmdaindex=param_index,
                                                                   dist=dist,
                                                                   user_obs=obs_df[win],
                                                                   user_obs_indices=gw_obs_indices,
                                                                   pump_ens=pumping_ens,
                                                                   annual_pump=annual_pump,
                                                                   listdaily_pump=listdaily_pump,
                                                                   init_state=[t_ic, h_ic])

            sub_m.append(temp_sub_m)

        # Get the approximated parameters
        averages = np.average(sub_m[win]["mprior"][-1], axis=1)
        # For each parameter ensemble, save CI
        well_CI = []
        n_param = len(averages)
        for i in range(n_param):

            CI = stats.t.interval(0.90, df=ne-1,
                                  loc=averages[i],
                                  scale=stats.sem(
                                      sub_m[win]["mprior"][-1, i, :]))
            well_CI.append(np.insert(CI, 1, averages[i]))

        ens_CIresults["CI"].append(well_CI)
        ens_mresults["mprior"].append(averages)

    # Saving
    # Path to save models
    path = os.path.abspath("models//ESMDA//BKK//na8ne1000//")

    # Write
    for na_i in range(na+1):
        fig_name1 = wellnestlist[
            0] + "_Dpred_" + p_multop[1] + "_na" + str(na_i) + ".csv"
        full_figpath = os.path.join(path, fig_name1)
        np.savetxt(full_figpath, sub_m["dpred"][na_i], delimiter=",")

        fig_name1 = wellnestlist[
            0] + "_Mprior_" + p_multop[1] + "_na" + str(na_i) + ".csv"
        full_figpath = os.path.join(path, fig_name1)
        np.savetxt(full_figpath, sub_m["mprior"][na_i], delimiter=",")
