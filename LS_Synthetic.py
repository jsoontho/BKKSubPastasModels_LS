# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:29:28 2024

@author: jtsoonthornran

SYNTHETIC CASE:

Code to run LEAST SQUARES with state vector including pastas parameters,
pumping time series, and subsidence multiplier parameters

Synthetic case where true pumping is arbitrarily lower for 1970-1990
and best info is basin-wide pumping

Code for plotting is separate

For well nest with one well (ie one time series)
"""

###############################################################################
# import statements
###############################################################################

import os
import pandas as pd
import numpy as np
import warnings
import lmfit
import time

# Bangkok Subsidence Model Package
import bkk_sub_gw

# Ignoring Pastas warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

# Changing current directory to locaiton of python script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# %%##########################################################################
# Generate ensemble of pumping time series settings
###############################################################################

st = time.time()
st2 = time.process_time()


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

    # Ann pump log
    ann_log = ann_pump.copy()
    ann_log.Pump = np.log(ann_log.Pump)
    ann_log.Std = abs(.1 * ann_log.Pump)

    # For each ensemble member (each consisting of a time series)
    for i in range(n):

        temp_list = []

        if option[0] == "normal":

            # For each t
            for t in range(n_pump):

                # Choosing between normal distribution
                # (mean, std)
                temp = np.random.normal(ann_log.iloc[t, 1],
                                        ann_log.iloc[t, 2])

                # Make sure temp is > 0 or < 500
                if np.exp(temp) < 0:

                    temp_list.append(np.log(0))
                if np.exp(temp) > 500:

                    temp_list.append(np.log(500))
                else:

                    temp_list.append(temp)
        # EXP
        temp_list = np.exp(temp_list)
        mat = pd.DataFrame(temp_list, index=ann_pump.index,
                           columns=[i])
        df = pd.concat([df, mat], join="outer",
                       keys=["Date", "Date"], axis=1)
        df.columns = df.columns.droplevel()

    return df


# For each well nest
wellnestlist = ["LCBKK018"]

pumpexperiment = "pumpingcase1"

if pumpexperiment == "pumpingcase1":
    # Folder to save/import graph and model
    path = os.path.abspath("models//ESMDA//bangkok-based//")

# Cyclical pump
elif pumpexperiment == "cyclical":
    # Folder to save/import graph and model
    path = os.path.abspath("models//cyclical//")


# Reading
fig_name1 = wellnestlist[0] + "_GWObs.csv"
full_figpath = os.path.join(path, fig_name1)
gw_obs = pd.read_csv(full_figpath, delim_whitespace="\t")
gw_obs.Date = pd.to_datetime(gw_obs.Date, format="%Y-%m-%d")
gw_obs.index = gw_obs.Date

fig_name1 = wellnestlist[0] + "_SubObs.csv"
full_figpath = os.path.join(path, fig_name1)
sub_obs = pd.read_csv(full_figpath, delim_whitespace="\t")
sub_obs.Date = pd.to_datetime(sub_obs.Date, format="%Y-%m-%d")
sub_obs.index = sub_obs.Date
# Saving truth
truth = []
fig_name1 = wellnestlist[0] + "_SubTruth.csv"
full_figpath = os.path.join(path, fig_name1)
truth_temp = pd.read_csv(full_figpath, delim_whitespace="\t")
truth_temp.Date = pd.to_datetime(truth_temp.Date, format="%Y-%m-%d")
truth_temp.index = truth_temp.Date
truth.append(truth_temp)

fig_name1 = wellnestlist[0] + "_GWTruth.csv"
full_figpath = os.path.join(path, fig_name1)
truth_temp = pd.read_csv(full_figpath, delim_whitespace="\t")
truth_temp.Date = pd.to_datetime(truth_temp.Date, format="%Y-%m-%d")
truth_temp.index = truth_temp.Date
truth.append(truth_temp)

fig_name1 = wellnestlist[0] + "_PumpTruth.csv"
full_figpath = os.path.join(path, fig_name1)
pumptrue = pd.read_csv(full_figpath, delim_whitespace="\t")
pumptrue.Date = pd.to_datetime(pumptrue.Date, format="%Y-%m-%d")
pumptrue.index = pumptrue.Date

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
# Only until 2024
annual_pump2 = annual_pump.copy()
annual_pump = annual_pump[annual_pump.index <= "2023"]
listdaily_pump = listdaily_pump[listdaily_pump.index <= "2023"]

# %% LS Implementation
importing = 0
saving = 1

# True parameters
Atrue = -.1
ntrue = 2.5
atrue = 50
dtrue = 2

# Calibration period
calitime_min = "1978"
calitime_max = "2005"

# Creating observations for LS
dobs = pd.Series(np.empty(1, dtype=object))
tempdata = sub_obs[sub_obs.index <= "1996"].iloc[:, -1]
dobs = pd.concat([dobs, tempdata])
# Saving gw dates
gw_obs_indices = [tempdata.index]
tempdata = gw_obs[gw_obs.index <= calitime_max].iloc[:, -1]
tempdata = gw_obs[gw_obs.index >= calitime_min].iloc[:, -1]
dobs = pd.concat([dobs, tempdata])
dobs = dobs[1:]

# Reading in groundwater data
# Total path
tot_path = os.path.abspath("inputs")
full_path = os.path.join(tot_path, wellnestlist[0] + ".xlsx")
data = pd.read_excel(full_path, skiprows=3)

# Well names
well_names = data.columns[-(len(data.columns)-2):]
num_wells = len(well_names)

for well_i in range(num_wells):

    tempdata = gw_obs.Date[well_i*int(
        len(gw_obs)/num_wells):(well_i+1)*int(len(gw_obs)/num_wells)]
    tempdata = tempdata[tempdata <= calitime_max]
    tempdata = tempdata[tempdata >= calitime_min]
    gw_obs_indices.append(tempdata)

Wellnest_name = wellnestlist[0]

tmin = "1978"
tmax = "2020"

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

    if pumpexperiment == "pumpingcase1":
        # Folder to save/import graph and model
        mpath = os.path.abspath("models//bangkok-based//")

    # Cyclical pump
    elif pumpexperiment == "cyclical":
        # Folder to save/import graph and model
        mpath = os.path.abspath("models//cyclical//")

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

# Ensemble size (Geir used for the figures 1e7)
# (Reduce to speed up; costly is the pd-estimation for plotting,
# not ES-MDA)
ne = int(1)
na = 8  # Number of assimilation steps
obs_error = 3.5

# Index of interested parameters
param_index = np.array([0, 1, 2, 3])

par_error = [30, 100, 100, 30]
dist = ["norm", "norm", "norm", "norm"]

esmdaflag = "ls_gwparam_subSELECT_pump"

# Generate pumping ensemble
option = ["normal", .99]

# Pumping error in prior
pump_err = .5
annual_pump["Std"] = annual_pump['Pump'] * pump_err
pumping_ens = generate_pumping_ens(annual_pump, ne, option)

lambda_ = 15
# Number of ensembles
n = 1
for n_ens in range(n):
    ls_sub_m, \
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
                                                           na=na,
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
                                                           lambda_=lambda_)

lmfit.fit_report(ls_sub_m)
ls_sub_m.params.pretty_print()

# save fit report to a file:
with open(mpath + "//" + wellnestlist[0] + '_LSreg' +
          str(lambda_) + '_modelresult.txt', 'w') as fh:
    fh.write(lmfit.fit_report(ls_sub_m))

et = time.time()
et2 = time.process_time()
elapse_time = et - st
elapse_time2 = et2 - st2
print('Execution time:', elapse_time/(60*60), 'hours')
print('CPU Execution time:', elapse_time2/(60*60), 'hours')
