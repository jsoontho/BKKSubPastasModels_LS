# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:47:51 2024

@author: jtsoonthornran

REAL CASE:

Code to plot ESMDA results with state vector including pastas parameters,
pumping time series, and subsidence multiplier parameters

Best info is basin-wide pumping

Code for running ESMDA is in "ESMDA_SyntheticOneWell"

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
from matplotlib.lines import Line2D
import pastas as ps
from sklearn.metrics import mean_squared_error
from matplotlib.pyplot import cm
import matplotlib.dates as mdates
from itertools import compress
from matplotlib.patches import ConnectionPatch
import sys
import warnings

# Bangkok Subsidence Model Package
import bkk_sub_gw

# Ignoring Pastas warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

# Changing current directory to locaiton of python script
os.chdir(os.path.dirname(os.path.abspath(__file__)))


# %% Plotting

wellnestlist = ["LCBKK018"]

Wellnest_name = wellnestlist[0]

# Folder to save/import graph and model
modelpath = os.path.abspath("models//")

# Total path
tot_path = os.path.abspath("inputs")

# Saving sub and gw obs
dobs = pd.Series(np.empty(1, dtype=object))
gw_obs_indices = []

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

sub_obs = bench.iloc[:, 0].dropna()

# Reading in groundwater data
full_path = os.path.join(tot_path, Wellnest_name + ".xlsx")
data = pd.read_excel(full_path, skiprows=3)

# Importing model
# Model files
modelfiles = os.listdir(modelpath)

# Groundwater observations
# For all wells in well nest
models = []
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
        models.append(model)

        gw_obs = model.observations()
        # Saving groundwater observations and indices
        gw_obs_indices.append(model.observations().index)
        dobs = pd.concat([dobs, model.observations()])
    except:

        sys.exit("Model doesn't exist")

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

p_multop = [True, "SsK"]

par_error = [4, 2, 4]

# Ensemble size (Geir used for the figures 1e7)
# (Reduce to speed up; costly is the pd-estimation for plotting,
# not ES-MDA)
ne = int(250)
na = 8  # Number of assimilation steps
na_win = 1
obs_error = 3.5
par_error.extend([30, 100, 100, 30])
dist = ["norm", "norm", "norm", "norm"]

# Index of interested parameters
param_index = np.array([0, 1, 2, 3])

n_pump2 = len(annual_pump2)

n_param = len(param_index)

# Number of wells
num_wells = len(data.columns[-(len(data.columns)-2):])
# WEll names
well_names = [data.columns[-(len(data.columns)-2):][0]]
# Number of pumping
# Number of pumping, pastas, sub parameters
n_pump = len(annual_pump)

# Number of subsidence parameters
if p_multop[1] == "Sskv" or p_multop[1] == "K" or p_multop[1] == "Sske":
    n_sub = 1

elif p_multop[1] == "Ss" or p_multop[1] == "SsK":
    n_sub = 2

elif p_multop[1] == "all":
    n_sub = 3

# Creating blank dictionary
sub_m = {"dpred": np.zeros((na+1, len(sub_obs)+len(gw_obs),
                            ne)),
         "mprior": np.zeros((na+1, n_pump+n_param+n_sub,
                            ne))}

# Path to save models
modelpath = os.path.abspath("models//ESMDA//BKK//")
# Saving ESMDA results
gw_obs_indices = [sub_obs.index, gw_obs.index]
for na_i in range(na+1):

    fig_name1 = wellnestlist[0] + "_Dpred_" + p_multop[1] + "_na" + str(na_i) + ".csv"
    full_figpath = os.path.join(modelpath, fig_name1)
    sub_m["dpred"][na_i] = pd.read_csv(full_figpath, delimiter=",", header=None,
                                       dtype=str)

    sub_m["dpred"][na_i] = sub_m["dpred"][na_i].astype("float64")

    fig_name1 = wellnestlist[
        0] + "_Mprior_" + p_multop[1] + "_na" + str(na_i) + ".csv"
    full_figpath = os.path.join(modelpath, fig_name1)

    sub_m["mprior"][na_i] = pd.read_csv(full_figpath, delimiter=",", header=None,
                                        dtype=str)

    sub_m["mprior"][na_i] = sub_m["mprior"][na_i].astype("float64")

ens_mresults = {}
ens_mresults["mprior"] = []
averages = np.average(sub_m["mprior"][-1, :, :], axis=1)
ens_mresults["mprior"].append(averages)

# Reading in groundwater data
full_path = os.path.join(tot_path, Wellnest_name + ".xlsx")
data = pd.read_excel(full_path, skiprows=3)

models_plot = []
time_mins_plot = []
time_maxs_plot = []
pastas_param = []

# Pumping index
pump_index0 = n_param*num_wells
pump_index1 = n_param*num_wells + n_pump

# For all wells in well nest
for well_i, wells in enumerate(well_names):

    # Name of well as a string
    well_name = wells
    #######################################################################

    mean_dict = {}
    meanpastas_dict = {}
    meansub_dict = {}

    for num_ens, ens in enumerate(ens_mresults[
            str(list(ens_mresults.keys())[-1])]):
        if num_ens == 0:

            # Initializing pumping, pastas, sub dictionaries
            for pump_i in range(pump_index0, pump_index1):
                mean_dict[str(pump_i)] = []

            for pastas_i in range(n_param*num_wells):
                meanpastas_dict[str(pastas_i)] = []

            for sub_i in range(pump_index1, pump_index1+n_sub):
                meansub_dict[str(sub_i)] = []

        for param_i, param in enumerate(ens):

            # If pastas parameters
            if param_i < (num_wells)*n_param:
                # meanpastas_dict[str(param_i)].append(np.mean(param,
                #                                              axis=0))
                meanpastas_dict[str(param_i)].append(param)
            # If pumping
            elif param_i >= (num_wells)*n_param and param_i < (pump_index1):

                # mean_dict[str(param_i-(num_wells)*4)].append(np.mean(param,
                #                                              axis=0))
                mean_dict[str(param_i)].append(param)

            # Sub
            elif param_i >= (pump_index1) and param_i < pump_index1+n_sub:
                meansub_dict[str(param_i)].append(param)

    # Mean time series to mean complete time series
    # Adds new pumping
    # Taking mean of the ensemble
    pump_mean = [sum(element) / len(element)
                 for element in mean_dict.values()]

    mean_ = pd.Series(pump_mean)
    mean_.index = annual_pump.index

    # Isolating pumping data
    pump_df = pd.DataFrame(mean_, index=annual_pump.index,
                           columns=["0"])
    pump_df.index = annual_pump.index
    df = pd.DataFrame(index=listdaily_pump.index)
    df = pd.concat([df, pump_df], join="outer",
                   keys=["Date", "Date"], axis=1)
    df.columns = df.columns.droplevel()

    # Interpolating pumping data
    pump_interp = df.interpolate(method="cubic")
    pump_interp = pump_interp.dropna()

    model_plot = models[well_i]
    model_plot.del_stressmodel("well")  # Deletes previous pumping
    EstTotPump_ = ps.StressModel(pump_interp,
                                 rfunc=ps.Gamma(), name="well",
                                 settings="well", up=False)
    model_plot.add_stressmodel(EstTotPump_)
    # Assigns parameters to previous optimal parameters and SD
    model_plot.parameters["optimal"][param_index] = [sum(element) / len(element)
                                                     for element in
                                                     (meanpastas_dict.values())]
    model_plot.parameters["initial"][param_index] = [sum(element) / len(element)
                                                     for element in
                                                     (meanpastas_dict.values())]

    models_plot.append(model_plot)

    ptime_min = "1978"
    ptime_max = "2015"

    # Saving time_mins and time_maxs
    time_mins_plot.append(ptime_min)
    time_maxs_plot.append(ptime_max)

    # Parameter save
    pastas_param.append(model_plot.parameters["optimal"].rename(well_name))

# parameter save for well nest
# pastaspath = os.path.join(os.path.abspath("inputs"),
#                           Wellnest_name +
#                           "_ESMDA_PastasParam_na4ne500pyESMDA2.xlsx")
# pastas_param = pd.DataFrame(pastas_param)
# pastas_param.to_excel(pastaspath)
figpath = os.path.join(os.path.abspath("figures"),
                       "ESMDA\\All_Pump_ANAD")

plt.rc("font", size=12)  # controls default text size
plt.rc("axes", titlesize=5)  # fontsize of the title
plt.rc("axes", labelsize=6)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=6)  # fontsize of the x tick labels
plt.rc("ytick", labelsize=6)  # fontsize of the y tick labels
plt.rc("legend", fontsize=6)  # fontsize of the legend

bkk_sub_gw.bkk_plotting.Pastas_results(models_plot, Wellnest_name,
                                       well_names, time_mins_plot,
                                       time_maxs_plot, figpath, save=1)

# %% LS Import

# save fit report to a file:
with open(os.path.abspath(
        "models//ESMDA//BKK//" +
        wellnestlist[0] + "_LS_REAL.txt"),
        'r') as fh:
    temp = fh.readlines()
    temp = [x.replace("\n", "") for x in temp]
    temp = temp[
        temp.index("[[Variables]]") + 1:temp.index(
            "[[Correlations]] (unreported correlations are < 0.100)")]
fh.close()

# Saving variables and values
ls_sub_m = {}

# Saving values from string
for temp_i in range(len(temp)):
    variable_name = temp[temp_i].split(":")[0].strip()

    # Constant d lines are too long
    if "constant_d" in variable_name:
        ls_sub_m[variable_name] = float(temp[temp_i][16:28])
    else:
        ls_sub_m[variable_name] = float(temp[temp_i][15:27])

# Saving model plots for least squares
model_plotls = []

# For all wells in well nest
for well_i, wells in enumerate(well_names):

    # Name of well as a string
    well_name = wells
    #######################################################################

    # Saving pumping
    pump_mean = []
    for pump_i in range(n_pump):
        pump_mean.append(ls_sub_m["pump"+str(pump_i)])

    mean_ = pd.Series(pump_mean)
    mean_.index = annual_pump.index

    # Isolating pumping data
    pump_df = pd.DataFrame(mean_, index=annual_pump.index,
                           columns=["0"])
    pump_df.index = annual_pump.index
    df = pd.DataFrame(index=listdaily_pump.index)
    df = pd.concat([df, pump_df], join="outer",
                   keys=["Date", "Date"], axis=1)
    df.columns = df.columns.droplevel()

    # Interpolating pumping data
    pump_interp = df.interpolate(method="cubic")
    pump_interp = pump_interp.dropna()
    pump_interp = pump_interp.rename(columns={"0": well_name})
    model_plot = models_plot[well_i]
    model_plot.del_stressmodel("well")  # Deletes previous pumping
    EstTotPump_ = ps.StressModel(pump_interp,
                                 rfunc=ps.Gamma(), name="well",
                                 settings="well", up=False)
    model_plot.add_stressmodel(EstTotPump_)
    # Assigns parameters to previous optimal parameters and SD
    for param_i in param_index:
        model_plot.parameters["optimal"][param_i] = ls_sub_m[
            model_plot.parameters[
                "optimal"].index[param_i]+str(well_i)]
        model_plot.parameters["initial"][param_i] = ls_sub_m[
            model_plot.parameters[
                "optimal"].index[param_i]+str(well_i)]

    model_plotls.append(model_plot)

# %% BEST SUB TRY

tmin = "1978"
tmax = "2020"

# Reading in thickness and storage data
path = os.path.join(os.path.abspath("inputs"),
                    "SUBParametersPriortoManual.xlsx")
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

# Random multipler for each well nest
for wellnest in wellnestlist:

    Sskv_data.loc[wellnest] *= np.exp(ls_sub_m[p_multop[1]+"0"])

    Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * .15
    Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

    K_data.loc[wellnest] *= np.exp(ls_sub_m[p_multop[1]+"1"])

# Mode can be "raw" as in raw groundwater data vs "Pastas" for importing Pastas
# simulated groundwater in the aquifers
mode = "Pastas"

modelpath = os.path.abspath("models//")

# Pumping flag, for PASTAS, if changing pumping scenario
pumpflag = 1
ppath = os.path.join(os.path.abspath("inputs"), "BasinPumping.xlsx")
psheet = "EstTotalPump_54-60_Int50"

# Convergence criteria
CC = 1 * 10**-5

# Number of nodes in clay
node_num = 10

# Using available heads as proxy for missing
proxyflag = 1

pump_interp.columns = well_names

# Calculates subsidence
all_results, sub_total, subv_total = bkk_sub_gw.\
    bkk_sub.bkk_subsidence(wellnestlist,
                           mode, tmin,
                           tmax,
                           Thick_data,
                           K_data,
                           Sskv_data,
                           Sske_data,
                           CC=CC,
                           Nz=node_num,
                           ic_run=True,
                           proxyflag=proxyflag,
                           pumpflag=pumpflag,
                           pump_path=ppath,
                           pump_sheet=psheet,
                           pump_series=pump_interp,
                           model_path=modelpath, califlag=0,
                           esmdaflag=0, user_models=model_plotls)

# Post process data
sub_total, subv_total, ann_sub, \
    avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                             sub_total,
                                             subv_total,
                                             all_results)

bestsubtry = ann_sub

# %%###########################################################################
# Plotting settings
###############################################################################

plt.rc("font", size=28)  # controls default text size
plt.rc("axes", titlesize=24)  # fontsize of the title
plt.rc("axes", labelsize=18)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=16)  # fontsize of the x tick labels
plt.rc("ytick", labelsize=16)  # fontsize of the y tick labels
plt.rc("legend", fontsize=14)  # fontsize of the legend

# Head plot
fig, axs = plt.subplots(2, 2, figsize=(16, 10), sharey="row")
# For each model
for mod_i in range(len(model_plotls)):

    well_name = well_names[mod_i]
    # Getting obs, index and plotting obs and truth
    obs = model_plotls[mod_i].observations()
    index = model_plotls[mod_i].observations().index
    axs[0, 0].plot(obs, "o", label="Observations", color="black",
                   markersize=4)
    axs[0, 0].annotate(well_name,
                       xy=(-.05, 1.05), xycoords="axes fraction",
                       fontsize=20, horizontalalignment="center",
                       weight="bold",
                       bbox=dict(boxstyle="round", fc="0.8"),
                       verticalalignment="baseline")

    # For each alpha (ESMDA)
    for n in range(na+1):

        params = {
            'color': 'blue' if n == na else 'C3',  # Last blue, rest red
            'lw': 1 if 1 <= n <= na-1 else 4,      # First/last thick
            'alpha': 1 if n in [0, na] else n/na,  # start faint
            'label': ['Prior', *((na-2)*('',)), 'ESMDA Steps', 'Posterior'][n],
        }

        # Plotting ESMDA
        axs[0, 0].plot(gw_obs_indices[1],
                       np.mean(sub_m["dpred"][n], axis=1)[
                           len(gw_obs_indices[0]):],
                       **params)

    # add rsq to simulation
    # COmparing ESMDA with OBS
    rmse = mean_squared_error(np.mean(sub_m["dpred"][-1], axis=1)[
        len(gw_obs_indices[0]):], gw_obs, squared=False)

    axs[0, 0].set_ylabel("Head (m)")
    axs[0, 0].set_ylim([-30, 2])

    # GROUNDWATER ENSEMBLE PLOTTING
    # Index
    obs_index = gw_obs_indices[1]

    # For each ensemble member
    for n_ens in range(ne):

        axs[0, 1].plot(obs_index,
                       sub_m['dpred'][
                           -1, len(gw_obs_indices[0]):, n_ens], color="peru",
                       alpha=1, linewidth=.5)

    # Observations and esmda results
    axs[0, 1].plot(model_plotls[0].observations().index,
                   model_plotls[0].observations(), "ok", label="Observations",
                   markersize=4)
    axs[0, 1].plot(obs_index,
                   np.mean(sub_m["dpred"][
                       -1, len(gw_obs_indices[0]):], axis=1), "-b", linewidth=4,
                   label="Posterior")

    # Labeling
    axs[0, 0].annotate("(a)",
                       xy=(.1, .01), xycoords="axes fraction",
                       fontsize=20, horizontalalignment="right",
                       weight="bold",
                       verticalalignment="bottom")
    axs[0, 1].annotate("(b)",
                       xy=(.1, .01), xycoords="axes fraction",
                       fontsize=20, horizontalalignment="right",
                       weight="bold",
                       verticalalignment="bottom")
    axs[0, 1].annotate("RMSE: " + "{:.1f}".format(rmse) + " m",
                       xy=(.99, .01), xycoords="axes fraction",
                       fontsize=14, horizontalalignment="right",
                       weight="bold",
                       verticalalignment="bottom")
    axs[0, 1].set_ylim([-30, 2])

    axs[0, 1].plot(index,
                   model_plotls[mod_i].simulate()[index], "--",
                   color="fuchsia", linewidth=4,
                   label="Least Squares")

    # access legend objects automatically created from data
    handles, labels = axs[1, 0].get_legend_handles_labels()

    # Adding lines for ensemble members
    line = Line2D([0], [0], label='Ensemble Members', color='peru',
                  alpha=1)

    # plt.ylim([-30, 0])
    # add manual symbols to auto legend
    handles.extend([line])

# Subsidence plot
# For each model
for mod_i in range(len(model_plotls)):

    # Obs and index and plotting obs and truth
    obs = -sub_obs
    index = gw_obs_indices[0]
    axs[1, 0].plot(index, obs, "o", label="Observations", color="black",
                   markersize=8, zorder=10)

    if mod_i == 0:

        index0 = 0
        nd_temp_1 = len(gw_obs_indices[0])

    # New axis
    ax2 = axs[1, 0].twinx()
    # For each alpha (ESMDA)
    for n in range(na+1):

        params = {
            'color': 'blue' if n == na else 'C3',  # Last blue, rest red
            'lw': 1 if 1 <= n <= na-1 else 4,      # First/last thick
            'alpha': 1 if n in [0, na] else n/na,  # start faint
            'label': ['Prior', *((na-2)*('',)), 'ESMDA Steps', 'Posterior'][n],
        }

        if n <= 1:
            ax2.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
            ax2.plot(index,
                     -np.mean(sub_m["dpred"][n], axis=1)[
                         index0:nd_temp_1],
                     **params, zorder=10)

        else:
            # Plotting mean of ESMDA ensemble
            axs[1, 0].plot(index,
                           -np.mean(sub_m["dpred"][n], axis=1)[
                               index0:nd_temp_1],
                           **params)

    # add rsq to simulation
    # Comparing obs to ESMDA posterior mean
    rmse = mean_squared_error(np.mean(sub_m["dpred"][-1], axis=1)[
        index0:nd_temp_1], -obs, squared=False)

    # Annotations
    axs[1, 0].annotate("(c)",
                       xy=(.1, .01), xycoords="axes fraction",
                       fontsize=20, horizontalalignment="right",
                       weight="bold",
                       verticalalignment="bottom")
    axs[1, 1].annotate("(d)",
                       xy=(.1, .01), xycoords="axes fraction",
                       fontsize=20, horizontalalignment="right",
                       weight="bold",
                       verticalalignment="bottom")
    axs[1, 1].annotate("RMSE: " + "{:.1f}".format(rmse) + " cm/yr",
                       xy=(.99, .06), xycoords="axes fraction",
                       fontsize=14, horizontalalignment="right",
                       weight="bold",
                       verticalalignment="top")
    axs[1, 0].annotate(Wellnest_name,
                       xy=(-.05, 1.1), xycoords="axes fraction",
                       fontsize=20, horizontalalignment="center",
                       weight="bold",
                       bbox=dict(boxstyle="round", fc="0.8"),
                       verticalalignment="baseline")

    # X axis labeling
    axs[1, 0].xaxis.set_major_locator(mdates.YearLocator(2))
    axs[1, 0].xaxis.set_minor_locator(mdates.YearLocator())

    # access legend objects automatically created from data
    handles, labels = axs[1, 0].get_legend_handles_labels()

    # Adding lines for ensemble members
    line = Line2D([0], [0], label='Prior', color='tab:red',
                  alpha=1, linewidth=4)

    # add manual symbols to auto legend
    handles.insert(2, line)
    axs[1, 0].set_zorder(2)
    axs[1, 0].set_frame_on(False)
    ax2.set_zorder(1)
    axs[1, 0].legend(handles=handles).set_zorder(1)
    axs[1, 0].set_ylabel("Rates (cm/yr)")
    axs[1, 1].set_xlabel("Date")
    ax2.set_ylabel("Prior Rates (cm/yr)", rotation=270, labelpad=20)

    # SUBSIDENCE ENSEMBLE PLOTTING
    # Index
    obs_index = gw_obs_indices[0]

    # FOr each ensemble member
    for n_ens in range(ne):

        axs[1, 1].plot(obs_index,
                       -sub_m['dpred'][-1, :len(obs_index), n_ens], color="peru",
                       alpha=1, linewidth=.5)

    # Plotting observation, truth, esmda, least squares
    axs[1, 1].plot(obs_index,
                   -sub_obs, "ok", label="Observations",
                   zorder=6, markersize=8)
    axs[1, 1].plot(obs_index,
                   -np.mean(sub_m[
                       "dpred"][-1, :len(obs_index)], axis=1), "-b", linewidth=4,
                   label="Posterior")
    axs[1, 1].plot(bestsubtry[0][1].iloc[:, 2][index].index,
                   -bestsubtry[0][1].iloc[:, 2][index]*100, "--",
                   color="fuchsia",
                   label="Least Squares", linewidth=4)
    axs[1, 1].legend()
    axs[1, 1].xaxis.set_major_locator(mdates.YearLocator(2))
    axs[1, 1].xaxis.set_minor_locator(mdates.YearLocator())

    # access legend objects automatically created from data
    handles, labels = axs[1, 1].get_legend_handles_labels()

    # Adding lines for ensemble members
    line = Line2D([0], [0], label='Ensemble Members', color='peru',
                  alpha=1)

    # plt.ylim([-30, 0])
    # add manual symbols to auto legend
    handles.extend([line])

    axs[1, 1].legend(handles=handles)
    fig.subplots_adjust(hspace=0.35)
    # Saving
    fig_name1 = "SUB_Head_ESMDA_na.png"
    full_figpath = os.path.join(figpath, fig_name1)
    plt.savefig(full_figpath, bbox_inches="tight", format="png")

    fig_name1 = "SUB_Head_ESMDA_na.eps"
    full_figpath = os.path.join(figpath, fig_name1)
    plt.savefig(full_figpath, bbox_inches="tight", format="eps")


# Pastas parameter plots
# For each model
def kde(data, points):
    return stats.gaussian_kde(data).evaluate(points)


# Number of pastas parameters
n_param = len(param_index)

fig, axs = plt.subplots(2, 3, figsize=(16, 10))
# All Pastas parameters
if n_param == 4:

    # For each pastas models
    for mod_i in range(len(model_plotls)):

        # For each parameter
        for param_i in range(n_param):

            # A
            if param_i == 0:
                row = 0
                col = 0
                axs[row, col].set_title('A', fontsize=20)
                axs[row, col].axvline(ls_sub_m["well_A" + str(mod_i)],
                                      linestyle="--", color="fuchsia", linewidth=4,
                                      label="Least Squares")
                axs[row, col].annotate("(a)",
                                       xy=(.1, 1.01), xycoords="axes fraction",
                                       fontsize=20, horizontalalignment="right",
                                       weight="bold",
                                       verticalalignment="bottom")
                axs[row, col].set(ylabel='Frequency')
            # n
            elif param_i == 1:
                row = 0
                col = 1
                axs[row, col].set_title('n', fontsize=20)
                axs[row, col].axvline(ls_sub_m["well_n" + str(mod_i)],
                                      linestyle="--", color="fuchsia", linewidth=4,
                                      label="Least Squares")
                axs[row, col].annotate("(b)",
                                       xy=(.1, 1.01), xycoords="axes fraction",
                                       fontsize=20, horizontalalignment="right",
                                       weight="bold",
                                       verticalalignment="bottom")

            # a
            elif param_i == 2:
                row = 1
                col = 0
                axs[row, col].set_title('a', fontsize=20)
                axs[row, col].axvline(ls_sub_m["well_a" + str(mod_i)],
                                      linestyle="--", color="fuchsia", linewidth=4,
                                      label="Least Squares")
                axs[row, col].annotate("(d)",
                                       xy=(.1, 1.01), xycoords="axes fraction",
                                       fontsize=20, horizontalalignment="right",
                                       weight="bold",
                                       verticalalignment="bottom")
                axs[row, col].set(ylabel='Frequency')

            # d
            else:
                row = 1
                col = 1
                axs[row, col].set_title('d', fontsize=20)
                axs[row, col].axvline(ls_sub_m["constant_d" + str(mod_i)],
                                      linestyle="--", color="fuchsia", linewidth=4,
                                      label="Least Squares")
                axs[row, col].set_xlabel('Values', fontsize=20)
                axs[row, col].annotate("(e)",
                                       xy=(.1, 1.01), xycoords="axes fraction",
                                       fontsize=20, horizontalalignment="right",
                                       weight="bold",
                                       verticalalignment="bottom")

            # # Plot step (ESMDA)
            for i in range(na+1):
                params = {
                    'color': 'C0' if i == na else 'C3',  # Last blue, rest red
                    'lw': 2 if i in [0, na] else 1,      # First/last thick
                    'alpha': 1 if i in [0, na] else i/na,  # start faint
                    'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
                }
                paramshist = {
                    'color': 'blue' if i == na else 'C3',  # Last blue, rest red
                    'alpha': 1 if i in [0, na] else i/na,  # start faint
                    'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
                }

                # Histogram of parameter ensemble
                axs[row, col].hist(
                    sub_m['mprior'][i, mod_i*n_param+param_i, :], **paramshist)

# For each sub param
for param_i in range(n_sub):
    # # Plot step (ESMDA)
    for i in range(na+1):
        params = {
            'color': 'C0' if i == na else 'C3',  # Last blue, rest red
            'lw': 2 if i in [0, na] else 1,      # First/last thick
            'alpha': 1 if i in [0, na] else i/na,  # start faint
            'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
        }
        paramshist = {
            'color': 'blue' if i == na else 'C3',  # Last blue, rest red
            'alpha': 1 if i in [0, na] else i/na,  # start faint
            'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
        }

        # Histogram of ensemble
        axs[param_i, 2].hist(
            sub_m['mprior'][i, mod_i*n_param+param_i+n_pump, :], **paramshist)

    # Truth val, title, save title
    if p_multop[1] == "Sskv" or p_multop[1] == "SsK":

        # Sskv
        if param_i == 0:

            # plt.title('Sskv', fontsize=14)
            axs[param_i, 2].set_title('Log Sskv Multiplier', fontsize=20)

            axs[param_i, 2].axvline(ls_sub_m["SsK0"], linestyle="--",
                                    color="fuchsia", linewidth=4)

            axs[param_i, 2].annotate("(c)",
                                     xy=(.1, 1.01), xycoords="axes fraction",
                                     fontsize=20, horizontalalignment="right",
                                     weight="bold",
                                     verticalalignment="bottom")

        # Sske
        elif param_i == 1:

            axs[param_i, 2].set_title('Log K Multiplier', fontsize=20)

            axs[param_i, 2].axvline(ls_sub_m["SsK1"], linestyle="--",
                                    color="fuchsia", linewidth=4,
                                    label="Least Squares")
            axs[param_i, 2].annotate("(f)",
                                     xy=(.1, 1.01), xycoords="axes fraction",
                                     fontsize=20, horizontalalignment="right",
                                     weight="bold",
                                     verticalalignment="bottom")

            axs[param_i, 2].legend()

fig_name1 = "GWSub_ESMDA_1param_LOG" + p_multop[1] + ".png"
fig_name2 = "GWSub_ESMDA_1param_LOG" + p_multop[1] + ".eps"
fig.set_rasterized(True)
# Saving figure
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

full_figpath = os.path.join(figpath, fig_name2)
plt.savefig(full_figpath, bbox_inches="tight", format="eps")

# Paper pumping

fig, axs = plt.subplots(2, 1, figsize=(16, 10))

# Pumping subset
obs_index = gw_obs_indices[1]

pump_index2 = annual_pump[
    np.logical_and(
        annual_pump.index >= obs_index[0], annual_pump.index <= obs_index[-1])]
pump_indices2 = np.logical_and(
    annual_pump.index >= obs_index[0], annual_pump.index <= obs_index[-1])

# Pumping plot
for mod_i in range(len(model_plotls)):

    # # Plot step (ESMDA)
    for i in range(na+1):
        params = {
            'color': 'tab:blue' if i == na else 'tab:red',  # Last blue, rest red
            'lw': .5 if i in [0, na] else 1,      # First/last thick
            'alpha': .1 if i in [0, na] else i/na,  # start faint
        }
        if 0 < i < na:
            continue
        par_domain = np.linspace(-.2, 0, 1001)
        axs[0].plot(pump_index2.index,
                    sub_m['mprior'][i,
                                    num_wells *
                                    n_param:num_wells *
                                    n_param+n_pump][
                                        pump_indices2],
                    **params)
    axs[0].plot(pump_index2.index,
                np.mean(sub_m['mprior'][-1,
                                        num_wells*n_param:num_wells*n_param+n_pump],
                        axis=1)[pump_indices2], color="blue",
                alpha=1, linewidth=4, zorder=10,
                marker="o")
    axs[0].plot(pump_index2.index,
                list(compress(pump_mean, pump_indices2)), color="fuchsia",
                alpha=1, linewidth=4,
                linestyle="--")
    axs[0].set_ylabel("Pumping Rate * 10$^4$ m$^3$/day")
    axs[0].set_xlabel("Date")

# access legend objects automatically created from data
handles, labels = axs[0].get_legend_handles_labels()

# Adding lines for ensemble members
linetruth = Line2D([0], [0], label='Truth', color='c',
                   alpha=1, linewidth=4, marker='o')
linepostmean = Line2D([0], [0], label='Posterior Mean', color='blue',
                      alpha=1, linewidth=4, marker='o')
lineleast = Line2D([0], [0], label='Least Squares', color='fuchsia',
                   alpha=1, linewidth=4, linestyle="--")
line = Line2D([0], [0], label='Prior Ensemble', color='red',
              alpha=.8)
line2 = Line2D([0], [0], label='Posterior Ensemble', color='blue',
               alpha=.8)

# plt.ylim([-30, 0])
# add manual symbols to auto legend
handles.extend([linetruth, linepostmean, lineleast, line, line2])

axs[0].legend(handles=handles)

pump_2020Index = annual_pump.index <= "2020"
# Gw obs range
gw_min_index = gw_obs_indices[1][0]
gw_max_index = gw_obs_indices[1][-1]

for mod_i in range(len(model_plotls)):

    # Posterior mean
    axs[1].plot(annual_pump.index[pump_2020Index],
                np.mean(sub_m['mprior'][-1,
                                        num_wells*n_param:num_wells*n_param+n_pump],
                        axis=1)[pump_2020Index], color="blue",
                alpha=1, linewidth=4,
                marker="o")

    # Observations
    axs[1].plot(annual_pump.index, annual_pump.Pump, "-", color="tab:red", linewidth=3)
    # LS
    axs[1].plot(pump_interp.index, pump_interp.iloc[:, -1], "--",
                color="fuchsia", linewidth=3)
    axs[1].set_ylabel("Pumping Rate * 10$^4$ m$^3$/day")
    # Labeling
    axs[0].annotate("(a)",
                    xy=(0.04, .01), xycoords="axes fraction",
                    fontsize=20, horizontalalignment="right",
                    weight="bold",
                    verticalalignment="bottom")
    axs[1].annotate("(b)",
                    xy=(0.04, .01), xycoords="axes fraction",
                    fontsize=20, horizontalalignment="right",
                    weight="bold",
                    verticalalignment="bottom")
    axs[1].set_xlabel("Date")
    axs[1].legend()

    # Shade obs period
    axs[1].fill_betweenx([-50, 500],
                         gw_min_index, gw_max_index,
                         color='lightgray')

    # Adding connection
    conn = ConnectionPatch(xyA=(0, 0), coordsA='axes fraction', axesA=axs[0],
                           xyB=(gw_min_index, 500), coordsB='data', axesB=axs[1],
                           color='lightgray', linewidth=3
                           )
    conn2 = ConnectionPatch(xyA=(1, 0), coordsA='axes fraction', axesA=axs[0],
                            xyB=(gw_max_index, 500), coordsB='data', axesB=axs[1],
                            color='lightgray', linewidth=3
                            )
    axs[1].add_artist(conn)
    axs[1].add_artist(conn2)
    conn.set_in_layout(False)  # remove from layout calculations
    conn2.set_in_layout(False)  # remove from layout calculations

# access legend objects automatically created from data
handles, labels = axs[1].get_legend_handles_labels()

# Adding lines for ensemble members
linetruth = Line2D([0], [0], label='Truth', color='c',
                   alpha=1, linewidth=4, marker='o')
linepriormean = Line2D([0], [0], label='Prior Mean', color='red',
                       alpha=1, linewidth=4, marker='o')
linepostmean = Line2D([0], [0], label='Posterior Mean', color='blue',
                      alpha=1, linewidth=4, marker='o')
lineleast = Line2D([0], [0], label='Least Squares', color='fuchsia',
                   alpha=1, linewidth=4, marker='o')

# add manual symbols to auto legend
handles.extend([linetruth, linepriormean, linepostmean, lineleast])

axs[1].legend(handles=handles)

# Saving
fig_name1 = "Pumping_ESMDA_Paper.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

fig_name1 = "Pumping_ESMDA_Paper.eps"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="eps")
sys.exit()
# %%
# Head plot
plt.figure()
# For each model
for mod_i in range(len(models_plot)):

    plt.figure()
    obs = models_plot[mod_i].observations()
    index = models_plot[mod_i].observations().index
    plt.plot(obs, "-o", label="Obs", color="black")

    if mod_i == 0:

        index0 = len(gw_obs_indices[0])
        nd_temp_1 = index0 + len(obs)

    else:

        index0 += len(models_plot[mod_i-1].observations())
        nd_temp_1 += len(obs)

    # For each alpha
    for n in range(na+1):

        params = {
            'color': 'C0' if n == na else 'C3',  # Last blue, rest red
            'lw': 2 if n in [1, na-1] else 6,      # First/last thick
            'alpha': 1 if n in [0, na] else n/na,  # start faint
            'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][n],
        }
        plt.plot(index,
                 np.mean(sub_m["dpred"][n], axis=1)[
                     index0:nd_temp_1],
                 **params)

        # add rsq to simulation
        rmse = mean_squared_error(np.mean(sub_m["dpred"][n], axis=1)[
            index0:nd_temp_1], obs, squared=False)
        plt.title("ESMDA: Head (m) for Well"
                  "\nRMSE (m): " + f"{rmse:.2f}")

    plt.legend()
    # Saving figure
    fig_name1 = "Head_ESMDA_1param_" + p_multop[1] + "Individ.png"
    full_figpath = os.path.join(figpath, fig_name1)
    plt.savefig(full_figpath, bbox_inches="tight", format="png")

    # Saving figure
    # fig_name1 = "Prediction Interval_" + paramname + "Individ.png"
    # full_figpath = os.path.join(figpath, fig_name1)
    # plt.savefig(full_figpath, bbox_inches="tight", format="png")

# Subsidence plot
plt.figure()
# For each model
for mod_i in range(len(models_plot)):

    plt.figure()
    index = gw_obs_indices[0]
    obs = -dobs.iloc[:len(index)]
    plt.plot(obs, "-o", label="Obs", color="black")

    if mod_i == 0:

        index0 = 0
        nd_temp_1 = len(gw_obs_indices[0])

    # For each alpha
    for n in range(na+1):

        params = {
            'color': 'C0' if n == na else 'C3',  # Last blue, rest red
            'lw': 2 if n in [1, na-1] else 6,      # First/last thick
            'alpha': 1 if n in [0, na] else n/na,  # start faint
            'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][n],
        }
        plt.plot(index,
                 -np.mean(sub_m["dpred"][n], axis=1)[
                     index0:nd_temp_1],
                 **params)

        # add rsq to simulation
        rmse = mean_squared_error(np.mean(sub_m["dpred"][n], axis=1)[
            index0:nd_temp_1], obs, squared=False)
        plt.title("ESMDA: Head (m) for Well"
                  "\nRMSE (m): " + f"{rmse:.2f}")

    plt.legend()
    fig_name1 = "SUB_ESMDA_na.png"
    full_figpath = os.path.join(figpath, fig_name1)
    plt.savefig(full_figpath, bbox_inches="tight", format="png")

n_param = len(param_index)

# For n, d
if n_param == 2:

    if [0, 3] in param_index:
        for mod_i in range(len(models_plot)):
            fig, axs = plt.subplots(1, 2)
            plt.suptitle("ESMDA: Pastas Parameters for Well")

            for param_i in range(n_param):
                # # Plot step
                for i in range(na+1):
                    params = {
                        'color': 'C0' if i == na else 'C3',  # Last blue, rest red
                        'lw': 2 if i in [0, na] else 1,      # First/last thick
                        'alpha': 1 if i in [0, na] else i/na,  # start faint
                        'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
                    }
                    paramshist = {
                        'color': 'C0' if i == na else 'C3',  # Last blue, rest red
                        'alpha': 1 if i in [0, na] else i/na,  # start faint
                        'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
                    }

                    par_domain = np.linspace(0.01, .3, 1001)
                    axs[param_i].hist(sub_m['mprior'][
                        i, mod_i*n_param+param_i, :], **paramshist)
                    if param_i == 0:
                        axs[param_i].set_title('A', fontsize=14)
                        fig_name1 = "Head_ESMDA_1param_A_Individ.png"

                    else:
                        axs[param_i].set_title('d', fontsize=14)
                        fig_name1 = "Head_ESMDA_1param_d_Individ.png"
            print("Analysis: " +
                  f"{np.mean(sub_m['mprior'][-1, mod_i*n_param+param_i, :]): .2f}")
            for ax in axs.flat:
                # ax.set(xlabel='x-label', ylabel='y-label')
                ax.set(ylabel='Frequency')

            plt.legend()
            # Saving figure
            full_figpath = os.path.join(figpath, fig_name1)
            plt.savefig(full_figpath, bbox_inches="tight", format="png")

# Subsidence parameter plotting
for mod_i in range(len(models_plot)):

    if p_multop[1] == "Sskv":
        sub_param = 1
    elif p_multop[1] == "all":
        sub_param = 3
    for param_i in range(n_sub):
        fig, axs = plt.subplots()
        if p_multop[1] == "Sskv":
            plt.suptitle("ESMDA: Sskv Parameters for Well")

        elif p_multop[1] == "all":
            plt.suptitle("ESMDA: Sub Parameters for Well")

        if param_i == 0:
            plt.title("Sskv", fontsize=14)
            fig_name1 = "Sub_ESMDA_1param_Sskv_Individ.png"

        elif param_i == 1:

            plt.title("Sske", fontsize=14)
            fig_name1 = "Sub_ESMDA_1param_Sske_Individ.png"

        elif param_i == 2:

            plt.title("K", fontsize=14)
            fig_name1 = "Sub_ESMDA_1param_K_Individ.png"

        # # Plot step
        for i in range(na+1):
            params = {
                'color': 'C0' if i == na else 'C3',  # Last blue, rest red
                'lw': 2 if i in [0, na] else 1,      # First/last thick
                'alpha': 1 if i in [0, na] else i/na,  # start faint
                'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
            }
            paramshist = {
                'color': 'C0' if i == na else 'C3',  # Last blue, rest red
                'alpha': 1 if i in [0, na] else i/na,  # start faint
                'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
            }
            plt.hist(
                np.exp(sub_m['mprior'][i,
                                       mod_i * n_param +
                                       n_param + n_pump + param_i, :]),
                **paramshist)
            plt.xlim([0, 100])
        full_figpath = os.path.join(figpath, fig_name1)
        plt.savefig(full_figpath, bbox_inches="tight", format="png")

# Subsidence parameter plotting LOG
for mod_i in range(len(models_plot)):
    if p_multop[1] == "Sskv":
        sub_param = 1
    elif p_multop[1] == "all":
        sub_param = 3
    for param_i in range(n_sub):
        fig, axs = plt.subplots()
        if p_multop[1] == "Sskv":
            plt.suptitle("ESMDA: Sskv Parameters for Well")

        elif p_multop[1] == "all":
            plt.suptitle("ESMDA: Sub Parameters for Well")

        if param_i == 0:

            plt.title("Sskv", fontsize=14)
            fig_name1 = "Sub_ESMDA_1param_LOGSskv_Individ.png"

        elif param_i == 1:

            plt.title("Sske", fontsize=14)
            fig_name1 = "Sub_ESMDA_1param_LOGSske_Individ.png"

        elif param_i == 2:

            plt.title("K", fontsize=14)
            fig_name1 = "Sub_ESMDA_1param_LOGK_Individ.png"

        # # Plot step
        for i in range(na+1):
            params = {
                'color': 'C0' if i == na else 'C3',  # Last blue, rest red
                'lw': 2 if i in [0, na] else 1,      # First/last thick
                'alpha': 1 if i in [0, na] else i/na,  # start faint
                'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
            }
            paramshist = {
                'color': 'C0' if i == na else 'C3',  # Last blue, rest red
                'alpha': 1 if i in [0, na] else i/na,  # start faint
                'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
            }

            plt.hist(sub_m['mprior'][i,
                                     mod_i * n_param + n_param +
                                     n_pump+param_i, :],
                     **paramshist)
            plt.title('Sskv', fontsize=14)

        full_figpath = os.path.join(figpath, fig_name1)
        plt.savefig(full_figpath, bbox_inches="tight", format="png")

# Plotting observations with simulation
# SUBSIDNECE
plt.figure()
color = cm.YlOrBr(np.linspace(0, 1, ne))

obs_index = gw_obs_indices[0]
for n_ens in range(ne):

    plt.plot(obs_index,
             -sub_m['dpred'][-1, :len(obs_index), n_ens], color=color[n_ens],
             alpha=0.2)

plt.plot(obs_index,
         -dobs.iloc[:len(obs_index)], "-ok", label="Observations", zorder=10)
plt.plot(obs_index,
         -np.mean(sub_m["dpred"][-1, :len(obs_index)], axis=1), "-b", linewidth=8,
         label="ESMDA Analysis")
plt.plot(bestsubtry[0][1].iloc[:, 2][obs_index].index,
         -bestsubtry[0][1].iloc[:, 2][obs_index]*100, "--",
         color="fuchsia",
         label="Least Squares", linewidth=6)
plt.plot
plt.legend()
fig_name1 = "Sub_ESMDA_evol.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

# GW
plt.figure()
color = cm.YlOrBr(np.linspace(0, 1, ne))

obs_index = gw_obs_indices[1]
for n_ens in range(ne):

    plt.plot(obs_index,
             sub_m['dpred'][-1, len(gw_obs_indices[0]):, n_ens], color=color[n_ens],
             alpha=0.2)
plt.plot(model.observations().index,
         model.observations(), "ok", label="Observations")
plt.plot(obs_index,
         np.mean(
             sub_m["dpred"][-1, len(gw_obs_indices[0]):], axis=1), "-b", linewidth=8,
         label="ESMDA Analysis")

# LS simulation plotting
ls_sim = model.simulate(tmin=str(obs_index[0].year),
                        tmax=str(obs_index[-1].year+1))
plt.plot(obs_index, model_plotls[0].simulate()[obs_index],
         "--", color="fuchsia",
         linewidth=4, label="Least Squares", zorder=10)

plt.xlabel("Date")
plt.ylabel("Head (m)")
plt.title("Spaghetti Plot for Heads (m)")

# access legend objects automatically created from data
handles, labels = plt.gca().get_legend_handles_labels()

# Adding lines for ensemble members
line = Line2D([0], [0], label='Ensemble Members', color='tab:orange',
              alpha=.8)

# plt.ylim([-30, 0])
# add manual symbols to auto legend
handles.extend([line])

plt.legend(handles=handles)
plt.show()

fig_name1 = "Head_ESMDA_ne.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

index = -1
# Pumping plot
for mod_i in range(len(models_plot)):
    fig, axs = plt.subplots(figsize=(10, 4))
    plt.suptitle("ESMDA: Pumping for Well")

    if index == 0:
        color = "tab:red"
    elif index == -1:
        color = "tab:blue"
    axs.plot(annual_pump.index,
             np.mean(sub_m['mprior'][index,
                                     num_wells*n_param:num_wells*n_param+n_pump],
                     axis=1), color="tab:blue",
             alpha=1, linewidth=3, label="Posterior Mean",
             marker="o")
    axs.plot(annual_pump.index, annual_pump.Pump, "--m", linewidth=3,
             label="Prior Mean", alpha=.4,
             marker="o")
    axs.set_ylabel("Pumping Rate * 10000 m$^3$/day")
    axs.set_xlabel("Date")
    plt.legend()

    # New figure with alphas
    plt.figure()
    # # Plot step
    for i in range(na+1):
        params = {
            'color': 'C0' if i == na else 'C3',  # Last blue, rest red
            'lw': 2 if i in [0, na] else 1,      # First/last thick
            'alpha': 1 if i in [0, na] else i/na,  # start faint
        }
        paramshist = {
            'color': 'C0' if i == na else 'C3',  # Last blue, rest red
            'alpha': .5,  # start faint
            'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
        }

        par_domain = np.linspace(-.2, 0, 1001)
        plt.plot(annual_pump.index,
                 sub_m['mprior'][i,
                                 num_wells*n_param:num_wells*n_param+n_pump],
                 **params)
    plt.plot(annual_pump.index,
             np.mean(sub_m['mprior'][index,
                                     num_wells*n_param:num_wells*n_param+n_pump],
                     axis=1), color="blue",
             alpha=1, linewidth=3, label="Posterior Mean", zorder=10,
             marker="o")
    plt.title('ESMDA Pumping')
    plt.ylabel("Pumping Rate * 10000 m$^3$/day")
    plt.xlabel("Date")

# access legend objects automatically created from data
handles, labels = plt.gca().get_legend_handles_labels()

# Adding lines for ensemble members
line = Line2D([0], [0], label='Prior', color='red',
              alpha=.8)
line2 = Line2D([0], [0], label='Posterior', color='blue',
               alpha=.8)

# plt.ylim([-30, 0])
# add manual symbols to auto legend
handles.extend([line, line2])

plt.legend(handles=handles)
plt.show()

fig_name1 = "Pumping_ESMDA_all.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")
