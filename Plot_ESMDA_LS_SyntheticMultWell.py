# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:47:51 2024

@author: jtsoonthornran

SYNTHETIC CASE:

Code to plot ESMDA results with state vector including pastas parameters,
pumping time series, and subsidence multiplier parameters

Synthetic case where true pumping is arbitrarily lower for 1970-1990
and best info is basin-wide pumping

Code for running ESMDA is in "ESMDA_SyntheticMultWell"

For well nest with multiple wells (ie multiple time series)
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
from matplotlib.lines import Line2D
import pastas as ps
from sklearn.metrics import mean_squared_error
from matplotlib.pyplot import cm
from itertools import compress
from matplotlib.patches import ConnectionPatch
import warnings

# Bangkok Subsidence Model Package
import bkk_sub_gw

# Ignoring Pastas warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

# Changing current directory to locaiton of python script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# True parameters for Pastas and subsidence mult (a, b, c)
Atrue = -.1
ntrue = 2.5
atrue = 50
dtrue = 2

a = 3.68
b = .15
c = 4.8

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
annual_pump = annual_pump[annual_pump.index <= "2023"]
listdaily_pump = listdaily_pump[listdaily_pump.index <= "2023"]

wellnestlist = Wellnest_name = ["LCBKK018"]

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

    modelpath = os.path.abspath("models//ESMDA//all//random_pump_NA")

# Path to save models
tot_path = os.path.abspath("inputs")
figpath = os.path.join(os.path.abspath("figures"),
                       "ESMDA_Syn\\Synthetic\\all\\PastasNA+pump\\randompumping\\")

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

# If running initial condiiton
ic_run = True

# CALIBRATION SETTINGS
return_sub = False
p_multop = [True, "SsK"]

par_error = [4, 2, 4]

# If mode is Sskv
if p_multop[1] == "Sskv":
    Sske_data.loc[Wellnest_name][::2] = Sskv_data.loc[Wellnest_name][::2] * b
    Sske_data.loc[Wellnest_name][1::2] = Sske_data.loc[Wellnest_name][0::2] / 10

    K_data.loc[Wellnest_name] *= c

# If mode is Sskv and K
elif p_multop[1] == "SsK":

    Sske_data.loc[Wellnest_name][::2] = Sskv_data.loc[Wellnest_name][::2] * b
    Sske_data.loc[Wellnest_name][1::2] = Sske_data.loc[Wellnest_name][0::2] / 10

# Ensemble size (Geir used for the figures 1e7)
# (Reduce to speed up; costly is the pd-estimation for plotting,
# not ES-MDA)
ne = int(250)
na = 8  # Number of assimilation steps
obs_error = 3.5
par_error.extend([30, 100, 100, 30])
dist = ["norm", "norm", "norm", "norm"]

# Esmda mode
esmdaflag = "my_gwparam_subSELECT_pump"

# Index of interested parameters
param_index = np.array([1, 2])

# Reading in groundwater data
full_path = os.path.join(tot_path, Wellnest_name[0] + ".xlsx")
data = pd.read_excel(full_path, skiprows=3)

# Number of pumping time series
n_pump = len(annual_pump)
n_pump2 = 107

# Number of well/time series models
num_wells = len(data.columns)-2

# Number of pastas parameters
n_param = len(param_index)

# Number of subsidence parameters
if p_multop[1] == "Sskv" or p_multop[1] == "K" or p_multop[1] == "Sske":
    n_sub = 1

elif p_multop[1] == "Ss" or p_multop[1] == "SsK":
    n_sub = 2

elif p_multop[1] == "all":
    n_sub = 3

# Getting well names
well_names = data.columns[-(len(data.columns)-2):].values
well_names = well_names.tolist()

# Reading
# GW obs
fig_name1 = wellnestlist[0] + "_GWObs.csv"
full_figpath = os.path.join(modelpath, fig_name1)
gw_obs = pd.read_csv(full_figpath, delim_whitespace="\t")
gw_obs.Date = pd.to_datetime(gw_obs.Date, format="%Y-%m-%d")

# Sub obs
fig_name1 = wellnestlist[0] + "_SubObs.csv"
full_figpath = os.path.join(modelpath, fig_name1)
sub_obs = pd.read_csv(full_figpath, delim_whitespace="\t")
sub_obs.Date = pd.to_datetime(sub_obs.Date, format="%Y-%m-%d")

# Saving gw dates
gw_obs_indices = [sub_obs.Date]

for well_i in range(num_wells):

    gw_obs_indices.append(gw_obs.Date[well_i*int(
        len(gw_obs)/num_wells):(well_i+1)*int(len(gw_obs)/num_wells)])

# Saving truth
truth = []
# Sub truth
fig_name1 = wellnestlist[0] + "_SubTruth.csv"
full_figpath = os.path.join(modelpath, fig_name1)
truth_temp = pd.read_csv(full_figpath, delim_whitespace="\t")
truth_temp.Date = pd.to_datetime(truth_temp.Date, format="%Y-%m-%d")
truth_temp.index = truth_temp.Date
truth.append(truth_temp)

# GW Truth
fig_name1 = wellnestlist[0] + "_GWTruth.csv"
full_figpath = os.path.join(modelpath, fig_name1)
truth_temp = pd.read_csv(full_figpath, delim_whitespace="\t")
truth_temp.Date = pd.to_datetime(truth_temp.Date, format="%Y-%m-%d")
truth_temp.index = truth_temp.Date
truth.append(truth_temp)

# Pumping truth
fig_name1 = wellnestlist[0] + "_PumpTruth.csv"
full_figpath = os.path.join(modelpath, fig_name1)
pumptrue = pd.read_csv(full_figpath, delim_whitespace="\t")
pumptrue.Date = pd.to_datetime(pumptrue.Date, format="%Y-%m-%d")
pumptrue.index = pumptrue.Date

# Creating blank dictionary
sub_m = {"dpred": np.zeros((na+1, len(sub_obs)+len(gw_obs),
                            ne)),
         "mprior": np.zeros((na+1, n_pump+num_wells*n_param+n_sub,
                            ne))}

# Saving ESMDA results
for na_i in range(na+1):

    fig_name1 = wellnestlist[0] + "_Dpred_" + p_multop[1] + "_na" + str(na_i) + ".csv"
    full_figpath = os.path.join(modelpath, fig_name1)
    sub_m["dpred"][na_i] = pd.read_csv(full_figpath, delimiter=",", header=None,
                                       dtype=str)

    sub_m["dpred"][na_i] = sub_m["dpred"][na_i].astype("float64")

    fig_name1 = wellnestlist[0] + "_Mprior_" + p_multop[1] + "_na" + str(na_i) + ".csv"
    full_figpath = os.path.join(modelpath, fig_name1)

    sub_m["mprior"][na_i] = pd.read_csv(full_figpath, delimiter=",", header=None,
                                        dtype=str)

    sub_m["mprior"][na_i] = sub_m["mprior"][na_i].astype("float64")
ens_mresults = {}
ens_mresults["mprior"] = []
averages = np.average(sub_m["mprior"][-1, :, :], axis=1)
ens_mresults["mprior"].append(averages)

models_plot = []
time_mins_plot = []
time_maxs_plot = []
pastas_param = []

# Pumping index
pump_index0 = n_param*num_wells
pump_index1 = n_param*num_wells + n_pump

# Random number generator
rng = np.random.default_rng()

# Keeping best try Pastas models
besttry_Pastasmodels = []

# For all wells in well nest
for num_well, wells in enumerate(data.columns[-(len(data.columns)-2):]):

    # Name of well as a string
    well_name = wells

    #######################################################################
    # Assuming importing Pastas Model
    #######################################################################

    # Model files
    modelfiles = os.listdir(modelpath)

    # Load existing model
    wellmodel = [s for s in modelfiles
                 if np.logical_and(Wellnest_name[0] in s, well_name in s)][0]
    ml = ps.io.load(modelpath + "/" + wellmodel)

    besttry_Pastasmodels.append(ml)

# Preallocation
mean_dict = {}
meanpastas_dict = {}
meansub_dict = {}

# ESMDA Mean results
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
            meanpastas_dict[str(param_i)].append(param)
        # If pumping
        elif param_i >= (num_wells)*n_param and param_i < (pump_index1):

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

# For all wells in well nest
for well_i, wells in enumerate(well_names):

    # Name of well as a string
    well_name = wells
    #######################################################################

    # Getting modela nd adding pumping
    model_plot = besttry_Pastasmodels[well_i]
    model_plot.del_stressmodel("well")  # Deletes previous pumping
    EstTotPump_ = ps.StressModel(pump_interp,
                                 rfunc=ps.Gamma(), name="well",
                                 settings="well", up=False)
    model_plot.add_stressmodel(EstTotPump_)

    # Param names
    param_names = model_plot.parameters.index.values

    # Assigns parameters to previous optimal parameters and SD
    for param_i, param in enumerate(param_index):
        model_plot.set_parameter(name=param_names[param_i],
                                 initial=[sum(element) / len(element)
                                          for element in
                                          (meanpastas_dict.values())][
                                              well_i*n_param+param_i],
                                 optimal=[sum(element) / len(element)
                                          for element in
                                          (meanpastas_dict.values())][
                                              well_i*n_param+param_i])
    # Other parameters
    # Other index
    range_ = set(range(0, len(model_plot.get_parameters())))
    other_i = np.array(list(
        set(param_index).symmetric_difference(range_)))

    if len(other_i) != 0:

        # Setting to optimal parameters
        # model_plot.parameters["optimal"][other_i] = ml.get_parameters()[other_i]
        # model_plot.parameters["initial"][other_i] = ml.get_parameters()[other_i]

        # Setting to truth parameters
        model_plot.set_parameter(name=param_names[0],
                                 initial=Atrue,
                                 optimal=Atrue)

        model_plot.set_parameter(name=param_names[-1],
                                 initial=dtrue,
                                 optimal=dtrue)

    # Saves model
    models_plot.append(model_plot)

    ptime_min = "1978"
    ptime_max = "2015"

    # Saving time_mins and time_maxs
    time_mins_plot.append(ptime_min)
    time_maxs_plot.append(ptime_max)

    # Parameter save
    pastas_param.append(model_plot.parameters["optimal"].rename(well_name))

plt.rc("font", size=12)  # controls default text size
plt.rc("axes", titlesize=5)  # fontsize of the title
plt.rc("axes", labelsize=6)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=6)  # fontsize of the x tick labels
plt.rc("ytick", labelsize=6)  # fontsize of the y tick labels
plt.rc("legend", fontsize=6)  # fontsize of the legend

bkk_sub_gw.bkk_plotting.Pastas_results(models_plot, Wellnest_name[0],
                                       well_names, time_mins_plot,
                                       time_maxs_plot, figpath, save=1)

# %% LS Import

# save fit report to a file:
with open(os.path.abspath(
        "models//ESMDA//all//na8ne250_pumpingcase_NA//" +
        wellnestlist[0] + "_LS_modelresult.txt"),
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
        ls_sub_m[variable_name] = float(temp[temp_i][15:28])
    else:
        ls_sub_m[variable_name] = float(temp[temp_i][14:25])

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
    model_plot = besttry_Pastasmodels[well_i]
    model_plot.del_stressmodel("well")  # Deletes previous pumping
    EstTotPump_ = ps.StressModel(pump_interp,
                                 rfunc=ps.Gamma(), name="well",
                                 settings="well", up=False)
    model_plot.add_stressmodel(EstTotPump_)

    # Assigns parameters to previous optimal parameters and SD
    for param_i in param_index:

        model_plot.set_parameter(name=param_names[param_i],
                                 initial=ls_sub_m[
                                     model_plot.parameters[
                                         "optimal"].index[param_i]+str(well_i)],
                                 optimal=ls_sub_m[
                                     model_plot.parameters[
                                         "optimal"].index[param_i]+str(well_i)])

    if len(other_i) != 0:

        # Setting to optimal parameters
        # model_plot.parameters["optimal"][other_i] = ml.get_parameters()[other_i]
        # model_plot.parameters["initial"][other_i] = ml.get_parameters()[other_i]

        # Setting to truth parameters
        model_plot.set_parameter(name=param_names[0],
                                 initial=Atrue,
                                 optimal=Atrue)

        model_plot.set_parameter(name=param_names[-1],
                                 initial=dtrue,
                                 optimal=dtrue)

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

    Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
    Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

    K_data.loc[wellnest] *= np.exp(ls_sub_m[p_multop[1]+"1"])

# Mode can be "raw" as in raw groundwater data vs "Pastas" for importing Pastas
# simulated groundwater in the aquifers
mode = "Pastas"

modelpath = os.path.abspath("models//ESMDA//all//na8ne250_pumpingcase_NA")

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

# Repeat interp for each well
pump_interp = pd.concat(
    [pump_interp.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

# Reorder well list to shallow to deep aquifers
# BK, PD, NL, NB
well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
              if y in x]

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
# For each model
fig, axs = plt.subplots(num_wells+1, 2, figsize=(16, 16), sharey="row")
for mod_i in range(len(besttry_Pastasmodels)):

    well_name = well_names[mod_i]

    # Coordinate specific
    if mod_i == 0:
        index0 = len(gw_obs_indices[0])
        index1 = len(gw_obs_indices[1]) + index0

        # Axes coord
        row = mod_i
        col = 0

        axs[row, col].set(xticklabels=[])

        # Annotations
        axs[row, col].annotate("(a)",
                               xy=(.1, .01), xycoords="axes fraction",
                               fontsize=20, horizontalalignment="right",
                               weight="bold",
                               verticalalignment="bottom")

    else:
        index0 += len(gw_obs_indices[mod_i])
        index1 += len(gw_obs_indices[mod_i+1])

        # Axes coord
        # If mod_i == 1
        if mod_i == 1:

            row = mod_i
            col = 0
            axs[row, col].set(xticklabels=[])

            axs[row, col].annotate("(c)",
                                   xy=(.1, .01), xycoords="axes fraction",
                                   fontsize=20, horizontalalignment="right",
                                   weight="bold",
                                   verticalalignment="bottom")

        elif mod_i == 2:
            row = mod_i
            col = 0

            axs[row, col].annotate("(e)",
                                   xy=(.1, .01), xycoords="axes fraction",
                                   fontsize=20, horizontalalignment="right",
                                   weight="bold",
                                   verticalalignment="bottom")

        elif mod_i == 3:
            row = mod_i
            col = 0

            axs[row, col].annotate("(g)",
                                   xy=(.1, .01), xycoords="axes fraction",
                                   fontsize=20, horizontalalignment="right",
                                   weight="bold",
                                   verticalalignment="bottom")

    # PLOTTING
    # Getting obs, index and plotting obs and truth
    obs = besttry_Pastasmodels[mod_i].observations()
    index = besttry_Pastasmodels[mod_i].observations().index
    ls_sim = besttry_Pastasmodels[mod_i].simulate(tmin=str(index[0].year),
                                                  tmax=str(index[-1].year+1))
    axs[row, col].plot(obs, "o", label="Observations", color="black",
                       markersize=4)
    axs[row, col].plot(index, truth[1].loc[index].iloc[:, -1],
                       "-*c", label="Truth", linewidth=4)

    # For each alpha ESMDA
    for n in range(na+1):

        params = {
            'color': 'blue' if n == na else 'red',  # Last blue, rest red
            'lw': 2 if 0 <= n <= na-1 else 4,      # First/last thick
            'alpha': 1 if n in [0, na] else n/na,  # start faint
            'label': ['Prior', 'ESMDA Steps', *((na-2)*('',)), 'Posterior'][n],
        }

        axs[row, col].plot(gw_obs_indices[mod_i+1],
                           np.mean(sub_m["dpred"][n], axis=1)[
                               index0:index1],
                           **params)

    # Annotating
    axs[row, col].annotate(well_name,
                           xy=(-.1, 1.2), xycoords="axes fraction",
                           fontsize=20, horizontalalignment="center",
                           weight="bold",
                           bbox=dict(boxstyle="round", fc="0.8"),
                           verticalalignment="baseline")
    axs[row, col].set_ylabel("Head\n(m)")

    # add rsq to simulation
    rmse = mean_squared_error(np.mean(sub_m["dpred"][-1], axis=1)[
        index0:index1], gw_obs.iloc[(
            index0-len(gw_obs_indices[0])):(index1-len(gw_obs_indices[0])), -1],
        squared=False)

    # Plotting another subplot (GW)
    for n_ens in range(ne):

        axs[mod_i, 1].plot(index,
                           sub_m['dpred'][-1, index0:index1, n_ens], color="peru",
                           alpha=1, linewidth=1)

    # Axes coord
    if mod_i == 0:

        # Axes coord
        row = mod_i
        col = 1

        axs[row, col].set(xticklabels=[])

        axs[row, col].annotate("(b)",
                               xy=(.1, .01), xycoords="axes fraction",
                               fontsize=20, horizontalalignment="right",
                               weight="bold",
                               verticalalignment="bottom")

    elif mod_i == 1:

        row = mod_i
        col = 1
        axs[row, col].set(xticklabels=[])

        axs[row, col].annotate("(d)",
                               xy=(.1, .01), xycoords="axes fraction",
                               fontsize=20, horizontalalignment="right",
                               weight="bold",
                               verticalalignment="bottom")
    elif mod_i == 2:
        row = mod_i
        col = 1
        axs[row, col].set(xticklabels=[])

        axs[row, col].annotate("(f)",
                               xy=(.1, .01), xycoords="axes fraction",
                               fontsize=20, horizontalalignment="right",
                               weight="bold",
                               verticalalignment="bottom")

    elif mod_i == 3:
        row = mod_i
        col = 1

        axs[row, col].annotate("(h)",
                               xy=(.1, .01), xycoords="axes fraction",
                               fontsize=20, horizontalalignment="right",
                               weight="bold",
                               verticalalignment="bottom")

    axs[row, col].plot(index, model_plotls[mod_i].simulate()[index],
                       "--", color="fuchsia",
                       linewidth=4, label="Least Squares", zorder=10)
    axs[row, col].plot(obs, "o", label="Observations", color="black",
                       markersize=4)
    axs[row, col].plot(index, truth[1].loc[index].iloc[:, -1],
                       "-*c", label="Truth", linewidth=4)
    axs[row, col].plot(gw_obs_indices[mod_i+1],
                       np.mean(sub_m["dpred"][-1], axis=1)[
                           index0:index1],
                       **params)

    # Annotating
    axs[row, col].annotate(well_name,
                           xy=(1.1, 1.2), xycoords="axes fraction",
                           fontsize=20, horizontalalignment="center",
                           weight="bold",
                           bbox=dict(boxstyle="round", fc="0.8"),
                           verticalalignment="baseline")

    axs[row, col].annotate("RMSE: " + "{:.1f}".format(rmse) + " m",
                           xy=(.99, .01), xycoords="axes fraction",
                           fontsize=14, horizontalalignment="right",
                           weight="bold",
                           verticalalignment="bottom")

    # access legend objects automatically created from data
    handles, labels = plt.gca().get_legend_handles_labels()

    # Adding lines for ensemble members
    line = Line2D([0], [0], label='Ensemble Members', color='tab:orange',
                  alpha=.8)

    # plt.ylim([-30, 0])
    # add manual symbols to auto legend
    handles.extend([line])

    plt.legend(handles=handles)
    plt.legend()

    # Print observed range
    obsgw_range = obs.max() - obs.min()
    print("Observed range of head is : (m) " + f'{obsgw_range:.2f}')
    print("RMSE (gwtruth): (m)" + str(mean_squared_error(
        np.mean(sub_m["dpred"][-1], axis=1)[
            index0:index1], truth[1].iloc[
                index0:index1, -1],
        squared=False)))

obssub_range = sub_obs.iloc[:, -1].max() - sub_obs.iloc[:, -1].min()
print("Observed range of sub is : (cm/yr) " + f'{obssub_range:.2f}')
print("RMSE (subtruth): (m) " + str(
    mean_squared_error(np.mean(sub_m["dpred"][-1], axis=1)[
        0:len(gw_obs_indices[0])], truth[0].iloc[:, -1],
        squared=False)))

# Subsidence plot
if num_wells > 3:

    row = num_wells
    col = 0

elif num_wells == 3:

    row = 3
    col = 0

elif num_wells == 1:

    row = 1
    col = 0

elif num_wells == 2:

    row = 2
    col = 0


# Plotting
obs = -sub_obs.iloc[:, -1]
index = gw_obs_indices[0]
axs[row, col].plot(index, obs, "o", label="Observations", color="black",
                   markersize=8, zorder=2)
axs[row, col].plot(index, -truth[0].dropna().iloc[:, -1],
                   "-*c", label="Truth", linewidth=6)

axs[row, col].annotate("(i)",
                       xy=(.1, .01), xycoords="axes fraction",
                       fontsize=20, horizontalalignment="right",
                       weight="bold",
                       verticalalignment="bottom")

# New axis
ax2 = axs[row, col].twinx()

# For each alpha
for n in range(na+1):

    params = {
        'color': 'blue' if n == na else 'red',  # Last blue, rest red
        'lw': 2 if 0 <= n <= na-1 else 4,      # First/last thick
        'alpha': 1 if n in [0, na] else n/na,  # start faint
        'label': ['Prior', 'ESMDA Steps', *((na-2)*('',)), 'Posterior'][n],
    }

    if n == 0:

        ax2.plot(index, -np.mean(sub_m["dpred"][n], axis=1)[
            0:len(gw_obs_indices[0])],
            **params)
    else:
        axs[row, col].plot(index,
                           -np.mean(sub_m["dpred"][n], axis=1)[
                               0:len(gw_obs_indices[0])],
                           **params)

# Annotating and labeling
axs[row, col].set_ylabel("Rates\n(cm/yr)")
axs[row, col].set_xlabel("Date")

# add rsq to simulation
rmse = mean_squared_error(np.mean(sub_m["dpred"][-1], axis=1)[
    0:len(gw_obs_indices[0])], sub_obs.iloc[:, -1], squared=False)

ax2.set_ylabel("Prior Rates", rotation=270, labelpad=15)

axs[row, col].set_zorder(1)  # default zorder is 0 for ax1 and ax2
axs[row, col].patch.set_visible(False)

axs[row, col].annotate(Wellnest_name[0][2:],
                       xy=(-.13, 1.2), xycoords="axes fraction",
                       fontsize=20, horizontalalignment="center",
                       weight="bold",
                       bbox=dict(boxstyle="round", fc="0.8"),
                       verticalalignment="baseline")

# access legend objects automatically created from data
handles, labels = axs[row, col].get_legend_handles_labels()

# Adding lines for ensemble members
line = Line2D([0], [0], label='Prior', color='tab:red',
              alpha=1, linewidth=2)

# add manual symbols to auto legend
handles.insert(3, line)

# Plotting
# FOr each ensemble member
for n_ens in range(ne):

    axs[row, 1].plot(index,
                     -sub_m['dpred'][-1, :len(gw_obs_indices[0]), n_ens], color="peru",
                     alpha=1, linewidth=1)

# Plot truth, least squares, ESMDA posterior, observations
axs[row, 1].plot(index, -truth[0].dropna().iloc[:, -1],
                 "-*c", label="Truth", linewidth=4)
axs[row, 1].plot(bestsubtry[0][1].iloc[:, 2][index].index,
                 -bestsubtry[0][1].iloc[:, 2][index]*100, "--",
                 color="fuchsia",
                 label="Least Squares", linewidth=4)
axs[row, 1].plot(index,
                 -np.mean(sub_m["dpred"][
                     -1, :len(gw_obs_indices[0])], axis=1),
                 "-b", linewidth=4,
                 label="ESMDA Analysis")
axs[row, 1].plot(index, obs, "o", label="Observations", color="black",
                 markersize=8, zorder=10)

# Annotating and formatting
axs[row, 1].annotate("RMSE: " + "{:.1f}".format(rmse) + " cm/year",
                     xy=(.99, .01), xycoords="axes fraction",
                     fontsize=14, horizontalalignment="right",
                     weight="bold",
                     verticalalignment="bottom")
axs[row, 1].annotate(Wellnest_name[0][2:],
                     xy=(1.13, 1.2), xycoords="axes fraction",
                     fontsize=20, horizontalalignment="center",
                     weight="bold",
                     bbox=dict(boxstyle="round", fc="0.8"),
                     verticalalignment="baseline")
axs[row, 1].annotate("(j)",
                     xy=(.1, .01), xycoords="axes fraction",
                     fontsize=20, horizontalalignment="right",
                     weight="bold",
                     verticalalignment="bottom")
# FORMATTING
axs[row, 1].set_ylim([-10, 40])
# Shrink current axis by 20%
box = axs[row, 1].get_position()

# Put a legend to the right of the current axis
axs[row, 1].legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.25))

fig.subplots_adjust(hspace=0.5)
fig.set_rasterized(True)

# Saving
fig_name1 = Wellnest_name[0] + "_SUB_GW_ESMDA_na.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

fig_name1 = Wellnest_name[0] + "_SUB_GW_ESMDA_na.eps"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="eps")

# %%


# Pastas parameter plots
# For each model
def kde(data, points):
    return stats.gaussian_kde(data).evaluate(points)


# Number of pastas parameters
n_param = len(param_index)

fig, axs = plt.subplots(2, 3, figsize=(16, 10))
# All Pastas parameters
if n_param == 4:

    # Plotting only first well
    mod_i = 0

    # For each parameter
    for param_i in range(n_param):

        # A
        if param_i == 0:
            row = 0
            col = 0
            axs[row, col].set_title('A', fontsize=20)
            axs[row, col].axvline(ls_sub_m["well_A"+str(mod_i)],
                                  linestyle="--", color="fuchsia", linewidth=4,
                                  label="Least Squares")
            truthval = Atrue
            axs[row, col].set(ylabel='Frequency')
            axs[row, col].annotate("(a)",
                                   xy=(.1, 1.01), xycoords="axes fraction",
                                   fontsize=20, horizontalalignment="right",
                                   weight="bold",
                                   verticalalignment="bottom")
        # n
        elif param_i == 1:
            row = 0
            col = 1
            axs[row, col].set_title('n', fontsize=20)
            axs[row, col].axvline(ls_sub_m["well_n"+str(mod_i)],
                                  linestyle="--", color="fuchsia", linewidth=4,
                                  label="Least Squares")
            axs[row, col].annotate("(b)",
                                   xy=(.1, 1.01), xycoords="axes fraction",
                                   fontsize=20, horizontalalignment="right",
                                   weight="bold",
                                   verticalalignment="bottom")
            truthval = ntrue

        # a
        elif param_i == 2:
            row = 1
            col = 0
            axs[row, col].set_title('a', fontsize=20)
            axs[row, col].axvline(ls_sub_m["well_a"+str(mod_i)],
                                  linestyle="--", color="fuchsia", linewidth=4,
                                  label="Least Squares")
            truthval = atrue
            axs[row, col].set(ylabel='Frequency')
            axs[row, col].annotate("(d)",
                                   xy=(.1, 1.01), xycoords="axes fraction",
                                   fontsize=20, horizontalalignment="right",
                                   weight="bold",
                                   verticalalignment="bottom")

        # d
        else:
            row = 1
            col = 1
            axs[row, col].set_title('d', fontsize=20)
            axs[row, col].axvline(ls_sub_m["constant_d"+str(mod_i)],
                                  linestyle="--", color="fuchsia", linewidth=4,
                                  label="Least Squares")
            axs[row, col].set_xlabel('Values', fontsize=20)
            axs[row, col].annotate("(e)",
                                   xy=(.1, 1.01), xycoords="axes fraction",
                                   fontsize=20, horizontalalignment="right",
                                   weight="bold",
                                   verticalalignment="bottom")
            truthval = dtrue

        # # Plot step (ESMDA)
        for i in range(na+1):
            params = {
                'color': 'blue' if i == na else 'red',  # Last blue, rest red
                'lw': 2 if i in [0, na] else 1,      # First/last thick
                'alpha': 1 if i in [0, na] else i/na,  # start faint
                'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
            }
            paramshist = {
                'color': 'blue' if i == na else 'red',  # Last blue, rest red
                'alpha': 1 if i in [0, na] else i/na,  # start faint
                'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
            }

            # Histogram of parameter ensemble
            axs[row, col].hist(
                sub_m['mprior'][i, mod_i*n_param+param_i, :], **paramshist)

        # Plotting truth
        axs[row, col].axvline(truthval, color="c", linewidth=5,
                              label="Truth")

# All Pastas parameters
if n_param == 2:

    # Plotting only first well
    mod_i = 0

    # For each parameter
    for param_i in range(n_param):

        # A
        if param_i == 0:
            row = 0
            col = 1
            axs[row, col].set_title('n', fontsize=20)
            axs[row, col].axvline(ls_sub_m["well_n"+str(mod_i)],
                                  linestyle="--", color="fuchsia", linewidth=4,
                                  label="Least Squares")
            axs[row, col].annotate("(b)",
                                   xy=(.1, 1.01), xycoords="axes fraction",
                                   fontsize=20, horizontalalignment="right",
                                   weight="bold",
                                   verticalalignment="bottom")
            truthval = ntrue
        # n
        elif param_i == 1:
            row = 1
            col = 1
            axs[row, col].set_title('a', fontsize=20)
            axs[row, col].axvline(ls_sub_m["well_a"+str(mod_i)],
                                  linestyle="--", color="fuchsia", linewidth=4,
                                  label="Least Squares")
            truthval = atrue
            axs[row, col].set(ylabel='Frequency')
            axs[row, col].annotate("(d)",
                                   xy=(.1, 1.01), xycoords="axes fraction",
                                   fontsize=20, horizontalalignment="right",
                                   weight="bold",
                                   verticalalignment="bottom")

        # # Plot step (ESMDA)
        for i in range(na+1):
            params = {
                'color': 'blue' if i == na else 'red',  # Last blue, rest red
                'lw': 2 if i in [0, na] else 1,      # First/last thick
                'alpha': 1 if i in [0, na] else i/na,  # start faint
                'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
            }
            paramshist = {
                'color': 'blue' if i == na else 'red',  # Last blue, rest red
                'alpha': 1 if i in [0, na] else i/na,  # start faint
                'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
            }

            # Histogram of parameter ensemble
            axs[row, col].hist(
                sub_m['mprior'][i, mod_i*n_param+param_i, :], **paramshist)

        # Plotting truth
        axs[row, col].axvline(truthval, color="c", linewidth=5,
                              label="Truth")

# For each sub param
for param_i in range(n_sub):
    # # Plot step (ESMDA)
    for i in range(na+1):
        params = {
            'color': 'blue' if i == na else 'red',  # Last blue, rest red
            'lw': 2 if i in [0, na] else 1,      # First/last thick
            'alpha': 1 if i in [0, na] else i/na,  # start faint
            'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
        }
        paramshist = {
            'color': 'blue' if i == na else 'red',  # Last blue, rest red
            'alpha': 1 if i in [0, na] else i/na,  # start faint
            'label': ['Initial', *((na-2)*('',)), 'MDA steps', 'MDA'][i],
        }

        # Histogram of ensemble
        axs[param_i, 2].hist(
            sub_m['mprior'][i, num_wells*n_param+n_pump+param_i, :], **paramshist)

    # Truth val, title, save title
    if p_multop[1] == "Sskv" or p_multop[1] == "SsK":

        # Sskv
        if param_i == 0:

            # plt.title('Sskv', fontsize=14)
            truthval = a
            axs[param_i, 2].set_title('Log Sskv Multiplier', fontsize=20)

            axs[param_i, 2].axvline(np.log(truthval), label="Truth",
                                    color="c", linewidth=4)
            axs[param_i, 2].axvline(ls_sub_m["SsK0"], linestyle="--",
                                    color="fuchsia", linewidth=4,
                                    label="Least Squares")
            axs[param_i, 2].set_xlim([-2, 2])
            axs[param_i, 2].legend()
            axs[param_i, 2].annotate("(c)",
                                     xy=(.1, 1.01), xycoords="axes fraction",
                                     fontsize=20, horizontalalignment="right",
                                     weight="bold",
                                     verticalalignment="bottom")

        # Sske
        elif param_i == 1:
            # plt.suptitle("ESMDA: LOG K Parameters for Well")

            # plt.title('K', fontsize=14)
            truthval = c
            axs[param_i, 2].set_title('Log K Multiplier', fontsize=20)

            axs[param_i, 2].axvline(np.log(truthval), color="c", linewidth=4,
                                    label="Truth")
            axs[param_i, 2].axvline(ls_sub_m["SsK1"], linestyle="--",
                                    color="fuchsia", linewidth=4,
                                    label="Least Squares")
            axs[param_i, 2].annotate("(f)",
                                     xy=(.1, 1.01), xycoords="axes fraction",
                                     fontsize=20, horizontalalignment="right",
                                     weight="bold",
                                     verticalalignment="bottom")

            axs[param_i, 2].set_xlim([-10, 10])

fig_name1 = wellnestlist[0] + "_GWSub_ESMDA_1param_LOG" + p_multop[1] + ".png"
fig_name2 = wellnestlist[0] + "_GWSub_ESMDA_1param_LOG" + p_multop[1] + ".eps"
fig.set_rasterized(True)
# Saving figure
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

full_figpath = os.path.join(figpath, fig_name2)
plt.savefig(full_figpath, bbox_inches="tight", format="eps")


# Plotting observations with simulation
# SUBSIDNECE
plt.figure(figsize=(10, 4))
color = cm.YlOrBr(np.linspace(0, 1, ne))

obs_index = gw_obs_indices[0]
for n_ens in range(ne):

    plt.plot(obs_index,
             -sub_m['dpred'][-1, :len(obs_index), n_ens], color=color[n_ens],
             alpha=0.2)

plt.plot(obs_index,
         -sub_obs.iloc[:, 1], "-ok", label="Observations", zorder=10)
plt.plot(obs_index,
         -truth[0].dropna().iloc[:, -1], "-c", label="Truth", linewidth=5)
plt.plot(obs_index,
         -np.mean(sub_m["dpred"][-1, :len(obs_index)], axis=1), "-b", linewidth=8,
         label="ESMDA Analysis")
plt.plot(bestsubtry[0][1].iloc[:, 2][index].index,
         -bestsubtry[0][1].iloc[:, 2][index]*100, "--", color="fuchsia",
         label="Least Squares", linewidth=6)

plt.ylabel("Rates\n(cm/yr)")
plt.legend()
fig_name1 = Wellnest_name[0] + "Sub_ESMDA_evol.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")


# %% PDF

# GW Domain
sub_index = gw_obs_indices[0]

curr_index = len(sub_index)
curr_index1 = curr_index
# For each model
for mod_i in range(num_wells):
    obs_index = gw_obs_indices[mod_i+1]
    curr_index1 += len(obs_index)
    gw_quant05 = []
    gw_quant95 = []

    # For each time
    for i in range(len(obs_index)):
        i += len(sub_index)
        if mod_i >= 1:
            i += len(obs_index) * mod_i

        min_ = min(sub_m["dpred"][-1][i, :])
        max_ = min(sub_m["dpred"][-1][i, :])
        domain = np.linspace(min_,
                             max_,
                             ne)

        # Pdf
        pdf = kde(sub_m["dpred"][-1][i, :], domain)

        # Quantiles
        a, b = np.quantile(sub_m["dpred"][-1][i, :], [0.05, .95])
        gw_quant05.append(a)
        gw_quant95.append(b)

    plt.figure(figsize=(10, 4))
    plt.plot(obs_index,
             np.mean(sub_m[
                 "dpred"][-1, curr_index:curr_index1], axis=1), "-b", linewidth=8,
             label="ESMDA Analysis")
    plt.plot(obs_index, truth[1].loc[obs_index].iloc[:, -1],
             "-c", linewidth=5, label="Truth")
    plt.plot(obs_index, gw_quant05, "--", color="tab:purple",
             linewidth=6,
             label="5% Quantile")
    plt.plot(obs_index, gw_quant95, "--", color="tab:purple",
             linewidth=6,
             label="95% Quantile")
    plt.ylabel("Head (m)")
    plt.xlabel("Date")
    plt.legend()

    plt.plot(besttry_Pastasmodels[mod_i].observations().index,
             besttry_Pastasmodels[mod_i].observations(), "ok", label="Observations")

    fig_name1 = "Well" + str(mod_i) + "_HeadQuantile.png"
    full_figpath = os.path.join(figpath, fig_name1)
    plt.savefig(full_figpath, bbox_inches="tight", format="png")

    curr_index += len(obs_index)

sub_quant05 = []
sub_quant95 = []

# For each time
for i in range(len(sub_index)):

    min_ = min(sub_m["dpred"][-1][i, :])
    max_ = min(sub_m["dpred"][-1][i, :])
    domain = np.linspace(min_,
                         max_,
                         ne)

    # Pdf
    pdf = kde(sub_m["dpred"][-1][i, :], domain)

    # Quantiles
    a, b = np.quantile(sub_m["dpred"][-1][i, :], [0.05, .95])
    sub_quant05.append(-1*a)
    sub_quant95.append(-1*b)

plt.figure()
plt.plot(sub_index,
         -np.mean(sub_m["dpred"][-1, :len(gw_obs_indices[0])], axis=1), "-b", linewidth=8,
         label="ESMDA Analysis")
obs = -sub_obs.dropna().iloc[:, -1]
plt.plot(sub_index, obs, "-o", label="Obs", color="black")
plt.plot(sub_index, -truth[0].dropna().iloc[:, -1], "-*c", label="Truth", linewidth=6)

plt.plot(sub_index, sub_quant05, "--", color="tab:purple",
         linewidth=6,
         label="5% Quantile")
plt.plot(sub_index, sub_quant95, "--", color="tab:purple",
         linewidth=6,
         label="95% Quantile")
plt.ylabel("Subsidence Rates (cm/yr)")
plt.xlabel("Date")
plt.legend()

fig_name1 = "SubQuantile.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

# %% Paper pumping

fig, axs = plt.subplots(2, 1, figsize=(16, 10))

# Pumping subset
obs_index = gw_obs_indices[1]
pump_index = pumptrue[
    np.logical_and(
        pumptrue.index >= obs_index[0], pumptrue.index <= obs_index.iloc[-1])]
pump_indices = np.logical_and(
    pumptrue.index >= obs_index[0], pumptrue.index <= obs_index.iloc[-1])

pump_index2 = annual_pump[
    np.logical_and(
        annual_pump.index >= obs_index[0], annual_pump.index <= obs_index.iloc[-1])]
pump_indices2 = np.logical_and(
    annual_pump.index >= obs_index[0], annual_pump.index <= obs_index.iloc[-1])

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
    axs[0].plot(pump_index.index,
                sub_m['mprior'][i,
                                num_wells *
                                n_param:num_wells *
                                n_param+n_pump][
                                    pump_indices],
                **params)
axs[0].plot(pump_index.index, pumptrue[pump_indices].Pump,
            color="c", linewidth=4,
            marker="o")
axs[0].plot(pump_index.index,
            np.mean(sub_m['mprior'][-1,
                                    num_wells*n_param:num_wells*n_param+n_pump],
                    axis=1)[pump_indices], color="blue",
            alpha=1, linewidth=4, zorder=10,
            marker="o")
axs[0].plot(pump_index.index,
            list(compress(pump_mean, pump_indices)), color="fuchsia",
            alpha=1, linewidth=4,
            linestyle="--")
axs[0].set_ylabel("Pumping Rate * 10$^4$ (m$^3$/day)")
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

pump_2020Index = pumptrue.index <= "2020"
# Gw obs range
gw_min_index = gw_obs_indices[1][0]
gw_max_index = gw_obs_indices[1].iloc[-1]

for mod_i in range(len(besttry_Pastasmodels)):

    # Posterior mean
    axs[1].plot(pumptrue.index[pump_2020Index],
                np.mean(sub_m['mprior'][-1,
                                        num_wells*n_param:num_wells*n_param+n_pump],
                        axis=1)[pump_2020Index], color="blue",
                alpha=1, linewidth=4,
                marker="o")

    # Truth
    axs[1].plot(pumptrue.index[pump_2020Index], pumptrue.Pump[pump_2020Index],
                color="c", linewidth=4,
                marker="o")
    # Observations
    axs[1].plot(annual_pump.index, annual_pump.Pump, "-", color="tab:red", linewidth=3)
    # LS
    axs[1].plot(pump_interp.index, pump_interp.iloc[:, -1], "--",
                color="fuchsia", linewidth=3)
    axs[1].set_ylabel("Pumping Rate * 10$^4$ (m$^3$/day)")
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
fig_name1 = wellnestlist[0] + "_Pumping_ESMDA_Paper.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

fig_name1 = wellnestlist[0] + "_Pumping_ESMDA_Paper.eps"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="eps")

# %% Pumping during obs

fig, axs = plt.subplots(1, 1, figsize=(16, 5))

# Pumping subset
obs_index = gw_obs_indices[1]
pump_index = pumptrue[
    np.logical_and(
        pumptrue.index >= obs_index[0], pumptrue.index <= obs_index.iloc[-1])]
pump_indices = np.logical_and(
    pumptrue.index >= obs_index[0], pumptrue.index <= obs_index.iloc[-1])

pump_index2 = annual_pump[
    np.logical_and(
        annual_pump.index >= obs_index[0], annual_pump.index <= obs_index.iloc[-1])]
pump_indices2 = np.logical_and(
    annual_pump.index >= obs_index[0], annual_pump.index <= obs_index.iloc[-1])

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
    axs.plot(pump_index.index,
             sub_m['mprior'][i,
                             num_wells * n_param:num_wells *
                             n_param+n_pump][
                                 pump_indices],
             **params)
axs.plot(pump_index.index, pumptrue[pump_indices].Pump,
         color="c", linewidth=4,
         marker="o")
axs.plot(pump_index.index,
         np.mean(sub_m['mprior'][-1,
                                 num_wells*n_param:num_wells*n_param+n_pump],
                 axis=1)[pump_indices], color="blue",
         alpha=1, linewidth=4, zorder=10,
         marker="o")
axs.plot(pump_index.index,
         np.mean(sub_m['mprior'][0,
                                 num_wells*n_param:num_wells*n_param+n_pump],
                 axis=1)[pump_indices], color="red",
         alpha=1, linewidth=4, zorder=1,
         marker="o")
axs.plot(pump_index.index,
         list(compress(pump_mean, pump_indices)), color="fuchsia",
         alpha=1, linewidth=4,
         linestyle="--")
axs.set_ylabel("Pumping Rate * 10$^4$ (m$^3$/day)")
axs.set_xlabel("Date")

# access legend objects automatically created from data
handles, labels = axs.get_legend_handles_labels()

# Adding lines for ensemble members
linetruth = Line2D([0], [0], label='Truth', color='c',
                   alpha=1, linewidth=4, marker='o')
linepostmean = Line2D([0], [0], label='Posterior Mean', color='blue',
                      alpha=1, linewidth=4, marker='o')
linepriormean = Line2D([0], [0], label='Prior Mean', color='red',
                       alpha=1, linewidth=4, marker='o')
lineleast = Line2D([0], [0], label='Least Squares', color='fuchsia',
                   alpha=1, linewidth=4, linestyle="--")
line = Line2D([0], [0], label='Prior Ensemble', color='red',
              alpha=.8)
line2 = Line2D([0], [0], label='Posterior Ensemble', color='blue',
               alpha=.8)

# plt.ylim([-30, 0])
# add manual symbols to auto legend
handles.extend([linetruth, linepriormean, linepostmean, lineleast, line, line2])

axs.legend(handles=handles)
fig.set_rasterized(True)

# Saving
fig_name1 = wellnestlist[0] + "_ObsPumping_ESMDA_Paper.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

fig_name1 = wellnestlist[0] + "_ObsPumping_ESMDA_Paper.eps"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="eps")
