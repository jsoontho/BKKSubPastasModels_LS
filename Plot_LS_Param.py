"""
Created on Thu Sep 12 10:47:1 2024

@author: jtsoonthornran

SYNTHETIC CASE:

Code to plot results with state vector including pastas parameters,
pumping time series, and subsidence multiplier parameters

Synthetic case where true pumping is arbitrarily lower for 1970-1990
and best info is basin-wide pumping

Code for running LS is in "LS"

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
from matplotlib.dates import DateFormatter
import warnings

# Bangkok Subsidence Model Package
import bkk_sub_gw

# Ignoring Pastas warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

# Changing current directory to locaiton of python script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# True parameters for Pastas and subsidence mult (a, b, c)
Atrue = -.1
ntrue = 1.2
atrue = 50
dtrue = 2

a = 3.68
b = .15
c = 4.8

pumpexperiment = "simpleexp"

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

    # pumping case 1: true 1980-1990 arbitrarily lower
    if pumpexperiment == "pumpingcase1":
        # Folder to save/import graph and model
        modelpath = os.path.abspath("models//bangkok-based//")

    # Cyclical pump
    elif pumpexperiment == "cyclical":
        # Folder to save/import graph and model
        modelpath = os.path.abspath("models//cyclical//")

    # pumping case 1: true 1980-1990 arbitrarily lower
    elif pumpexperiment == "simpleexp":
        # Folder to save/import graph and model
        modelpath = os.path.abspath(
            "models//cowboyhat_SUB//")

# Path to save models
tot_path = os.path.abspath("inputs")

# FIgure path
if pumpexperiment == "simpleexp":
    # Folder to save/import graph and model
    figpath = os.path.abspath(
        "figures//cowboyhat_SUB//")

elif pumpexperiment == "pumpingcase1":
    # Folder to save/import graph and model
    figpath = os.path.abspath(
        "figures//bangkok-based//")

elif pumpexperiment == "cyclical":
    # Folder to save/import graph and model
    figpath = os.path.abspath(
        "figures//cyclical//")

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

# Calibration period
calitime_min = "1978"
calitime_max = "2005"

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
param_index = np.array([0, 1, 2, 3])

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
gw_obs.index = gw_obs.Date
fullgwobs = gw_obs.copy()

gw_obs = gw_obs[gw_obs.index >= calitime_min]
gw_obs = gw_obs[gw_obs.index <= calitime_max]

# Sub obs
fig_name1 = wellnestlist[0] + "_SubObs.csv"
full_figpath = os.path.join(modelpath, fig_name1)
sub_obs = pd.read_csv(full_figpath, delim_whitespace="\t")
sub_obs.Date = pd.to_datetime(sub_obs.Date, format="%Y-%m-%d")
sub_obs.index = sub_obs.Date
fullsubobs = sub_obs.copy()

sub_obs = sub_obs[sub_obs.index <= "1996"]

# Saving gw dates
gw_obs_indices = [sub_obs.Date]
# gw_obs_indices = []
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
pumptrue = pumptrue.rename(columns={0: "Pump"})

models_plot = []
time_mins_plot = []
time_maxs_plot = []
pastas_param = []

# Pumping index
# pump_index0 = n_param*num_wells
# pump_index1 = n_param*num_wells + n_pump

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

    param_names = ml.parameters.index.values

pumptrue = pumptrue.rename(columns={"0": well_name})
# %% LS Import

if pumpexperiment == "simpleexp":
    # Folder to save/import graph and model
    lspath = os.path.abspath("models//cowboyhat_SUB//")

elif pumpexperiment == "cyclical":
    # Folder to save/import graph and model
    lspath = os.path.abspath("models//cyclical//")

elif pumpexperiment == "pumpingcase1":
    # Folder to save/import graph and model
    lspath = os.path.abspath("models//bangkok-based//")

# save fit report to a file:
with open(os.path.abspath(
        lspath + "//" +
        wellnestlist[0] + "_LSreg0_modelresult.txt"),
        'r') as fh:
    temp = fh.readlines()
    temp = [x.replace("\n", "") for x in temp]
    temp = temp[
        temp.index("[[Variables]]") + 1:temp.index(
            "[[Correlations]] (unreported correlations are < 0.100)")]
    # temp = temp[
    #     temp.index("[[Variables]]") + 1:]
fh.close()

# Saving variables and values
ls_sub_m = {}
ls_error = []

# Saving values from string
for temp_i in range(len(temp)):
    variable_name = temp[temp_i].split(":")[0].strip()
    # Saving error
    ls_error.append(float(temp[temp_i][33:43]))
    # Constant d lines are too long
    if "constant_d" in variable_name:
        ls_sub_m[variable_name] = float(temp[temp_i][16:28])

    else:
        ls_sub_m[variable_name] = float(temp[temp_i][14:27])

# Saving model plots for least squares
model_plotls = []

# For all wells in well nest
for well_i, wells in enumerate(well_names):

    # Name of well as a string
    well_name = wells
    #######################################################################

    # Saving pumping
    # pump_mean = []
    # for pump_i in range(n_pump):
    #     pump_mean.append(ls_sub_m["pump"+str(pump_i)])

    # mean_ = pd.Series(np.exp(pump_mean))
    # mean_.index = annual_pump.index

    # # Isolating pumping data
    # pump_df = pd.DataFrame(mean_, index=annual_pump.index,
    #                        columns=["0"])
    # pump_df.index = annual_pump.index
    # df = pd.DataFrame(index=listdaily_pump.index)
    # df = pd.concat([df, pump_df], join="outer",
    #                keys=["Date", "Date"], axis=1)
    # df.columns = df.columns.droplevel()

    # # Interpolating pumping data
    # pump_interp = df.interpolate(method="cubic")
    # pump_interp = pump_interp.dropna()
    # pump_interp = pump_interp.rename(columns={"0": well_name})
    model_plot = besttry_Pastasmodels[well_i].copy()
    # model_plot.del_stressmodel("well")  # Deletes previous pumping
    # EstTotPump_ = ps.StressModel(pump_interp,
    #                              rfunc=ps.Gamma(), name="well",
    #                              settings="well", up=False)
    # model_plot.add_stressmodel(EstTotPump_)

    # Assigns parameters to previous optimal parameters and SD
    for param_i in param_index:

        model_plot.set_parameter(name=param_names[param_i],
                                 initial=ls_sub_m[
                                     model_plot.parameters[
                                         "optimal"].index[param_i]+str(well_i)],
                                 optimal=ls_sub_m[
                                     model_plot.parameters[
                                         "optimal"].index[param_i]+str(well_i)])

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

if pumpexperiment == "pumpingcase1":
    # Folder to save/import graph and model
    modelpath = os.path.abspath("models//bangkok-based//")

# Cyclical pump
elif pumpexperiment == "cyclical":
    # Folder to save/import graph and model
    modelpath = os.path.abspath("models//cyclical//")

# Simple pump
elif pumpexperiment == "simpleexp":
    # Folder to save/import graph and model
    modelpath = os.path.abspath("models//cowboyhat_SUB//")

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
# pump_interp = pd.concat(
#     [pump_interp.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

# Reorder well list to shallow to deep aquifers
# BK, PD, NL, NB
well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
              if y in x]

# pump_interp.columns = well_names

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
                           pump_series=pumptrue,
                           model_path=modelpath, califlag=0,
                           esmdaflag=0, user_models=model_plotls)

# Post process data
sub_total, subv_total, ann_sub, \
    avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                             sub_total,
                                             subv_total,
                                             all_results)

bestsubtry = ann_sub.copy()


# %%###########################################################################
# Init run
###############################################################################

# Initial ls
# Saving variables and values
ls_init = {}

# Saving values from string
for temp_i in range(len(temp)):
    variable_name = temp[temp_i].split(":")[0].strip()
    ls_init[variable_name] = float(
        temp[temp_i].split("=")[1].strip()[:-1])

# Saving model plots for least squares
model_lsinit = []

# For all wells in well nest
for well_i, wells in enumerate(well_names):

    # Name of well as a string
    well_name = wells
    #######################################################################

    # Saving pumping
    # pump_mean_init = []
    # for pump_i in range(n_pump):
    #     pump_mean_init.append(ls_init["pump"+str(pump_i)])

    # mean_ = pd.Series(np.exp(pump_mean_init))
    # mean_.index = annual_pump.index

    # # Isolating pumping data
    # pump_df = pd.DataFrame(mean_, index=annual_pump.index,
    #                        columns=["0"])
    # pump_df.index = annual_pump.index
    # df = pd.DataFrame(index=listdaily_pump.index)
    # df = pd.concat([df, pump_df], join="outer",
    #                keys=["Date", "Date"], axis=1)
    # df.columns = df.columns.droplevel()

    # # Interpolating pumping data
    # pump_interp_init = df.interpolate(method="cubic")
    # pump_interp_init = pump_interp_init.dropna()
    # pump_interp_init = pump_interp_init.rename(columns={"0": well_name})
    model_plot_init = besttry_Pastasmodels[well_i].copy()
    # model_plot_init.del_stressmodel("well")  # Deletes previous pumping
    # EstTotPump_ = ps.StressModel(pump_interp_init,
    #                              rfunc=ps.Gamma(), name="well",
    #                              settings="well", up=False)
    # model_plot_init.add_stressmodel(EstTotPump_)

    # Assigns parameters to previous optimal parameters and SD
    for param_i in param_index:

        model_plot_init.set_parameter(name=param_names[param_i],
                                      initial=ls_init[
                                          model_plot.parameters[
                                              "optimal"].index[param_i]+str(well_i)],
                                      optimal=ls_init[
                                          model_plot.parameters[
                                              "optimal"].index[param_i]+str(well_i)])

    model_lsinit.append(model_plot_init)

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

    Sskv_data.loc[wellnest] *= np.exp(ls_init[p_multop[1]+"0"])

    Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
    Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

    K_data.loc[wellnest] *= np.exp(ls_init[p_multop[1]+"1"])

# Mode can be "raw" as in raw groundwater data vs "Pastas" for importing Pastas
# simulated groundwater in the aquifers
mode = "Pastas"

if pumpexperiment == "pumpingcase1":
    # Folder to save/import graph and model
    modelpath = os.path.abspath("models//bangkok-based//")

# Cyclical pump
elif pumpexperiment == "cyclical":
    # Folder to save/import graph and model
    modelpath = os.path.abspath("models//cyclical//")

# Simple pump
elif pumpexperiment == "simpleexp":
    # Folder to save/import graph and model
    modelpath = os.path.abspath("models//cowboyhat_SUB//")

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
# pump_interp_init = pd.concat(
#     [pump_interp_init.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

# Reorder well list to shallow to deep aquifers
# BK, PD, NL, NB
well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
              if y in x]

# pump_interp_init.columns = well_names

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
                           # pump_series=pump_interp_init,
                           model_path=modelpath, califlag=0,
                           esmdaflag=0, user_models=model_lsinit)

# Post process data
sub_total, subv_total, ann_sub, \
    avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                             sub_total,
                                             subv_total,
                                             all_results)

initsubtry = ann_sub.copy()

# %% 30 regularized least squares

# # save fit report to a file:
# with open(os.path.abspath(
#         lspath + "//" +
#         wellnestlist[0] + "_LSreg30.0_modelresult.txt"),
#         'r') as fh:
#     temp = fh.readlines()
#     temp = [x.replace("\n", "") for x in temp]
#     temp = temp[
#         temp.index("[[Variables]]") + 1:temp.index(
#             "[[Correlations]] (unreported correlations are < 0.100)")]
#     # temp = temp[
#     #     temp.index("[[Variables]]") + 1:]
# fh.close()

# # Saving variables and values
# ls_sub30 = {}

# # Saving values from string
# for temp_i in range(len(temp)):
#     variable_name = temp[temp_i].split(":")[0].strip()

#     # Constant d lines are too long
#     if "constant_d" in variable_name:
#         ls_sub30[variable_name] = float(temp[temp_i][16:28])
#     else:
#         ls_sub30[variable_name] = float(temp[temp_i][14:25])

# # Saving model plots for least squares
# model_plotls30 = []

# # For all wells in well nest
# for well_i, wells in enumerate(well_names):

#     # Name of well as a string
#     well_name = wells
#     #######################################################################

#     # Saving pumping
#     pump_mean30 = []
#     for pump_i in range(n_pump):
#         pump_mean30.append(ls_sub30["pump"+str(pump_i)])

#     mean_ = pd.Series(np.exp(pump_mean30))
#     mean_.index = annual_pump.index

#     # Isolating pumping data
#     pump_df = pd.DataFrame(mean_, index=annual_pump.index,
#                            columns=["0"])
#     pump_df.index = annual_pump.index
#     df = pd.DataFrame(index=listdaily_pump.index)
#     df = pd.concat([df, pump_df], join="outer",
#                    keys=["Date", "Date"], axis=1)
#     df.columns = df.columns.droplevel()

#     # Interpolating pumping data
#     pump_interp30 = df.interpolate(method="cubic")
#     pump_interp30 = pump_interp30.dropna()
#     pump_interp30 = pump_interp30.rename(columns={"0": well_name})
#     model_plot30 = besttry_Pastasmodels[well_i].copy()
#     model_plot30.del_stressmodel("well")  # Deletes previous pumping
#     EstTotPump_ = ps.StressModel(pump_interp30,
#                                  rfunc=ps.Gamma(), name="well",
#                                  settings="well", up=False)
#     model_plot30.add_stressmodel(EstTotPump_)

#     # Assigns parameters to previous optimal parameters and SD
#     for param_i in param_index:

#         model_plot30.set_parameter(name=param_names[param_i],
#                                    initial=ls_sub30[
#                                        model_plot30.parameters[
#                                            "optimal"].index[param_i]+str(well_i)],
#                                    optimal=ls_sub30[
#                                        model_plot30.parameters[
#                                            "optimal"].index[param_i]+str(well_i)])

#     model_plotls30.append(model_plot30)

# # %% BEST SUB TRY

# tmin = "1978"
# tmax = "2020"

# # Reading in thickness and storage data
# path = os.path.join(os.path.abspath("inputs"),
#                     "SUBParametersPriortoManual.xlsx")
# Thick_data = pd.read_excel(path, sheet_name="Thickness",
#                            index_col=0)  # Thickness
# Sskv_data = pd.read_excel(path,
#                           sheet_name="Sskv",
#                           index_col=0)  # Sskv
# Sske_data = pd.read_excel(path,
#                           sheet_name="Sske",
#                           index_col=0)  # Ssk
# K_data = pd.read_excel(path,
#                        sheet_name="K",
#                        index_col=0)  # K

# # Random multipler for each well nest
# for wellnest in wellnestlist:

#     Sskv_data.loc[wellnest] *= np.exp(ls_sub30[p_multop[1]+"0"])

#     Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
#     Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

#     K_data.loc[wellnest] *= np.exp(ls_sub30[p_multop[1]+"1"])

# # Mode can be "raw" as in raw groundwater data vs "Pastas" for importing Pastas
# # simulated groundwater in the aquifers
# mode = "Pastas"

# if pumpexperiment == "pumpingcase1":
#     # Folder to save/import graph and model
#     modelpath = os.path.abspath("models//bangkok-based//")

# # Cyclical pump
# elif pumpexperiment == "cyclical":
#     # Folder to save/import graph and model
#     modelpath = os.path.abspath("models//cyclical//")

# # Simple pump
# elif pumpexperiment == "simpleexp":
#     # Folder to save/import graph and model
#     modelpath = os.path.abspath("models//cowboyhat_SUB//")

# # Pumping flag, for PASTAS, if changing pumping scenario
# pumpflag = 1
# ppath = os.path.join(os.path.abspath("inputs"), "BasinPumping.xlsx")
# psheet = "EstTotalPump_54-60_Int50"

# # Convergence criteria
# CC = 1 * 10**-5

# # Number of nodes in clay
# node_num = 10

# # Using available heads as proxy for missing
# proxyflag = 1

# # Repeat interp for each well
# pump_interp30 = pd.concat(
#     [pump_interp30.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

# # Reorder well list to shallow to deep aquifers
# # BK, PD, NL, NB
# well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
#               if y in x]

# pump_interp30.columns = well_names

# # Calculates subsidence
# all_results, sub_total, subv_total = bkk_sub_gw.\
#     bkk_sub.bkk_subsidence(wellnestlist,
#                            mode, tmin,
#                            tmax,
#                            Thick_data,
#                            K_data,
#                            Sskv_data,
#                            Sske_data,
#                            CC=CC,
#                            Nz=node_num,
#                            ic_run=True,
#                            proxyflag=proxyflag,
#                            pumpflag=pumpflag,
#                            pump_path=ppath,
#                            pump_sheet=psheet,
#                            pump_series=pump_interp30,
#                            model_path=modelpath, califlag=0,
#                            esmdaflag=0, user_models=model_plotls30)

# # Post process data
# sub_total, subv_total, ann_sub, \
#     avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
#                                              sub_total,
#                                              subv_total,
#                                              all_results)

# bestsubtry30 = ann_sub.copy()

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
fig, axs = plt.subplots(num_wells+1, 1, figsize=(16, 16), sharey="row")
for mod_i in range(len(besttry_Pastasmodels)):

    well_name = well_names[mod_i]

    # PLOTTING
    # Getting obs, index and plotting obs and truth
    obs = besttry_Pastasmodels[mod_i].observations()[:"2004"]
    index = besttry_Pastasmodels[mod_i].observations()[:"2004"].index

    # Axes coord
    if mod_i == 0:

        # Axes coord
        row = mod_i
        col = 1
        axs[row].set_xlabel("Date")
        axs[row].annotate("(a)",
                          xy=(.05, .01), xycoords="axes fraction",
                          fontsize=20, horizontalalignment="right",
                          weight="bold",
                          verticalalignment="bottom")

    elif mod_i == 1:

        row = mod_i
        col = 1
        axs[row].set(xticklabels=[])

        axs[row].annotate("(b)",
                          xy=(.1, .01), xycoords="axes fraction",
                          fontsize=20, horizontalalignment="right",
                          weight="bold",
                          verticalalignment="bottom")
    elif mod_i == 2:
        row = mod_i
        col = 1
        axs[row].set(xticklabels=[])

        axs[row].annotate("(c)",
                          xy=(.1, .01), xycoords="axes fraction",
                          fontsize=20, horizontalalignment="right",
                          weight="bold",
                          verticalalignment="bottom")

    elif mod_i == 3:
        row = mod_i
        col = 1

        axs[row].annotate("(d)",
                          xy=(.1, .01), xycoords="axes fraction",
                          fontsize=20, horizontalalignment="right",
                          weight="bold",
                          verticalalignment="bottom")

    axs[row].plot(index, model_plotls[mod_i].simulate()[index],
                  "--", color="fuchsia",
                  linewidth=4, label="Ordinary Least Squares", zorder=10)
    # axs[row].plot(index, model_plotls30[mod_i].simulate()[index],
    #               "-", color="blue",
    #               linewidth=4, label="Regularized Least Squares", zorder=1)
    axs[row].plot(obs, "o", label="Observations", color="black",
                  markersize=4)
    axs[row].plot(index, truth[1].loc[index].iloc[:, -1],
                  "-*c", label="Truth", linewidth=4)
    axs[row].plot(index, model_lsinit[mod_i].simulate()[index],
                  "-", color="red", alpha=.3,
                  linewidth=4, label="Initial", zorder=10)

    # Annotating
    # axs[row].annotate(well_name,
    #                   xy=(.01, 1.1), xycoords="axes fraction",
    #                   fontsize=20, horizontalalignment="center",
    #                   weight="bold",
    #                   bbox=dict(boxstyle="round", fc="0.8"),
    #                   verticalalignment="baseline")

    # add rsq to simulation
    # rmse_data = gw_obs.merge(
    #     model_plotls[mod_i].simulate(),
    #     left_on=gw_obs.index, right_on=model_plotls[mod_i].simulate().index,
    #     how="inner")
    # rmse = mean_squared_error(rmse_data["0"], rmse_data.Simulation,
    #                           squared=False)
    # print("RMSE (OLSvsGWObs): " + "{:.1f}".format(rmse) + " m")
    # rmse_data = gw_obs.merge(
    #     model_plotls[mod_i].simulate(),
    #     left_on=gw_obs.index, right_on=model_plotls30[mod_i].simulate().index,
    #     how="inner")
    # rmse = mean_squared_error(rmse_data["0"], rmse_data.Simulation,
    #                           squared=False)
    # print("RMSE (RLSvsGWObs): " + "{:.1f}".format(rmse) + " m")
    axs[row].set_ylabel("Head (m)")

    # Print observed range
    # obsgw_range = obs.max() - obs.min()
    # rmse_data = truth[0].merge(
    #     model_plotls[mod_i].simulate(),
    #     left_on=truth[0].index, right_on=model_plotls[mod_i].simulate().index,
    #     how="inner")
    # rmse_data.index = rmse_data.Date
    # rmse_data = rmse_data[:"2005"]
    # print("Observed range of head is : (m) " + f'{obsgw_range:.2f}')
    # print("RMSE (OLSvsGWtruth): (m)" + str(mean_squared_error(
    #     rmse_data["0"], rmse_data.Simulation,
    #     squared=False)))

    # rmse_data = truth[1].merge(
    #     model_plotls[mod_i].simulate(),
    #     left_on=truth[1].index, right_on=model_plotls30[mod_i].simulate().index,
    #     how="inner")
    # rmse_data.index = rmse_data.Date
    # rmse_data = rmse_data[:"2005"]
    # print("Observed range of head is : (m) " + f'{obsgw_range:.2f}')
    # print("RMSE (RLSvsGWtruth): (m)" + str(mean_squared_error(
    #     rmse_data["0"], rmse_data.Simulation,
    #     squared=False)))

# obssub_range = sub_obs.iloc[:, -1].max() - sub_obs.iloc[:, -1].min()
# print("Observed range of sub is : (cm/yr) " + f'{obssub_range:.2f}')

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
axs[row].plot(index, obs, "o", label="Observations", color="black",
              markersize=8, zorder=2)
axs[row].plot(index, -truth[0].dropna().loc[index].iloc[:, -1],
              "-*c", label="Truth", linewidth=6)

axs[row].annotate("(" + chr(ord('`')+num_wells+1) + ")",
                  xy=(.05, .01), xycoords="axes fraction",
                  fontsize=20, horizontalalignment="right",
                  weight="bold",
                  verticalalignment="bottom")

# Annotating and labeling
axs[row].set_ylabel("Subsidence Rates (cm/yr)")
axs[row].set_xlabel("Date")

# add rsq to simulation
rmse = mean_squared_error(
    bestsubtry[0][1].iloc[:, 2][index]*100, sub_obs.iloc[:, -1], squared=False)

print("RMSE (OLSvsSubObs): " + "{:.1f}".format(rmse) + " cm/yr",)

# add rsq to simulation
# rmse = mean_squared_error(
#     bestsubtry30[0][1].iloc[:, 2][index]*100, sub_obs.iloc[:, -1], squared=False)

# print("RMSE (RLSvsSubObs): " + "{:.1f}".format(rmse) + " cm/yr",)

axs[row].set_zorder(1)  # default zorder is 0 for ax1 and ax2
axs[row].patch.set_visible(False)

# axs[row].annotate(Wellnest_name[0][2:],
#                   xy=(.01, 1.2), xycoords="axes fraction",
#                   fontsize=20, horizontalalignment="center",
#                   weight="bold",
#                   bbox=dict(boxstyle="round", fc="0.8"),
#                   verticalalignment="baseline")

# access legend objects automatically created from data
handles, labels = axs[row].get_legend_handles_labels()

# Adding lines for ensemble members
line = Line2D([0], [0], label='Initial', color='tab:red',
              alpha=1, linewidth=4)
linels = Line2D([0], [0], label="Ordinary Least Squares", color='fuchsia',
                linestyle="--",
                alpha=1, linewidth=4)
linels30 = Line2D([0], [0], label='Regularized Least Squares', color='blue',
                  linestyle="-",
                  alpha=1, linewidth=4)

# add manual symbols to auto legend
handles.insert(1, line)
handles.insert(2, linels)
handles.insert(3, linels30)

# Plot truth, least squares, observations
axs[row].plot(index, -truth[0].dropna().loc[index].iloc[:, -1],
              "-*c", label="Truth", linewidth=4)
axs[row].plot(bestsubtry[0][1].iloc[:, 2][index].index,
              -bestsubtry[0][1].iloc[:, 2][index]*100, "--",
              color="fuchsia",
              label="Ordinary Least Squares", linewidth=4)
# axs[row].plot(bestsubtry30[0][1].iloc[:, 2][index].index,
#               -bestsubtry30[0][1].iloc[:, 2][index]*100, "-",
#               color="blue", zorder=1,
#               label="Regularized Least Squares", linewidth=4)
axs[row].plot(initsubtry[0][1].iloc[:, 2][index].index,
              -initsubtry[0][1].iloc[:, 2][index]*100, "-",
              color="red", alpha=.3,
              label="Initial", linewidth=4)
axs[row].plot(index, obs, "o", label="Observations", color="black",
              markersize=8, zorder=10)

# FORMATTING
# axs[row].set_ylim([-5, 20])
# Shrink current axis by 20%
box = axs[row].get_position()

# Put a legend to the right of the current axis
# axs[row].legend(handles=handles, loc='center left', bbox_to_anchor=(.7, .3))

# fig.subplots_adjust(hspace=0.5)
fig.set_rasterized(True)
# axs[row].xaxis.set_major_locator(YearLocator(30))
axs[row].xaxis.set_major_formatter(DateFormatter("%Y"))

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


plt.rc("xtick", labelsize=13)  # fontsize of the x tick labels


# Number of pastas parameters
n_param = len(param_index)
labels = ["Initial", "Ordinary\nLeast\nSquares", "Truth"]
bar_colors = ['red', 'fuchsia', 'cyan']
bar_labels = ["Initial", "Ordinary Least Squares", "Truth"]
fig, axs = plt.subplots(2, 3, figsize=(20, 10))
# fig, axs = plt.subplots(2, 2, figsize=(20, 10))
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
            axs[row, col].set_title('A', fontsize=23, weight="bold")
            values = [ls_init["well_A"+str(mod_i)],
                      ls_sub_m["well_A"+str(mod_i)],
                      # ls_sub30["well_A"+str(mod_i)],
                      Atrue]
            axs[row, col].bar(labels, values, label=bar_labels, color=bar_colors)
            axs[row, col].errorbar(labels, values,
                                   yerr=[np.nan, ls_error[param_i], np.nan],
                                   fmt="o", color="k")
            axs[row, col].set(ylabel="Value")
            axs[row, col].annotate("(a)",
                                   xy=(.1, 1.01), xycoords="axes fraction",
                                   fontsize=20, horizontalalignment="right",
                                   weight="bold",
                                   verticalalignment="bottom")
            truthval = Atrue
            # mean_ = ls_sub30["well_A"+str(mod_i)]
            # percerr = abs(abs(mean_ - truthval)/truthval) * 100
            # print(mean_, truthval)
            # print("Param " + str(param_i) + ": " + f'{percerr:.2f}')

        # n
        elif param_i == 1:
            row = 0
            col = 1
            axs[row, col].set_title('n', fontsize=23, weight="bold")
            values = [ls_init["well_n"+str(mod_i)],
                      ls_sub_m["well_n"+str(mod_i)],
                      # ls_sub30["well_n"+str(mod_i)],
                      ntrue]
            axs[row, col].bar(labels, values, label=bar_labels, color=bar_colors)
            axs[row, col].errorbar(labels, values,
                                   yerr=[np.nan, ls_error[param_i], np.nan],
                                   fmt="o", color="k")
            axs[row, col].set(ylabel="Value")
            axs[row, col].annotate("(b)",
                                   xy=(.1, 1.01), xycoords="axes fraction",
                                   fontsize=20, horizontalalignment="right",
                                   weight="bold",
                                   verticalalignment="bottom")
            truthval = ntrue
            # mean_ = ls_sub30["well_n"+str(mod_i)]
            # percerr = abs(abs(mean_ - truthval)/truthval) * 100
            # print(mean_, truthval)
            # print("Param " + str(param_i) + ": " + f'{percerr:.2f}')

        # a
        elif param_i == 2:
            row = 1
            col = 0
            axs[row, col].set_title('a', fontsize=23, weight="bold")
            values = [ls_init["well_a"+str(mod_i)],
                      ls_sub_m["well_a"+str(mod_i)],
                      # ls_sub30["well_a"+str(mod_i)],
                      atrue]
            axs[row, col].bar(labels, values, label=bar_labels, color=bar_colors)
            axs[row, col].errorbar(labels, values,
                                   yerr=[np.nan, ls_error[param_i], np.nan],
                                   fmt="o", color="k")
            axs[row, col].set(ylabel="Value")
            axs[row, col].annotate("(d)",
                                   xy=(.1, 1.01), xycoords="axes fraction",
                                   fontsize=20, horizontalalignment="right",
                                   weight="bold",
                                   verticalalignment="bottom")
            truthval = atrue
            # Print percent error
            # mean_ = ls_sub30["well_a"+str(mod_i)]
            # percerr = abs(abs(mean_ - truthval)/truthval) * 100
            # print(mean_, truthval)
            # print("Param " + str(param_i) + ": " + f'{percerr:.2f}')

        # d
        else:
            row = 1
            col = 1
            axs[row, col].set_title('d', fontsize=23, weight="bold")
            values = [ls_init["constant_d"+str(mod_i)],
                      ls_sub_m["constant_d"+str(mod_i)],
                      # ls_sub30["constant_d"+str(mod_i)],
                      dtrue]
            axs[row, col].bar(labels, values, label=bar_labels, color=bar_colors)
            axs[row, col].errorbar(labels, values,
                                   yerr=[np.nan, ls_error[param_i], np.nan],
                                   fmt="o", color="k")
            axs[row, col].axhline(0, color="k", linewidth=1)
            axs[row, col].set(ylabel="Value")
            axs[row, col].annotate("(e)",
                                   xy=(.1, 1.01), xycoords="axes fraction",
                                   fontsize=20, horizontalalignment="right",
                                   weight="bold",
                                   verticalalignment="bottom")
            truthval = dtrue
            # Print percent error
            # mean_ = ls_sub30["constant_d"+str(mod_i)]
            # percerr = abs(abs(mean_ - truthval)/truthval) * 100
            # print(mean_, truthval)
            # print("Param " + str(param_i) + ": " + f'{percerr:.2f}')

# For each sub param
for param_i in range(n_sub):

    # Truth val, title, save title
    if p_multop[1] == "Sskv" or p_multop[1] == "SsK":

        # Sskv
        if param_i == 0:

            # plt.title('Sskv', fontsize=14)
            truthval = a
            axs[param_i, 2].set_title(
                'S$_{skv}$ Multiplier', fontsize=23, weight="bold")
            values = [np.exp(ls_init["SsK0"]),
                      np.exp(ls_sub_m["SsK0"]),
                      # np.exp(ls_sub30["SsK0"]),
                      truthval]
            axs[param_i, 2].axhline(0, color="k", linewidth=1)
            axs[param_i, 2].set(ylabel="Value")
            axs[param_i, 2].bar(labels, values, label=bar_labels, color=bar_colors)
            axs[param_i, 2].errorbar(labels, values,
                                     yerr=[np.nan, ls_error[param_i], np.nan],
                                     fmt="o", color="k")
            axs[param_i, 2].annotate("(c)",
                                     xy=(.1, 1.01), xycoords="axes fraction",
                                     fontsize=20, horizontalalignment="right",
                                     weight="bold",
                                     verticalalignment="bottom")
            axs[param_i, 2].annotate(f'{np.exp(ls_init["SsK0"]):.2f}',
                                     (0.07, 0.05), xycoords='axes fraction',
                                     fontsize=22)
            # mean_ = np.exp(ls_sub30["SsK0"])
            # percerr = abs(abs(mean_ - truthval)/truthval) * 100
            # print(mean_, truthval)
            # print("Param " + str(param_i) + ": " + f'{percerr:.2f}')

        # Sske
        elif param_i == 1:
            # plt.suptitle("ESMDA: LOG K Parameters for Well")

            # plt.title('K', fontsize=14)
            truthval = c
            axs[param_i, 2].set_title('K Multiplier', fontsize=23, weight="bold")

            values = [np.exp(ls_init["SsK1"]),
                      np.exp(ls_sub_m["SsK1"]),
                      # np.exp(ls_sub30["SsK1"]),
                      truthval]
            axs[param_i, 2].set(ylabel="Value")
            bars = axs[param_i, 2].bar(labels, values, label=bar_labels, color=bar_colors)
            axs[param_i, 2].errorbar(labels, values,
                                     yerr=[np.nan, ls_error[param_i], np.nan],
                                     fmt="o", color="k")

            # axs[param_i, 2].annotate(f'{np.exp(ls_sub30["SsK1"]):.2f}',
            #                          (0.54, 0.05), xycoords='axes fraction',
            #                          fontsize=22)
            axs[param_i, 2].annotate("(f)",
                                     xy=(.1, 1.01), xycoords="axes fraction",
                                     fontsize=20, horizontalalignment="right",
                                     weight="bold",
                                     verticalalignment="bottom")
            # mean_ = np.exp(ls_sub30["SsK1"])
            # percerr = abs(abs(mean_ - truthval)/truthval) * 100
            # print(mean_, truthval)
            # print("Param " + str(param_i) + ": " + f'{percerr:.2f}')

plt.show()
plt.subplots_adjust(top=0.99, bottom=0.01, hspace=.5)
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

# plt.plot(obs_index,
#          -sub_obs.iloc[:, 1], "-ok", label="Observations", zorder=10)
# plt.plot(obs_index,
#          -truth[0].dropna().loc[obs_index].iloc[:, -1], "-c", label="Truth", linewidth=5)
# plt.plot(bestsubtry[0][1].iloc[:, 2][index].index,
#          -bestsubtry[0][1].iloc[:, 2][index]*100, "--", color="fuchsia",
#          label="Ordinary Least Squares", linewidth=6)
# plt.plot(bestsubtry30[0][1].iloc[:, 2][index].index,
#          -bestsubtry30[0][1].iloc[:, 2][index]*100, "--", color="blue",
#          label="Regularized Least Squares", linewidth=6)

# plt.ylabel("Rates\n(cm/yr)")
# plt.legend()
# fig_name1 = Wellnest_name[0] + "Sub_ESMDA_evol.png"
# full_figpath = os.path.join(figpath, fig_name1)
# plt.savefig(full_figpath, bbox_inches="tight", format="png")


# %% Pumping during obs

# fig, axs = plt.subplots(1, 1, figsize=(16, 5))

# # Pumping subset
# obs_index = gw_obs_indices[1]
# pump_index = pumptrue[
#     np.logical_and(
#         pumptrue.index >= obs_index[0], pumptrue.index <= obs_index.iloc[-1])]
# pump_indices = np.logical_and(
#     pumptrue.index >= obs_index[0], pumptrue.index <= obs_index.iloc[-1])

# pump_index2 = annual_pump[
#     np.logical_and(
#         annual_pump.index >= obs_index[0], annual_pump.index <= obs_index.iloc[-1])]
# pump_indices2 = np.logical_and(
#     annual_pump.index >= obs_index[0], annual_pump.index <= obs_index.iloc[-1])

# axs.plot(pump_index.index, pumptrue[pump_indices].iloc[:, -1],
#          color="c", linewidth=4,
#          marker="o")
# axs.plot(pump_index2.index,
#          list(compress(np.exp(pump_mean), pump_indices2)), color="fuchsia",
#          alpha=1, linewidth=4,
#          linestyle="--")
# # axs.plot(pump_index2.index,
# #          list(compress(np.exp(pump_mean30), pump_indices2)), color="blue",
# #          alpha=1, linewidth=4,
# #          linestyle="-")
# axs.plot(pump_index2.index,
#          list(compress(np.exp(pump_mean_init), pump_indices2)), color="red",
#          alpha=.3, linewidth=4,
#          linestyle="-")
# axs.set_ylabel("Pumping Rate * 10$^4$ (m$^3$/day)")
# axs.set_xlabel("Date")

# # access legend objects automatically created from data
# handles, labels = axs.get_legend_handles_labels()

# # Adding lines for ensemble members
# linetruth = Line2D([0], [0], label='Truth', color='c',
#                    alpha=1, linewidth=4, marker='o')
# linepriormean = Line2D([0], [0], label='Initial', color='red',
#                        alpha=1, linewidth=4, marker='o')
# lineleast = Line2D([0], [0], label="Ordinary Least Squares", color='fuchsia',
#                    alpha=1, linewidth=4, linestyle="--")
# lineleast30 = Line2D([0], [0], label='Regularized Least Squares', color='blue',
#                      alpha=1, linewidth=4, linestyle="-")

# # plt.ylim([-30, 0])
# # add manual symbols to auto legend
# handles.extend([linepriormean, lineleast, lineleast30, linetruth])

# axs.legend(handles=handles)
# fig.set_rasterized(True)

# # Saving
# fig_name1 = wellnestlist[0] + "_ObsPumping_ESMDA_Paper.png"
# full_figpath = os.path.join(figpath, fig_name1)
# plt.savefig(full_figpath, bbox_inches="tight", format="png")

# fig_name1 = wellnestlist[0] + "_ObsPumping_ESMDA_Paper.eps"
# full_figpath = os.path.join(figpath, fig_name1)
# plt.savefig(full_figpath, bbox_inches="tight", format="eps")

# # %%###########################################################################
# # Plotting settings
# ###############################################################################

# plt.rc("font", size=28)  # controls default text size
# plt.rc("axes", titlesize=24)  # fontsize of the title
# plt.rc("axes", labelsize=18)  # fontsize of the x and y labels
# plt.rc("xtick", labelsize=16)  # fontsize of the x tick labels
# plt.rc("ytick", labelsize=16)  # fontsize of the y tick labels
# plt.rc("legend", fontsize=14)  # fontsize of the legend


# # %% Pumping during obs

# fig, axs = plt.subplots(1, 1, figsize=(16, 5))

# # Pumping subset
# obs_index = gw_obs_indices[1]
# pump_index = pumptrue[
#     np.logical_and(
#         pumptrue.index >= obs_index[0], pumptrue.index <= obs_index.iloc[-1])]
# pump_indices = np.logical_and(
#     pumptrue.index >= obs_index[0], pumptrue.index <= obs_index.iloc[-1])

# pump_index2 = annual_pump[
#     np.logical_and(
#         annual_pump.index >= obs_index[0], annual_pump.index <= obs_index.iloc[-1])]
# pump_indices2 = np.logical_and(
#     annual_pump.index >= obs_index[0], annual_pump.index <= obs_index.iloc[-1])

# axs.plot(pump_index.index, -dtrue+(-Atrue)*pumptrue[pump_indices].iloc[:, -1],
#          color="c", linewidth=4,
#          marker="o")

# # A and d least squares
# Als = -ls_sub_m[
#     model_plot.parameters[
#         "optimal"].index[0]+str(well_i)]
# dls = ls_sub_m[
#     model_plot.parameters[
#         "optimal"].index[3]+str(well_i)]

# # A and d initial
# Als_init = -ls_init[
#     model_plot.parameters[
#         "optimal"].index[0]+str(well_i)]
# dls_init = ls_init[
#     model_plot.parameters[
#         "optimal"].index[3]+str(well_i)]

# # A and d regularized least squares
# # Als30 = -ls_sub30[
# #     model_plot.parameters[
# #         "optimal"].index[0]+str(well_i)]
# # dls30 = ls_sub30[
# #     model_plot.parameters[
# #         "optimal"].index[3]+str(well_i)]

# # Plotting
# axs.plot(pump_index2.index,
#          -dls+Als*pd.Series(list(compress(
#              np.exp(pump_mean), pump_indices2))), color="fuchsia",
#          alpha=1, linewidth=4,
#          linestyle="--")
# # axs.plot(pump_index2.index,
# #          -dls30+Als30*pd.Series(list(compress(
# #              np.exp(pump_mean30), pump_indices2))), color="blue",
# #          alpha=1, linewidth=4,
# #          linestyle="-")
# axs.plot(pump_index2.index,
#          -dls_init+Als_init*pd.Series(list(compress(
#              np.exp(pump_mean_init), pump_indices2))), color="red",
#          alpha=.3, linewidth=4,
#          linestyle="-")
# axs.set_ylabel("A*Pump-d [m]")
# axs.set_xlabel("Date")

# # access legend objects automatically created from data
# handles, labels = axs.get_legend_handles_labels()

# # Adding lines for ensemble members
# linetruth = Line2D([0], [0], label='Truth', color='c',
#                    alpha=1, linewidth=4, marker='o')
# linepriormean = Line2D([0], [0], label='Initial', color='red',
#                        alpha=1, linewidth=4, marker='o')
# lineleast = Line2D([0], [0], label="Ordinary Least Squares", color='fuchsia',
#                    alpha=1, linewidth=4, linestyle="--")
# lineleast30 = Line2D([0], [0], label='Regularized Least Squares', color='blue',
#                      alpha=1, linewidth=4, linestyle="-")

# # plt.ylim([-30, 0])
# # add manual symbols to auto legend
# handles.extend([linepriormean, lineleast, lineleast30, linetruth])

# axs.legend(handles=handles)
# fig.set_rasterized(True)

# # Saving
# fig_name1 = wellnestlist[0] + "_-AxPump.png"
# full_figpath = os.path.join(figpath, fig_name1)
# plt.savefig(full_figpath, bbox_inches="tight", format="png")

# fig_name1 = wellnestlist[0] + "_-AxPump.eps"
# full_figpath = os.path.join(figpath, fig_name1)
# plt.savefig(full_figpath, bbox_inches="tight", format="eps")

# # %% Norm Pumping during obs

# fig, axs = plt.subplots(1, 1, figsize=(16, 5))

# # Pumping subset
# obs_index = gw_obs_indices[1]

# pump_index2 = annual_pump[
#     np.logical_and(
#         annual_pump.index >= obs_index[0], annual_pump.index <= obs_index.iloc[-1])]
# pump_indices2 = np.logical_and(
#     annual_pump.index >= obs_index[0], annual_pump.index <= obs_index.iloc[-1])

# data = pumptrue[pump_indices].iloc[:, -1]
# min_data = np.min(data)
# max_data = np.max(data)

# norm_data = (data - min_data)/(max_data - min_data)
# axs.plot(pump_index.index, norm_data,
#          color="c", linewidth=4,
#          marker="o")

# # Ordinary least squares
# data = pd.Series(list(compress(np.exp(pump_mean), pump_indices2)))
# min_data = np.min(data)
# max_data = np.max(data)

# norm_data = (data - min_data)/(max_data - min_data)
# axs.plot(annual_pump.index[
#     pump_indices2], norm_data, color="fuchsia", alpha=1, linewidth=4, linestyle="--")

# # Reg least squares
# # data = pd.Series(list(compress(np.exp(pump_mean30), pump_indices2)))
# # min_data = np.min(data)
# # max_data = np.max(data)

# # norm_data = (data - min_data)/(max_data - min_data)
# # axs.plot(annual_pump.index[
# #     pump_indices2], norm_data, color="blue", alpha=1, linewidth=4, linestyle="-")

# # Initial values
# data = pd.Series(list(compress(np.exp(pump_mean_init), pump_indices2)))
# min_data = np.min(data)
# max_data = np.max(data)

# norm_data = (data - min_data)/(max_data - min_data)
# axs.plot(annual_pump.index[
#     pump_indices2], norm_data, color="red", alpha=.3, linewidth=4, linestyle="-")
# plt.show()
# axs.set_ylabel("Normalized Pumping")
# axs.set_xlabel("Date")

# # access legend objects automatically created from data
# handles, labels = axs.get_legend_handles_labels()

# # Adding lines for ensemble members
# linepriormean = Line2D([0], [0], label='Initial', color='red',
#                        alpha=.3, linewidth=4)
# lineleast = Line2D([0], [0], label="Ordinary Least Squares", color='fuchsia',
#                    alpha=1, linewidth=4, linestyle="--")
# lineRLS = Line2D([0], [0], label="Regularized Least Squares", color='blue',
#                  alpha=1, linewidth=4, linestyle="-")
# linetruth = Line2D([0], [0], label='Truth', color='cyan',
#                    alpha=.8, linewidth=4)

# # plt.ylim([-30, 0])
# # add manual symbols to auto legend
# handles.extend([linepriormean, lineleast, lineRLS, linetruth])

# axs.legend(handles=handles)
# fig.set_rasterized(True)

# # Saving
# fig_name1 = wellnestlist[0] + "_normPump.png"
# full_figpath = os.path.join(figpath, fig_name1)
# plt.savefig(full_figpath, bbox_inches="tight", format="png")

# fig_name1 = wellnestlist[0] + "_normPump.eps"
# full_figpath = os.path.join(figpath, fig_name1)
# plt.savefig(full_figpath, bbox_inches="tight", format="eps")
