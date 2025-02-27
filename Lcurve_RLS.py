# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 13:54:55 2025

@author: jtsoonthornran

L shaped curve for lambda determination (regularized least squares)

"""

# ##############################################################################

###############################################################################
# import statements
###############################################################################

import os
import pandas as pd
import random
import numpy as np
import matplotlib.pyplot as plt
import pastas as ps
import warnings

# Bangkok Subsidence Model Package
import bkk_sub_gw

# Ignoring Pastas warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

# Changing current directory to locaiton of python script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

pumpexperiment = "pumpingcase1"

# True parameters for Pastas and subsidence mult (a, b, c)
Atrue = -.1
ntrue = 1.2
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
            "models//cowboyhat//")

# Path to save models
tot_path = os.path.abspath("inputs")

# FIgure path
if pumpexperiment == "simpleexp":
    # Folder to save/import graph and model
    figpath = os.path.abspath(
        "figures//cowboyhat//")

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

    param_names = ml.parameters.index.values

# %% LS Import

if pumpexperiment == "simpleexp":
    # Folder to save/import graph and model
    lspath = os.path.abspath("models//cowboyhat//")

elif pumpexperiment == "cyclical":
    # Folder to save/import graph and model
    lspath = os.path.abspath("models//cyclical//")

elif pumpexperiment == "pumpingcase1":
    # Folder to save/import graph and model
    lspath = os.path.abspath("models//bangkok-based//")

lsreg_gw = []
lsreg_sub = []
lsreg_pump = []

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

# Saving values from string
for temp_i in range(len(temp)):
    variable_name = temp[temp_i].split(":")[0].strip()

    # Constant d lines are too long
    if "constant_d" in variable_name:
        ls_sub_m[variable_name] = float(temp[temp_i][16:28])
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

    mean_ = pd.Series(np.exp(pump_mean))
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
    model_plot = besttry_Pastasmodels[well_i].copy()
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
    modelpath = os.path.abspath("models//cowboyhat//")

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

bestsubtry = ann_sub.copy()

lsreg_sub.append(bestsubtry)
lsreg_gw.append(model_plotls)
lsreg_pump.append(pd.Series(np.exp(pump_mean)))

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
    pump_mean_init = []
    for pump_i in range(n_pump):
        pump_mean_init.append(ls_init["pump"+str(pump_i)])

    mean_ = pd.Series(np.exp(pump_mean_init))
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
    pump_interp_init = df.interpolate(method="cubic")
    pump_interp_init = pump_interp_init.dropna()
    pump_interp_init = pump_interp_init.rename(columns={"0": well_name})
    model_plot_init = besttry_Pastasmodels[well_i].copy()
    model_plot_init.del_stressmodel("well")  # Deletes previous pumping
    EstTotPump_ = ps.StressModel(pump_interp_init,
                                 rfunc=ps.Gamma(), name="well",
                                 settings="well", up=False)
    model_plot_init.add_stressmodel(EstTotPump_)

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
    modelpath = os.path.abspath("models//cowboyhat//")

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
pump_interp_init = pd.concat(
    [pump_interp_init.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

# Reorder well list to shallow to deep aquifers
# BK, PD, NL, NB
well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
              if y in x]

pump_interp_init.columns = well_names

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
                           pump_series=pump_interp_init,
                           model_path=modelpath, califlag=0,
                           esmdaflag=0, user_models=model_lsinit)

# Post process data
sub_total, subv_total, ann_sub, \
    avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                             sub_total,
                                             subv_total,
                                             all_results)

initsubtry = ann_sub.copy()

if pumpexperiment == "simpleexp":
    # %% 10 regularized least squares

    # save fit report to a file:
    with open(os.path.abspath(
            lspath + "//" +
            wellnestlist[0] + "_LSreg10.0_modelresult.txt"),
            'r') as fh:
        temp = fh.readlines()
        temp = [x.replace("\n", "") for x in temp]
        # temp = temp[
        #     temp.index("[[Variables]]") + 1:temp.index(
        #         "[[Correlations]] (unreported correlations are < 0.100)")]
        temp = temp[
            temp.index("[[Variables]]") + 1:]
    fh.close()

    # Saving variables and values
    ls_sub10 = {}

    # Saving values from string
    for temp_i in range(len(temp)):
        variable_name = temp[temp_i].split(":")[0].strip()

        # Constant d lines are too long
        if "constant_d" in variable_name:
            ls_sub10[variable_name] = float(temp[temp_i][16:28])
        else:
            ls_sub10[variable_name] = float(temp[temp_i][14:25])

    # Saving model plots for least squares
    model_plotls10 = []

    # For all wells in well nest
    for well_i, wells in enumerate(well_names):

        # Name of well as a string
        well_name = wells
        #######################################################################

        # Saving pumping
        pump_mean10 = []
        for pump_i in range(n_pump):
            pump_mean10.append(ls_sub10["pump"+str(pump_i)])

        mean_ = pd.Series(np.exp(pump_mean10))
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
        pump_interp10 = df.interpolate(method="cubic")
        pump_interp10 = pump_interp10.dropna()
        pump_interp10 = pump_interp10.rename(columns={"0": well_name})
        model_plot10 = besttry_Pastasmodels[well_i].copy()
        model_plot10.del_stressmodel("well")  # Deletes previous pumping
        EstTotPump_ = ps.StressModel(pump_interp10,
                                     rfunc=ps.Gamma(), name="well",
                                     settings="well", up=False)
        model_plot10.add_stressmodel(EstTotPump_)

        # Assigns parameters to previous optimal parameters and SD
        for param_i in param_index:

            model_plot10.set_parameter(name=param_names[param_i],
                                       initial=ls_sub10[
                                           model_plot10.parameters[
                                               "optimal"].index[param_i]+str(well_i)],
                                       optimal=ls_sub10[
                                           model_plot10.parameters[
                                               "optimal"].index[param_i]+str(well_i)])

        model_plotls10.append(model_plot10)

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

        Sskv_data.loc[wellnest] *= np.exp(ls_sub10[p_multop[1]+"0"])

        Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
        Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

        K_data.loc[wellnest] *= np.exp(ls_sub10[p_multop[1]+"1"])

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
        modelpath = os.path.abspath("models//cowboyhat//")

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
    pump_interp10 = pd.concat(
        [pump_interp10.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

    # Reorder well list to shallow to deep aquifers
    # BK, PD, NL, NB
    well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
                  if y in x]

    pump_interp10.columns = well_names

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
                               pump_series=pump_interp10,
                               model_path=modelpath, califlag=0,
                               esmdaflag=0, user_models=model_plotls10)

    # Post process data
    sub_total, subv_total, ann_sub, \
        avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                                 sub_total,
                                                 subv_total,
                                                 all_results)

    bestsubtry10 = ann_sub.copy()

    lsreg_sub.append(bestsubtry10)
    lsreg_gw.append(model_plotls10)
    lsreg_pump.append(pd.Series(np.exp(pump_mean10)))

    # %% 20 regularized least squares

    # save fit report to a file:
    with open(os.path.abspath(
            lspath + "//" +
            wellnestlist[0] + "_LSreg20.0_modelresult.txt"),
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
    ls_sub20 = {}

    # Saving values from string
    for temp_i in range(len(temp)):
        variable_name = temp[temp_i].split(":")[0].strip()

        # Constant d lines are too long
        if "constant_d" in variable_name:
            ls_sub20[variable_name] = float(temp[temp_i][16:28])
        else:
            ls_sub20[variable_name] = float(temp[temp_i][14:25])

    # Saving model plots for least squares
    model_plotls20 = []

    # For all wells in well nest
    for well_i, wells in enumerate(well_names):

        # Name of well as a string
        well_name = wells
        #######################################################################

        # Saving pumping
        pump_mean20 = []
        for pump_i in range(n_pump):
            pump_mean20.append(ls_sub20["pump"+str(pump_i)])

        mean_ = pd.Series(np.exp(pump_mean20))
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
        pump_interp20 = df.interpolate(method="cubic")
        pump_interp20 = pump_interp20.dropna()
        pump_interp20 = pump_interp20.rename(columns={"0": well_name})
        model_plot20 = besttry_Pastasmodels[well_i].copy()
        model_plot20.del_stressmodel("well")  # Deletes previous pumping
        EstTotPump_ = ps.StressModel(pump_interp20,
                                     rfunc=ps.Gamma(), name="well",
                                     settings="well", up=False)
        model_plot20.add_stressmodel(EstTotPump_)

        # Assigns parameters to previous optimal parameters and SD
        for param_i in param_index:

            model_plot20.set_parameter(name=param_names[param_i],
                                       initial=ls_sub20[
                                           model_plot20.parameters[
                                               "optimal"].index[param_i]+str(well_i)],
                                       optimal=ls_sub20[
                                           model_plot20.parameters[
                                               "optimal"].index[param_i]+str(well_i)])

        model_plotls20.append(model_plot20)

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

        Sskv_data.loc[wellnest] *= np.exp(ls_sub20[p_multop[1]+"0"])

        Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
        Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

        K_data.loc[wellnest] *= np.exp(ls_sub20[p_multop[1]+"1"])

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
        modelpath = os.path.abspath("models//cowboyhat//")

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
    pump_interp20 = pd.concat(
        [pump_interp20.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

    # Reorder well list to shallow to deep aquifers
    # BK, PD, NL, NB
    well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
                  if y in x]

    pump_interp20.columns = well_names

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
                               pump_series=pump_interp20,
                               model_path=modelpath, califlag=0,
                               esmdaflag=0, user_models=model_plotls20)

    # Post process data
    sub_total, subv_total, ann_sub, \
        avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                                 sub_total,
                                                 subv_total,
                                                 all_results)

    bestsubtry20 = ann_sub.copy()

    lsreg_sub.append(bestsubtry20)
    lsreg_gw.append(model_plotls20)
    lsreg_pump.append(pd.Series(np.exp(pump_mean20)))

    # %% 25 regularized least squares

    # save fit report to a file:
    with open(os.path.abspath(
            lspath + "//" +
            wellnestlist[0] + "_LSreg25.0_modelresult.txt"),
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
    ls_sub25 = {}

    # Saving values from string
    for temp_i in range(len(temp)):
        variable_name = temp[temp_i].split(":")[0].strip()

        # Constant d lines are too long
        if "constant_d" in variable_name:
            ls_sub25[variable_name] = float(temp[temp_i][16:28])
        else:
            ls_sub25[variable_name] = float(temp[temp_i][14:25])

    # Saving model plots for least squares
    model_plotls25 = []

    # For all wells in well nest
    for well_i, wells in enumerate(well_names):

        # Name of well as a string
        well_name = wells
        #######################################################################

        # Saving pumping
        pump_mean25 = []
        for pump_i in range(n_pump):
            pump_mean25.append(ls_sub25["pump"+str(pump_i)])

        mean_ = pd.Series(np.exp(pump_mean25))
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
        pump_interp25 = df.interpolate(method="cubic")
        pump_interp25 = pump_interp25.dropna()
        pump_interp25 = pump_interp25.rename(columns={"0": well_name})
        model_plot25 = besttry_Pastasmodels[well_i].copy()
        model_plot25.del_stressmodel("well")  # Deletes previous pumping
        EstTotPump_ = ps.StressModel(pump_interp25,
                                     rfunc=ps.Gamma(), name="well",
                                     settings="well", up=False)
        model_plot25.add_stressmodel(EstTotPump_)

        # Assigns parameters to previous optimal parameters and SD
        for param_i in param_index:

            model_plot25.set_parameter(name=param_names[param_i],
                                       initial=ls_sub25[
                                           model_plot25.parameters[
                                               "optimal"].index[param_i]+str(well_i)],
                                       optimal=ls_sub25[
                                           model_plot25.parameters[
                                               "optimal"].index[param_i]+str(well_i)])

        model_plotls25.append(model_plot25)

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

        Sskv_data.loc[wellnest] *= np.exp(ls_sub25[p_multop[1]+"0"])

        Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
        Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

        K_data.loc[wellnest] *= np.exp(ls_sub25[p_multop[1]+"1"])

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
        modelpath = os.path.abspath("models//cowboyhat//")

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
    pump_interp25 = pd.concat(
        [pump_interp25.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

    # Reorder well list to shallow to deep aquifers
    # BK, PD, NL, NB
    well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
                  if y in x]

    pump_interp25.columns = well_names

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
                               pump_series=pump_interp25,
                               model_path=modelpath, califlag=0,
                               esmdaflag=0, user_models=model_plotls25)

    # Post process data
    sub_total, subv_total, ann_sub, \
        avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                                 sub_total,
                                                 subv_total,
                                                 all_results)

    bestsubtry25 = ann_sub.copy()

    lsreg_sub.append(bestsubtry25)
    lsreg_gw.append(model_plotls25)
    lsreg_pump.append(pd.Series(np.exp(pump_mean25)))

    # %% 30 regularized least squares

    # save fit report to a file:
    with open(os.path.abspath(
            lspath + "//" +
            wellnestlist[0] + "_LSreg30.0_modelresult.txt"),
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
    ls_sub30 = {}

    # Saving values from string
    for temp_i in range(len(temp)):
        variable_name = temp[temp_i].split(":")[0].strip()

        # Constant d lines are too long
        if "constant_d" in variable_name:
            ls_sub30[variable_name] = float(temp[temp_i][16:28])
        else:
            ls_sub30[variable_name] = float(temp[temp_i][14:25])

    # Saving model plots for least squares
    model_plotls30 = []

    # For all wells in well nest
    for well_i, wells in enumerate(well_names):

        # Name of well as a string
        well_name = wells
        #######################################################################

        # Saving pumping
        pump_mean30 = []
        for pump_i in range(n_pump):
            pump_mean30.append(ls_sub30["pump"+str(pump_i)])

        mean_ = pd.Series(np.exp(pump_mean30))
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
        pump_interp30 = df.interpolate(method="cubic")
        pump_interp30 = pump_interp30.dropna()
        pump_interp30 = pump_interp30.rename(columns={"0": well_name})
        model_plot30 = besttry_Pastasmodels[well_i].copy()
        model_plot30.del_stressmodel("well")  # Deletes previous pumping
        EstTotPump_ = ps.StressModel(pump_interp30,
                                     rfunc=ps.Gamma(), name="well",
                                     settings="well", up=False)
        model_plot30.add_stressmodel(EstTotPump_)

        # Assigns parameters to previous optimal parameters and SD
        for param_i in param_index:

            model_plot30.set_parameter(name=param_names[param_i],
                                       initial=ls_sub30[
                                           model_plot30.parameters[
                                               "optimal"].index[param_i]+str(well_i)],
                                       optimal=ls_sub30[
                                           model_plot30.parameters[
                                               "optimal"].index[param_i]+str(well_i)])

        model_plotls30.append(model_plot30)

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

        Sskv_data.loc[wellnest] *= np.exp(ls_sub30[p_multop[1]+"0"])

        Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
        Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

        K_data.loc[wellnest] *= np.exp(ls_sub30[p_multop[1]+"1"])

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
        modelpath = os.path.abspath("models//cowboyhat//")

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
    pump_interp30 = pd.concat(
        [pump_interp30.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

    # Reorder well list to shallow to deep aquifers
    # BK, PD, NL, NB
    well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
                  if y in x]

    pump_interp30.columns = well_names

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
                               pump_series=pump_interp30,
                               model_path=modelpath, califlag=0,
                               esmdaflag=0, user_models=model_plotls30)

    # Post process data
    sub_total, subv_total, ann_sub, \
        avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                                 sub_total,
                                                 subv_total,
                                                 all_results)

    bestsubtry30 = ann_sub.copy()

    lsreg_sub.append(bestsubtry30)
    lsreg_gw.append(model_plotls30)
    lsreg_pump.append(pd.Series(np.exp(pump_mean30)))

    # %% 31 regularized least squares

    # save fit report to a file:
    with open(os.path.abspath(
            lspath + "//" +
            wellnestlist[0] + "_LSreg31.0_modelresult.txt"),
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
    ls_sub31 = {}

    # Saving values from string
    for temp_i in range(len(temp)):
        variable_name = temp[temp_i].split(":")[0].strip()

        # Constant d lines are too long
        if "constant_d" in variable_name:
            ls_sub31[variable_name] = float(temp[temp_i][16:28])
        else:
            ls_sub31[variable_name] = float(temp[temp_i][14:25])

    # Saving model plots for least squares
    model_plotls31 = []

    # For all wells in well nest
    for well_i, wells in enumerate(well_names):

        # Name of well as a string
        well_name = wells
        #######################################################################

        # Saving pumping
        pump_mean31 = []
        for pump_i in range(n_pump):
            pump_mean31.append(ls_sub31["pump"+str(pump_i)])

        mean_ = pd.Series(np.exp(pump_mean31))
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
        pump_interp31 = df.interpolate(method="cubic")
        pump_interp31 = pump_interp31.dropna()
        pump_interp31 = pump_interp31.rename(columns={"0": well_name})
        model_plot31 = besttry_Pastasmodels[well_i].copy()
        model_plot31.del_stressmodel("well")  # Deletes previous pumping
        EstTotPump_ = ps.StressModel(pump_interp31,
                                     rfunc=ps.Gamma(), name="well",
                                     settings="well", up=False)
        model_plot31.add_stressmodel(EstTotPump_)

        # Assigns parameters to previous optimal parameters and SD
        for param_i in param_index:

            model_plot31.set_parameter(name=param_names[param_i],
                                       initial=ls_sub31[
                                           model_plot31.parameters[
                                               "optimal"].index[param_i]+str(well_i)],
                                       optimal=ls_sub31[
                                           model_plot31.parameters[
                                               "optimal"].index[param_i]+str(well_i)])

        model_plotls31.append(model_plot31)

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

        Sskv_data.loc[wellnest] *= np.exp(ls_sub31[p_multop[1]+"0"])

        Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
        Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

        K_data.loc[wellnest] *= np.exp(ls_sub31[p_multop[1]+"1"])

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
        modelpath = os.path.abspath("models//cowboyhat//")

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
    pump_interp31 = pd.concat(
        [pump_interp31.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

    # Reorder well list to shallow to deep aquifers
    # BK, PD, NL, NB
    well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
                  if y in x]

    pump_interp31.columns = well_names

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
                               pump_series=pump_interp31,
                               model_path=modelpath, califlag=0,
                               esmdaflag=0, user_models=model_plotls31)

    # Post process data
    sub_total, subv_total, ann_sub, \
        avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                                 sub_total,
                                                 subv_total,
                                                 all_results)

    bestsubtry31 = ann_sub.copy()

    lsreg_sub.append(bestsubtry31)
    lsreg_gw.append(model_plotls31)
    lsreg_pump.append(pd.Series(np.exp(pump_mean31)))

    # %% 35 regularized least squares

    # save fit report to a file:
    with open(os.path.abspath(
            lspath + "//" +
            wellnestlist[0] + "_LSreg35.0_modelresult.txt"),
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
    ls_sub35 = {}

    # Saving values from string
    for temp_i in range(len(temp)):
        variable_name = temp[temp_i].split(":")[0].strip()

        # Constant d lines are too long
        if "constant_d" in variable_name:
            ls_sub35[variable_name] = float(temp[temp_i][16:28])
        else:
            ls_sub35[variable_name] = float(temp[temp_i][14:25])

    # Saving model plots for least squares
    model_plotls35 = []

    # For all wells in well nest
    for well_i, wells in enumerate(well_names):

        # Name of well as a string
        well_name = wells
        #######################################################################

        # Saving pumping
        pump_mean35 = []
        for pump_i in range(n_pump):
            pump_mean35.append(ls_sub35["pump"+str(pump_i)])

        mean_ = pd.Series(np.exp(pump_mean35))
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
        pump_interp35 = df.interpolate(method="cubic")
        pump_interp35 = pump_interp35.dropna()
        pump_interp35 = pump_interp35.rename(columns={"0": well_name})
        model_plot35 = besttry_Pastasmodels[well_i].copy()
        model_plot35.del_stressmodel("well")  # Deletes previous pumping
        EstTotPump_ = ps.StressModel(pump_interp35,
                                     rfunc=ps.Gamma(), name="well",
                                     settings="well", up=False)
        model_plot35.add_stressmodel(EstTotPump_)

        # Assigns parameters to previous optimal parameters and SD
        for param_i in param_index:

            model_plot35.set_parameter(name=param_names[param_i],
                                       initial=ls_sub35[
                                           model_plot35.parameters[
                                               "optimal"].index[param_i]+str(well_i)],
                                       optimal=ls_sub35[
                                           model_plot35.parameters[
                                               "optimal"].index[param_i]+str(well_i)])

        model_plotls35.append(model_plot35)

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

        Sskv_data.loc[wellnest] *= np.exp(ls_sub35[p_multop[1]+"0"])

        Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
        Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

        K_data.loc[wellnest] *= np.exp(ls_sub35[p_multop[1]+"1"])

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
        modelpath = os.path.abspath("models//cowboyhat//")

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
    pump_interp35 = pd.concat(
        [pump_interp35.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

    # Reorder well list to shallow to deep aquifers
    # BK, PD, NL, NB
    well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
                  if y in x]

    pump_interp35.columns = well_names

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
                               pump_series=pump_interp35,
                               model_path=modelpath, califlag=0,
                               esmdaflag=0, user_models=model_plotls35)

    # Post process data
    sub_total, subv_total, ann_sub, \
        avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                                 sub_total,
                                                 subv_total,
                                                 all_results)

    bestsubtry35 = ann_sub.copy()

    lsreg_sub.append(bestsubtry35)
    lsreg_gw.append(model_plotls35)
    lsreg_pump.append(pd.Series(np.exp(pump_mean35)))

    # %% 50 regularized least squares

    # save fit report to a file:
    with open(os.path.abspath(
            lspath + "//" +
            wellnestlist[0] + "_LSreg50.0_modelresult.txt"),
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
    ls_sub50 = {}

    # Saving values from string
    for temp_i in range(len(temp)):
        variable_name = temp[temp_i].split(":")[0].strip()

        # Constant d lines are too long
        if "constant_d" in variable_name:
            ls_sub50[variable_name] = float(temp[temp_i][16:28])
        else:
            ls_sub50[variable_name] = float(temp[temp_i][14:25])

    # Saving model plots for least squares
    model_plotls50 = []

    # For all wells in well nest
    for well_i, wells in enumerate(well_names):

        # Name of well as a string
        well_name = wells
        #######################################################################

        # Saving pumping
        pump_mean50 = []
        for pump_i in range(n_pump):
            pump_mean50.append(ls_sub50["pump"+str(pump_i)])

        mean_ = pd.Series(np.exp(pump_mean50))
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
        pump_interp50 = df.interpolate(method="cubic")
        pump_interp50 = pump_interp50.dropna()
        pump_interp50 = pump_interp50.rename(columns={"0": well_name})
        model_plot50 = besttry_Pastasmodels[well_i].copy()
        model_plot50.del_stressmodel("well")  # Deletes previous pumping
        EstTotPump_ = ps.StressModel(pump_interp50,
                                     rfunc=ps.Gamma(), name="well",
                                     settings="well", up=False)
        model_plot50.add_stressmodel(EstTotPump_)

        # Assigns parameters to previous optimal parameters and SD
        for param_i in param_index:

            model_plot50.set_parameter(name=param_names[param_i],
                                       initial=ls_sub50[
                                           model_plot50.parameters[
                                               "optimal"].index[param_i]+str(well_i)],
                                       optimal=ls_sub50[
                                           model_plot50.parameters[
                                               "optimal"].index[param_i]+str(well_i)])

        model_plotls50.append(model_plot50)

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

        Sskv_data.loc[wellnest] *= np.exp(ls_sub50[p_multop[1]+"0"])

        Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
        Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

        K_data.loc[wellnest] *= np.exp(ls_sub50[p_multop[1]+"1"])

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
        modelpath = os.path.abspath("models//cowboyhat//")

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
    pump_interp50 = pd.concat(
        [pump_interp50.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

    # Reorder well list to shallow to deep aquifers
    # BK, PD, NL, NB
    well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
                  if y in x]

    pump_interp50.columns = well_names

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
                               pump_series=pump_interp50,
                               model_path=modelpath, califlag=0,
                               esmdaflag=0, user_models=model_plotls50)

    # Post process data
    sub_total, subv_total, ann_sub, \
        avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                                 sub_total,
                                                 subv_total,
                                                 all_results)

    bestsubtry50 = ann_sub.copy()

    lsreg_sub.append(bestsubtry50)
    lsreg_gw.append(model_plotls50)
    lsreg_pump.append(pd.Series(np.exp(pump_mean50)))

    def generate_simppumping_ens(ann_pump, n):
        """Generates ensemble of time series of groundwater pumping.

        Input:
        ann_pump - annual pumping rates (initial, mean) + std
        n - number of ensemble members

        Returns:
        interp_pump - list_dates with interpolated pumping m3/day
        """

        # ann_dates has Pump2 for each t (mean) and std which will be
        # used to generate normal dist at each t
        # New list for all values about to be randomly chosen and
        # interpolated
        ann_pump.index = pd.to_datetime(ann_pump.index)
        df = pd.DataFrame(index=ann_pump.index)

        for i in range(n):

            y = np.ones(len(ann_pump.index)) * random.randint(0, 500)

            # Dataframe
            mat = pd.DataFrame(y, index=ann_pump.index,
                               columns=[i])
            df = pd.concat([df, mat], join="outer",
                           keys=["Date", "Date"], axis=1)
            df.columns = df.columns.droplevel()

        return df


elif pumpexperiment == "pumpingcase1":
    # %% 15 regularized least squares

    # save fit report to a file:
    with open(os.path.abspath(
            lspath + "//" +
            wellnestlist[0] + "_LSreg15.0_modelresult.txt"),
            'r') as fh:
        temp = fh.readlines()
        temp = [x.replace("\n", "") for x in temp]
        # temp = temp[
        #     temp.index("[[Variables]]") + 1:temp.index(
        #         "[[Correlations]] (unreported correlations are < 0.100)")]
        temp = temp[
            temp.index("[[Variables]]") + 1:]
    fh.close()

    # Saving variables and values
    ls_sub15 = {}

    # Saving values from string
    for temp_i in range(len(temp)):
        variable_name = temp[temp_i].split(":")[0].strip()

        # Constant d lines are too long
        if "constant_d" in variable_name:
            ls_sub15[variable_name] = float(temp[temp_i][16:28])
        else:
            ls_sub15[variable_name] = float(temp[temp_i][14:25])

    # Saving model plots for least squares
    model_plotls15 = []

    # For all wells in well nest
    for well_i, wells in enumerate(well_names):

        # Name of well as a string
        well_name = wells
        #######################################################################

        # Saving pumping
        pump_mean15 = []
        for pump_i in range(n_pump):
            pump_mean15.append(ls_sub15["pump"+str(pump_i)])

        mean_ = pd.Series(np.exp(pump_mean15))
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
        pump_interp15 = df.interpolate(method="cubic")
        pump_interp15 = pump_interp15.dropna()
        pump_interp15 = pump_interp15.rename(columns={"0": well_name})
        model_plot15 = besttry_Pastasmodels[well_i].copy()
        model_plot15.del_stressmodel("well")  # Deletes previous pumping
        EstTotPump_ = ps.StressModel(pump_interp15,
                                     rfunc=ps.Gamma(), name="well",
                                     settings="well", up=False)
        model_plot15.add_stressmodel(EstTotPump_)

        # Assigns parameters to previous optimal parameters and SD
        for param_i in param_index:

            model_plot15.set_parameter(name=param_names[param_i],
                                       initial=ls_sub15[
                                           model_plot15.parameters[
                                               "optimal"].index[param_i]+str(well_i)],
                                       optimal=ls_sub15[
                                           model_plot15.parameters[
                                               "optimal"].index[param_i]+str(well_i)])

        model_plotls15.append(model_plot15)

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

        Sskv_data.loc[wellnest] *= np.exp(ls_sub15[p_multop[1]+"0"])

        Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
        Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

        K_data.loc[wellnest] *= np.exp(ls_sub15[p_multop[1]+"1"])

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
        modelpath = os.path.abspath("models//cowboyhat//")

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
    pump_interp15 = pd.concat(
        [pump_interp15.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

    # Reorder well list to shallow to deep aquifers
    # BK, PD, NL, NB
    well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
                  if y in x]

    pump_interp15.columns = well_names

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
                               pump_series=pump_interp15,
                               model_path=modelpath, califlag=0,
                               esmdaflag=0, user_models=model_plotls15)

    # Post process data
    sub_total, subv_total, ann_sub, \
        avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                                 sub_total,
                                                 subv_total,
                                                 all_results)

    bestsubtry15 = ann_sub.copy()

    lsreg_sub.append(bestsubtry15)
    lsreg_gw.append(model_plotls15)
    lsreg_pump.append(pd.Series(np.exp(pump_mean15)))

    # %% 60 regularized least squares

    # save fit report to a file:
    with open(os.path.abspath(
            lspath + "//" +
            wellnestlist[0] + "_LSreg60.0_modelresult.txt"),
            'r') as fh:
        temp = fh.readlines()
        temp = [x.replace("\n", "") for x in temp]
        # temp = temp[
        #     temp.index("[[Variables]]") + 1:temp.index(
        #         "[[Correlations]] (unreported correlations are < 0.100)")]
        temp = temp[
            temp.index("[[Variables]]") + 1:]
    fh.close()

    # Saving variables and values
    ls_sub60 = {}

    # Saving values from string
    for temp_i in range(len(temp)):
        variable_name = temp[temp_i].split(":")[0].strip()

        # Constant d lines are too long
        if "constant_d" in variable_name:
            ls_sub60[variable_name] = float(temp[temp_i][16:28])
        else:
            ls_sub60[variable_name] = float(temp[temp_i][14:25])

    # Saving model plots for least squares
    model_plotls60 = []

    # For all wells in well nest
    for well_i, wells in enumerate(well_names):

        # Name of well as a string
        well_name = wells
        #######################################################################

        # Saving pumping
        pump_mean60 = []
        for pump_i in range(n_pump):
            pump_mean60.append(ls_sub60["pump"+str(pump_i)])

        mean_ = pd.Series(np.exp(pump_mean60))
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
        pump_interp60 = df.interpolate(method="cubic")
        pump_interp60 = pump_interp60.dropna()
        pump_interp60 = pump_interp60.rename(columns={"0": well_name})
        model_plot60 = besttry_Pastasmodels[well_i].copy()
        model_plot60.del_stressmodel("well")  # Deletes previous pumping
        EstTotPump_ = ps.StressModel(pump_interp60,
                                     rfunc=ps.Gamma(), name="well",
                                     settings="well", up=False)
        model_plot60.add_stressmodel(EstTotPump_)

        # Assigns parameters to previous optimal parameters and SD
        for param_i in param_index:

            model_plot60.set_parameter(name=param_names[param_i],
                                       initial=ls_sub60[
                                           model_plot60.parameters[
                                               "optimal"].index[param_i]+str(well_i)],
                                       optimal=ls_sub60[
                                           model_plot60.parameters[
                                               "optimal"].index[param_i]+str(well_i)])

        model_plotls60.append(model_plot60)

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

        Sskv_data.loc[wellnest] *= np.exp(ls_sub60[p_multop[1]+"0"])

        Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
        Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

        K_data.loc[wellnest] *= np.exp(ls_sub60[p_multop[1]+"1"])

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
        modelpath = os.path.abspath("models//cowboyhat//")

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
    pump_interp60 = pd.concat(
        [pump_interp60.iloc[:, -1]]*num_wells, ignore_index=True, axis=1)

    # Reorder well list to shallow to deep aquifers
    # BK, PD, NL, NB
    well_names = [x for y in ["BK", "PD", "NL", "NB"] for x in well_names
                  if y in x]

    pump_interp60.columns = well_names

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
                               pump_series=pump_interp60,
                               model_path=modelpath, califlag=0,
                               esmdaflag=0, user_models=model_plotls60)

    # Post process data
    sub_total, subv_total, ann_sub, \
        avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                                 sub_total,
                                                 subv_total,
                                                 all_results)

    bestsubtry60 = ann_sub.copy()

    lsreg_sub.append(bestsubtry60)
    lsreg_gw.append(model_plotls60)
    lsreg_pump.append(pd.Series(np.exp(pump_mean60)))

# %% Objective Function calculation

if pumpexperiment == "simpleexp":
    lsreg = [0, 10, 20, 25, 30, 31, 35, 50]
elif pumpexperiment == "pumpingcase1":
    lsreg = [0, 15, 60]
index = gw_obs_indices[0]

mobjfc = []
robjfc = []

# For each model
for mod_i in range(len(besttry_Pastasmodels)):

    # For each regularization
    for reg_i, reg in enumerate(lsreg[1:]):
        temp_gwdata = gw_obs.merge(
            model_plotls[mod_i].simulate(),
            left_on=gw_obs.index, right_on=lsreg_gw[reg_i+1][mod_i].simulate().index,
            how="inner")
        y_gw = temp_gwdata.Simulation
        d_gw = temp_gwdata["0"]

        y_sub = lsreg_sub[reg_i+1][0][1].iloc[:, 2][index]*100
        d_sub = sub_obs.iloc[:, -1]

        if pumpexperiment == "simpleexp":
            d_pump = generate_simppumping_ens(annual_pump, 1)
            y_pump = lsreg_pump[reg_i+1]

        elif pumpexperiment == "pumpingcase1":
            d_pump = annual_pump.Pump.values
            y_pump = lsreg_pump[reg_i+1].values

        y = pd.concat([y_gw, y_sub], axis=0, ignore_index=True)
        d = pd.concat([d_gw, d_sub], axis=0, ignore_index=True)

        mobjfc.append(sum((y-d)**2)/np.var(y))
        robjfc.append(sum((y_pump-d_pump)**2)/np.var(y_pump))

# %% Plotting L curve

plt.rc("font", size=28)  # controls default text size
plt.rc("axes", titlesize=24)  # fontsize of the title
plt.rc("axes", labelsize=18)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=16)  # fontsize of the x tick labels
plt.rc("ytick", labelsize=16)  # fontsize of the y tick labels
plt.rc("legend", fontsize=14)  # fontsize of the legend

fig, ax = plt.subplots()
plt.scatter(robjfc, mobjfc)
plt.ylabel("Regularized Objective Function")
plt.xlabel("Measurement Objective Function")
for i, txt in enumerate(lsreg[1:]):
    ax.annotate(txt/100, (robjfc[i], mobjfc[i]))
