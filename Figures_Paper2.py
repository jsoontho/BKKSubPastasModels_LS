# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 16:30:29 2024

Creating figures for the second paper
@author: jtsoonthornran
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import warnings
from sklearn.metrics import mean_squared_error
import random
import os
import pastas as ps
import sys
import datetime as dt
from mpl_toolkits.basemap import Basemap
import importlib
import math

# Bangkok Subsidence Model Package
import bkk_sub_gw

# Ignoring Pastas warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

# Changing current directory to locaiton of python script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# %%##########################################################################
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

        # Taking from normal distriubtion
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

        # Formatting
        mat = pd.DataFrame(temp_list, index=ann_pump.index,
                           columns=[i])
        df = pd.concat([df, mat], join="outer",
                       keys=["Date", "Date"], axis=1)
        df.columns = df.columns.droplevel()

    return df


# %%###########################################################################
# Plotting settings
###############################################################################

plt.rc("font", size=28)  # controls default text size
plt.rc("axes", titlesize=24)  # fontsize of the title
plt.rc("axes", labelsize=18)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=16)  # fontsize of the x tick labels
plt.rc("ytick", labelsize=16)  # fontsize of the y tick labels
plt.rc("legend", fontsize=14)  # fontsize of the legend

# True parameters
Atrue = -.1
ntrue = 2.5
atrue = 50
dtrue = 2

# Number of ensemble members
ne = 250
par_error = [30, 100, 100, 30]
dist = ["norm", "norm", "norm", "norm"]

wellnestlist = ["LCBKK018"]
tot_path = os.path.abspath("inputs")
modelpath = os.path.abspath("models//synthetic//perfectmodel//")
files = wellnestlist

for Wellnest_name in files:

    # Reading in groundwater data
    full_path = os.path.join(tot_path, Wellnest_name + ".xlsx")
    data = pd.read_excel(full_path, skiprows=3)

    num_wells = len(data.columns)-2

    besttry_Pastasmodels = []

    # For all wells in well nest
    for num_well, wells in enumerate(data.columns[-(len(data.columns)-2):]):

        # Name of well as a string
        well_name = wells

        # if np.logical_or("PD" in well_name, "NL" in well_name):
        #     continue
        #######################################################################
        # Assuming importing Pastas Model
        #######################################################################

        # Model files
        modelfiles = os.listdir(modelpath)

        # Load existing model
        wellmodel = [s for s in modelfiles
                     if np.logical_and(Wellnest_name in s, well_name in s)][0]
        ml = ps.io.load(modelpath + "/" + wellmodel)

        besttry_Pastasmodels.append(ml)

# Get initial parameters
init = ml.get_init_parameters()
# opt = model.parameters

# Params
params = pd.DataFrame()

# Stats module
mod = importlib.import_module("scipy.stats")

# Random values
# Set the initial parameters to a normal distribution
for name_i, name in enumerate(ml.parameters.index):

    loc = init.initial[name]  # Pastas initial
    scale = abs(par_error[name_i]/100 * loc)

    # Sampled
    sampled = []

    minimum = ml.parameters["pmin"][name]
    maximum = ml.parameters["pmax"][name]
    if math.isnan(minimum):
        minimum = -10000

    if math.isnan(maximum):
        maximum = 10000

    # If the samples don't equal the number of ensemble members wanted
    while len(sampled) != ne:

        if dist[name_i] == "lognorm":
            data = np.log(getattr(mod, dist[name_i]).rvs(
                s=scale, loc=0, scale=np.exp(loc), size=ne-len(sampled)))
        elif dist[name_i] == "uniform":
            data = getattr(mod, dist[name_i]).rvs(loc=minimum,
                                                  scale=maximum-minimum,
                                                  size=ne-len(sampled))
        else:
            data = getattr(mod, dist[name_i]).rvs(loc=loc, scale=scale,
                                                  size=ne-len(sampled))

        sampled.extend(data[(minimum <= data) & (data <= maximum)])

    params[name] = sampled

# %% Plotting

# Importing pumping
# Annual pumping data (mean), std
pumppath = os.path.join(os.path.abspath("inputs"),
                        "BasinPumping_Annual_ESMDA.xlsx")
pumpsheet = "EstTotalPump_54-60"
annual_pump = pd.read_excel(pumppath, sheet_name=pumpsheet,
                            index_col=0, parse_dates=["Date"])

# Only until 2024
annual_pump = annual_pump[annual_pump.index <= "2023"]

# Importing subsidence boundaries
# PARAMETER BOUNDARIES!
# parambound_path = os.path.join(os.path.abspath("inputs"),
#                                "SUBParametersCali.xlsx")
parambound_path = os.path.join(os.path.abspath("inputs"),
                               "SUBParametersPriortoManual.xlsx")

parambound = pd.read_excel(parambound_path,
                           sheet_name="bounds_mult",
                           index_col=0)
parambound = pd.DataFrame(parambound)
parambound.iloc[:, 0:2] = np.log(parambound.iloc[:, 0:2])

m_bounds = list(zip(parambound.loc[
    parambound.Wellnest == wellnestlist[0]].iloc[0:3, 0],
    parambound.loc[
    parambound.Wellnest == wellnestlist[0]].iloc[0:3, 1]))

nm = 2
p_mult = []
[p_mult.append(random.uniform(m_bounds[x][0], m_bounds[x][1]))
 for x in np.arange(0, nm+1, nm)]

# Prior: Let's start average of the bounds
mprior_mean = p_mult

mprior_mean = np.array(mprior_mean)

# mprior_mean = np.ones(nm) * random.uniform(m_bounds[x][0], m_bounds[x][1])

# Parameter error
par_error = [4, 4]
par_std = par_error

# Random number generator
rng = np.random.default_rng()

# Prior distribution of parameters
param_init = rng.normal(loc=mprior_mean, scale=par_std,
                        size=(ne, nm)).T

param_init = pd.DataFrame(param_init)

# Pumping error in prior
pump_err = .5
# Generate pumping ensemble
option = ["normal", .99]
annual_pump["Std"] = annual_pump['Pump'] * pump_err
pumping_ens = generate_pumping_ens(annual_pump, ne, option)

# Parameters to be plotted
mydata = params.iloc[:, -1]
fig, axs = plt.subplots(3, 1, figsize=(15, 10))
axs[0].hist(mydata,
            weights=np.zeros_like(mydata) + 1. / mydata.size,
            color="tab:red", edgecolor="k", linewidth=3,
            alpha=.9)
# Truth
# axs[0].axvline(dtrue, ymin=0, ymax=1, label=("Truth (" + str(dtrue) + " m)"),
#                linewidth=6,
#                color="c")
axs[0].axvline(np.mean(mydata), ymin=0, ymax=1, label=(
    "Prior Ens. Mean (" + f'{np.mean(mydata):.1f}' + " m)"), linewidth=6,
    color="maroon", alpha=1)
[axs[0].axvline(_x, ymin=0, ymax=.05, linewidth=1, color='lightcoral') for _x in mydata]
axs[0].set_xlabel("$\it{d}$ Values (m)")
axs[0].set_ylabel("Relative\nFrequency (-)")
axs[0].legend()

# Plotting pumping
axs[1].plot(pumping_ens, alpha=.05, color="tab:red")
axs[1].plot(annual_pump.Pump,
            label="Prior Ens. Mean\n(Basin-wide pumping)", color="maroon", linewidth=6)
axs[1].set_xlabel("Date")
axs[1].set_ylabel("Pumping Rate * 10$^4$\n(m$^3$/day)")
axs[1].legend()

# Plotting subsidence parameter
axs[2].hist(param_init.iloc[0, :],
            weights=np.zeros_like(param_init.iloc[0, :]) +
            1. / param_init.iloc[0, :].size,
            color="tab:red", edgecolor="k", linewidth=3,
            alpha=.9)
[axs[2].axvline(_x,
                ymin=0,
                ymax=.05, linewidth=1,
                color='lightcoral') for _x in param_init.iloc[0, :]]

axs[2].set_xlabel("Log $\it{Sskv}$ Multiplier Values")
axs[2].set_ylabel("Relative\nFrequency (-)")
axs[2].axvline(np.mean(param_init.iloc[0, :]), ymin=0, ymax=1, label=(
    "Prior Ens. Mean (" + f'{np.mean(param_init.iloc[0, :]):.1f}' + ")"), linewidth=6,
    color="maroon", alpha=1)
axs[2].legend()

# Labeling
axs[0].annotate("(a)",
                xy=(0.04, .85), xycoords="axes fraction",
                fontsize=20, horizontalalignment="right",
                weight="bold",
                verticalalignment="bottom")
axs[1].annotate("(b)",
                xy=(0.04, .85), xycoords="axes fraction",
                fontsize=20, horizontalalignment="right",
                weight="bold",
                verticalalignment="bottom")
axs[2].annotate("(c)",
                xy=(0.04, .85), xycoords="axes fraction",
                fontsize=20, horizontalalignment="right",
                weight="bold",
                verticalalignment="bottom")

fig.tight_layout()
plt.subplots_adjust(hspace=.3)

# Saving figure
# Path to save models
figpath = os.path.abspath("figures")
fig.set_rasterized(True)
fig_name1 = "Constantd_Pump_subPrior.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

fig_name1 = "Constantd_Pump_subPrior.eps"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="eps")

# %% Plotting inelastic/elastic values

# Reading in thickness and storage data
path = os.path.join(os.path.abspath("inputs"),
                    "SUBParameters.xlsx")
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
# Getting a list of all the wells
# Total path
tot_path = os.path.abspath("inputs")

files = os.listdir(tot_path)
files = [i.replace(".xlsx", "") for i in files
         if i.startswith("LC") and "_" not in i]

# Subsets of the res_tmax data, and list of wellnest names
well_data = []
wellnest_list = []

for wellnest in files:
    well_data.extend([wellnest, wellnest, K_data.loc[wellnest][::2][0]])
    wellnest_list.append(wellnest)

# Importing spatial coordinates
full_path = os.path.join(tot_path, "GroundwaterWellLocs.xls")
gwwell_locs = pd.read_excel(full_path)

# Locations of wellnests; removing duplicates
gwwell_locs = gwwell_locs.drop_duplicates("WellNest_Name", keep="first")

# Preallocation
# Empty dictionary
d_dict = {}

# Preallocation
# Saving relevant xs, ys, and tmax
xs = []
ys = []
cs = []

# Unique x, ys only
# Dissuades if multiple well nests at the same location
unique = []

# Getting rid of repeating wells and data points
# zip joins x and y coordinates in pairs
for x, y in zip(gwwell_locs.Long, gwwell_locs.Lat):

    # Check if x, y is unique
    if (x, y) not in unique:

        # Saves location for plotting
        unique.append((x, y))

        # Label is well nest name
        label = gwwell_locs.loc[
            gwwell_locs.Long == x]["WellNest_Name"].tolist()

        # Specific well nest does not have a well in the aquifer
        if label[0] not in wellnest_list:
            continue

        # If well nest has a well in the aquifer
        else:

            # Saving data
            xs.append(x)
            ys.append(y)
            cs.append(well_data[well_data.index(label[0])+2])

# Plot settings
# Initializing figure
fig, ax = plt.subplots(figsize=(3.2, 2.2), dpi=400)

data_lim = [1E-7, 1E-3]

plt.set_cmap("plasma")  # Color map colors

# Plots
map = Basemap(llcrnrlon=100.3, llcrnrlat=13.4, urcrnrlon=100.8, urcrnrlat=14,
              resolution="h", ellps="WGS84", lat_0=13.6, lon_0=100.4)
bkk_sub_gw.bkk_plotting.draw_basemap(map, xs, ys, cs, fig=fig, ax=ax,
                                     datalim=data_lim, mode="Param_val", save=0,
                                     perc=0,
                                     figpath=figpath)

# Saving Figures
fig_name1 = "VscK_map.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

fig_name1 = "VscK_map.eps"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="eps")

# %% Plotting RMSE over for every well nest

# Getting a list of all the wells
# Total path
tot_path = os.path.abspath("inputs")

files = os.listdir(tot_path)
files = [i.replace(".xlsx", "") for i in files
         if i.startswith("LC") and "_" not in i]

# Subsets of the res_tmax data, and list of wellnest names
well_data = []
wellnest_list = []

for wellnest in files:

    # Folder to save/import graph and model
    modelpath = os.path.abspath("models//BKK//")

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

    subdata = pd.read_excel(loc, sheet_name=wellnest + "_Leveling",
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
    bench[bench == 0] = np.nan  # 0's to nan

    # Gets the last date of each year
    lastdate = bench.groupby(pd.DatetimeIndex(bench["date"]).year,
                             as_index=False).agg(
                                 {"date": max}).reset_index(drop=True)
    bench = bench.loc[lastdate.date]

    # Subsidence observations
    dobs = pd.concat([dobs, bench.iloc[:, 0].dropna()])
    # Saving obs indices
    gw_obs_indices.append(bench.iloc[:, 0].dropna().index)

    sub_obs = bench.iloc[:, 0].dropna()

    # Reading in groundwater data
    full_path = os.path.join(tot_path, wellnest + ".xlsx")
    data = pd.read_excel(full_path, skiprows=3)

    # Importing model
    # Model files
    modelfiles = os.listdir(modelpath)

    # Groundwater observations
    # For all wells in well nest
    models = []
    gw_obs = pd.Series(np.empty(1, dtype=object))
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
                         if np.logical_and(wellnest in s,
                                           well_name in s)][0]
            model = ps.io.load(modelpath + "/" + wellmodel)
            models.append(model)

            gw_obs = pd.concat([gw_obs, model.observations()])

            # Saving groundwater observations and indices
            gw_obs_indices.append(model.observations().index)
            dobs = pd.concat([dobs, model.observations()])
        except:

            sys.exit("Model doesn't exist")
    gw_obs = gw_obs[1:]
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
    well_names = data.columns[-(len(data.columns)-2):].values.tolist()

    # Number of pumping
    # Number of pumping, pastas, sub parameters
    n_pump = len(annual_pump)

    # Path to save models
    modelpath = os.path.abspath("models//BKK//ESMDA//na8ne250")

    # Number of subsidence parameters
    if p_multop[1] == "Sskv" or p_multop[1] == "K" or p_multop[1] == "Sske":
        n_sub = 1

    elif p_multop[1] == "Ss" or p_multop[1] == "SsK":
        n_sub = 2

    elif p_multop[1] == "all":
        n_sub = 3

    if wellnest == "LCBKK013":
        # Creating blank dictionary
        sub_m = {"dpred": np.zeros((na+1, 709,
                                    500)),
                 "mprior": np.zeros((na+1, n_pump+num_wells*n_param+n_sub,
                                    500))}
    else:
        # Creating blank dictionary
        sub_m = {"dpred": np.zeros((na+1, len(sub_obs)+len(gw_obs),
                                    ne)),
                 "mprior": np.zeros((na+1, n_pump+num_wells*n_param+n_sub,
                                    ne))}

    # Saving ESMDA results
    for na_i in range(na+1):

        fig_name1 = wellnest + "_Dpred_" + p_multop[1] + "_na" + str(na_i) + ".csv"
        full_figpath = os.path.join(modelpath, fig_name1)
        sub_m["dpred"][na_i] = pd.read_csv(full_figpath, delimiter=",", header=None,
                                           dtype=str)

        sub_m["dpred"][na_i] = sub_m["dpred"][na_i].astype("float64")

        fig_name1 = wellnest + "_Mprior_" + p_multop[1] + "_na" + str(na_i) + ".csv"
        full_figpath = os.path.join(modelpath, fig_name1)

        sub_m["mprior"][na_i] = pd.read_csv(full_figpath, delimiter=",", header=None,
                                            dtype=str)

        sub_m["mprior"][na_i] = sub_m["mprior"][na_i].astype("float64")

    ens_mresults = {}
    ens_mresults["mprior"] = []
    averages = np.average(sub_m["mprior"][-1, :, :], axis=1)
    ens_mresults["mprior"].append(averages)

    # Calculating rmse for sub
    # add rsq to simulation
    rmse = mean_squared_error(np.mean(sub_m["dpred"][-1], axis=1)[
        0:len(gw_obs_indices[0])], sub_obs, squared=False)

    well_data.extend([wellnest, wellnest, rmse])
    wellnest_list.append(wellnest)

# Importing spatial coordinates
full_path = os.path.join(tot_path, "GroundwaterWellLocs.xls")
gwwell_locs = pd.read_excel(full_path)

# Locations of wellnests; removing duplicates
gwwell_locs = gwwell_locs.drop_duplicates("WellNest_Name", keep="first")

# Preallocation
# Empty dictionary
d_dict = {}

# Preallocation
# Saving relevant xs, ys, and tmax
xs = []
ys = []
cs = []

# Unique x, ys only
# Dissuades if multiple well nests at the same location
unique = []

# Getting rid of repeating wells and data points
# zip joins x and y coordinates in pairs
for x, y in zip(gwwell_locs.Long, gwwell_locs.Lat):

    # Check if x, y is unique
    if (x, y) not in unique:

        # Saves location for plotting
        unique.append((x, y))

        # Label is well nest name
        label = gwwell_locs.loc[
            gwwell_locs.Long == x]["WellNest_Name"].tolist()

        # Specific well nest does not have a well in the aquifer
        if label[0] not in wellnest_list:
            continue

        # If well nest has a well in the aquifer
        else:

            # Saving data
            xs.append(x)
            ys.append(y)
            cs.append(well_data[well_data.index(label[0])+2])

print(np.mean(cs))
# Plot settings
# Initializing figure
fig, ax = plt.subplots(figsize=(3.2, 2.2), dpi=400)

data_lim = [.4, 3.9]

plt.set_cmap("coolwarm")  # Color map colors

# Plots
figpath = os.path.abspath("figures")
map = Basemap(llcrnrlon=100.3, llcrnrlat=13.4, urcrnrlon=100.8, urcrnrlat=14,
              resolution="h", ellps="WGS84", lat_0=13.6, lon_0=100.4)
bkk_sub_gw.bkk_plotting.draw_basemap(map, xs, ys, cs, fig=fig, ax=ax,
                                     datalim=data_lim, mode="Sub_RMSE_paper2", save=0,
                                     perc=0,
                                     figpath=figpath)

# Saving Figures
fig_name1 = "RMSE_sub_map.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

fig_name1 = "RMSE_sub_map.eps"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="eps")

# %% Plotting pumping over for every well nest

# Getting a list of all the wells
# Total path
tot_path = os.path.abspath("inputs")

files = os.listdir(tot_path)
files = [i.replace(".xlsx", "") for i in files
         if i.startswith("LC") and "_" not in i]

# Subsets of the res_tmax data, and list of wellnest names
well_data = []
wellnest_list = []

for wellnest in files:

    # Folder to save/import graph and model
    modelpath = os.path.abspath("models//BKK//")

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

    subdata = pd.read_excel(loc, sheet_name=wellnest + "_Leveling",
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
    bench[bench == 0] = np.nan  # 0's to nan

    # Gets the last date of each year
    lastdate = bench.groupby(pd.DatetimeIndex(bench["date"]).year,
                             as_index=False).agg(
                                 {"date": max}).reset_index(drop=True)
    bench = bench.loc[lastdate.date]

    # Subsidence observations
    dobs = pd.concat([dobs, bench.iloc[:, 0].dropna()])
    # Saving obs indices
    gw_obs_indices.append(bench.iloc[:, 0].dropna().index)

    sub_obs = bench.iloc[:, 0].dropna()

    # Reading in groundwater data
    full_path = os.path.join(tot_path, wellnest + ".xlsx")
    data = pd.read_excel(full_path, skiprows=3)

    # Importing model
    # Model files
    modelfiles = os.listdir(modelpath)

    # Groundwater observations
    # For all wells in well nest
    models = []
    gw_obs = pd.Series(np.empty(1, dtype=object))
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
                         if np.logical_and(wellnest in s,
                                           well_name in s)][0]
            model = ps.io.load(modelpath + "/" + wellmodel)
            models.append(model)

            gw_obs = pd.concat([gw_obs, model.observations()])

            # Saving groundwater observations and indices
            gw_obs_indices.append(model.observations().index)
            dobs = pd.concat([dobs, model.observations()])
        except:

            sys.exit("Model doesn't exist")
    gw_obs = gw_obs[1:]
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
    well_names = data.columns[-(len(data.columns)-2):].values.tolist()

    # Number of pumping
    # Number of pumping, pastas, sub parameters
    n_pump = len(annual_pump)

    # Path to save models
    modelpath = os.path.abspath("models//BKK//ESMDA//na8ne250//")

    # Number of subsidence parameters
    if p_multop[1] == "Sskv" or p_multop[1] == "K" or p_multop[1] == "Sske":
        n_sub = 1

    elif p_multop[1] == "Ss" or p_multop[1] == "SsK":
        n_sub = 2

    elif p_multop[1] == "all":
        n_sub = 3

    if wellnest == "LCBKK013":
        # Creating blank dictionary
        sub_m = {"dpred": np.zeros((na+1, 709,
                                    500)),
                 "mprior": np.zeros((na+1, n_pump+num_wells*n_param+n_sub,
                                    500))}
    else:
        # Creating blank dictionary
        sub_m = {"dpred": np.zeros((na+1, len(sub_obs)+len(gw_obs),
                                    ne)),
                 "mprior": np.zeros((na+1, n_pump+num_wells*n_param+n_sub,
                                    ne))}

    # Saving ESMDA results
    for na_i in range(na+1):

        fig_name1 = wellnest + "_Dpred_" + p_multop[1] + "_na" + str(na_i) + ".csv"
        full_figpath = os.path.join(modelpath, fig_name1)
        sub_m["dpred"][na_i] = pd.read_csv(full_figpath, delimiter=",", header=None,
                                           dtype=str)

        sub_m["dpred"][na_i] = sub_m["dpred"][na_i].astype("float64")

        fig_name1 = wellnest + "_Mprior_" + p_multop[1] + "_na" + str(na_i) + ".csv"
        full_figpath = os.path.join(modelpath, fig_name1)

        sub_m["mprior"][na_i] = pd.read_csv(full_figpath, delimiter=",", header=None,
                                            dtype=str)

        sub_m["mprior"][na_i] = sub_m["mprior"][na_i].astype("float64")

    ens_mresults = {}
    ens_mresults["mprior"] = []
    averages = np.average(sub_m["mprior"][-1, :, :], axis=1)
    ens_mresults["mprior"].append(averages)

    # Pumping index
    pump_index0 = n_param*num_wells
    pump_index1 = n_param*num_wells + n_pump

    # For all wells in well nest
    for well_i, wells in enumerate(well_names):

        if well_i > 0:
            continue
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

    well_data.extend([wellnest, wellnest, np.mean(pump_interp.loc["1997"])])
    wellnest_list.append(wellnest)

# Importing spatial coordinates
full_path = os.path.join(tot_path, "GroundwaterWellLocs.xls")
gwwell_locs = pd.read_excel(full_path)

# Locations of wellnests; removing duplicates
gwwell_locs = gwwell_locs.drop_duplicates("WellNest_Name", keep="first")

# Preallocation
# Empty dictionary
d_dict = {}

# Preallocation
# Saving relevant xs, ys, and tmax
xs = []
ys = []
cs = []

# Unique x, ys only
# Dissuades if multiple well nests at the same location
unique = []

# Getting rid of repeating wells and data points
# zip joins x and y coordinates in pairs
for x, y in zip(gwwell_locs.Long, gwwell_locs.Lat):

    # Check if x, y is unique
    if (x, y) not in unique:

        # Saves location for plotting
        unique.append((x, y))

        # Label is well nest name
        label = gwwell_locs.loc[
            gwwell_locs.Long == x]["WellNest_Name"].tolist()

        # Specific well nest does not have a well in the aquifer
        if label[0] not in wellnest_list:
            continue

        # If well nest has a well in the aquifer
        else:

            # Saving data
            xs.append(x)
            ys.append(y)
            cs.append(well_data[well_data.index(label[0])+2])

# Plot settings
# Initializing figure
fig, ax = plt.subplots(figsize=(3.2, 2.2), dpi=400)

data_lim = [100, 300]

plt.set_cmap("Purples_r")  # Color map colors
figpath = os.path.abspath("figures")
# Plots
map = Basemap(llcrnrlon=100.3, llcrnrlat=13.4, urcrnrlon=100.8, urcrnrlat=14,
              resolution="h", ellps="WGS84", lat_0=13.6, lon_0=100.4)
bkk_sub_gw.bkk_plotting.draw_basemap(map, xs, ys, cs, fig=fig, ax=ax,
                                     datalim=data_lim, mode="pumping", save=0,
                                     perc=0,
                                     figpath=figpath)

# Saving Figures

fig_name1 = "Pumping1997_map.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

fig_name1 = "Pumping1997_map.eps"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="eps")

# %% GW Locs

bkk_sub_gw.bkk_plotting.gwlocs_map(figpath, save=1)

# %% RMSE na Sensitivity (BANGKOK)

wellnestlist = ["LCBKK005"]

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
bench[bench == 0] = np.nan  # 0's to nan

# Gets the last date of each year
lastdate = bench.groupby(pd.DatetimeIndex(bench["date"]).year,
                         as_index=False).agg(
                             {"date": max}).reset_index(drop=True)
bench = bench.loc[lastdate.date]

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
modelpath = os.path.abspath("models//BKK")
modelfiles = os.listdir(modelpath)

# Groundwater observations
# For all wells in well nest
models = []
gw_obs = pd.Series(np.empty(1, dtype=object))
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

        gw_obs = pd.concat([gw_obs, model.observations()])

        # Saving groundwater observations and indices
        gw_obs_indices.append(model.observations().index)
        dobs = pd.concat([dobs, model.observations()])
    except:

        sys.exit("Model doesn't exist")
gw_obs = gw_obs[1:]
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
well_names = data.columns[-(len(data.columns)-2):].values.tolist()

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

# For each path (na)
modelpaths = [os.path.abspath("models//BKK//ESMDA//na1ne250//"),
              os.path.abspath("models//BKK//ESMDA//na2ne250//"),
              os.path.abspath("models//BKK//ESMDA//na3ne250//"),
              os.path.abspath("models//BKK//ESMDA//na7ne250//"),
              os.path.abspath("models//BKK//ESMDA//na8ne250//"),
              os.path.abspath("models//BKK//ESMDA//na9ne250//"),
              os.path.abspath("models//BKK//ESMDA//na16ne250//"),]

nas = [1, 2, 3, 7, 8, 9, 16]

# Saving rmses
rmses = []

for modelpathi, modelpath in enumerate(modelpaths):

    # Creating blank dictionary
    sub_m = {"dpred": np.zeros((nas[modelpathi]+1, len(sub_obs)+len(gw_obs),
                                ne)),
             "mprior": np.zeros((nas[modelpathi]+1, n_pump+num_wells*n_param+n_sub,
                                ne))}

    # Saving ESMDA results
    for na_i in range(nas[modelpathi]+1):

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

    # Calculating rmse for sub
    # add rsq to simulation
    rmse_sub = mean_squared_error(np.mean(sub_m["dpred"][-1], axis=1)[
        0:len(gw_obs_indices[0])], sub_obs, squared=False)

    rmses.append(rmse_sub)

# Plotting rmse over time

plt.rc("font", size=28)  # controls default text size
plt.rc("axes", titlesize=24)  # fontsize of the title
plt.rc("axes", labelsize=18)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=16)  # fontsize of the x tick labels
plt.rc("ytick", labelsize=16)  # fontsize of the y tick labels
plt.rc("legend", fontsize=14)  # fontsize of the legend

plt.figure(figsize=(10, 4))

plt.plot(nas, rmses, "-o", linewidth=4, color="tab:brown")
plt.xlabel("$N_{mda}$ (total number of ESMDA iterations)")
plt.ylabel("Subsidence RMSE (cm/yr)")

# Saving
figpath = os.path.abspath("figures")
fig_name1 = wellnest + "_sub_nas_rmses.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

fig_name1 = wellnest + "_sub_nas_rmses.eps"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="eps")


# %% RMSE na Sensitivity (synthetic)

wellnestlist = ["LCBKK011"]

Wellnest_name = wellnestlist[0]

# For each path (na)
modelpaths = [os.path.abspath("models//synthetic//na1ne250//"),
              os.path.abspath("models//synthetic//na2ne250//"),
              os.path.abspath("models//synthetic//na3ne250//"),
              os.path.abspath("models//synthetic//na8ne250//"),
              os.path.abspath("models//synthetic//na9ne250//"),
              os.path.abspath("models//synthetic//na16ne250//"),]

nas = [1, 2, 3, 8, 9, 16]

# Saving rmses
rmses = []

for modelpathi, modelpath in enumerate(modelpaths):

    # Total path
    tot_path = os.path.abspath("inputs")
    # Reading in groundwater data
    full_path = os.path.join(tot_path, Wellnest_name + ".xlsx")
    data = pd.read_excel(full_path, skiprows=3)

    # Saving sub and gw obs
    dobs = pd.Series(np.empty(1, dtype=object))
    gw_obs_indices = []

    # Sub obs
    fig_name1 = Wellnest_name + "_SubObs.csv"
    full_figpath = os.path.join(modelpath, fig_name1)
    sub_obs = pd.read_csv(full_figpath, delim_whitespace="\t")
    sub_obs.Date = pd.to_datetime(sub_obs.Date, format="%Y-%m-%d")

    # Reading
    # GW obs
    fig_name1 = wellnestlist[0] + "_GWObs.csv"
    full_figpath = os.path.join(modelpath, fig_name1)
    gw_obs = pd.read_csv(full_figpath, delim_whitespace="\t")
    gw_obs.Date = pd.to_datetime(gw_obs.Date, format="%Y-%m-%d")

    # Saving gw dates
    gw_obs_indices = [sub_obs.Date]

    for well_i in range(num_wells):

        gw_obs_indices.append(gw_obs.Date[well_i*int(
            len(gw_obs)/num_wells):(well_i+1)*int(len(gw_obs)/num_wells)])

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
    param_index = np.array([1, 2])

    n_pump2 = len(annual_pump2)

    n_param = len(param_index)

    # Number of wells
    num_wells = len(data.columns[-(len(data.columns)-2):])
    # WEll names
    well_names = data.columns[-(len(data.columns)-2):].values.tolist()

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
    sub_m = {"dpred": np.zeros((nas[modelpathi]+1, len(sub_obs)+len(gw_obs),
                                ne)),
             "mprior": np.zeros((nas[modelpathi]+1, n_pump+num_wells*n_param+n_sub,
                                ne))}

    # Saving ESMDA results
    for na_i in range(nas[modelpathi]+1):

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

    # Calculating rmse for sub
    # add rsq to simulation
    rmse_sub = mean_squared_error(np.mean(sub_m["dpred"][-1], axis=1)[
        0:len(gw_obs_indices[0])], sub_obs.iloc[:, -1], squared=False)

    rmses.append(rmse_sub)

# Plotting rmse over time

plt.rc("font", size=28)  # controls default text size
plt.rc("axes", titlesize=24)  # fontsize of the title
plt.rc("axes", labelsize=18)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=16)  # fontsize of the x tick labels
plt.rc("ytick", labelsize=16)  # fontsize of the y tick labels
plt.rc("legend", fontsize=14)  # fontsize of the legend

plt.figure(figsize=(10, 4))

plt.plot(nas, rmses, "-o", linewidth=4, color="tab:brown")
plt.xlabel("$N_{mda}$ (total number of ESMDA iterations)")
plt.ylabel("Subsidence RMSE (cm/yr)")

# Saving
figpath = os.path.abspath("figures")
fig_name1 = wellnest + "_sub_nas_rmses.png"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="png")

fig_name1 = wellnest + "_sub_nas_rmses.eps"
full_figpath = os.path.join(figpath, fig_name1)
plt.savefig(full_figpath, bbox_inches="tight", format="eps")
