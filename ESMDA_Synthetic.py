# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:29:28 2024

@author: jtsoonthornran

SYNTHETIC CASE:

Code to run ESMDA with state vector including pastas parameters,
pumping time series, and subsidence multiplier parameters

Synthetic case where true pumping is arbitrarily lower for 1970-1990
and best info is basin-wide pumping

Code for plotting is separate

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
from scipy.special import gammainc, gammaincinv
import datetime as dt
import random
import pickle
import pastas as ps
import warnings

# Bangkok Subsidence Model Package
import bkk_sub_gw

# Ignoring Pastas warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

# Changing current directory to locaiton of python script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# %% ESMDA: parameter and control estimation
# Basin-wide pumping added in


############################################################################
# Generate ensemble of pumping time series settings
###############################################################################

def generate_pumping_ens(ann_pump, n, option):
    """Generates ensemble of time series of groundwater pumping.

    Input:
    ann_pump - annual pumping rates (initial, mean) + std
    n - number of ensemble members
    option - completely random within normal dist

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

        if option[0] == "normal":

            # For each t
            for t in range(n_pump):

                # Choosing between normal distribution
                # (mean, std)
                temp = np.random.normal(ann_pump.iloc[t, 1],
                                        ann_pump.iloc[t, 2])

                # Make sure temp is > 0
                if temp < 0:

                    temp_list.append(0)

                else:

                    temp_list.append(temp)

        # Dataframe
        mat = pd.DataFrame(temp_list, index=ann_pump.index,
                           columns=[i])
        df = pd.concat([df, mat], join="outer",
                       keys=["Date", "Date"], axis=1)
        df.columns = df.columns.droplevel()

    return df


############################################################################
# Pastas settings
###############################################################################

# If saving model, save_model = 1
save_model = 1

# Folder to save/import graph and model
modelpath = os.path.abspath("models//synthetic//na16ne250")

# Pumping response function
pump_rfunc = ps.Gamma()

# Calibration period
calitime_min = "1978"
calitime_max = "2005"

# Noise model
noise_TF = False

# Total path
tot_path = os.path.abspath("inputs")

# ###############################################################################
# # Creating Synthetic Observations
# ###############################################################################


def gamma_tmax(A, n, a, cutoff=0.999):
    return gammaincinv(n, cutoff) * a


def gamma_step(A, n, a, cutoff=0.999):
    tmax = gamma_tmax(A, n, a, cutoff)
    t = np.arange(0, tmax, 1)
    s = A * gammainc(n, t / a)
    return s


def gamma_block(A, n, a, cutoff=0.999):
    # returns the gamma block response starting at t=0 with intervals of delt = 1
    s = gamma_step(A, n, a, cutoff)
    return np.append(s[0], s[1:] - s[:-1])


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

# Making a full copy and keeping subset
annual_pump2 = annual_pump.copy()
annual_pump = annual_pump[annual_pump.index <= "2023"]
listdaily_pump = listdaily_pump[listdaily_pump.index <= "2023"]

# Make 1970-1990 arbitrary lower
adjustedann_pump = annual_pump.copy()
adjustedann_pump["Pump"].loc["1980":"1990"] /= 2

# Generating ensemble of time series
for i in range(1):

    # Setting pumping std to be less
    adjustedann_pump["Std"] = adjustedann_pump["Pump"] * .2
    pumptrue = adjustedann_pump.Pump

    # Plotting
    # plt.figure()
    # plt.plot(listdaily_pump.index, listdaily_pump.Pump, "-ok",
    #          linewidth=3, label="Initial Basin-wide (Best Guess)")
    # plt.plot(adjustedann_pump.index, adjustedann_pump.Pump, "-o",
    #          color="tab:orange",
    #          linewidth=3, label="Altered 1970-1980 (Truth)")
    # plt.xlabel("Date", fontsize=14)
    # plt.ylabel("Pumping rate * 10,000 (m$^{3}$/day)", fontsize=14)
    # plt.legend()

# Interpolate pumping
df = pd.DataFrame(index=listdaily_pump.index)
df = pd.concat([df, pumptrue], join="outer",
               keys=["Date", "Date"], axis=1)
df.columns = df.columns.droplevel()

# Interpolating pumping data
pumptrue_interp = df.interpolate(method="cubic")
pumptrue_interp = pumptrue_interp.dropna()
step = gamma_block(Atrue, ntrue, atrue)[1:]
lenstep = len(step)
h_pump = dtrue * np.ones(len(pumptrue_interp) + lenstep)

# Calculating head
for i in range(len(pumptrue_interp)):

    h_pump[i: i + lenstep] += pumptrue_interp.iloc[
        i].values[0] * step
head_pump = pd.DataFrame(
    index=pumptrue_interp.index,
    data=h_pump[: len(pumptrue_interp)],
)

# Plotting true head
# plt.figure(figsize=(12, 5))
# plt.plot(head_pump, "k.", label="True Head")
# plt.legend(loc=0)
# plt.ylabel("Head (m)")
# plt.xlabel("Time (years)")
# plt.title("True Head vs Head Observations")

random_seed = np.random.RandomState(15892)

# Saving model under fake wellnestname
wellnestname = "LCBKK011"
wellnestlist = [wellnestname]

# Reading in groundwater data
full_path = os.path.join(tot_path, wellnestname + ".xlsx")
data = pd.read_excel(full_path, skiprows=3)

# Well names
well_names = data.columns[-(len(data.columns)-2):]

# Saving gw_obs
gw_obs_list = []

# For all wells in well nest
for wells in well_names:

    # Name of well as a string
    well_name = wells

    # Adding nosie
    obs_std = 0.5
    noise = random_seed.normal(
        0, obs_std, len(head_pump)) * np.std(head_pump.values) * 0.5
    head_pump_noise = head_pump[0] + noise

    # Head with noise within certain years
    head_pump_noise = head_pump_noise[np.logical_and(
        head_pump_noise.index >= calitime_min, head_pump_noise.index <= calitime_max)]

    # Random heads for only sample size
    obs_sample = 300
    head_pump_noise = head_pump_noise.sample(obs_sample).sort_index()
    gw_obs_list.append(head_pump_noise)

    # Initializing Pastas model
    # BEST PASTAS MODEL
    # WITH OBS AND WRONG BASIN WIDE PUMPING
    ml3 = ps.Model(head_pump_noise)
    # Creating stress model
    ml3.del_noisemodel()
    EstTotPump_ = ps.StressModel(annual_pump.Pump, rfunc=pump_rfunc,
                                 name="well", settings="well",
                                 up=False)
    ml3.add_stressmodel(EstTotPump_)
    ml3.solve(tmin=calitime_min, tmax=calitime_max)

    # If saving model
    if save_model == 1:
        ml3.to_file(modelpath + "/" + wellnestname + "_" + well_name +
                    "_GW_" + calitime_min + "_" + calitime_max +
                    "_SYNBESTMODEL.pas")

modelpath = os.path.abspath("models//synthetic//na16ne250//perfect")

# For all wells in well nest
for wells in well_names:

    # Name of well as a string
    well_name = wells

    # Initializing Pastas model
    # Truth PASTAS MODEL
    # WITH OBS AND true BASIN WIDE PUMPING
    ml4 = ps.Model(head_pump)
    # Creating stress model
    ml4.del_noisemodel()
    EstTotPump_ = ps.StressModel(pumptrue_interp.Pump, rfunc=pump_rfunc,
                                 name="well", settings="well",
                                 up=False)
    ml4.add_stressmodel(EstTotPump_)
    ml4.solve(tmin=calitime_min, tmax=calitime_max)

    # If saving model
    if save_model == 1:
        ml4.to_file(modelpath + "/" + wellnestname + "_" + well_name +
                    "_GW_" + calitime_min + "_" + calitime_max +
                    "_SYNTRUTH.pas")

# Parameter Boundaries!
parambound_path = os.path.join(os.path.abspath("inputs"),
                               "SUBParametersPriortoManual.xlsx")

parambound = pd.read_excel(parambound_path,
                           sheet_name="bounds_mult",
                           index_col=0)
parambound = pd.DataFrame(parambound)

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

# Saving synthetic truth multipliers
syntruth = []

# True multipler for each well nest
for wellnest in wellnestlist:

    Sskv_data.loc[wellnest] *= a

    Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
    Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

    K_data.loc[wellnest] *= c

    syntruth.append([a, b, c])

# Mode can be "raw" as in raw groundwater data vs "Pastas" for importing Pastas
# simulated groundwater in the aquifers
mode = "Pastas"

# If mode is Pastas, need model path
if mode == "Pastas":

    mpath = os.path.abspath("models//synthetic//na16ne250//perfect")

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
                           model_path=mpath, califlag=0,
                           esmdaflag=0)

# Post process data
sub_total, subv_total, ann_sub, \
    avgsub = bkk_sub_gw.bkk_sub.bkk_postproc(wellnestlist,
                                             sub_total,
                                             subv_total,
                                             all_results)

# Saving annual subsidnece rates
annual_data = ann_sub

# True subsidence
syndata = pd.DataFrame()
# Observed subsdience
obssyndata = pd.DataFrame()

# Subsidence error
sub_error = 1

# Preallocation for truth
truth = []
obs_noise = []

# Creating observations for ESMDA
dobs = pd.Series(np.empty(1, dtype=object))
# Saving gw obs index
gw_obs_indices = []

wellnestlist = ["LCBKK011"]

for num_well, wellnest in enumerate(wellnestlist):

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

    # Gets the last date of each year
    lastdate = bench.groupby(pd.DatetimeIndex(bench["date"]).year,
                             as_index=False).agg(
                                 {"date": max}).reset_index(drop=True)
    bench = bench.loc[lastdate.date]

    # Getting rid of initial 0
    bench = bench.iloc[1:, :]
    syndata[wellnest] = bench.iloc[:, 0]

    plot_data = df.merge(annual_data[num_well][1]*100, left_on=df.date,
                         right_on=annual_data[num_well][1].index,
                         how="left")

    plot_data = plot_data.rename(columns={"key_0": "key0"})

    plot_data = plot_data.merge(bench, left_on=plot_data.key0,
                                right_on=bench.index,
                                how="left")
    # Renaming for other merge
    plot_data = plot_data.rename(columns={"key_0": "key1"})

    for idx, value in enumerate(syndata[wellnest]):

        if ~np.isnan(value):

            matchdate = plot_data.date_x == syndata.index[idx]

            if any(matchdate):

                syndata[wellnest][idx] = plot_data.AnnRates[
                    matchdate]

    obssyndata[wellnest] = bench.iloc[:, 0]
    for idx, value in enumerate(obssyndata[wellnest]):

        if ~np.isnan(value):

            matchdate = plot_data.date_x == obssyndata.index[idx]

            if any(matchdate):

                obssyndata[wellnest][idx] = plot_data.AnnRates[
                    matchdate] + random.uniform(0, sub_error *
                                                np.mean(plot_data.AnnRates))

    syndata["Year"] = syndata.index.year
    syndata = syndata.rename(columns={wellnest: 'Land_' + wellnest})

    obssyndata["Year"] = obssyndata.index.year
    obssyndata = obssyndata.rename(columns={wellnest: 'Land_' + wellnest})

    # Saving data
    truth.append(syndata)
    obs_noise.append(obssyndata)

# Adding subsidence data to observations
dobs = pd.concat([dobs, obssyndata.dropna().iloc[:, 0]])
gw_obs_indices.append(obssyndata.dropna().index)

# Saving truth
truth.append(head_pump)

# Saving obs
obs_noise.append(head_pump_noise)

# Adding groundwater observations to subsidence observations
for welli in range(len(data.columns[-(len(data.columns)-2):])):
    gw_obs_indices.append(gw_obs_list[welli].index)
    dobs = pd.concat([dobs, gw_obs_list[welli]])
dobs = dobs[1:]

gw_obs_df = pd.concat(gw_obs_list)

# Path to save models
savepath = os.path.abspath("models//synthetic//na16ne250")

# GW obs
fig_name1 = wellnestlist[0] + "_GWObs.csv"
full_figpath = os.path.join(savepath, fig_name1)
gw_obs_df.to_csv(full_figpath, sep="\t")
# GW Truth
fig_name1 = wellnestlist[0] + "_GWTruth.csv"
full_figpath = os.path.join(savepath, fig_name1)
head_pump.to_csv(full_figpath, sep="\t")
# Pumping truth
fig_name1 = wellnestlist[0] + "_PumpTruth.csv"
full_figpath = os.path.join(savepath, fig_name1)
pumptrue.to_csv(full_figpath, sep="\t")
# Sub obs
fig_name1 = wellnestlist[0] + "_SubObs.csv"
full_figpath = os.path.join(savepath, fig_name1)
obssyndata.dropna().iloc[:, 0].to_csv(full_figpath, index=True, sep='\t')
# Sub truth
fig_name1 = wellnestlist[0] + "_SubTruth.csv"
full_figpath = os.path.join(savepath, fig_name1)
syndata.dropna().iloc[:, 0].to_csv(full_figpath, index=True,
                                   sep='\t')

# %% ESMDA Implementation

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

    mpath = os.path.abspath("models//synthetic//na16ne250")

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

# If mode is Sskv
if p_multop[1] == "Sskv":
    Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
    Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

    K_data.loc[wellnest] *= c

# If mode is Sskv and K
elif p_multop[1] == "SsK":

    Sske_data.loc[wellnest][::2] = Sskv_data.loc[wellnest][::2] * b
    Sske_data.loc[wellnest][1::2] = Sske_data.loc[wellnest][0::2] / 10

# Ensemble size (Geir used for the figures 1e7)
# (Reduce to speed up; costly is the pd-estimation for plotting,
# not ES-MDA)
ne = int(250)
na = 16  # Number of assimilation steps
# obs_error = 3.5  # Synthetic originally
obs_error = 5.5
par_error.extend([30, 100, 100, 30])
dist = ["norm", "norm", "norm", "norm"]

# Esmda mode
esmdaflag = "my_gwparam_subSELECT_pump"

# Index of interested parameters
param_index = np.array([1, 2])

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

# Number of ensembles
n = 1
for n_ens in range(n):
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
                                                           listdaily_pump=listdaily_pump)

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

    # Write ESMDA results
    for na_i in range(na+1):
        fig_name1 = wellnestlist[0] + "_Dpred_" + p_multop[1] + "_na" + str(na_i) + ".csv"
        full_figpath = os.path.join(savepath, fig_name1)
        np.savetxt(full_figpath, sub_m["dpred"][na_i], delimiter=",")

        fig_name1 = wellnestlist[
            0] + "_Mprior_" + p_multop[1] + "_na" + str(na_i) + ".csv"
        full_figpath = os.path.join(savepath, fig_name1)
        np.savetxt(full_figpath, sub_m["mprior"][na_i], delimiter=",")
