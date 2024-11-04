# ##############################################################################
"""Plotting subsidence and groundwater model results.

Output:
Subsidence
1. Bar graphs of annual subsidence (cm) for each well nest during 1978-2020
(Shown in the main text and supplemental information)
2. Line graphs of annual subsidence (cm) for sensitivity analyses of each parameter
(Sskv, Sske, K, thickness) for one well nest (long run time so only calculating for
one well nest at a time) (Shown in supplemental information)
3. Line graphs of cumulative subsidence (cm) into the future depending on the
pumping scenario for each well nest during 1978-2060 (Shown in the main text and
supplemental information)
4. Spatial maps showing cumulative subsidence for each scenario for each
well nest
Groundwater
4. Groundwater model graphical results shown in the paper and supplemental
information. Simulated vs observed groundwater and inputs
5. Spatial maps that show the groundwater RMSE and t190 results for each well in
each well nest. Imports previously created Pastas models

Article Title: Hybrid data-driven, physics-based modeling of ground-
water and subsidence with application to Bangkok, Thailand

Jenny Soonthornrangsan 2023
TU Delft

"""
# ##############################################################################

###############################################################################
# import statements
###############################################################################

import os
import pandas as pd
import datetime as dt
from sklearn.metrics import mean_squared_error
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mycolorpy import colorlist as mcp
from matplotlib.ticker import (AutoMinorLocator)
import matplotlib.ticker as mticker
from matplotlib.collections import PatchCollection
from matplotlib.patches import Wedge
from statistics import median
import string
import matplotlib.lines as mlines
import matplotlib.ticker as ticker


# %%###########################################################################
# Plotting settings
###############################################################################

plt.rc("font", size=12)  # controls default text size
plt.rc("axes", titlesize=5)  # fontsize of the title
plt.rc("axes", labelsize=6)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=6)  # fontsize of the x tick labels
plt.rc("ytick", labelsize=6)  # fontsize of the y tick labels
plt.rc("legend", fontsize=8)  # fontsize of the legend


# %%###########################################################################
# Plotting results
##############################################################################

def nse(predictions, targets):
    return 1-(np.sum((targets-predictions)**2) /
              np.sum((targets-np.mean(targets))**2))


def InSAR_preproc(wellnest, tmin, tmax, ann=0):
    """Function to preprocess InSAR data
    Inputs:
    wellnest - wellnest name in strings
    tmin - model time min
    tmax - model time max
    ann - annual (1) or cumulative (0) data

    Outputs:
    ERS - dataframe with preprocessed ERS data for one well nest
    Sent - dataframe with preprocessed Sentinel data for one well nest
    merged - merged dataframe of ERS and Sentinel data
    """

    # If dealing with annual data
    if ann == 1:
        locERS = os.path.join(os.path.abspath("inputs"),
                              "ERS_InSAR_GWWells_AnnSub500m.xlsx")
        ERS_data = pd.read_excel(locERS, index_col=0)

        locSent = os.path.join(os.path.abspath("inputs"),
                               "Sentinel_InSAR_GWWells_AnnSub500m.xlsx")
        Sent_data = pd.read_excel(locSent, index_col=0)

        # Cumulative rates for that year
        ERS = ERS_data.groupby(pd.DatetimeIndex(ERS_data.index.year),
                               as_index=False).sum()
        ERS.index = np.unique(ERS_data.index.year)
        ERS.index = pd.to_datetime(ERS.index, format="%Y")
        ERS.index = ERS.index.shift(+1, freq="Y")

        # Offsetting InSAR dates by one year to add starting
        # value of 0 of previous year before starting date
        # Adds an extra year to the end
        ERS.loc[ERS.index[0] +
                pd.offsets.DateOffset(years=-1)] = 0
        ERS = ERS.sort_index()

        # Replacing no data at a well nest with nan
        for i in range(len(ERS.columns)):

            # If no data
            if (ERS.iloc[:, i] == 0).all():
                ERS.iloc[:, i] = ERS.iloc[:, i].replace(0, np.nan, inplace=False)

        # mm to cm; subsidence is pos
        ERS = ERS.divide(-10)

        Sent = Sent_data.groupby(pd.DatetimeIndex(Sent_data.index.year),
                                 as_index=False).sum()
        Sent.index = np.unique(Sent_data.index.year)
        Sent.index = pd.to_datetime(Sent.index, format="%Y")
        Sent.index = Sent.index.shift(+1, freq="Y")

        # Replacing no data at a well nest with nan
        for i in range(len(Sent.columns)):

            # If no data
            if (Sent.iloc[:, i] == 0).all():
                Sent.iloc[:, i] = Sent.iloc[:, i].replace(0, np.nan, inplace=False)

        # mm to cm; subsidence is pos
        Sent = Sent.divide(-10)

        # Merging both
        merged = pd.concat([ERS[wellnest], Sent[wellnest]])

        return merged, ERS[wellnest], Sent[wellnest]

    # If cumulative
    else:

        locERS = os.path.join(os.path.abspath("inputs"),
                              "ERS_InSAR_GWWells_CumSub500m.xlsx")
        ERS_data = pd.read_excel(locERS, index_col=0)

        locSent = os.path.join(os.path.abspath("inputs"),
                               "Sentinel_InSAR_GWWells_CumSub500m.xlsx")
        Sent_data = pd.read_excel(locSent, index_col=0)

        # Gets last value of each year cuase thats the cum sum so far
        N = int(np.ceil((ERS_data.index.max()-ERS_data.index.min()) /
                        np.timedelta64(1, "Y"))) + 2
        bins = [dt.datetime.strptime(str(ERS_data.index.year.max()) + str(1231),
                                     "%Y%m%d") -
                pd.offsets.DateOffset(years=y)
                for y in range(N)][::-1]

        # New annual dataframe
        ERS_annual = ERS_data.groupby(pd.cut(ERS_data.index, bins)).last()
        ERS_annual["year"] = ERS_data.index.year.unique()
        ERS_annual.index = pd.to_datetime(ERS_annual.year, format="%Y")
        ERS_annual.index = ERS_annual.index.shift(+1, freq="Y")

        # Gets last value of each year cuase thats the cum sum so far
        N = int(np.ceil((Sent_data.index.max()-Sent_data.index.min()) /
                        np.timedelta64(1, "Y"))) + 1
        bins = [dt.datetime.strptime(str(Sent_data.index.year.max()) + str(1231),
                                     "%Y%m%d") -
                pd.offsets.DateOffset(years=y)
                for y in range(N)][::-1]

        # New annual dataframe
        Sent_annual = Sent_data.groupby(pd.cut(Sent_data.index, bins)).last()
        Sent_annual["year"] = Sent_data.index.year.unique()
        Sent_annual.index = pd.to_datetime(Sent_annual.year, format="%Y")
        Sent_annual.index = Sent_annual.index.shift(+1, freq="Y")

        # Replacing no data at a well nest with nan
        for i in range(len(ERS_annual.columns)):

            # If no data
            if (ERS_annual.iloc[:, i] == 0).all():
                ERS_annual.iloc[:, i] = ERS_annual.iloc[:, i].replace(0, np.nan,
                                                                      inplace=False)

        # mm to cm; subsidence is pos
        ERS_annual = ERS_annual.divide(-10)

        # Replacing no data at a well nest with nan
        for i in range(len(Sent_annual.columns)):

            # If no data
            if (Sent_annual.iloc[:, i] == 0).all():
                Sent_annual.iloc[:, i] = Sent_annual.iloc[
                    :, i].replace(0, np.nan, inplace=False)

        # mm to cm; subsidence is pos
        Sent_annual = Sent_annual.divide(-10)

        # Merging both
        merged = pd.concat([ERS_annual[wellnest], Sent_annual[wellnest]])

        return merged, ERS_annual[wellnest], Sent_annual[wellnest]


# Function to preprocess GPS data
def GPS_preproc(wellnest, tmin, tmax):
    """Function to preprocess GPS data
    Inputs:
    wellnest - wellnest name in strings
    tmin - model time min
    tmax - model time max
    Outputs:
    GPSdata - dataframe with preprocessed GPS data
    GPSVar - dataframe with GPS variances
    Identifies which GPS station to use with which well
    """
    loc = os.path.join(os.path.abspath("inputs"),
                       "GPS_GWNest.xlsx")
    catalog = pd.read_excel(loc, sheet_name="WellGPS",
                            index_col=0)

    # If GPS station exists for this well nest
    try:
        sheet = catalog.Station.loc[wellnest]

    except:

        return None, None

    # Getting GPS time series
    GPSdata = pd.read_excel(loc, sheet_name=sheet,
                            index_col=1)
    GPSdata = pd.DataFrame(GPSdata)
    GPSdata.index = pd.to_datetime(GPSdata.index)

    # Getting rid of benchmarks outside time period
    GPSdata = GPSdata[(GPSdata.index.year >= int(tmin)) &
                      (GPSdata.index.year <= int(tmax))]

    # First data of GPS will be 0, new relative height calculated by subtracting
    # the reference value from all other values
    # Cumulative diff because each day is referenced to
    # the start
    GPSdata["NewRelHei"] = GPSdata.RelativeHeight_cm - GPSdata.RelativeHeight_cm[0]

    # Rates of each day
    GPSdata["DailyDiff"] = GPSdata["NewRelHei"].diff().fillna(GPSdata["NewRelHei"])

    # SD of diff calculated as sqrt(var 1 + var2)
    # https://study.com/skill/learn/how-to-calculate-the-standard-deviation-of-the-difference-of-two-random-variables-explanation.html
    GPSdata["DailyDiffVar"] = GPSdata["SD_cm"] * 0
    GPSdata["DailyDiffVar"][1:] = [x**2 + GPSdata["SD_cm"][index - 1]**2
                                   for index, x in enumerate(GPSdata["SD_cm"][1:])]

    # Annual sum of variances
    GPSVar_annual = GPSdata.DailyDiffVar.groupby(GPSdata.index.year).sum()

    # Annual sum of rates
    GPS_annual = GPSdata.DailyDiff.groupby(GPSdata.index.year).sum()
    # Converts year index to datetime
    GPS_annual.index = pd.to_datetime(GPS_annual.index, format="%Y")
    GPSVar_annual.index = pd.to_datetime(GPS_annual.index, format="%Y")

    # Offsetting GPS dates by one year to add starting
    # value of 0 of previous year before starting date
    # Adds an extra year to the end
    GPS_annual.loc[GPS_annual.index[-1] +
                   pd.offsets.DateOffset(years=1)] = 0
    GPSVar_annual.loc[GPSVar_annual.index[-1] +
                      pd.offsets.DateOffset(years=1)] = 0
    GPS_annual = GPS_annual.shift(1)  # Shifts all values down one year
    GPSVar_annual = GPSVar_annual.shift(1)   # Shifts all values down one year
    GPS_annual.iloc[0] = 0  # Sets first value as 0
    GPSVar_annual.iloc[0] = 0  # Sets first value as 0
    GPS_annual.index = GPS_annual.index.shift(-1, freq="D")
    GPSVar_annual.index = GPSVar_annual.index.shift(-1, freq="D")

    return GPS_annual, GPSVar_annual


def sub_bar(path, wellnestlist, all_results,
            sub_total, subv_total,
            annual_data, tmin=None, tmax=None, save=0,
            benchflag=0, insarflag=0, GPSflag=0):
    """Plot annual subsidence results.

    Bar graphs of annual subsidence (cm) for each well nest during 1978-2020
    (Shown in the main text and supplemental information)

    path - str: path to save figures
    wellnestlist - list of wellnests that were simualted
    all_results - lists of lists: wellnestname, well name, head matrix with
    head series of each node
    sub_total - lists of lists: wellnestname, well_name, total cum sub
    subv_total - lists of lists: wellnestname, well_name, inelastic cum sub
    annual_data - lists of lists: wellnestname, well_name, total cum sub for
    all four clay at a wellnest location
    save - if 1, save; if 0, don't save
    benchflag: no benchmark - if 0, no plot, if 1, benchmark, plot
    insarflag: no insar - if 0, no plot, if 1, insar, plot
    GPS flag: no GPS - if 0, no plot, if 1, GPS, plot
    Assume also that benchmark comparison starts at 0

    # ASSUMES FOUR WELLS IN WELLNEST
    """
    # saving rmse
    rmse = []

    # For each wellnest in list
    # num_well is the index, wellnest = name
    # Figures for each well nest
    for num_well, wellnest in enumerate(wellnestlist):

        if benchflag == 1:

            # # Averaged benchmarks
            # loc = os.path.join(os.path.abspath("inputs"),
            #                    "SurveyLevel_GW_AvgRates5m.xlsx")
            # try:
            #     subdata = pd.read_excel(loc, index_col=0)
            #     subdata = pd.DataFrame(subdata)
            #     subdata.index = pd.to_datetime(subdata.index)
            #     # Benchmarks should start at 0 at the first year.
            #     bench = subdata.loc[:, subdata.columns.str.contains(wellnest)]

            # Subsidence plotting
            # Getting benchmark time series
            loc = os.path.join(os.path.abspath("inputs"), "SurveyingLevels.xlsx")
            try:
                subdata = pd.read_excel(loc, sheet_name=wellnest + "_Leveling",
                                        index_col=3)
                subdata = pd.DataFrame(subdata)
                subdata.index = pd.to_datetime(subdata.index)

                # Getting rid of benchmarks outside time period
                subdata = subdata[(subdata.Year <= 2020)]

                # Benchmarks should start at 0 at the first year.
                bench = subdata.loc[:, subdata.columns.str.contains("Land")]
                bench = bench.fillna(0)

                if (bench.iloc[0] != 0).any():
                    bench.iloc[0] = 0

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

            except:

                bench = pd.DataFrame()

        # InSAR plotting
        if insarflag == 1:

            InSAR, _, _ = InSAR_preproc(wellnest, tmin, tmax, ann=1)

        # GPS plotting
        if GPSflag == 1:

            GPS, GPS_Var = GPS_preproc(wellnest, tmin, tmax)

        # BAR PLOT preparation
        daterange = pd.date_range(dt.datetime(1978, 12, 31), periods=43,
                                  freq="Y").tolist()
        df = pd.DataFrame(daterange, columns=["date"])

        x = np.arange(43)

        # If insar, bench, gps (various combos)
        if sum([insarflag, GPSflag, benchflag]) == 1:

            width = .5

        elif sum([insarflag, GPSflag, benchflag]) == 2:

            width = .333

        elif sum([insarflag, GPSflag, benchflag]) == 0:

            width = 1

        elif sum([insarflag, GPSflag, benchflag]) == 3:

            width = .25

        # Figure plotting model results against measurements
        # Converts to cm to match measurements
        # set fig size certain way if running batch well nests
        # Supplemental Information
        if len(wellnestlist) > 1:

            plt.figure(figsize=(6.75, 3.38), dpi=400)

        # Paper size
        else:

            plt.figure(figsize=(6.75, 2), dpi=400)

        # Bar graph
        # annual data in cm
        plot_data = df.merge(annual_data[num_well][1]*100, left_on=df.date,
                             right_on=annual_data[num_well][1].index,
                             how="left")
        # Renaming for second merge
        plot_data = plot_data.rename(columns={"key_0": "key0"})

        # Filling na with 0
        plot_data = plot_data.fillna(0)

        plt.bar(x,
                -plot_data.AnnRates,
                label="Analysis", width=width,
                linewidth=.5, edgecolor="k")

        # plt.plot(x,
        #          -plot_data.AnnRates, "-ob",
        #          label="Analysis",
        #          linewidth=.5, markersize=1)

        # Plotting benchmarks
        if benchflag == 1:

            if not bench.empty:

                # Measurements
                # Bar plot
                # Benchamrks already in cm
                plot_data = plot_data.merge(bench, left_on=plot_data.key0,
                                            right_on=bench.index,
                                            how="left")
                # Renaming for other merge
                plot_data = plot_data.rename(columns={"key_0": "key1"})

                # Filling na with 0
                plot_data = plot_data.fillna(0)

                plt.bar(x+width, -plot_data[
                    plot_data.columns[
                        plot_data.columns.str.contains("Land")].item()],
                        color="orange", linewidth=.5,
                        label="Observed", width=width, edgecolor="k")

                # plt.plot(x+width, -plot_data[
                #     plot_data.columns[
                #          plot_data.columns.str.contains("Land")].item()],
                #          "-o", markersize=1,
                #          color="orange", linewidth=.5,
                #          label="Observed")

                # plt.bar(x+width, -plot_data[
                #     plot_data.columns[
                #         plot_data.columns.str.contains(wellnest)].item()],
                #         color="orange", linewidth=.5,
                #         label="Observed", width=width, edgecolor="k")

            # Dropping NAs
            plot_data = plot_data.dropna()

            rms_data = plot_data[
                plot_data.columns[
                    plot_data.columns.str.contains("Land")].item()]

            # Calculating RMSE
            rms = mean_squared_error(rms_data[rms_data != 0],
                                     plot_data.AnnRates[rms_data != 0],
                                     squared=False)

            print(rms)

            # rms = mean_squared_error(plot_data[
            #     plot_data.columns[
            #         plot_data.columns.str.contains(wellnest)].item()],
            #     plot_data.AnnRates, squared=False)

            plt.annotate("RMSE: " + "{:.1f}".format(rms) + " cm/year",
                         xy=(.99, .97), xycoords="axes fraction",
                         fontsize=10, horizontalalignment="right",
                         verticalalignment="top")

            # saving rmse
            rmse.append(rms)

        # If insar data
        if insarflag == 1:

            # If data for either insars
            if ~InSAR.isna().all():

                # Measurements
                # Bar plot
                # Benchamrks already in cm
                plot_data = plot_data.merge(InSAR, left_on=plot_data.key0,
                                            right_on=InSAR.index,
                                            how="left")
                # Renaming for other merge
                plot_data = plot_data.rename(columns={"key_0": "key2"})

                # Filling na with 0
                plot_data = plot_data.fillna(0)

                # If insar, bench, gps (various combos), affecting where to
                # plot x and width
                if sum([insarflag, GPSflag, benchflag]) == 1:

                    actual_x = x + width

                elif sum([insarflag, GPSflag, benchflag]) == 2:

                    actual_x = x + 2 * width

                # in case
                else:

                    actual_x = x + width

                plt.bar(actual_x, plot_data[wellnest],
                        color="purple", linewidth=.5,
                        label="InSAR", width=width, edgecolor="k")

        # If GPS data
        if GPSflag == 1:

            # If data for either insars
            if not isinstance(GPS, type(None)):

                # Converting sub to positive
                GPS *= -1

                # Measurements
                # Bar plot
                # GPS already in cm
                plot_data = plot_data.merge(GPS, left_on=plot_data.key0,
                                            right_on=GPS.index,
                                            how="left")
                # Renaming for other merge
                plot_data = plot_data.rename(columns={"key_0": "key3"})

                # Filling na with 0
                plot_data = plot_data.fillna(0)

                # If insar, bench, gps (various combos), affecting where to
                # plot x and width
                if sum([insarflag, GPSflag, benchflag]) == 1:

                    actual_x = x + width

                elif sum([insarflag, GPSflag, benchflag]) == 2:

                    actual_x = x + 2*width

                elif sum([insarflag, GPSflag, benchflag]) == 3:

                    actual_x = x + 3*width

                # SD
                plot_data = plot_data.merge(np.sqrt(GPS_Var),
                                            left_on=plot_data.key0,
                                            right_on=GPS_Var.index,
                                            how="left")
                # Renaming for other merge
                plot_data = plot_data.rename(columns={"key_0": "key4"})

                # Filling na with 0
                plot_data = plot_data.fillna(0)

                plt.bar(actual_x, plot_data["DailyDiff"],
                        yerr=plot_data["DailyDiffVar"],
                        error_kw=dict(lw=.35, capsize=2, capthick=.35),
                        color="red", linewidth=.5,
                        label="GPS", width=width, edgecolor="k")

        # Plotting settings
        plt.legend(loc="lower left")
        ax = plt.gca()
        plt.draw()
        plt.axhline(y=0, color="k", linestyle="-", linewidth=1)
        ax.set_xticklabels(ax.get_xticks(), rotation=45)

        # If insar, bench, gps (various combos), affecting where to
        # plot x and width
        if sum([insarflag, GPSflag, benchflag]) == 1:

            actual_x = x + width

        elif sum([insarflag, GPSflag, benchflag]) == 2:

            actual_x = x + 2*width

        elif sum([insarflag, GPSflag, benchflag]) == 3:

            actual_x = x + 3*width

        else:
            actual_x = x + width
        plt.xticks(actual_x, ["1978", "", "1980", "", "1982",
                              "", "1984", "", "1986", "",
                              "1988", "", "1990", "", "1992",
                              "", "1994", "", "1996", "",
                              "1998", "", "2000", "", "2002",
                              "", "2004", "", "2006", "",
                              "2008", "", "2010", "", "2012",
                              "", "2014", "", "2016", "",
                              "2018", "", "2020"])

        # set y limits/title only if running batch well nests
        if len(wellnestlist) > 1:
            # plt.ylim((-2, 10))
            plt.title(wellnest)
        plt.ylabel("Annual Subsidence Rate (cm/yr)")
        plt.xlabel("Years")
        # plt.ylim((-.5, 5))
        # Setting fig size again
        # set fig size certain way if running batch well nests
        # Supplemental Information
        if len(wellnestlist) > 1:

            plt.gcf().set_size_inches(6.75, 3.38)

        # Paper size
        else:

            plt.gcf().set_size_inches(6.75, 2)

        # If saving figure
        if np.logical_and(save == 1, benchflag == 1):

            # set name of file certain way if running batch well nests
            if len(wellnestlist) > 1:

                # fig_name = wellnest + "_BenchvsImplicit_AnnSubTotal.eps"
                # full_figpath = os.path.join(path, fig_name)
                # plt.savefig(full_figpath, format="eps", bbox_inches="tight")

                fig_name = wellnest + "_BenchvsImplicit_AnnSubTotal.png"
                full_figpath = os.path.join(path, fig_name)
                plt.savefig(full_figpath, format="png", bbox_inches="tight")

            else:

                # fig_name = wellnest + "_BenchvsImplicit_AnnSubTotal_PAPER.eps"
                # full_figpath = os.path.join(path, fig_name)
                # plt.savefig(full_figpath, format="eps", bbox_inches="tight")

                fig_name = wellnest + "_BenchvsImplicit_AnnSubTotal_PAPER.png"
                full_figpath = os.path.join(path, fig_name)
                plt.savefig(full_figpath, format="png", bbox_inches="tight")

        else:
  
            fig_name = wellnest + "_AnnSubTotal.png"
            full_figpath = os.path.join(path, fig_name)
            plt.savefig(full_figpath, format="png", bbox_inches="tight")


def sub_plot(path, wellnestlist, all_results,
             sub_total, subv_total,
             annual_data, tmin=None, tmax=None, save=0,
             benchflag=0, insarflag=0, GPSflag=0):
    """Plot cumulative subsidence results.

    Line graphs of annual subsidence (cm) for each well nest during 1978-2020

    path - str: path to save figures
    wellnestlist - list of wellnests that were simualted
    all_results - lists of lists: wellnestname, well name, head matrix with
    head series of each node
    sub_total - lists of lists: wellnestname, well_name, total cum sub
    subv_total - lists of lists: wellnestname, well_name, inelastic cum sub
    annual_data - lists of lists: wellnestname, well_name, total cum sub for
    all four clay at a wellnest location
    save - if 1, save; if 0, don't save
    benchflag: no benchmark - if 0, no plot, if 1, benchmark, plot
    insarflag: no insar - if 0, no plot, if 1, insar, plot
    GPS flag: no GPS - if 0, no plot, if 1, GPS, plot
    Assume also that benchmark comparison starts at 0

    ASSUMES FOUR WELLS IN WELLNEST
    """

    # For each wellnest in list
    # num_well is the index, wellnest = name
    # Figures for each well nest
    for num_well, wellnest in enumerate(wellnestlist):

        if benchflag == 1:

            # # Averaged benchmarks
            # loc = os.path.join(os.path.abspath("inputs"),
            #                    "SurveyLevel_GW_AvgRates5m.xlsx")
            # try:
            #     subdata = pd.read_excel(loc, index_col=0)
            #     subdata = pd.DataFrame(subdata)
            #     subdata.index = pd.to_datetime(subdata.index)
            #     # Benchmarks should start at 0 at the first year.
            #     bench = subdata.loc[:, subdata.columns.str.contains(wellnest)]

            # Subsidence plotting
            # Getting benchmark time series
            loc = os.path.join(os.path.abspath("inputs"), "SurveyingLevels.xlsx")
            try:
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

            except:

                bench = pd.DataFrame()

        # InSAR plotting
        if insarflag == 1:

            InSAR, ERS, Sent = InSAR_preproc(wellnest, tmin, tmax, ann=0)

            # Subsidence is negative here
            InSAR *= -1
            ERS *= -1
            Sent *= -1

        # GPS plotting
        if GPSflag == 1:

            GPS, GPS_Var = GPS_preproc(wellnest, tmin, tmax)

        # Figure plotting model results against measurements
        # Converts to cm to match measurements
        plt.figure()
        plt.plot(annual_data[num_well][1].index,
                 annual_data[num_well][1].CumTotSum*100,
                 label="Simulated")

        # Plotting benchmarks
        if benchflag == 1:

            if not bench.empty:

                # First overlapping year
                overlapyr = next((i for i in bench.index if i in
                                  bench[~np.isnan(bench).any(axis=1)].index),
                                 None)

                # If overlapping year, but benchmark leveling does
                # not have data
                if np.logical_or(len(
                        bench[bench.index == overlapyr].values) == 0,
                        bench.empty):

                    # Adding starting simulated measurement to
                    # bench data
                    # To start at "0" together
                    benchdata = np.cumsum(bench) + 100 * annual_data[
                        num_well][1].CumTotSum[
                            annual_data[
                                num_well][1].index == bench.index[0]].values

                # If benchmark leveling and simulated data overlap and exists
                else:

                    # Adding starting benchmark measurement to
                    # bench data
                    # To start at "0" together
                    benchdata = np.cumsum(bench.iloc[:, 0]) + 100 * annual_data[
                        num_well][1].CumTotSum[
                            annual_data[
                                num_well][1].index == overlapyr].values[0]

                # Measurements
                plt.plot(bench.index, benchdata, "-o",
                         color="orange", label="Leveling")

        # Plottign GPS
        if GPSflag == 1:

            if not isinstance(GPS, type(None)):

                # If GPS and bench, check for first overlapping year
                if benchflag == 1:

                    # First overlapping year
                    overlapyr = next((i for i in GPS.index if i in
                                      bench[~bench.isnull().any(axis=1)].index),
                                     None)

                    # If overlapping year, but benchmark leveling does
                    # not have data
                    if np.logical_or(len(
                            bench[bench.index == overlapyr].values) == 0,
                            bench.empty):

                        # Adding starting simulated measurement to
                        # GPS data
                        # To start at "0" together
                        GPSdata = np.cumsum(GPS) + 100 * annual_data[
                            num_well][1].CumTotSum[
                                annual_data[
                                    num_well][1].index == GPS.index[0]].values

                    # If benchmark leveling and GPS data overlap and exists
                    else:

                        # Adding starting benchmark measurement to
                        # GPS data
                        # To start at "0" together

                        GPSdata = np.cumsum(GPS) + benchdata.loc[
                            bench.index == overlapyr].values[0]

                # Measurements
                plt.plot(GPSdata.index, GPSdata, "-o", color="red",
                         label="GPS")
                plt.errorbar(GPSdata.index, GPSdata, yerr=GPS_Var,
                             lw=.35, capsize=2, capthick=.35,
                             color="black")

        # Plottign INSAR
        if insarflag == 1:

            # If data for either insars
            if ~InSAR.isna().all():

                # If InSAR and bench, check for first overlapping year
                if benchflag == 1:

                    # First overlapping year
                    overlapyrERS = next((i for i in ERS.index if i in
                                        bench[~bench.isnull().any(axis=1)].index),
                                        None)
                    overlapyrSent = next((i for i in Sent.index if i in
                                          bench[~bench.isnull().any(axis=1)].index),
                                         None)

                    # If overlapping year, but benchmark leveling does
                    # not have data
                    if np.logical_or(len(
                            bench[bench.index == overlapyrERS].values) == 0,
                            bench.empty):

                        # Adding starting simulated measurement to
                        # GPS data
                        # To start at "0" together
                        ERSdata = ERS + 100 * annual_data[
                            num_well][1].CumTotSum[annual_data[
                                num_well][1].index == ERS.index[0]].values

                    # If benchmark leveling and InSARdata data overlap and exists
                    else:

                        # Adding starting benchmark measurement to
                        # InSARdata data
                        # To start at "0" together

                        ERSdata = ERS + benchdata.loc[
                            bench.index == overlapyrERS].values[0]

                    # If overlapping year, but benchmark leveling does
                    # not have data
                    if np.logical_or(len(
                            bench[bench.index == overlapyrSent].values) == 0,
                            bench.empty):

                        # Adding starting simulated measurement to
                        # GPS data
                        # To start at "0" together
                        Sentdata = Sent + 100 * annual_data[
                            num_well][1].CumTotSum[annual_data[
                                num_well][1].index == Sent.index[0]].values

                    # If benchmark leveling and InSARdata data overlap and exists
                    else:

                        # Adding starting benchmark measurement to
                        # InSARdata data
                        # To start at "0" together

                        Sentdata = Sent + benchdata.loc[
                            bench.index == overlapyrSent].values[0]

                # Measurements
                plt.plot(ERSdata.index, ERSdata, "-o", color="purple",
                         label="InSAR")
                # Measurements
                plt.plot(Sentdata.index, Sentdata, "-o", color="purple",
                         label="InSAR")

        plt.legend()
        plt.ylabel("Cumulative Subsidence (cm)")
        plt.xlabel("Years")
        plt.title(wellnest +
                  "\nObserved vs. Simulated Cumulative Subsidence")

        # If saving figure
        if save == 1:

            fig_name = wellnest + "_BenchvsImplicit_CumSubTotal.png"
            full_figpath = os.path.join(path, fig_name)
            plt.savefig(full_figpath, dpi=500, bbox_inches="tight")


def gwlocs_map(path, save=0):
    """Spatial mapping of groundwater well nest locations in BKK.

    path - path to save figures
    save - save or not save figures
    """
    # Importing spatial coordinates
    full_path = os.path.join(os.path.abspath("inputs"), "GroundwaterWellLocs.xls")
    gwwell_locs = pd.read_excel(full_path)

    # Locations of wellnests removing duplicates
    gwwell_locs = gwwell_locs.drop_duplicates("WellNest_Name", keep="first")

    # Unique well nests and locations only
    unique = []

    # Preallocation
    # Saving relevant xs, ys
    xs = []
    ys = []

    # Labels
    labels = []

    # list of wellnest names
    wellnestlist = ["LCBKK003",
                    "LCBKK005",
                    "LCBKK006",
                    "LCBKK007",
                    "LCBKK009",
                    "LCBKK011",
                    "LCBKK012",
                    "LCBKK013",
                    "LCBKK014",
                    "LCBKK015",
                    "LCBKK016",
                    "LCBKK018",
                    "LCBKK020",
                    "LCBKK021",
                    "LCBKK026",
                    "LCBKK027",
                    "LCBKK036",
                    "LCBKK038",
                    "LCBKK041",
                    "LCNBI003",
                    "LCNBI007",
                    "LCSPK007",
                    "LCSPK009"]

    # Getting rid of repeating wells and data points
    # zip joins x and y coordinates in pairs
    for x, y in zip(gwwell_locs.Long, gwwell_locs.Lat):

        # Check if x, y is unique
        if (x, y) not in unique:

            # Saves this location for plotting
            unique.append((x, y))

            # Label is well nest name
            label = gwwell_locs.loc[
                gwwell_locs.Long == x]["WellNest_Name"].tolist()

            # Specific well nest of interest
            if label[0] in wellnestlist:

                # Saving data
                xs.append(x)
                ys.append(y)
                labels.append(label[0])

            # If well nest irrelevant for this paper
            else:
                continue

    # Printing average subsidence rmse (cm/yr)
    # Initializing figure
    fig, ax = plt.subplots(figsize=(3.2, 2.2), dpi=400)
    datalim = None
    map = Basemap(llcrnrlon=100.3, llcrnrlat=13.4, urcrnrlon=100.8, urcrnrlat=14,
                  resolution="h", ellps="WGS84", lat_0=13.6, lon_0=100.4)
    draw_basemap(map, xs, ys, labels=labels, fig=fig, ax=ax,
                 datalim=datalim, mode="GW_WellNests", save=0,
                 figpath=path)

    # If saving figure
    if save == 1:

        fig_name1 = "Map_GWLocs.eps"
        full_figpath = os.path.join(path, fig_name1)
        plt.savefig(full_figpath, dpi=400, bbox_inches="tight", format="eps")

        fig_name1 = "Map_GWLocs.png"
        full_figpath = os.path.join(path, fig_name1)
        plt.savefig(full_figpath, dpi=400, bbox_inches="tight", format="png")


def sub_rmse_map(path, wellnestlist, all_results,
                 sub_total, subv_total,
                 annual_data, tmin=None, tmax=None, save=0):
    """Spatial mapping of simulated subsidence and observed.

    path - path to save figures
    wellnestlist - list of wellnests that were simualted
    all_results - lists of lists: wellnestname, well name, head matrix with
    head series of each node
    sub_total - lists of lists: wellnestname, well_name, total cum sub
    subv_total - lists of lists: wellnestname, well_name, inelastic cum sub
    annual_data - lists of lists: wellnestname, well_name, total cum sub for
    all four clay at a wellnest location
    save - if 1, save; if 0, don't save

    ASSUMES FOUR WELLS IN WELLNEST
    """
    # Importing spatial coordinates
    full_GWpath = os.path.join(os.path.abspath("inputs"),
                               "GroundwaterWellLocs.xls")
    gwwell_locs = pd.read_excel(full_GWpath)

    # Locations of wellnests; removing duplicates
    gwwell_locs = gwwell_locs.drop_duplicates("WellNest_Name", keep="first")

    # BAR PLOT preparation
    daterange = pd.date_range(dt.datetime(1978, 12, 31), periods=43,
                              freq="Y").tolist()
    df = pd.DataFrame(daterange, columns=["date"])

    # Saving relevant xs, ys, and cumsum
    xs = []
    ys = []
    cs_rmse = []

    # For each wellnest in list
    # num_well is the index, wellnest = name
    # Figures for each well nest
    for num_well, wellnest in enumerate(wellnestlist):

        # Benchmark
        loc = os.path.join(os.path.abspath("inputs"), "SurveyingLevels.xlsx")
        subdata = pd.read_excel(loc, sheet_name=wellnest + "_Leveling",
                                index_col=3)
        subdata = pd.DataFrame(subdata)
        subdata.index = pd.to_datetime(subdata.index)

        # Getting rid of benchmarks outside time period
        subdata = subdata[(subdata.Year <= 2020)]

        # Benchmarks should start at 0 at the first year.
        bench = subdata.loc[:, subdata.columns.str.contains("Land")]
        bench = bench.fillna(0)

        if (bench.iloc[0] != 0).any():
            bench.iloc[0] = 0

        # IMPORTANT INFO
        # For benchmark measurements, the first year is 0, the second year is
        # the compaction rate over that first year.
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
        # preparation
        daterange = pd.date_range(dt.datetime(1978, 12, 31), periods=43,
                                  freq="Y").tolist()
        df = pd.DataFrame(daterange, columns=["date"])

        # annual data in cm
        plot_data = df.merge(annual_data[num_well][1]*100, left_on=df.date,
                             right_on=annual_data[num_well][1].index,
                             how="left")
        # Renaming for second merge
        plot_data = plot_data.rename(columns={"key_0": "key0"})

        # Filling na with 0
        plot_data = plot_data.fillna(0)

        # Benchamrks already in cm
        plot_data = plot_data.merge(bench, left_on=plot_data.key0,
                                    right_on=bench.index,
                                    how="left")
        # Renaming for other merge
        plot_data = plot_data.rename(columns={"key_0": "key1"})

        # Filling na with 0
        plot_data = plot_data.fillna(0)

        plot_data = plot_data.dropna()

        landlevel = plot_data[
            plot_data.columns[
                plot_data.columns.str.contains(
                    "Land")].item()]

        # Calculating rmse
        rms = mean_squared_error(landlevel[landlevel != 0],
                                 plot_data.AnnRates[landlevel != 0],
                                 squared=False)
        cs_rmse.append(rms)

        x_ = gwwell_locs.Long[gwwell_locs.WellNest_Name == wellnest].item()
        y_ = gwwell_locs.Lat[gwwell_locs.WellNest_Name == wellnest].item()
        xs.append(x_)
        ys.append(y_)

    # Printing average subsidence rmse (cm/yr)
    # Initializing figure
    fig, ax = plt.subplots(figsize=(3.2, 2.2), dpi=400)
    datalim = [.6, 3.9]
    map = Basemap(llcrnrlon=100.3, llcrnrlat=13.4, urcrnrlon=100.8, urcrnrlat=14,
                  resolution="h", ellps="WGS84", lat_0=13.6, lon_0=100.4)
    draw_basemap(map, xs, ys, cs_rmse, fig=fig, ax=ax,
                 datalim=datalim, mode="Sub_RMSE", save=0,
                 time_min=tmin, time_max=tmax, figpath=path)
    print("Avg: " + str("%.2f" % np.average(cs_rmse)) + " cm/yr")
    print("Max: " + str("%.2f" % max(cs_rmse)) + " cm/yr")

    # If saving figure
    if save == 1:

        # fig_name1 = "Map_Sub_RMSE_" + tmin + "_" + tmax + "_50_2020.eps"
        # full_figpath = os.path.join(path, fig_name1)
        # plt.savefig(full_figpath, dpi=400, bbox_inches="tight", format="eps")

        fig_name1 = "Map_Sub_RMSE_" + tmin[0] + "_" + tmax[0] + "_50_2020.png"
        full_figpath = os.path.join(path, fig_name1)
        plt.savefig(full_figpath, dpi=400, bbox_inches="tight", format="png")


def ens_sub_plot(path, wellnestlist, all_results,
                 sub_total, subv_total,
                 annual_data, tmin=None, tmax=None, mode=None, n=None,
                 save=0,
                 benchflag=0, insarflag=0, GPSflag=0):
    """Sensitivity analysis on subsidence based on either Sskv, Sske, K, thickness.

    path - path to save figures
    wellnestlist - list of wellnests that were simualted
    all_results - lists of lists: wellnestname, well name, head matrix with
    head series of each node
    sub_total - lists of lists: wellnestname, well_name, total cum sub
    subv_total - lists of lists: wellnestname, well_name, inelastic cum sub
    annual_data - lists of lists: wellnestname, well_name, total cum sub for
    all four clay at a wellnest location
    mode - which parameter for ensemble (Sskv, Sske, K)
    n - number of ensemble members
    save - if 1, save; if 0, don't save
    benchflag: no benchmark - if 0, no plot, if 1, benchmark, plot
    insarflag: no insar - if 0, no plot, if 1, insar, plot
    GPS flag: no GPS - if 0, no plot, if 1, GPS, plot
    Assume also that benchmark comparison starts at 0

    ASSUMES FOUR WELLS IN WELLNEST
    """

    # For each wellnest in list
    # num_well is the index, wellnest = name
    # Figures for each well nest
    for num_well, wellnest in enumerate(wellnestlist):

        if benchflag == 1:

            # # Averaged benchmarks
            # loc = os.path.join(os.path.abspath("inputs"),
            #                    "SurveyLevel_GW_AvgRates.xlsx")
            # try:
            #     subdata = pd.read_excel(loc, index_col=0)
            #     subdata = pd.DataFrame(subdata)
            #     subdata.index = pd.to_datetime(subdata.index)
            #     # Benchmarks should start at 0 at the first year.
            #     bench = subdata.loc[:, subdata.columns.str.contains(wellnest)]

            # Subsidence plotting
            # Getting benchmark time series
            loc = os.path.join(os.path.abspath("inputs"), "SurveyingLevels.xlsx")
            try:
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

            except:

                bench = pd.DataFrame()

        # InSAR plotting
        if insarflag == 1:

            InSAR, ERS, Sent = InSAR_preproc(wellnest, tmin, tmax, ann=0)

            # Subsidence is negative here
            InSAR *= -1
            ERS *= -1
            Sent *= -1

        # GPS plotting
        if GPSflag == 1:

            GPS, GPS_Var = GPS_preproc(wellnest, tmin, tmax)

        # Figure plotting model results against measurements
        # Converts to cm to match measurements
        plt.figure(figsize=(6.75, 3.38), dpi=400)

        # For each ensemble member + calibrated param
        for i in range(n+1):

            # If first value, calibrated parameters
            # after i = 0, members
            if i == 0:
                plt.plot(annual_data[i][num_well][1].index,
                         annual_data[i][num_well][1].CumTotSum*100,
                         color="blue", label="Simulated")
            else:
                plt.plot(annual_data[i][num_well][1].index,
                         annual_data[i][num_well][1].CumTotSum*100,
                         color="grey")

        # Plotting benchmarks
        if benchflag == 1:

            if not bench.empty:

                # First overlapping year
                overlapyr = next((i for i in bench.index if i in
                                  bench[~np.isnan(bench).any(axis=1)].index),
                                 None)

                # If overlapping year, but benchmark leveling does
                # not have data
                if np.logical_or(len(
                        bench[bench.index == overlapyr].values) == 0,
                        bench.empty):

                    # Adding starting simulated measurement to
                    # bench data
                    # To start at "0" together
                    benchdata = np.cumsum(bench) + 100 * annual_data[0][
                        num_well][1].CumTotSum[
                            annual_data[0][
                                num_well][1].index == bench.index[0]].values

                # If benchmark leveling and simulated data overlap and exists
                else:

                    # Adding starting benchmark measurement to
                    # bench data
                    # To start at "0" together
                    benchdata = np.cumsum(bench.iloc[:, 0]) + 100 * annual_data[0][
                        num_well][1].CumTotSum[
                            annual_data[0][
                                num_well][1].index == overlapyr].values[0]

                # Measurements
                plt.plot(bench.index, benchdata, "-o",
                         color="orange", label="Leveling")

        # Plottign GPS
        if GPSflag == 1:

            if not isinstance(GPS, type(None)):

                # If GPS and bench, check for first overlapping year
                if benchflag == 1:

                    # First overlapping year
                    overlapyr = next((i for i in GPS.index if i in
                                      bench[~bench.isnull().any(axis=1)].index),
                                     None)

                    # If overlapping year, but benchmark leveling does
                    # not have data
                    if np.logical_or(len(
                            bench[bench.index == overlapyr].values) == 0,
                            bench.empty):

                        # Adding starting simulated measurement to
                        # GPS data
                        # To start at "0" together
                        GPSdata = np.cumsum(GPS) + 100 * annual_data[0][
                            num_well][1].CumTotSum[
                                annual_data[
                                    num_well][1].index == GPS.index[0]].values

                    # If benchmark leveling and GPS data overlap and exists
                    else:

                        # Adding starting benchmark measurement to
                        # GPS data
                        # To start at "0" together

                        GPSdata = np.cumsum(GPS) + benchdata.loc[
                            bench.index == overlapyr].values[0]

                # Measurements
                plt.plot(GPSdata.index, GPSdata, "-o", color="red",
                         label="GPS")
                plt.errorbar(GPSdata.index, GPSdata, yerr=GPS_Var,
                             lw=.35, capsize=2, capthick=.35,
                             color="black")

        # Plottign INSAR
        if insarflag == 1:

            # If data for either insars
            if ~InSAR.isna().all():

                # If InSAR and bench, check for first overlapping year
                if benchflag == 1:

                    # First overlapping year
                    overlapyrERS = next((i for i in ERS.index if i in
                                        bench[~bench.isnull().any(axis=1)].index),
                                        None)
                    overlapyrSent = next((i for i in Sent.index if i in
                                          bench[~bench.isnull().any(axis=1)].index),
                                         None)

                    # If overlapping year, but benchmark leveling does
                    # not have data
                    if np.logical_or(len(
                            bench[bench.index == overlapyrERS].values) == 0,
                            bench.empty):

                        # Adding starting simulated measurement to
                        # GPS data
                        # To start at "0" together
                        ERSdata = ERS + 100 * annual_data[0][
                            num_well][1].CumTotSum[annual_data[
                                num_well][1].index == ERS.index[0]].values

                    # If benchmark leveling and InSARdata data overlap and exists
                    else:

                        # Adding starting benchmark measurement to
                        # InSARdata data
                        # To start at "0" together

                        ERSdata = ERS + benchdata.loc[
                            bench.index == overlapyrERS].values[0]

                    # If overlapping year, but benchmark leveling does
                    # not have data
                    if np.logical_or(len(
                            bench[bench.index == overlapyrSent].values) == 0,
                            bench.empty):

                        # Adding starting simulated measurement to
                        # GPS data
                        # To start at "0" together
                        Sentdata = Sent + 100 * annual_data[0][
                            num_well][1].CumTotSum[annual_data[0][
                                num_well][1].index == Sent.index[0]].values

                    # If benchmark leveling and InSARdata data overlap and exists
                    else:

                        # Adding starting benchmark measurement to
                        # InSARdata data
                        # To start at "0" together

                        Sentdata = Sent + benchdata.loc[
                            bench.index == overlapyrSent].values[0]

                # Measurements
                plt.plot(ERSdata.index, ERSdata, "-o", color="purple",
                         label="InSAR")
                # Measurements
                plt.plot(Sentdata.index, Sentdata, "-o", color="purple",
                         label="InSAR")

        plt.legend()
        plt.ylabel("Cumulative Subsidence (cm)")
        plt.xlabel("Years")
        plt.title(wellnest +
                  "\nObserved vs. Simulated Cumulative Subsidence")

        # If saving figure
        if save == 1:

            fig_name = wellnest + "_BenchvsImplicit_CumSubTotal.png"
            full_figpath = os.path.join(path, fig_name)
            plt.savefig(full_figpath, dpi=500, bbox_inches="tight")


def sub_sens_line(path, wellnestlist, all_results,
                  sub_total, subv_total,
                  annual_data, tmin=None, tmax=None, mode=None, num=None,
                  save=0):
    """Sensitivity analysis on subsidence based on either Sskv, Sske, K, thickness.

    path - path to save figures
    wellnestlist - list of wellnests that were simualted
    all_results - lists of lists: wellnestname, well name, head matrix with
    head series of each node
    sub_total - lists of lists: wellnestname, well_name, total cum sub
    subv_total - lists of lists: wellnestname, well_name, inelastic cum sub
    annual_data - lists of lists: wellnestname, well_name, total cum sub for
    all four clay at a wellnest location
    mode - which parameter is being adjusted for sensitivity (Sskv, Sske, K,
    thickness)
    num - number of parameter increases in sensitivity
    save - if 1, save; if 0, don't save

    ASSUMES FOUR WELLS IN WELLNEST
    """
    plt.figure(figsize=(6.75, 3.38), dpi=400)
    color1 = mcp.gen_color(cmap="rainbow", n=num)

    # Coeff for sensitivity percentage and plotting colors
    coeff = 50
    color_coeff = 1

    # Saving rates
    lastrates = []
    change2020_2060 = []

    # For each sensitivity
    for i in range(num):

        # For each wellnest in list
        # num_well is the index, wellnest = name
        # Figures for each well nest
        for num_well, wellnest in enumerate(wellnestlist):

            # If plotting for everything except elastic specific storage for sand
            if mode != "Sske_sand":

                plt.plot(annual_data[i][num_well][1].index,
                         annual_data[i][num_well][1].CumTotSum*-100,
                         label=str(coeff) + "%", linewidth=1,
                         color=color1[i])

            # if plotting for elastic specific storage for sand
            else:

                # For 50-140% sensitivity
                if i != (num - 1):

                    plt.plot(annual_data[i][num_well][1].index,
                             annual_data[i][num_well][1].CumTotSum*-100,
                             label=str(coeff) + "%", linewidth=1,
                             color=color1[i])

                # For the plot for 150%, sensitivity analysis for if elastic
                # specific storage for sand is equal to elastic specific storage
                # for clay
                else:

                    plt.plot(annual_data[i][num_well][1].index,
                             annual_data[i][num_well][1].CumTotSum*-100,
                             label="Clay Sske Values", linewidth=1,
                             color=color1[i])

            lastrate = (annual_data[i][num_well][1].CumTotSum[-1] -
                        annual_data[i][num_well][1].CumTotSum[-2])*-1000  # mm
            lastrates.append(lastrate)

            # Get cumulative sum of highest year below tmax
            max_year = annual_data[i][num_well][1].year[
                annual_data[i][num_well][1].year < int(tmax)].max()

            # Saving difference in cum sum between a time period
            # Cumulative sub in cm
            change2020_2060.append((annual_data[i][num_well][1].CumTotSum[
                annual_data[i][num_well][1].year == max_year].item() -
                annual_data[i][num_well][1].CumTotSum[
                annual_data[i][num_well][1].year == int(tmin)].item()) *
                -100)

        coeff += 10
        color_coeff -= .1

    print("Cum sub in time period relative to total sub (%): " +
          str((annual_data[0][num_well][1].CumTotSum[-1] -
               annual_data[-1][num_well][1].CumTotSum[-1]) * -100 /
              (annual_data[5][num_well][1].CumTotSum[-1])))

    if "Sske" in mode:

        # Y location of text
        yloc = 170

    else:

        # Y location of text
        yloc = 70

    plt.annotate("2020-2060\nSubsidence\nChange (cm)",
                 xy=(360, 190),
                 xycoords="axes points",
                 color="k",
                 fontsize=6,
                 weight="bold")

    # Different positionings of text
    for i in range(len(change2020_2060)):

        plt.annotate("{:.1f}".format(change2020_2060[i]),
                     xy=(380, yloc),
                     xycoords="axes points",
                     color=color1[i], fontsize=5, weight="bold")

        if "Sske" in mode:

            # Y location of text
            yloc -= 10

        else:

            # Y location of text
            yloc += 10

    # Plotting settings
    plt.legend()
    plt.ylabel("Cumulative Subsidence (cm)")
    plt.xlabel("Years")
    plt.title(wellnest)
    plt.gcf().set_size_inches(6.75, 3.38)

    # Title and figure name changes based on mode
    # Inelastic specific storage
    if mode == "Sskv":

        plt.title("Inelastic Specific Storage\nSensitivity Analysis")

        fig_name1 = wellnest + "_CumSubTotal_SENS_" + \
            mode + ".eps"

        fig_name2 = wellnest + "_CumSubTotal_SENS_" + \
            mode + ".png"

    # Elastic specific storage for clay
    elif mode == "Sske_clay":

        plt.title("Clay Elastic Specific Storage\nSensitivity Analysis")

        fig_name1 = wellnest + "_CumSubTotal_SENS_" + \
            mode + ".eps"

        fig_name2 = wellnest + "_CumSubTotal_SENS_" + \
            mode + ".png"

    # Elastic specific storage for clay
    elif mode == "Sske_sand":

        plt.title("Sand Elastic Specific Storage\nSensitivity Analysis")

        fig_name1 = wellnest + "_CumSubTotal_SENS_" + \
            mode + ".eps"

        fig_name2 = wellnest + "_CumSubTotal_SENS_" + \
            mode + ".png"

    # Vertical hydraulic conductivity
    elif mode == "K":

        plt.title("Vertical Hydraulic Conductivity\nSensitivity Analysis")

        fig_name1 = wellnest + "_CumSubTotal_SENS_" + \
            mode + ".eps"

        fig_name2 = wellnest + "_CumSubTotal_SENS_" + \
            mode + ".png"

    # Thickness
    elif mode == "thick":

        plt.title("Thickness Sensitivity Analysis")

        fig_name1 = wellnest + "_CumSubTotal_SENS_" + \
            mode + ".eps"

        fig_name2 = wellnest + "_CumSubTotal_SENS_" + \
            mode + ".png"

    # Just in case
    else:
        fig_name1 = "False"
        fig_name2 = "False"

    # If saving figure
    if save == 1:

        # full_figpath = os.path.join(path, fig_name1)
        # plt.savefig(full_figpath, dpi=400, format="eps")

        full_figpath = os.path.join(path, fig_name2)
        plt.savefig(full_figpath, dpi=400, format="png")


def sub_forecast(path, wellnestlist, all_ann_subs, save=0):
    """Forecasts of subsidence based on five pumping scenarios.

    1. 500,000 m3/day
    2. 250,000 m3/day
    3. Delayed 250,000 m3/day
    4. 1,000,000 m3/day
    5. No pumping

    path - path to save figures
    wellnestlist - list of wellnests that were simualted
    all_ann_subs - lists of lists: wellnestname, dataframe with annual subsidence
    rates
    save - if 1, save; if 0, don't save

    ASSUMES FOUR WELLS IN WELLNEST
    """
    # Plotting settings
    plt.rc("font", size=5)  # controls default text size
    plt.rc("axes", titlesize=6)  # fontsize of the title
    plt.rc("axes", labelsize=6)  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=6)  # fontsize of the x tick labels
    plt.rc("ytick", labelsize=6)  # fontsize of the y tick labels
    plt.rc("legend", fontsize=6)  # fontsize of the legend

    # Saving each scenario last annual rate
    ann_2060_500 = []
    ann_2060_250 = []
    ann_2060_d250 = []
    ann_2060_1000 = []
    ann_2060_0 = []

    # For each well nest
    for num_well, wellnest in enumerate(wellnestlist):

        # Figure plotting model results forecast for each scenario
        fig, ax = plt.subplots(figsize=(3.2, 2.2), dpi=400)

        # -1000 is to convert to mm and negative because subsidence is positive
        # while uplift is negative
        # 500,000 m3/day scenario
        plt.plot(all_ann_subs[0][num_well][1].index,
                 all_ann_subs[0][num_well][1].CumTotSum*-100,
                 label="500,000 m$^3$/day", linewidth=1.5,
                 color="hotpink")
        lastrate = (all_ann_subs[0][num_well][1].CumTotSum[-1] -
                    all_ann_subs[0][num_well][1].CumTotSum[-2])*-1000  # mm
        ann_2060_500.append(lastrate)
        ax.annotate("mm/yr",
                    xy=(180, 120),
                    xycoords="axes points",
                    color="k",
                    weight="bold")
        ax.annotate("{:.1f}".format(lastrate),
                    xy=(185, 105),
                    xycoords="axes points",
                    color="hotpink")

        # 250,000 m3/day scenario
        plt.plot(all_ann_subs[1][num_well][1].index,
                 all_ann_subs[1][num_well][1].CumTotSum*-100,
                 label="250,000 m$^3$/day", linewidth=1.5,
                 color="tab:orange")
        lastrate = (all_ann_subs[1][num_well][1].CumTotSum[-1] -
                    all_ann_subs[1][num_well][1].CumTotSum[-2])*-1000   # mm
        ann_2060_250.append(lastrate)
        ax.annotate("{:.1f}".format(lastrate),
                    xy=(185, 95),
                    xycoords="axes points",
                    color="tab:orange")

        # 500,000 -> 250,000 m3/day scenario
        plt.plot(all_ann_subs[3][num_well][1].index,
                 all_ann_subs[3][num_well][1].CumTotSum*-100,
                 label="Delayed\n250,000 m$^3$/day", linewidth=1.5,
                 color="tab:green")
        lastrate = (all_ann_subs[3][num_well][1].CumTotSum[-1] -
                    all_ann_subs[3][num_well][1].CumTotSum[-2])*-1000  # mm
        ann_2060_d250.append(lastrate)
        ax.annotate("{:.1f}".format(lastrate),
                    xy=(185, 100),
                    xycoords="axes points",
                    color="tab:green")

        # 1,000,000 m3/day scenario
        plt.plot(all_ann_subs[2][num_well][1].index,
                 all_ann_subs[2][num_well][1].CumTotSum*-100,
                 label="1,000,000 m$^3$/day", linewidth=1.5,
                 color="tab:red")
        lastrate = (all_ann_subs[2][num_well][1].CumTotSum[-1] -
                    all_ann_subs[2][num_well][1].CumTotSum[-2])*-1000  # mm
        ann_2060_1000.append(lastrate)
        ax.annotate("{:.1f}".format(lastrate),
                    xy=(185, 110),
                    xycoords="axes points",
                    color="tab:red")

        # No pumping scenario
        plt.plot(all_ann_subs[4][num_well][1].index,
                 all_ann_subs[4][num_well][1].CumTotSum*-100,
                 label="No Pumping", linewidth=1.5,
                 color="tab:purple")
        lastrate = (all_ann_subs[4][num_well][1].CumTotSum[-1] -
                    all_ann_subs[4][num_well][1].CumTotSum[-2])*-1000  # mm
        ann_2060_0.append(lastrate)
        ax.annotate("{:.1f}".format(lastrate),
                    xy=(185, 90),
                    xycoords="axes points",
                    color="tab:purple")

        # Observed pumping
        plt.plot(all_ann_subs[4][num_well][1].index[:44],
                 all_ann_subs[4][num_well][1].CumTotSum.iloc[:44]*-100,  # mm
                 color="black", linewidth=1.5,
                 label="Observed Pumping")
        plt.legend()

        # Plotting settings
        plt.ylabel("Cumulative Subsidence (cm)")
        plt.xlabel("Years")
        plt.title(wellnest)
        fig.set_size_inches(3.2, 2.2)
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        plt.grid(True, linestyle=(0, (1, 10)), which="minor")
        plt.grid(True, linestyle="dashed", which="major")

        # Annotating specific points
        index_1990 = all_ann_subs[4][num_well][1].year == 1990
        index_2000 = all_ann_subs[4][num_well][1].year == 2000
        cum_value_1990 = all_ann_subs[4][num_well][1].CumTotSum[index_1990]*-100
        cum_value_2000 = all_ann_subs[4][num_well][1].CumTotSum[index_2000]*-100
        ann_value_1990 = all_ann_subs[4][num_well][1].AnnRates[index_1990]*-1000
        ann_value_2000 = all_ann_subs[4][num_well][1].AnnRates[index_2000]*-1000
        plt.scatter(cum_value_1990.index, cum_value_1990[0], color="cyan")
        plt.scatter(cum_value_2000.index, cum_value_2000[0], color="cyan")
        # annotation
        plt.text(cum_value_1990.index, cum_value_1990[0] - 4, "1990: " +
                 f"{ann_value_2000[0]:.1f}" + " mm/yr")
        # annotation
        plt.text(cum_value_2000.index, cum_value_2000[0] - 4, "2000: " +
                 f"{ann_value_1990[0]:.1f}" + " mm/yr")
        plt.show()

        # Saving figure
        if save == 1:
            # fig_name = wellnest + "_CumSubForecast_ALLPUMP.eps"
            # full_figpath = os.path.join(path, fig_name)
            # plt.savefig(full_figpath, bbox_inches="tight",
            #             format="eps")

            fig_name = wellnest + "_CumSubForecast_ALLPUMP.png"
            full_figpath = os.path.join(path, fig_name)
            plt.savefig(full_figpath, bbox_inches="tight",
                        format="png")

    # Printing statistics
    print("\n500,000 scenario min, avg, max, med 2060 rate mm/yr: " +
          f"{min(ann_2060_500):.4f}" + ", " +
          f"{np.average(ann_2060_500):.4f}" + ", " +
          f"{max(ann_2060_500):.4f}" + ", " +
          f"{median(ann_2060_500):.4f}")

    print("\n250,000 scenario min, avg, max, med 2060 rate mm/yr: " +
          f"{min(ann_2060_250):.4f}" + ", " +
          f"{np.average(ann_2060_250):.4f}" + ", " +
          f"{max(ann_2060_250):.4f}" + ", " +
          f"{median(ann_2060_250):.4f}")

    print("\nDelayed250,000 scenario min, avg, max, med 2060 rate mm/yr: " +
          f"{min(ann_2060_d250):.4f}" + ", " +
          f"{np.average(ann_2060_d250):.4f}" + ", " +
          f"{max(ann_2060_d250):.4f}" + ", " +
          f"{median(ann_2060_d250):.4f}")

    print("\n1,000,000 scenario min, avg, max, med 2060 rate mm/yr: " +
          f"{min(ann_2060_1000):.4f}" + ", " +
          f"{np.average(ann_2060_1000):.4f}" + ", " +
          f"{max(ann_2060_1000):.4f}" + ", " +
          f"{median(ann_2060_1000):.4f}")

    print("\nNo pumping scenario min, avg, max, med 2060 rate mm/yr: " +
          f"{min(ann_2060_0):.4f}" + ", " +
          f"{np.average(ann_2060_0):.4f}" + ", " +
          f"{max(ann_2060_0):.4f}" + ", " +
          f"{median(ann_2060_0):.4f}")


def draw_basemap(map, xs, ys, cs=None, fig=None, ax=None,
                 datalim=None, mode=None, save=0, aq=None, labels=None,
                 perc=None, time_min=None, time_max=None, figpath=None, crit=None):
    """Drawing basemap for BKK.

    mode - Mode can be RMSE_full (Pastas), step_full (Pastas t90)
    sub_RMSE (subsidence RMSE), sub_forecast (subsidence forecast)
    Map contains the basemap
    xs - x locations (longitude) in list
    ys - y locations (latitude) in list
    cs - data to plot in list
    labels - list of names for labelling purposes
    datalim - datalimits
    fig - figure object
    ax - axis object
    aq - specific aquifer
    perc - percentage of rmse for example to put into title
    time_min - time minimum as string
    time_max - time maximum as string
    fig_path - figure to save path
    """
    # Plotting map
    # Land as green and ocean as blue
    # Drawing coastline
    map.drawcoastlines(zorder=2, linewidth=1)
    map.drawmapboundary(fill_color="#c1d4ec")
    # Continents
    map.fillcontinents(color="#4d9c83")

    # Adding Thailand province boundaries
    map.readshapefile(os.path.join(os.path.abspath("inputs/GIS"), "provinces"),
                      name="provinces", drawbounds=True, zorder=1, linewidth=.5)

    # Drawing rivers
    map.drawrivers(color="teal", linewidth=1)

    # draw parallels and meridians
    map.drawparallels(np.arange(12.5, 14, .5), labels=[1, 0, 0, 0],
                      fontsize=4)
    map.drawmeridians(np.arange(99.5, 101.5, .25), labels=[0, 0, 1, 0],
                      fontsize=4)

    # Drawing Pastas/subsidence datapoints
    x, y = map(xs, ys)

    # FOR WEDGES: 0 is at 3 pm, 90 is at noon
    # RMSE mode for all wells and all well nests
    # wells as wedges
    if mode == "RMSE_full":

        # Angle of wedges
        theta1 = 90
        theta2 = 180
        r = .018  # radius

        # Patches
        patches = []

        # All cs's
        cs_all = []

        # For each well
        for item in cs.items():

            # For each location
            for j in range(len(item[1].x)):

                # Creating wedge
                wedge = Wedge((item[1].x[j], item[1].y[j]),
                              r, theta1, theta2)
                patches.append(wedge)

            # Saving cs's
            cs_all.extend(item[1].cs)

            # Updating theta
            theta1 -= 90
            theta2 -= 90

        # Adding collection
        p = PatchCollection(patches, zorder=10,
                            edgecolor="k",
                            linewidth=.5)
        p.set_array(cs_all)  # Colorbar
        ax.add_collection(p)

        # New ax with dimensions of the colorbar
        cbar_ax = fig.add_axes([0.292, 0.05, 0.44, 0.03])

        # Colorbar
        cb = fig.colorbar(p, ax=ax, cax=cbar_ax, orientation="horizontal",
                          pad=0.05, ticks=mticker.MultipleLocator(.5))
        cb.ax.tick_params(labelsize=5)
        cb.set_label("RMSE (m)", fontsize=5)
        plt.set_cmap("coolwarm")
        cb.mappable.set_clim(vmin=datalim[0],
                             vmax=datalim[1])
        cb.solids.set_rasterized(False)

        # Legend objects
        class WedgeObject(object):
            pass

        class WedgeObjectHandler(object):

            def legend_artist(self, legend, orig_handle, fontsize, handlebox):
                x0, y0 = handlebox.xdescent, handlebox.ydescent
                width, height = handlebox.width, handlebox.height
                hw = x0+.45*width
                hh = y0+.5*height
                r2 = 5
                colors = ["#42BCFF", "#FF8300", "#D30808", "#009914"]
                lgd_patches = [Wedge((hw, hh), r2, 90, 180, color=colors[0],
                                     label="BK"),
                               Wedge((hw, hh), r2, 0, 90, color=colors[1],
                                     label="PD"),
                               Wedge((hw, hh), r2, 180, 270, color=colors[2],
                                     label="NB"),
                               Wedge((hw, hh), r2, 270, 360, color=colors[3],
                                     label="NL")]

                lgd_elements = PatchCollection(lgd_patches,
                                               match_original=True,
                                               edgecolor="k",
                                               linewidth=.5)

                handlebox.add_artist(lgd_elements)
                return lgd_elements

        ax.legend([WedgeObject()], ["BK PD\nNB NL"],
                  handler_map={WedgeObject: WedgeObjectHandler()},
                  fontsize=5, loc="upper left")

        # Title
        avgRMSE = pd.concat([cs["NB"].cs, cs["NL"].cs, cs["PD"].cs,
                             cs["BK"].cs], ignore_index=True)
        print(str("%.2f" % np.average(cs["NB"].cs)) + "m RMSE for NB")
        print(str("%.2f" % np.average(cs["NL"].cs)) + "m RMSE for NL")
        print(str("%.2f" % np.average(cs["PD"].cs)) + "m RMSE for PD")
        print(str("%.2f" % np.average(cs["BK"].cs)) + "m RMSE for BK")
        print("Average RMSE for all four aquifers: " +
              str("%.2f" % np.average(avgRMSE)) + "m")

    # t90 mode for all wells and all well nests
    # wells as wedges
    elif mode == "step_full":

        # Angle of wedges
        theta1 = 90
        theta2 = 180
        r = .018  # radius

        # Patches
        patches = []

        # All cs's
        cs_all = []

        # For each well
        for item in cs.items():

            # For each location
            for j in range(len(item[1].x)):

                # Creating wedge
                wedge = Wedge((item[1].x[j], item[1].y[j]),
                              r, theta1, theta2)
                patches.append(wedge)

            # Saving cs's
            cs_all.extend(item[1].cs)

            # Updating theta
            theta1 -= 90
            theta2 -= 90

        # Adding collection
        p = PatchCollection(patches, zorder=10,
                            edgecolor="k",
                            linewidth=.5)
        p.set_array(cs_all)  # Colorbar
        ax.add_collection(p)

        # Colorbar
        # New ax with dimensions of the colorbar
        cbar_ax = fig.add_axes([0.292, 0.05, 0.44, 0.03])

        # Colorbar
        cb = fig.colorbar(p, ax=ax, cax=cbar_ax, orientation="horizontal",
                          pad=0.05,
                          ticks=mticker.MultipleLocator(5))
        cb.ax.tick_params(labelsize=5)
        cb.set_label("Years", fontsize=5)
        cb.mappable.set_clim(vmin=datalim[0],
                             vmax=datalim[1])
        plt.set_cmap("plasma")
        cb.solids.set_rasterized(False)

        class Wedge_obj(object):
            pass

        class WedgeHandler(object):

            def legend_artist(self, legend, orig_handle, fontsize, handlebox):
                x0, y0 = handlebox.xdescent, handlebox.ydescent
                width, height = handlebox.width, handlebox.height
                hw = x0+.45*width
                hh = y0+.5*height
                r2 = 5
                colors = ["#42BCFF", "#FF8300", "#D30808", "#009914"]
                lgd_patches = [Wedge((hw, hh), r2, 90, 180, color=colors[0],
                                     label="BK"),
                               Wedge((hw, hh), r2, 0, 90, color=colors[1],
                                     label="PD"),
                               Wedge((hw, hh), r2, 180, 270, color=colors[2],
                                     label="NB"),
                               Wedge((hw, hh), r2, 270, 360, color=colors[3],
                                     label="NL")]

                lgd_elements = PatchCollection(lgd_patches,
                                               match_original=True,
                                               edgecolor="k",
                                               linewidth=.5)

                handlebox.add_artist(lgd_elements)
                return lgd_elements

        ax.legend([Wedge_obj()], ["BK PD\nNB NL"],
                  handler_map={Wedge_obj: WedgeHandler()},
                  fontsize=5, loc="upper left")

    elif mode == "Sub_RMSE_paper1":

        x = np.array(x)
        y = np.array(y)
        cs = np.array(cs)

        # Grouping by regular (0), noisy (1), bad (2)
        cluster = np.array([0, 2, 2, 1, 0,
                            1, 2, 0, 0, 0,
                            0, 1, 1, 0, 0,
                            1, 0, 0, 0, 0,
                            1, 0, 0])
        # EGU limited number of wells
        cluster = np.array([0, 2, 2, 1,
                            1, 2, 0, 0, 0,
                            0, 1, 1, 0, 0,
                            1, 0, 0])

        # For the colorbar of each three
        mini, maxi = np.min(cs), np.max(cs)
        norm = plt.Normalize(mini, maxi)

        map.scatter(x[cluster == 2], y[cluster == 2], zorder=3, marker="^",
                    c=cs[cluster == 2], label="Bad Fits",
                    cmap="RdYlBu_r", norm=norm, s=25,
                    edgecolor="k", linewidth=.75)
        map.scatter(x[cluster == 1], y[cluster == 1],
                    c=cs[cluster == 1], marker="s",
                    label="Noisy Observations",
                    cmap="RdYlBu_r", norm=norm, s=25,
                    edgecolor="k", linewidth=.75)
        map.scatter(x[cluster == 0], y[cluster == 0], norm=norm, s=30,
                    c=cs[cluster == 0], zorder=3,
                    marker="o", edgecolor="k",
                    cmap="RdYlBu_r", linewidth=.75)

        black_cirle = mlines.Line2D([], [], color="k", marker="o",
                                    linestyle="None",
                                    markersize=2, label="Good Fits")
        black_triangle = mlines.Line2D([], [], color="k", marker="^",
                                       linestyle="None",
                                       markersize=2, label="Bad Fits")
        black_square = mlines.Line2D([], [], color="k", marker="s",
                                     linestyle="None",
                                     markersize=2, label="Noisy Observations")
        plt.legend(handles=[black_cirle, black_triangle, black_square],
                   loc="lower right", prop={"size": 4})

        # New ax with dimensions of the colorbar
        cbar_ax = fig.add_axes([0.292, 0.05, 0.44, 0.03])
        cb = plt.colorbar(orientation="horizontal", cax=cbar_ax, pad=0.05,
                          ticks=mticker.MultipleLocator(.5))
        cb.mappable.set_clim(vmin=datalim[0],
                             vmax=datalim[1])
        cb.ax.tick_params(labelsize=5)
        cb.set_label("RMSE (cm/year)", fontsize=5)
        cb.solids.set_rasterized(False)
        plt.show()

    elif mode == "Sub_RMSE_paper2":

        x = np.array(x)
        y = np.array(y)
        cs = np.array(cs)

        # For the colorbar of each three
        mini, maxi = np.min(cs), np.max(cs)
        norm = plt.Normalize(mini, maxi)
        plt.annotate("(a)",
                     xy=(-.01, .99), xycoords="axes fraction",
                     fontsize=8, horizontalalignment="right",
                     weight="bold",
                     verticalalignment="bottom")
        map.scatter(x, y, norm=norm, s=30,
                    c=cs, zorder=3,
                    marker="o", edgecolor="k",
                    cmap="RdYlBu_r", linewidth=.75)

        # New ax with dimensions of the colorbar
        cbar_ax = fig.add_axes([0.292, 0.05, 0.44, 0.03])
        cb = plt.colorbar(orientation="horizontal", cax=cbar_ax, pad=0.05,
                          ticks=mticker.MultipleLocator(.5))
        cb.mappable.set_clim(vmin=datalim[0],
                             vmax=datalim[1])
        cb.ax.tick_params(labelsize=5)
        cb.set_label("RMSE (cm/year)", fontsize=5)
        cb.solids.set_rasterized(False)
        plt.show()

    # Inelastic or elastic values
    elif mode == "Param_val":

        x = np.array(x)
        y = np.array(y)
        cs = np.array(cs)

        map.scatter(x, y, s=30,
                    c=cs, zorder=3,
                    marker="o", edgecolor="k",
                    cmap="viridis_r", linewidth=.75)

        # New ax with dimensions of the colorbar
        cbar_ax = fig.add_axes([0.292, 0.05, 0.44, 0.03])
        cbformat = ticker.ScalarFormatter()
        cbformat.set_scientific('%.2e')
        cbformat.set_powerlimits((-4,12))
        cbformat.set_useMathText(True)
        cb = plt.colorbar(orientation="horizontal", cax=cbar_ax, pad=0.05,
                          ticks=mticker.MultipleLocator(.00001),
                          format=cbformat)
        # cb.mappable.set_clim(vmin=datalim[0],
        #                      vmax=datalim[1])
        cb.ax.tick_params(labelsize=5)
        cb.set_label("K Values [m/day]", fontsize=5)
        cb.solids.set_rasterized(False)
        plt.show()

    # Inelastic or elastic values
    elif mode == "pumping":

        x = np.array(x)
        y = np.array(y)
        cs = np.array(cs)
        plt.annotate("(b)",
                     xy=(-.01, .99), xycoords="axes fraction",
                     fontsize=8, horizontalalignment="right",
                     weight="bold",
                     verticalalignment="bottom")
        map.scatter(x, y, s=30,
                    c=cs, zorder=3,
                    marker="o", edgecolor="k",
                    cmap="Purples_r", linewidth=.75)

        # New ax with dimensions of the colorbar
        cbar_ax = fig.add_axes([0.292, 0.05, 0.44, 0.03])
        cbformat = ticker.ScalarFormatter()
        cbformat.set_scientific('%.2e')
        cbformat.set_powerlimits((-4,12))
        cbformat.set_useMathText(True)
        cb = plt.colorbar(orientation="horizontal", cax=cbar_ax, pad=0.05,
                          ticks=mticker.MultipleLocator(20),
                          format=cbformat)
        # cb.mappable.set_clim(vmin=datalim[0],
        #                      vmax=datalim[1])
        cb.ax.tick_params(labelsize=5)
        cb.set_label("Average pumping rate ($m^3$/day)", fontsize=5)
        cb.solids.set_rasterized(False)
        plt.show()

    # Plotting GW well locations
    elif mode == "GW_WellNests":

        map.scatter(x, y, s=15,
                    zorder=3,
                    marker="o", edgecolor="k",
                    linewidth=.75, color="mediumorchid")

        # Labelling well nests
        for i, txt in enumerate(labels):

            # Different labeling locations for these well nests
            if txt in ["LCBKK038", "LCBKK007", "LCBKK003", "LCBKK041",
                       "LCBKK005", "LCBKK021"]:
                ax.annotate(txt[2:], (x[i], y[i]), xycoords="data",
                            fontsize=3.5, color="w",
                            xytext=(-8, 1.5), textcoords="offset points",
                            weight="bold")

            # Other well nests standard labeling
            else:
                ax.annotate(txt[2:], (x[i], y[i]), xycoords="data",
                            fontsize=3.5, color="w",
                            xytext=(-8, 3.5), textcoords="offset points",
                            weight="bold")

            # Different color for well nest LCBKK013 in paper
            if txt == "LCBKK005":

                map.scatter(x[i], y[i], s=15,
                            zorder=3,
                            marker="o", edgecolor="k",
                            linewidth=.75, color="yellow")

        plt.show()

    # Forecasting subsidence for all wells
    elif mode == "Sub_Forecast_Map":

        # Angle of wedges
        theta1 = 90
        theta2 = 162
        r = .018  # radius

        # All cs's
        cs_all = []

        # Patches
        patches = []

        # For each pumping scenario
        for item in cs.items():

            # For each location/well nest
            for j in range(len(item[1].x)):

                # Creating wedge
                wedge = Wedge((item[1].x[j], item[1].y[j]),
                              r, theta1, theta2)
                patches.append(wedge)

            cs_all.extend(item[1].cs)

            # Updating theta
            theta1 -= 72
            theta2 -= 72

        # Adding collection
        p = PatchCollection(patches, zorder=10,
                            edgecolor="k",
                            linewidth=.5, match_original=True)
        p.set_array(cs_all)  # Colorbar
        ax.add_collection(p)

        # Colorbar
        # New ax with dimensions of the colorbar
        cbar_ax = fig.add_axes([0.292, 0.05, 0.44, 0.03])
        plt.set_cmap("BrBG")
        cb = fig.colorbar(p, ax=ax, cax=cbar_ax, orientation="horizontal",
                          pad=0.05,
                          ticks=mticker.MultipleLocator(5))
        cb.ax.tick_params(labelsize=5)
        cb.set_label("Cumulative Subsidence (cm)", fontsize=5)
        plt.set_cmap("BrBG")
        cb.mappable.set_clim(vmin=datalim[0],
                             vmax=datalim[1])
        cb.solids.set_rasterized(False)

        # Legend objects
        class WedgeObject(object):
            pass

        class WedgeObjectHandler(object):

            def legend_artist(self, legend, orig_handle, fontsize, handlebox):
                x0, y0 = handlebox.xdescent, handlebox.ydescent
                width, height = handlebox.width, handlebox.height
                hw = x0+.45*width
                hh = y0+.5*height
                r2 = 5
                colors = ["hotpink", "tab:orange", "tab:green", "tab:red",
                          "tab:purple"]
                lgd_patches = [Wedge((hw, hh), r2, 90, 162, color=colors[0],
                                     label="500,000"),
                               Wedge((hw, hh), r2, 18, 90, color=colors[1],
                                     label="250,000"),
                               Wedge((hw, hh), r2, 162, 234, color=colors[-1],
                                     label="No Pumping"),
                               Wedge((hw, hh), r2, 234, 306, color=colors[2],
                                     label="Delayed 250,000"),
                               Wedge((hw, hh), r2, 306, 18, color=colors[-2],
                                     label="1,000,000")]

                lgd_elements = PatchCollection(lgd_patches,
                                               match_original=True,
                                               edgecolor="k",
                                               linewidth=.5)

                handlebox.add_artist(lgd_elements)
                return lgd_elements

        ax.legend([WedgeObject()],
                  ["        500,000    250,000\nNo Pumping          1,000,000\n" +
                   "         Delayed 250,000"],
                  handler_map={WedgeObject: WedgeObjectHandler()},
                  fontsize=4, loc="lower left")

        plt.show()

        return

    # Saaving graphs
    if save == 1:
        fig_name1 = aq + "_" + mode + "_" + time_min + "_" + time_max + "_maps.eps"
        full_figpath = os.path.join(figpath, fig_name1)
        plt.savefig(full_figpath, dpi=400, bbox_inches="tight", format="eps")


def sub_forecast_map(path, wellnestlist, all_ann_subs,
                     tmin=None, tmax=None, save=0):
    """Plot subsidence forecast maps that are in main paper.

    all_results - lists of lists: wellnestname, well name, head matrix with
    head series of each node
    all_ann_subs - list of list of list of subsidence results for each well nest
    for each pumping scenario
    [0] 500,000
    [1] 250,000
    [2] 1,000,000
    [3] delayed 250,000
    [4] no pumping
    tmin - time min
    tmax - time max
    save - if 1, save; if 0, don't save

    Assumes four wells in well nest
    """
    # Importing spatial coordinates
    full_GWpath = os.path.join(os.path.abspath("inputs"),
                               "GroundwaterWellLocs.xls")
    gwwell_locs = pd.read_excel(full_GWpath)

    # Locations of wellnests; removing duplicates
    gwwell_locs = gwwell_locs.drop_duplicates("WellNest_Name", keep="first")

    # Preallocation
    # Empty dictionary
    d_dict = {}

    # For each pumping scenario
    # num_scenario is the pumping scenario index,scen_res = pumping scenario result
    for num_scenario, scen_res in enumerate(all_ann_subs):

        # Preallocation
        # Saving relevant xs, ys, and cumsum
        xs = []
        ys = []
        cs = []

        # For each well nest
        for i in range(len(scen_res)):

            # Get cumulative sum of highest year below tmax
            max_year = scen_res[i][1].year[
                scen_res[i][1].year < int(tmax)].max()

            # Saving difference in cum sum between a time period
            # Cumulative sub in cm
            cs.append((scen_res[i][1].CumTotSum[
                scen_res[i][1].year == max_year].item() -
                scen_res[i][1].CumTotSum[
                    scen_res[i][1].year == int(tmin)].item()) * -100)

            x_ = gwwell_locs.Long[
                gwwell_locs.WellNest_Name == scen_res[i][0]].item()
            y_ = gwwell_locs.Lat[
                gwwell_locs.WellNest_Name == scen_res[i][0]].item()
            xs.append(x_)
            ys.append(y_)

        # Creates a dictionary with location and cum sub
        d_dict[num_scenario] = pd.DataFrame({"x": xs, "y": ys, "cs": cs})

    # Initializing figure
    fig, ax = plt.subplots(figsize=(3.2, 2.2), dpi=400)
    datalim = [-5, 35]
    plt.set_cmap("BrBG")
    map = Basemap(llcrnrlon=100.3, llcrnrlat=13.4, urcrnrlon=100.8, urcrnrlat=14,
                  resolution="h", ellps="WGS84", lat_0=13.6, lon_0=100.4)
    draw_basemap(map, xs, ys, d_dict, fig=fig, ax=ax,
                 datalim=datalim, mode="Sub_Forecast_Map", save=0,
                 time_min=tmin, time_max=tmax, figpath=path)

    # If saving figure
    if save == 1:

        # fig_name1 = "Map_CumSub_" + tmin + "_" + tmax + "_ALLPump.eps"
        # full_figpath = os.path.join(path, fig_name1)
        # plt.savefig(full_figpath, dpi=400, bbox_inches="tight", format="eps")

        fig_name1 = "Map_CumSub_" + tmin + "_" + tmax + "_ALLPump.png"
        full_figpath = os.path.join(path, fig_name1)
        plt.savefig(full_figpath, dpi=400, bbox_inches="tight", format="png")


###############################################################################
# Plotting settings
###############################################################################

plt.rc("font", size=10)  # controls default text size
plt.rc("axes", titlesize=10)  # fontsize of the title
plt.rc("axes", labelsize=6)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=6)  # fontsize of the x tick labels
plt.rc("ytick", labelsize=6)  # fontsize of the y tick labels
plt.rc("legend", fontsize=6)  # fontsize of the legend


def Pastas_results(models, Wellnest_name, well_names,
                   time_mins, time_maxs, figpath, save):
    """Plot Pastas graphs that are in main paper.

    Wellnest_name - Name of well nest as a string
    well_names - list of name of wells in well nest
    models - listPastas models that was previously created
    time_mins - list of time minimum as string
    time_maxs - list of time maximum as string
    figpath - string path to save figure
    save - 1 to save figure, 0 to not save
    """
    # For each well in well names, determine order
    # Which well to start with because want BK, PD, NL, NB order
    order = [well_names.index(x) for y in ["BK", "PD", "NL", "NB"]
             for x in well_names if y in x]

    # Creating figure
    fig = plt.figure(figsize=(5, 5), dpi=400)
    # fig.suptitle(Wellnest_name)  # Figure title

    # Adding subfigures
    gs = fig.add_gridspec(ncols=1, nrows=len(order)+1,
                          width_ratios=[.5])

    gs.update(wspace=0.025, hspace=0.4)

    # Axes
    axes = []

    # For each well/model
    for idx, i in enumerate(order):

        # Current well
        model = models[i]
        time_min = min(time_mins)
        time_max = max(time_maxs)
        well_name = well_names[i]

        # Setting colors
        if "BK" in well_name:
            color = "blue"
        elif "PD" in well_name:
            color = "orange"
        elif "NL" in well_name:
            color = "green"
        elif "NB" in well_name:
            color = "red"

        # Obs, simulation, residuals, stress time series
        # Observation time series
        o = model.observations(tmin=time_min, tmax=time_max)
        o_nu = model.oseries.series.drop(o.index)
        o_nu = o_nu[time_min: time_max]

        # Simulation time series
        sim = model.simulate(tmin=time_min, tmax=time_max,
                             warmup=365*30, return_warmup=False)

        # Setting y limits
        ylims = [(min([sim.min(), o[time_min:time_max].min()]),
                  max([sim.max(), o[time_min:time_max].max()]))]
        # ylims = [-30, -10]

        # Frame
        if idx > 0:
            ax1 = fig.add_subplot(gs[idx], sharex=axes[0])
        else:
            ax1 = fig.add_subplot(gs[idx])
        axes.append(ax1)

        # First subplot
        # observation plot
        o.plot(ax=ax1, marker=".",
               color="k", label="Observed",
               x_compat=True, markersize=1, linewidth=.5)
        if not o_nu.empty:
            # plot parts of the oseries that are not used in grey
            o_nu.plot(ax=ax1, linestyle="-", marker=".", color="0.5",
                      label="", linewidth=.2,
                      x_compat=True, zorder=-1,
                      markersize=1)

        # add rsq to simulation
        rmse = model.stats.rmse(tmin=time_min, tmax=time_max)

        # Simulation plot
        sim.plot(ax=ax1, x_compat=True,
                 label=f"Simulated (RMSE = {rmse:.2} m)",
                 linewidth=1.5, color=color)

        # Plot 1 settings
        ax1.set_ylabel("Head\n(m)", labelpad=0)
        ax1.legend(loc=(0, 1), ncol=3, frameon=False, numpoints=3,
                   fontsize=8)
        ax1.set_ylim(ylims[0])
        ax1.tick_params(axis="both", which="major", pad=0)
        # ax1.set_title(well_name, fontsize=6)
        ax1.set(xlabel=None)
        ax1.xaxis.set_ticks_position("none")
        ax1.annotate(well_name,
                     xy=(-.1, 1.1), xycoords="axes fraction",
                     fontsize=8, horizontalalignment="center",
                     weight="bold",
                     bbox=dict(boxstyle="round", fc="0.8"),
                     verticalalignment="baseline")

        # Nonpaper version
        # Stress time series
        ax0 = ax1.twinx()
        stress = np.multiply(model.get_stress("well", tmin=time_min,
                                              tmax=time_max),
                             1 * 10**4)

        stress.plot(ax=ax0, linestyle="-", color=color, alpha=.25,
                    label="Observed Pumping")

        # Plot 1 settings
        ax0.set_ylabel("Rate\n(m$^3$/day)", labelpad=0)
        # ax0.set_xlabel("Year")
        # ax0.legend(loc=(.1, 1), ncol=3, frameon=False,
        #            fontsize=8)
        ax0.tick_params(axis="both", which="major", pad=0)
        ax0.set_ylim([2.5 * 10**5, 5 * 10**6])

        # xlim sets minorticks back after plots:
        ax0.minorticks_off()

        ax0.set_xlim(time_min, time_max)

        # Add a row for each stressmodel
        rmin = 0  # tmin of the response
        rmax = 0  # tmax of the response
        axb = None

        # plot the step response (90% of response)
        response = model._get_response(block_or_step="step",
                                       name="well", dt=0.5,
                                       add_0=True) * 50
        tmax = model.get_response_tmax("well", cutoff=0.90)
        response = response[:round(tmax*0.9)]

        rmax = max(rmax, response.index.max())

        # Inset graph settings
        left, bottom, width, height = [0.82, 0.3, 0.15, .4]
        axb = ax1.inset_axes([left, bottom, width, height])
        response.plot(ax=axb)
        title = "Step response"
        axb.tick_params(axis="both", which="major", pad=0)
        axb.set_xlabel("Days", fontsize=5, labelpad=-5)
        axb.set_ylabel("Head", fontsize=5, labelpad=-5)
        axb.xaxis.set_label_coords(.23, -.4)
        axb.yaxis.set_label_coords(-.38, .2)
        axb.set_title(title, fontsize=5, pad=0)
        axb.tick_params(labelsize=5)
        axb.set_xlim(rmin, rmax)

        # If last well
        # if idx == len(order) - 1: ## Correct code
        if idx == len(order):
            ax0 = fig.add_subplot(gs[idx+1], sharex=ax1)
            # Stress time series
            stress = np.multiply(model.get_stress("well", tmin=time_min,
                                                  tmax=time_max),
                                 1 * 10**4)

            stress.plot(ax=ax0, linestyle="-", color="tab:purple",
                        label="Observed Pumping")

            # Plot 1 settings
            ax0.set_ylabel("Rate\n(m$^3$/day)", labelpad=0)
            ax0.set_xlabel("Year")
            ax0.legend(loc=(.1, 1), ncol=3, frameon=False,
                       fontsize=8)
            ax0.tick_params(axis="both", which="major", pad=0)
            ax0.set_ylim([2.5 * 10**5, 3 * 10**6])

            # xlim sets minorticks back after plots:
            ax0.minorticks_off()

            ax0.set_xlim(time_min, time_max)

            # Grids
            for index, ax in enumerate(fig.axes):
                ax.grid(True)

                # Graph labels
                ax.text(0.01, .1, "(" + string.ascii_lowercase[index] + ")",
                        transform=ax.transAxes,
                        size=10, weight="bold")

            # No grids for inset graph
            axb.grid(False)

        if isinstance(fig, plt.Figure):
            fig.tight_layout(pad=0)  # Before making the table

        plt.subplots_adjust(right=0.95)

    fig.set_size_inches(5, 5)

    # To save figure
    if save == 1:
        # Fig name
        fig_name3 = Wellnest_name + "_GW_" + \
            time_min + "_" + time_max + "_PAPER.png"
        # Fig path
        full_figpath = os.path.join(figpath, fig_name3)
        # Save fig
        plt.savefig(full_figpath, dpi=400, bbox_inches="tight",
                    format="png")

        # # Fig name
        # fig_name3 = Wellnest_name + "_GW_" + \
        #     time_min + "_" + time_max + "_PAPER.eps"
        # # Fig path
        # full_figpath = os.path.join(figpath, fig_name3)
        # # Save fig
        # plt.savefig(full_figpath, dpi=400, bbox_inches="tight",
        #             format="eps")

def esmda_plot(mode, data):
    """Plot esmda evaluation plots

    mode - "Kalman"
    data - data to be plotted
    """

    if mode == "Kalman":

        plt.title("Kalman Gain")