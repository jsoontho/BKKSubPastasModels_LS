"""

Script for Thailand GW data pre-processing and post-plotting

Thai GW data from http://tgms.dgr.go.th/

Article Title: Hybrid data-driven, physics-based modeling of ground-
water and subsidence with application to Bangkok, Thailand

Jenny Soonthornrangsan 2023
TU Delft

"""
###############################################################################
# import statements
###############################################################################

# import statements
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
import datetime as dt

# Closing all figures
plt.close("all")


###############################################################################
###############################################################################
###############################################################################

# String formatting
# Function to preprocess groundwater data from http://tgms.dgr.go.th/#/home
# Assuming excel format
def GW_Data_Process(GW_Data, well_name=None):
    # # Inputs:
    # GW data assuming as dataframe and first five lines irrevelevant with - to be
    # replaced with nan
    # Thai years need to be converted to english years
    # Well name - if None, returning all wells with same dates
    # # Outputs:
    # Data frame with dates and groundwater data (HEAD REL to 0 NOT DEPTH TO WATER)
    # Subset of data for only the well from well_name - If well name not None

    # Ignoring first five lines
    # Replacing no date (-) with nans
    data = GW_Data.iloc[0:len(GW_Data)-5, :]
    data = data.replace("-", np.nan)

    # Creating data frame
    well_data = data.iloc[:, 1:]
    df_data = well_data.rename(columns={"วันที่": "Date"})
    df_data.Date = df_data["Date"].astype(str)

    # Reformating date from thai years to english years
    len_head = len(df_data.iloc[:, 1])
    date_list = []
    df_data["Year"] = np.nan
    df_data["Month"] = np.nan
    df_data["Day"] = np.nan

    # For each data point
    for i in range(len_head):

        # Current date
        date = df_data.Date.loc[i]

        # If leap day
        if date[0:2] == "29" and date[3:5] == "02":

            # Thai years - 543 = english years
            df_data.Year.loc[i] = int(date[6:10]) - 543
            df_data.Month.loc[i] = int(date[3:5])
            df_data.Day.loc[i] = int(date[0:2])

            # Converting to date time format
            date_list.append(dt.datetime.strptime(str(df_data.Day.loc[i]).replace(
                ".0", "") + "/" +
                str(df_data.Month.loc[i]).replace(".0", "") +
                "/" + str(df_data.Year.loc[i]).replace(".0", ""),
                "%d/%m/%Y").date())

        # If not a leap date
        else:

            # Converting to date time format
            date_list.append(dt.datetime.strptime(date, "%d/%m/%Y %M:%S").date())
            df_data.Year.loc[i] = (date_list[i].year) - 543
            df_data.Month.loc[i] = (date_list[i].month)
            df_data.Day.loc[i] = (date_list[i].day)

    # Saving new english date
    df_data["EngDate"] = pd.to_datetime(df_data[["Year", "Month", "Day"]])

    # If individual well name given
    if well_name is not None:

        # Subset of dataframe to get ready for Pastas
        # DTW data converted to approximate head using negative sign
        head_subsetdata = {"Date": df_data.EngDate,
                           "Head": - df_data[well_name].astype("float")}
        Head_Data = pd.DataFrame(head_subsetdata, columns=["Date", "Head"])
        Head_Data.index = pd.to_datetime(Head_Data.Date)

    # If individual head date not given
    else:
        Head_Data = None

    # Cleaning all data up
    # Returning head data not depth to water!
    if len(well_data.columns.values[1:]) > 1:
        all_data = pd.concat([df_data["EngDate"],
                              -well_data.iloc[:, 1:].astype("float")],
                             axis=1,
                             keys=["Date"] + well_data.columns.values[1:])
    else:
        all_data = pd.concat([df_data["EngDate"],
                              -well_data.iloc[:, 1:].astype("float")],
                             axis=1)

    if (np.size(-well_data.iloc[:, 1:], axis=1)) > 1:
        all_data = all_data.droplevel(level=0, axis=1)

    return all_data, Head_Data


###############################################################################
###############################################################################
###############################################################################

# Checks for outliers
def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    """
    if np.size(points) == 1:
        points = points[:, None]
    median = np.median(points)
    diff = (points - median)**2
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh
