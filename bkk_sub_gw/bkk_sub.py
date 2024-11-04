# ##############################################################################
"""Calculate subsidence in BKK at wellnests with 8 aquifers but simulates top four.

BK, PD, NL, NB
all are confined and overlain by clay layer
Implicit method according to USGS SUB package Hoffman report pg. 14

Article Title: Hybrid data-driven, physics-based modeling of ground-
water and subsidence with application to Bangkok, Thailand

Jenny Soonthornrangsan 2023
TU Delft

"""
# ##############################################################################

###############################################################################
# import statements
###############################################################################

import re
import os
import pandas as pd
import numpy as np
import scipy.linalg as lin
import sys
import functools
import pastas as ps
import datetime as dt
import random
import lmfit
import importlib

# Importing script for pre-processing Thai GW data
import main_functions as mfs
import warnings

# Ignoring Pastas warnings
pd.options.mode.chained_assignment = None  # default='warn'
warnings.simplefilter(action="ignore", category=FutureWarning)


# %%###########################################################################
# Preprocessing GW well nest data
###############################################################################

def bkk_wellnest_preproc(wellnestname, tmin, tmax, proxyflag):
    """Take well nest name, load data, and clean it up.

    Returns data within tmin and tmax
    By cleaning, gets rid of thai characters, interpolates dates and head
    Keeps matching dates and data only between wells in the nest
    wellnestname-- (str) name of well nest
    tmin, tmax -- (str) minimum and maximum year, if min year = 1900 but
    data starts at 1960, will return data from 1960 onwards
    proxyflag - 1 if using available heads as proxy for missing

    Returns: well_data_dates - dataframe with only matching dates and data
    between wells
    well_data - dataframe with matching dates, some wells have missing data
    """
    # Reading in GW data
    # Path to GW data
    try:
        full_path = os.path.join(os.path.abspath("inputs"), wellnestname + ".xlsx")
        data = pd.read_excel(full_path, skiprows=3)

    # If well nest does not exist
    except:
        raise ValueError("\nWell nest or file for well nest does not exist.")

    # List of wells in well nest
    welllist = data.columns[-(len(data.columns) - 2):]

    # Reorder well list to shallow to deep aquifers
    # BK, PD, NL, NB
    welllist = [x for y in ["BK", "PD", "NL", "NB"] for x in welllist
                if y in x]

    # Returns all data, and specific well data if specified.
    # GW Head not DTW
    all_head_data, gw_well_head = mfs.GW_Data_Process(data)
    all_head_data.index = all_head_data.EngDate

    # Stores interpolated data
    interp_welldata = []

    # For each well in welllist
    for i in welllist:

        # Gets rid of NA, then resamples daily, then does cubic interpolation
        # Interpolating 'inside'
        interp_welldata.append(
            all_head_data.loc[:, i].dropna().
            resample("D").interpolate("linear"))

    lenlist = len(welllist)

    # If using available heads as proxy for missing heads
    if proxyflag == 1:

        well_data = []

        # For those missing wells, checks which well is missing
        if lenlist < 4:

            # If only three files
            # If first file in sorted list is PD, missing BK
            if lenlist == 3:

                if not any("BK" in substring for substring in welllist):

                    # print("\nPROXY\n")
                    # Using PD as BK Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy BK")
                    well_data.append(temp)

                    # Adding rest of the wells
                    [well_data.append(i) for i in interp_welldata]

            # If only two files
            # If first file in sorted list is PD and next is NB
            # missing BK and NL
            elif lenlist == 2:

                if np.logical_and(not any("BK" in substring
                                          for substring in welllist),
                                  not any("NL" in substring
                                          for substring in welllist)):

                    # Using PD as BK Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy BK")
                    well_data.append(temp)

                    # PD
                    temp = interp_welldata[0]
                    well_data.append(temp)

                    # Using NB as NL Proxy
                    temp = interp_welldata[1]
                    temp = temp.rename("Proxy NL")
                    well_data.append(temp)

                    # NB
                    temp = interp_welldata[1]
                    well_data.append(temp)

            # If only two files
            # If first file in sorted list is PD, next and NL
            # missing BK and NB

                elif np.logical_and(
                        not any("BK" in substring for substring in welllist),
                        not any("NB" in substring for substring in welllist)):

                    # Using PD as BK Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy BK")
                    well_data.append(temp)

                    # PD
                    temp = interp_welldata[0]
                    well_data.append(temp)

                    # NL
                    temp = interp_welldata[1]
                    well_data.append(temp)

                    # Using NL as NB proxy
                    temp = interp_welldata[1]
                    temp = temp.rename("Proxy NB")
                    well_data.append(temp)

            # If only two files
            # If first file in sorted list is NL, next and NB
            # missing BK and PD

                elif np.logical_and(
                        not any("BK" in substring for substring in welllist),
                        not any("PD" in substring for substring in welllist)):

                    # Using NL as BK Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy BK")
                    well_data.append(temp)

                    # Using NL as PD Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy PD")
                    well_data.append(temp)

                    # NL
                    temp = interp_welldata[0]
                    well_data.append(temp)

                    # Using NL as NB proxy
                    temp = interp_welldata[1]
                    well_data.append(temp)

            # If only 1 file
            # missing others
            if lenlist == 1:

                # If only has BK
                if any("BK" in substring for substring in welllist):

                    # BK
                    temp = interp_welldata[0]
                    well_data.append(temp)

                    # Using BK as PD Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy PD")
                    well_data.append(temp)

                    # Using BK as NL Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy NL")
                    well_data.append(temp)

                    # Using BK as NB Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy NB")
                    well_data.append(temp)

                # If only has PD
                elif any("PD" in substring for substring in welllist):

                    # Using PD as BK Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy BK")
                    well_data.append(temp)

                    # PD
                    temp = interp_welldata[0]
                    well_data.append(temp)

                    # Using PD as NL Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy NL")
                    well_data.append(temp)

                    # Using PD as NB Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy NB")
                    well_data.append(temp)

                # If only has NL
                elif any("NL" in substring for substring in welllist):

                    # Using NL as BK Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy BK")
                    well_data.append(temp)

                    # Using NL as PD Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy PD")
                    well_data.append(temp)

                    # Using NL
                    temp = interp_welldata[0]
                    well_data.append(temp)

                    # Using NL as NB Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy NB")
                    well_data.append(temp)

                # If only has NB
                elif any("NB" in substring for substring in welllist):

                    # Using NB as BK Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy BK")
                    well_data.append(temp)

                    # Using NB as PD Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy PD")
                    well_data.append(temp)

                    # Using NB as NL Proxy
                    temp = interp_welldata[0]
                    temp = temp.rename("Proxy NL")
                    well_data.append(temp)

                    # NB
                    temp = interp_welldata[0]
                    well_data.append(temp)

        # No missing wells
        else:

            well_data = interp_welldata

    # If using only available heads
    else:

        well_data = interp_welldata

        # Needs all four wells if proxyflag is not on
        if lenlist < 4:

            sys.exit("Needs all four wells if proxyflag is not on")

    # Well data with matching dates only
    well_data_dates = functools.reduce(lambda left, right:
                                       pd.merge(left, right, on=["EngDate"],
                                                how="inner"), well_data)
    well_data_dates = well_data_dates[(well_data_dates.index.year >= int(tmin)) &
                                      well_data_dates.index.year <= int(tmax)]

    # All well data
    all_well_data = functools.reduce(lambda left, right:
                                     pd.merge(left, right, on=["EngDate"],
                                              how="outer", sort=True),
                                     well_data)

    return well_data_dates, all_well_data


# %%###########################################################################
# Groundwater model for clay: Using Matrixes to solve for head
###############################################################################

class SolveFDM:
    """Solve for head for time step n using the finite difference implicit method.

    Solving for x (head) in Ax = b (See USGS SUB report by Hoffman pg 14)
    """

    def __init__(self, Nz, n, dz, Kv, Sskv, Sske,
                 dt, precon, CC, toplay):

        self.Nz = Nz  # Number of nodes
        self.n = n  # Current time step
        self.dz = dz  # Cell diff
        self.Kv = Kv  # Vertical hydraulic conductivity
        self.Sskv = Sskv  # Specific storage (inelastic)
        self.Sske = Sske  # Specific storage (elastic)
        self.dt = dt  # Time diff
        self.precon = precon  # Preconsolidated head
        self.CC = CC  # Convergence criteria
        self.toplay = toplay  # If top layer or not

    # Checking if elastic or inelastic for each clay node
    def ElasticInelastic(self, h_matr):
        """Check if elastic or inelastic for each clay node.

        Input: h_matrx - head matrix
        Output:
        Ss - array the size of precon with either Sskv or Sske for each node
        precons - array with the preconsolidated head for each node
        """
        # Creating Ss array
        Ss = np.zeros(np.shape(self.precon))

        # Creating array for updated preconsolidated head
        precons = self.precon.copy()

        # For each inner node, ignoring aquifer nodes
        for i in range(1, self.Nz+1):

            # If it is the first time step
            if self.n == 1:

                # If current head is less than or equal (INELASTIC)
                if h_matr[i, self.n-1] <= self.precon[i]:

                    # Saves inelastic storage term for this node
                    Ss[i] = self.Sskv

                    # Sets new preconsolidation head
                    precons[i] = h_matr[i, self.n-1]

                # ELASTIC
                else:

                    # Saves elastic storage term for this node
                    Ss[i] = self.Sske

            # All other time steps
            else:

                # Difference in head for the previouse two time steps
                dh = h_matr[i, self.n-1] - h_matr[i, self.n-2]

                # If current head is less than or equal (INELASTIC) and if
                # slope is negative
                if np.logical_and(h_matr[i, self.n-1] <= self.precon[i], dh < 0):

                    # Inelastic storage term saved
                    Ss[i] = self.Sskv

                    # Sets new preconsolidation head
                    precons[i] = h_matr[i, self.n-1]

                # ELASTIC
                else:
                    Ss[i] = self.Sske

        # Returning
        return Ss, precons

    # Building A matrix
    def buildCoeffMatrix(self, h_matr):
        """Build A matrix.

        Input: h_matrix - head matrix
        Output:
        A - A matrix
        Ss - array with either elastic/inelastic storage term for each node
        precon - updated preconsolidated head for each node
        """
        Ss, precon = self.ElasticInelastic(h_matr)

        # Preallocation
        Adiag_val = np.ones(self.Nz)

        # For each main diagonal except the first and last inner node
        # IMPORTANT: Ss array includes cell for aquifer on top and bottom
        # Diag only for inner nodes -> thus, inner node 2 has index of 1 in the
        # diagonal while index of 2 in the Ss (i+1)
        for i in range(1, self.Nz-1):
            Adiag_val[i] = (-2 * self.Kv / self.dz) - (self.dz / self.dt * Ss[i+1])

        # First value and last value of the main diagonal
        # Inner nodes that border the aquifer
        # First inner node near top aquifer
        # If not top clay layer, the top is an aquifer
        if self.toplay is False:
            Adiag_val[0] = (-3 * self.Kv / self.dz) - (self.dz / self.dt * Ss[1])

        # If top clay layer, the top is a noflow boundary
        else:
            Adiag_val[0] = -(self.Kv / self.dz) - (self.dz / (self.dt) * Ss[1])

        # Last inner node near bottom aquifer
        Adiag_val[-1] = (-3 * self.Kv / self.dz) - (self.dz / self.dt * Ss[-2])

        # Creating A matrix
        Aupper = np.diag(np.ones(self.Nz-1) * self.Kv / self.dz, 1)  # Upper diag
        Alower = np.diag(np.ones(self.Nz-1) * self.Kv / self.dz, -1)  # Lower diag
        Adiag = np.diag(Adiag_val)  # Main diagonal
        A = Alower + Aupper + Adiag

        # Returning
        return A, Ss, precon

    # Building b matrix
    def buildRHSVector(self, h_matr, Ss, precon):
        """Right hand side vector (b).

        Input:
        h_matr - head matrix
        Ss - array of either elastic or inelastic storage for each node
        precon - array of updated preconsolidated head for each node
        """
        # Preallocation
        b = np.ones(self.Nz)

        # For each inner node that is not the first/last
        # IMPORTANT: Ss/h_matr/precon array includes cell for aquifer on top
        # and bottom; b only for inner nodes -> thus, inner node 2 has index
        # of 1 in b while index of 2 in the other arrays (i+1)
        for i in range(1, self.Nz-1):
            b[i] = (self.dz/self.dt) * (-Ss[i+1] * self.precon[i+1] +
                                        self.Sske * (self.precon[i+1] -
                                                     h_matr[i+1, self.n-1]))

        # If not top clay layer, the top is an aquifer
        if self.toplay is False:
            # First inner node near top aquifer
            b[0] = (self.dz/self.dt) * (-Ss[1] * self.precon[1] +
                                        self.Sske * (self.precon[1] -
                                                     h_matr[1, self.n-1])) - \
                2 * self.Kv / self.dz * h_matr[0, self.n]

        # If top clay layer, the top is a noflow boundary
        else:
            # First inner node near top aquifer
            b[0] = (self.dz/self.dt) * (-Ss[1] * self.precon[1] +
                                        self.Sske * (self.precon[1] -
                                                     h_matr[1, self.n-1]))

        # Last inner node near bottom aquifer
        b[-1] = (self.dz/self.dt) * (-Ss[-2] * self.precon[-2] +
                                     self.Sske * (self.precon[-2] -
                                                  h_matr[-2, self.n-1])) - \
            2 * self.Kv / self.dz * h_matr[-1, self.n]

        # Returning
        return b

    def solveLinearSystem(self, A, b):
        """Solve linear system of matrices."""
        h = lin.solve(A, b)
        return h

    # Iterates through until all cells meet the convergence criteria for a time
    # step n
    def iterate(self, h_matr, precons_head):
        """Iterate until all cells meet the convergence criteria for a timestep n.

        Input:
        h_matr - head matrix
        precons_head - current preconsolidated head for each node at the
        start of the time step
        Output:
        h_matr - head matrix updated with new heads in n time step after
        iterating
        precons_head - updated preconsolidated head for each node at the end
        of the time step
        """
        # Preallocation for the head diff in each cell
        Cell_change = np.ones(self.Nz)

        # Sets the starting new heads
        h_new = h_matr[1:self.Nz+1, self.n].copy()

        # While head diff of cells is greater than convergence criteria
        # Iterates
        while np.sum(Cell_change > self.CC) > 0:

            # Remembers old head
            old_head = h_new.copy()

            # Creates new class with updated precons_head
            fdm = SolveFDM(self.Nz, self.n, self.dz, self.Kv, self.Sskv, self.Sske,
                           self.dt, precons_head, self.CC,
                           self.toplay)

            # Builds A matrix and updates Ss and preconsolidated head
            A, Ss, precons_head = fdm.buildCoeffMatrix(h_matr)

            # Builds right hand side array
            b = fdm.buildRHSVector(h_matr, Ss, precons_head)

            # Solves for head using A and RHS matrix b
            h_new = fdm.solveLinearSystem(A, b)

            # Checks for the difference between iterations
            Cell_change = np.abs(np.subtract(h_new, old_head))

        # Adds the aquifer heads at the top and bottom
        h_new = np.insert(h_new, 0, h_matr[0, self.n])
        h_new = np.append(h_new, h_matr[-1, self.n])

        # Saves new head array in the current time step in the head matrix
        h_matr[:, self.n] = h_new

        # Returning updated head matrix and preconsolidated head after iterate
        return h_matr, precons_head


# %%###########################################################################
# Subsidence model: calculates compaction for each layer
##############################################################################


# Calculates compaction based on groundwater levels in aquifer above and
# below a clay layer
def calc_deformation(timet, headt, headb, Kv, Sskv, Sske, Sske_sandt,
                     Sske_sandb, claythick, nclay, sandthickt, sandthickb,
                     Nt, CC, Nz=None, ic=None):
    """Calculate deformation for a single clay layer of user defined thickness.

    Use whatever units for time and length as desired, but they need to stay
    consistent
    Inputs:
    timet - a vector of same lenght as head with the times that head
    measurements are taken. Numeric (years or days, typically)
    headt - a vector of same length as time. Head of top aquifer
    headb - a vector of same length as time. Head of bottom aquifer
    Kv - vertical hydraulic conductivity
    Sske - Skeletal specific storage (elastic)
    Sskv - skeletalt specific storage (inelastic)
    Sske_sandt - Skeletal specific storage (elastic) of aq on top
    Sske_sandb - Skeletal specific storage (elastic) of aq on bottom
    claythick - thickness of single clay layer modeled
    nclay - number of clay layers
    sandthickt -  thickness of sand in top aquifer
    sandthickb -  thickness of sand in bottom aquifer
    Nz - number of layers in z direction, within the clay layer modeled.
    Nt - number of time steps
    CC - convergence criteria
    ic - if providing initial condition of clay, ndarray given

    Outputs:
    t - interpolated time
    deformation - cumulative sum of deformation of total clay layer (m)
    boundaryt - interpolated head at the top boundary
    boundaryb - interpolated head at the bottom boundary
    deformation_v - cumulative sum of inelastic deformation of total clay (m)
    h - aquifer heads row 0 and -1, rest are clay nodes head
    """
    # Storage coefficients of aquifers top and bottom
    Ske_sandt = Sske_sandt * sandthickt
    Ske_sandb = Sske_sandb * sandthickb

    # Interpolated time
    # The first time step (0) doesn't count so needs Nt + 1
    t = np.linspace(timet[0], timet[-1], int(Nt+1))

    # Interpolated head at the boundaries (aquifers)
    # Time and head have the same time steps. Interpolating head for t
    # Has the same number of steps as t
    # If not the top clay layer with no top aquifer
    if isinstance(headt, pd.Series):

        # Interpolating top and bottom aquifer heads
        boundaryt = np.interp(t, timet, headt)  # Linear
        boundaryb = np.interp(t, timet, headb)  # Linear

        # Initial conditions of head grid
        if isinstance(ic, np.ndarray):

            # Setting the head matrix to the initial condition
            h = np.tile(ic, (Nt+1, 1))
            h = h.transpose()

        h[0, :] = boundaryt
        h[-1, :] = boundaryb

        # It is not the top clay layer
        toplay = False

    # If top clay layer
    else:
        boundaryb = np.interp(t, timet, headb)  # Linear
        boundaryt = np.zeros(np.shape(boundaryb))

        # Initial conditions of head grid
        if isinstance(ic, np.ndarray):

            # Setting the head matrix to the initial condition
            h = np.tile(ic, (Nt+1, 1))
            h = h.transpose()

        else:

            # Initial conditions of head grid
            h = (headb[0])*np.ones((Nz+2, Nt+1))

        h[-1, :] = boundaryb

        # It is the top clay layer
        toplay = True

    # Initial precons head made of initial head in each node
    # Preconsolidated head set to head from first time step
    precons_head = h[:, 0].copy()

    # Preallocation for total/inelastic deformation
    deformation = np.zeros(np.shape(h))
    deformation_v = np.zeros(np.shape(h))

    # Length of z
    dz = claythick / (Nz)

    # For each time step
    # Starting at 1, because 0 doesn't count as a time step
    for n in range(1, Nt+1):

        # Difference in time
        dt2 = t[n] - t[n-1]

        # Finite difference implicit method solving for head at current time
        # step. Uses matrix. Iterative because Sskv and Sske can change
        fdm = SolveFDM(Nz, n, dz, Kv, Sskv, Sske, dt2,
                       precons_head, CC, toplay=toplay)
        h, precons_head = fdm.iterate(h, precons_head)

        # New head that is calculated is already saved to head matrix
        h_new = h[:, n]

        # Compute compaction
        for i in range(1, Nz+1):

            # Diff between new and old
            dh = (h_new[i] - h[i, n-1])

            # If head drops below preconsolidation head and slope is neg
            # INELASTIC
            if np.logical_and(h_new[i] <= precons_head[i], dh < 0):

                # Calculating deformation
                defm = dh * Sskv * dz

                # Adds new deformation to the min from before for this
                # node from time steps before
                deformation_v[i, n] = defm + np.min(deformation_v[i, 0:(n)])

            # ELASTIC
            else:

                # Calculating deformation
                defm = dh * Sske * dz

                # Next time step deformation equals this current time step
                deformation_v[i, n] = deformation_v[i, n-1]

            # Total deformation updated
            deformation[i, n] = defm + deformation[i, n-1]

    # Deformation multipled by number of clay layers
    deformation = np.sum(deformation, axis=0) * nclay
    deformation_v = np.sum(deformation_v, axis=0) * nclay

    # If not the top clay layer
    if isinstance(headt, pd.Series):

        # Interpolated head minus initial
        boundary0_t = boundaryt - boundaryt[0]

    # If the top clay layer
    else:
        boundary0_t = 0

        h[0, :] = h[1, :]

        # Top row equals the second row
        boundaryt = h[1, :]

    # Interpolated head minus initial
    boundary0_b = boundaryb - boundaryb[0]

    # Sand deformation (top and bottom aquifer)
    # All elastic
    sanddef = boundary0_t * Ske_sandt + boundary0_b * Ske_sandb

    # If adding deformation from sand
    deformation = deformation + sanddef

    # Returning
    return (t, deformation, boundaryt, boundaryb, deformation_v, h)


# %%###########################################################################
# Pastas model forecasting: simulates groundwater using different
# pumping scenarios
##############################################################################

# Future simulating Pastas with different pumping scenarios
def pastas_setparam(model, pump_series=None, pump_path=None, pump_sheet=None,
                    initoptiparam=None, well_name=None):
    """Pastas model that is already created.

    pump_series - if pump series given already
    pump_path - path to pumping excel sheet
    pump_sheet - sheet of specific pumping scenario
    initoptiparam - optimal parameters provided or not
    well_name - name of well
    """

    model.del_noisemodel()
    stdparam = model.parameters["stderr"]

    # Saves optimal parameters and SD
    if initoptiparam is None:

        optiparam = model.parameters["optimal"]

        # If no pump series provided
        if pump_series is None:

            # If pump path provided
            if pump_path is not None:
                model.del_stressmodel("well")  # Deletes previous pumping

                # Adds new pumping
                EstTotPump = pd.read_excel(pump_path,
                                           sheet_name=pump_sheet, index_col=0,
                                           parse_dates=["Date"])

                EstTotPump_ = ps.StressModel(EstTotPump.Pump,
                                             rfunc=ps.Gamma(), name="well",
                                             settings="well", up=False)
                model.add_stressmodel(EstTotPump_, replace=True)

        # If pumping series provided
        else:

            model.del_stressmodel("well")  # Deletes previous pumping

            EstTotPump_ = ps.StressModel(pump_series, rfunc=ps.Gamma(), name="well",
                                         settings="well", up=False)
            model.add_stressmodel(EstTotPump_)

        for name_i, name in enumerate(model.parameters.index):
            model.set_parameter(name=name, optimal=optiparam[name_i],
                                initial=optiparam[name_i])
        model.parameters["stderr"] = stdparam

    else:

        # If no pump series provided
        if pump_series is None:

            # If pump path provided
            if pump_path is not None:

                model.del_stressmodel("well")  # Deletes previous pumping

                # Adds new pumping
                EstTotPump = pd.read_excel(pump_path,
                                           sheet_name=pump_sheet, index_col=0,
                                           parse_dates=["Date"])
                EstTotPump_ = ps.StressModel(EstTotPump.Pump,
                                             rfunc=ps.Gamma(), name="well",
                                             settings="well", up=False)
                model.add_stressmodel(EstTotPump_)

        # If pumping series provided
        else:

            model.del_stressmodel("well")  # Deletes previous pumping

            EstTotPump_ = ps.StressModel(pump_series, rfunc=ps.Gamma(), name="well",
                                         settings="well", up=False)

            model.add_stressmodel(EstTotPump_)

        for name_i, name in enumerate(model.parameters.index):
            model.set_parameter(name=name, optimal=initoptiparam[name_i],
                                initial=initoptiparam[name_i])
        model.parameters["stderr"] = stdparam

    # Returns model
    return model


# %%###########################################################################
# Loads Pastas models
##############################################################################

def load_Pastas_models(Pastasfiles, model_path, SS_data):
    """ Loads Pastas models and saves models for use later (so that models can be
    loaded only once)

    Pastasfiles- list of Pastas file names to load
    model_path - path of where files are located
    SS_data - steady state heads to set the pmin and pmax of constant d

    Return:
        models - list of Pastas models
        well_names - list of well names (str)
        modelinitparam - initial parameters
    """
    # Number of files
    num_models = len(Pastasfiles)

    # saves models
    models = []

    # saves well names for each model
    well_names = []

    # Pastas optimized parameters saved
    pastas_optparam = []

    # For each model
    for num in range(num_models):

        # Loads model
        model = ps.io.load(model_path + "/" + Pastasfiles[num])
        s = Pastasfiles[num]
        result = re.search("_(.*)_GW", s)  # Gets well name
        well_name = result.group(1)
        # Saving optimized parameters
        # If A, initial
        temp_param = [model.parameters["initial"][0]]
        temp_param.extend(model.parameters["optimal"][1:].values)
        pastas_optparam.append(temp_param)
        # Well nest name
        Wellnest_name = Pastasfiles[num][0:8]

        # Getting steady state heads according to aquifer
        if "BK" in well_name:

            initial_d = SS_data.loc[Wellnest_name, "BK"]

        elif "PD" in well_name:

            initial_d = SS_data.loc[Wellnest_name, "PD"]

        elif "NL" in well_name:

            initial_d = SS_data.loc[Wellnest_name, "NL"]

        elif "NB" in well_name:

            initial_d = SS_data.loc[Wellnest_name, "NB"]

        # Setting d parameter to SS heads and to vary +/- initial
        # estimates
        model.set_parameter(name="constant_d",
                            initial=initial_d,
                            pmin=initial_d-10,
                            pmax=initial_d+10,
                            vary=True)

        # Saves
        models.append(model)
        well_names.append(well_name)

    return models, well_names, np.array(pastas_optparam)


def load_Pastas(Pastasfiles, lenfiles, proxyflag, models, well_names,
                model_path, pumpflag, tmin, tmax, pump_series=None,
                pump_path=None, pump_sheet=None, initoptiparam=None):
    """Loads Pastas models

    Pastasfiles - list of Pastas file names
    models - lsit of Pastas model instances
    well_names - list of well names (str)
    lenfiles - how many files there are
    proxyflag - 1 if using available heads as proxy for missing heads
    pump_path - path to pumping excel sheet
    pump_sheet - sheet of specific pumping scenario
    model_path - path to python models
    pumpflag - 1 if changing pumping scenario for Pastas
    tmin, tmax - (str) minimum and maximum year to calculate sub
    initoptiparam - optimal parameters provided

    Returns
    well_data_dates - Well data with matching dates only
    all_well4_data - all well data no matter the date
    """
    # well_data - where to save the data after loading Pasta files
    well_data = []

    # If using avaialbe heads as proxy for missing heads
    if proxyflag == 1:

        # For those missing wells, checks which well is missing
        if lenfiles < 4:

            # If only three files
            # If first file in sorted list is PD, missing BK
            if lenfiles == 3:

                # For each well model
                for num_model in range(len(well_names)):

                    # Loads model, PD as proxy for BK
                    model = models[num_model]

                    curr_well = well_names[num_model]

                    # If changing pumping scenario
                    if pumpflag == 1:

                        # If providing optimal parameters from ESMDA
                        if initoptiparam is None:

                            # If pumping time series given
                            if pump_series is not None:
                                model = pastas_setparam(model,
                                                        well_name=curr_well,
                                                        pump_series=pump_series[
                                                            curr_well])

                            # If path and sheet given instead
                            elif pump_path is not None:
                                # Updating model with new pumping scenario
                                model = pastas_setparam(model, pump_path=pump_path,
                                                        pump_sheet=pump_sheet)

                        else:

                            # If pumping time series given
                            if pump_path is None:
                                if pump_series is not None:
                                    model = pastas_setparam(model,
                                                            initoptiparam=initoptiparam.
                                                            loc[
                                                                curr_well],
                                                            pump_series=pump_series[
                                                                curr_well],
                                                            well_name=curr_well,)

                                else:
                                    model = pastas_setparam(model,
                                                            initoptiparam=initoptiparam.
                                                            loc[
                                                                curr_well],
                                                            well_name=curr_well,)

                            # If path and sheet given instead
                            else:

                                # Updating model with new pumping scenario
                                model = pastas_setparam(model,
                                                        initoptiparam=initoptiparam.
                                                        loc[
                                                            curr_well],
                                                        well_name=curr_well,
                                                        pump_path=pump_path,
                                                        pump_sheet=pump_sheet)

                    # If list starts with PD, means BK is missing
                    if "_PD" in Pastasfiles[0]:

                        # Identifies missing well and index
                        missing = "BK"

                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30, return_warmup=False)
                        temp = temp.rename("Proxy BK")  # Renames column
                        well_data.append(temp)  # Saves data

                        # Loads model, PD
                        temp = temp.rename(curr_well)  # Renames column
                        well_data.append(temp)  # Saves data

                    # Adds NL and NB simulations
                    else:

                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30, return_warmup=False)
                        temp = temp.rename(curr_well)  # Renames col
                        well_data.append(temp)  # Saves data

            # If only two files
            # If first file in sorted list is PD and next is NB
            # missing BK and NL
            if lenfiles == 2:

                # For each well model
                for num_model in range(len(well_names)):

                    # Loads model, PD as proxy for BK
                    model = models[num_model]

                    curr_well = well_names[num_model]

                    # If changing pumping scenario
                    if pumpflag == 1:

                        # If providing optimal parameters from ESMDA
                        if initoptiparam is None:

                            # If pumping time series given
                            if pump_series is not None:
                                model = pastas_setparam(model,
                                                        well_name=curr_well,
                                                        pump_series=pump_series[
                                                            curr_well])

                            # If path and sheet given instead
                            elif pump_path is not None:
                                # Updating model with new pumping scenario
                                model = pastas_setparam(model, pump_path=pump_path,
                                                        pump_sheet=pump_sheet)

                        else:

                            # If pumping time series given
                            if pump_path is None:
                                if pump_series is not None:
                                    model = pastas_setparam(model,
                                                            initoptiparam=initoptiparam.
                                                            loc[
                                                                curr_well],
                                                            pump_series=pump_series[
                                                                curr_well],
                                                            well_name=curr_well,)

                                else:
                                    model = pastas_setparam(model,
                                                            initoptiparam=initoptiparam.
                                                            loc[
                                                                curr_well],
                                                            well_name=curr_well,)

                            # If path and sheet given instead
                            else:

                                # Updating model with new pumping scenario
                                model = pastas_setparam(model,
                                                        initoptiparam=initoptiparam.
                                                        loc[
                                                            curr_well],
                                                        well_name=curr_well,
                                                        pump_path=pump_path,
                                                        pump_sheet=pump_sheet)

                    # If only PD and NB models
                    if np.logical_and("_PD" in Pastasfiles[0],
                                      "_NB" in Pastasfiles[1]):

                        # Identifies missing well and index
                        missing = "NL"

                        # If PD
                        if num_model == 0:

                            proxy_name = "Proxy BK"

                        # If NB
                        else:

                            proxy_name = "Proxy NL"

                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30, return_warmup=False)
                        temp = temp.rename(proxy_name)
                        well_data.append(temp)

                        # Loads model, PD
                        temp = temp.rename(curr_well)
                        well_data.append(temp)

                    # If only two files
                    # If first file in sorted list is PD, next and NL
                    # missing BK and NB
                    if np.logical_and("_PD" in Pastasfiles[0],
                                      "_NL" in Pastasfiles[1]):

                        # Identifies missing well and index
                        missing = "NB"

                        # If PD
                        if num_model == 0:

                            proxy_name = "Proxy BK"
                            temp = model.simulate(tmin="1950", tmax=tmax,
                                                  warmup=365*30,
                                                  return_warmup=False)
                            temp = temp.rename(proxy_name)
                            well_data.append(temp)

                            # Loads model, PD
                            temp = temp.rename(curr_well)
                            well_data.append(temp)

                        # If NB
                        else:

                            proxy_name = "Proxy NB"
                            temp = model.simulate(tmin="1950", tmax=tmax,
                                                  warmup=365*30,
                                                  return_warmup=False)
                            temp = temp.rename(curr_well)
                            well_data.append(temp)

                            # Loads model, NL as proxy for NB
                            temp = temp.rename(proxy_name)
                            well_data.append(temp)

                # If only two files
                # If first file in sorted list is NL, next and NB
                # missing BK and PD

                    if np.logical_and("_NL" in Pastasfiles[0],
                                      "_NB" in Pastasfiles[1]):

                        # Identifies missing well and index
                        missing = "PD"

                        # If NL
                        if num_model == 0:

                            temp = model.simulate(tmin="1950", tmax=tmax,
                                                  warmup=365*30,
                                                  return_warmup=False)
                            temp = temp.rename("Proxy BK")
                            well_data.append(temp)

                            # Loads model, NL as proxy for PD, BK
                            temp = temp.rename("Proxy PD")
                            well_data.append(temp)

                            # Loads model, NL
                            temp = temp.rename(curr_well)
                            well_data.append(temp)

                        # NB
                        else:

                            temp = model.simulate(tmin="1950", tmax=tmax,
                                                  warmup=365*30,
                                                  return_warmup=False)
                            temp = temp.rename(curr_well)
                            well_data.append(temp)

            # If only 1 file
            # missing others
            if lenfiles == 1:

                missing = "others"

                # For each well model
                for num_model in range(len(well_names)):

                    # Loads model, PD as proxy for BK
                    model = models[num_model]

                    curr_well = well_names[num_model]

                    # If changing pumping scenario
                    if pumpflag == 1:

                        # If providing optimal parameters from ESMDA
                        if initoptiparam is None:

                            # If pumping time series given
                            if pump_series is not None:
                                model = pastas_setparam(model,
                                                        well_name=curr_well,
                                                        pump_series=pump_series[
                                                            curr_well])

                            # If path and sheet given instead
                            elif pump_path is not None:
                                # Updating model with new pumping scenario
                                model = pastas_setparam(model, pump_path=pump_path,
                                                        pump_sheet=pump_sheet)

                        else:

                            # If pumping time series given
                            if pump_path is None:
                                if pump_series is not None:
                                    model = pastas_setparam(model,
                                                            initoptiparam=initoptiparam.
                                                            loc[
                                                                curr_well],
                                                            pump_series=pump_series[
                                                                curr_well],
                                                            well_name=curr_well,)

                                else:
                                    model = pastas_setparam(model,
                                                            initoptiparam=initoptiparam.
                                                            loc[
                                                                curr_well],
                                                            well_name=curr_well,)

                            # If path and sheet given instead
                            else:

                                # Updating model with new pumping scenario
                                model = pastas_setparam(model,
                                                        initoptiparam=initoptiparam.
                                                        loc[
                                                            curr_well],
                                                        well_name=curr_well,
                                                        pump_path=pump_path,
                                                        pump_sheet=pump_sheet)

                    if "BK" in curr_well:
                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30,
                                              return_warmup=False)
                        temp = temp.rename(curr_well)
                        well_data.append(temp)

                        # Others proxy
                        temp = temp.rename("Proxy PD")
                        well_data.append(temp)

                        # Others proxy
                        temp = temp.rename("Proxy NL")
                        well_data.append(temp)

                        # Others proxy
                        temp = temp.rename("Proxy NB")
                        well_data.append(temp)

                    elif "PD" in curr_well:

                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30,
                                              return_warmup=False)
                        temp = temp.rename("Proxy BK")
                        well_data.append(temp)

                        # Only head value as proxy for others that are
                        # missing
                        temp = temp.rename(curr_well)
                        well_data.append(temp)

                        # Proxy
                        temp = temp.rename("Proxy NL")
                        well_data.append(temp)

                        # Proxy
                        temp = temp.rename("Proxy NB")
                        well_data.append(temp)

                    elif "NL" in curr_well:

                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30,
                                              return_warmup=False)
                        temp = temp.rename("Proxy BK")
                        well_data.append(temp)

                        # Proxy
                        temp = temp.rename("Proxy PD")
                        well_data.append(temp)

                        # Only head value as proxy for others that are
                        # missing
                        temp = temp.rename(curr_well)
                        well_data.append(temp)

                        # Proxy
                        temp = temp.rename("Proxy NB")
                        well_data.append(temp)

                    elif "NB" in curr_well:

                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30,
                                              return_warmup=False)
                        temp = temp.rename("Proxy BK")
                        well_data.append(temp)

                        # Proxy
                        temp = temp.rename("Proxy PD")
                        well_data.append(temp)

                        # Proxy
                        temp = temp.rename("Proxy NL")
                        well_data.append(temp)

                        # Only head value as proxy for others that are
                        # missing
                        temp = temp.rename(curr_well)
                        well_data.append(temp)

        # No missing wells
        else:

            missing = None

    # If using only available heads
    else:
        missing = None

        # Needs all four wells if proxyflag is not on
        if lenfiles < 4:

            sys.exit("Needs all four wells if proxyflag is not on")

    # If there is no missing data
    if isinstance(missing, type(None)):

        # For each well model
        for num_model in range(len(well_names)):

            # Loads model, PD as proxy for BK
            model = models[num_model]

            curr_well = well_names[num_model]

            # If changing pumping scenario
            if pumpflag == 1:

                # If providing optimal parameters from ESMDA
                if initoptiparam is None:

                    # If pumping time series given
                    if pump_series is not None:
                        model = pastas_setparam(model,
                                                well_name=curr_well,
                                                pump_series=pump_series[
                                                    curr_well])

                    # If path and sheet given instead
                    elif pump_path is not None:
                        # Updating model with new pumping scenario
                        model = pastas_setparam(model, pump_path=pump_path,
                                                pump_sheet=pump_sheet)

                else:

                    # If pumping time series given
                    if pump_path is None:
                        if pump_series is not None:
                            model = pastas_setparam(model,
                                                    initoptiparam=initoptiparam.
                                                    loc[
                                                        curr_well],
                                                    pump_series=pump_series[
                                                        curr_well],
                                                    well_name=curr_well,)

                        else:
                            model = pastas_setparam(model,
                                                    initoptiparam=initoptiparam.
                                                    loc[
                                                        curr_well],
                                                    well_name=curr_well,)

                    # If path and sheet given instead
                    else:

                        # Updating model with new pumping scenario
                        model = pastas_setparam(model,
                                                initoptiparam=initoptiparam.loc[
                                                    curr_well],
                                                well_name=curr_well,
                                                pump_path=pump_path,
                                                pump_sheet=pump_sheet)

            temp = model.simulate(tmin="1950", tmax=tmax,
                                  warmup=365*30, return_warmup=False)
            temp = temp.rename(curr_well)
            well_data.append(temp)

    # Well data with matching dates only
    well_data_dates = functools.reduce(lambda left, right:
                                       pd.merge(left, right,
                                                left_index=True,
                                                right_index=True),
                                       well_data)
    well_data_dates = well_data_dates[
        (well_data_dates.index.year >= int(tmin)) &
        (well_data_dates.index.year <= int(tmax))]

    # All well data
    all_well4_data = functools.reduce(lambda left, right:
                                      pd.concat([left, right], axis=1,
                                                sort=True),
                                      well_data)

    # well_data_dates - Well data with matching dates only
    # all_well4_data - all well data no matter the date
    return well_data_dates, all_well4_data


# Get parameter sample for ESMDA

def get_parameter_sample(n, model, dist, std_mult):
    """


    Parameters
    ----------
    n : int
        number of ensemble members
    model : pastas model instance
    dist : kist of str
        distribution for each param
    std_mult: list of int
        multiplier to standard deviation for each param

    Returns
    -------
    df - dataframe of parameter ensemble [n x n_param] where n_param is
    the number of parameters in the model

    """

    # Get initial parameters
    init = model.get_init_parameters()
    # opt = model.parameters

    # Params
    params = pd.DataFrame()

    # Stats module
    mod = importlib.import_module("scipy.stats")

    # Random values
    # Set the initial parameters to a normal distribution
    for name_i, name in enumerate(model.parameters.index):

        loc = init.initial[name]  # Pastas initial
        scale = abs(std_mult[name_i]/100 * loc)

        # Sampled
        sampled = []

        minimum = model.parameters["pmin"][name]
        maximum = model.parameters["pmax"][name]

        # If the samples don't equal the number of ensemble members wanted
        while len(sampled) != n:

            if dist[name_i] == "lognorm":
                data = np.log(getattr(mod, dist[name_i]).rvs(
                    s=scale, loc=0, scale=np.exp(loc), size=n-len(sampled)))
            elif dist[name_i] == "uniform":
                data = getattr(mod, dist[name_i]).rvs(loc=minimum,
                                                      scale=maximum-minimum,
                                                      size=n-len(sampled))
            else:
                data = getattr(mod, dist[name_i]).rvs(loc=loc, scale=scale,
                                                      size=n-len(sampled))

            sampled.extend(data[(minimum <= data) & (data <= maximum)])

        params[name] = sampled

    return params


# %%###########################################################################
# Loads Pastas models for ESMDA for both pastas and subsidence
##############################################################################


def load_Pastas_ESMDA(Pastasfiles, lenfiles, proxyflag, models, well_names,
                      model_path, pumpflag, tmin, tmax, SS_data, init, std_mult=None,
                      ne=None,
                      pump_series=None,
                      pump_path=None, pump_sheet=None, initoptiparam=None,
                      dist="norm", mode=""):
    """Loads Pastas models

    Pastasfiles - list of Pastas file names
    lenfiles - how many files there are
    models - list of Pastas model instances
    well_names - list of well names (str)
    proxyflag - 1 if using available heads as proxy for missing heads
    pump_path - path to pumping excel sheet
    pump_sheet - sheet of specific pumping scenario
    model_path - path to python models
    SS_data - steady state head data
    pumpflag - 1 if changing pumping scenario for Pastas
    init - 0 if ESMDA run, 1 if initial run
    ne - number of ensemble members
    tmin, tmax - (str) minimum and maximum year to calculate sub
    initoptiparam - optimal parameters provided
    dist- distribution of prior
    mode - pumping or parameter or both estimation

    Returns
    models - with possibly new pumping and new parameters
    well_data_dates - Well data with matching dates only
    all_well4_data - all well data no matter the date
    obs4_data - all observation data no matter the date
    min_pastas - all min of pastas parameters
    max_pastas - all max of pastas parameters
    well_names - list of all well names (offical)
    """
    # well_data - where to save the data after loading Pasta files
    well_data = []

    # Observation data
    obs_data = []

    # Min of pastas parameters
    min_pastas = []
    # Max of pastas parameters
    max_pastas = []

    # Param ensemble
    param_ens = []

    # models
    newmodels = []

    # If using avaialbe heads as proxy for missing heads
    if proxyflag == 1:

        # For those missing wells, checks which well is missing
        if lenfiles < 4:

            # If only three files
            # If first file in sorted list is PD, missing BK
            if lenfiles == 3:
                # For each well model
                for num_model in range(len(well_names)):
                    # Loads model, PD as proxy for BK
                    model = models[num_model]
                    curr_well = well_names[num_model]

                    # If initial run
                    if init == 1:

                        obs_data.append(model.observations())  # Saves obs data

                        if mode != "my_pump" or mode != "my_sub":
                            # Saves min pastas values
                            min_ = model.parameters.loc[:, "pmin"].values
                            # rid of nan
                            # min_[np.isnan(min_)] = model.parameters[
                            #     "optimal"][np.isnan(min_)] - 10
                            min_pastas.append(min_)
                            # Saves max pastas values
                            max_ = model.parameters.loc[:, "pmax"].values
                            # max_[np.isnan(max_)] = model.parameters[
                            #     "optimal"][np.isnan(max_)] + 10
                            max_pastas.append(max_)

                            # Randomize initial starting points for parameter distribution
                            # startparam = []
                            # for n_min, n_max in zip(min_, max_):
                            #     startparam.append([
                            #         np.random.uniform(n_min, n_max, size=1)][0][0])

                            # Saves parameter ensemble
                            params = get_parameter_sample(ne, model, dist, std_mult)
                            param_ens.append(params)

                    # If changing pumping scenario
                    if pumpflag == 1:

                        # If providing optimal parameters from ESMDA
                        if initoptiparam is None:

                            # If pumping time series given
                            if pump_series is not None:
                                model = pastas_setparam(model,
                                                        well_name=curr_well,
                                                        pump_series=pump_series[
                                                            curr_well])

                            # If path and sheet given instead
                            elif pump_path is not None:

                                # Updating model with new pumping scenario
                                model = pastas_setparam(model, pump_path=pump_path,
                                                        pump_sheet=pump_sheet)

                        else:

                            # If pumping time series given
                            if pump_path is None:

                                if pump_series is not None:
                                    model = pastas_setparam(model,
                                                            initoptiparam=initoptiparam.
                                                            loc[
                                                                curr_well],
                                                            pump_series=pump_series[
                                                                curr_well],
                                                            well_name=curr_well,)

                                else:
                                    model = pastas_setparam(model,
                                                            initoptiparam=initoptiparam.
                                                            loc[
                                                                curr_well],
                                                            well_name=curr_well,)

                            # If path and sheet given instead
                            else:

                                # Updating model with new pumping scenario
                                model = pastas_setparam(model,
                                                        initoptiparam=initoptiparam.
                                                        loc[
                                                            curr_well],
                                                        well_name=curr_well,
                                                        pump_path=pump_path,
                                                        pump_sheet=pump_sheet)

                    # Saving new models if new
                    newmodels.append(model)

                    # If list starts with PD, means BK is missing
                    if "_PD" in Pastasfiles[0]:

                        # Identifies missing well and index
                        missing = "BK"
                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30, return_warmup=False)
                        temp = temp.rename("Proxy BK")  # Renames column
                        well_data.append(temp)  # Saves data

                        # Loads model, PD
                        temp = temp.rename(curr_well)  # Renames column
                        well_data.append(temp)  # Saves data

                    # Adds NL and NB simulations
                    else:

                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30, return_warmup=False)
                        temp = temp.rename(curr_well)  # Renames col
                        well_data.append(temp)  # Saves data

            # If only two files
            # If first file in sorted list is PD and next is NB
            # missing BK and NL
            if lenfiles == 2:

                # For each well model
                for num_model in range(len(well_names)):

                    # Loads model, PD as proxy for BK
                    model = models[num_model]

                    curr_well = well_names[num_model]

                    # If initial run
                    if init == 1:
                        obs_data.append(model.observations())  # Saves obs data

                        if mode != "my_pump" or mode != "my_sub":
                            # Saves min pastas values
                            min_ = model.parameters.loc[:, "pmin"].values
                            # rid of nan
                            # min_[np.isnan(min_)] = model.parameters[
                            #     "optimal"][np.isnan(min_)] - 10
                            min_pastas.append(min_)
                            # Saves max pastas values
                            max_ = model.parameters.loc[:, "pmax"].values
                            # max_[np.isnan(max_)] = model.parameters[
                            #     "optimal"][np.isnan(max_)] + 10
                            max_pastas.append(max_)

                            # Randomize initial starting points for parameter distribution
                            # startparam = []
                            # for n_min, n_max in zip(min_, max_):
                            #     startparam.append([
                            #         np.random.uniform(n_min, n_max, size=1)][0][0])

                            # Saves parameter ensemble
                            params = get_parameter_sample(ne, model, dist, std_mult)
                            param_ens.append(params)

                    # If changing pumping scenario
                    if pumpflag == 1:

                        # If providing optimal parameters from ESMDA
                        if initoptiparam is None:

                            # If pumping time series given
                            if pump_series is not None:
                                model = pastas_setparam(model,
                                                        well_name=curr_well,
                                                        pump_series=pump_series[
                                                            curr_well])

                            # If path and sheet given instead
                            elif pump_path is not None:
                                # Updating model with new pumping scenario
                                model = pastas_setparam(model, pump_path=pump_path,
                                                        pump_sheet=pump_sheet)

                        else:

                            # If pumping time series given
                            if pump_path is None:
                                if pump_series is not None:
                                    model = pastas_setparam(model,
                                                            initoptiparam=initoptiparam.
                                                            loc[
                                                                curr_well],
                                                            pump_series=pump_series[
                                                                curr_well],
                                                            well_name=curr_well,)

                                else:
                                    model = pastas_setparam(model,
                                                            initoptiparam=initoptiparam.
                                                            loc[
                                                                curr_well],
                                                            well_name=curr_well,)

                            # If path and sheet given instead
                            else:

                                # Updating model with new pumping scenario
                                model = pastas_setparam(model,
                                                        initoptiparam=initoptiparam.
                                                        loc[
                                                            curr_well],
                                                        well_name=curr_well,
                                                        pump_path=pump_path,
                                                        pump_sheet=pump_sheet)

                    # Saving new models if new
                    newmodels.append(model)

                    # If only PD and NB models
                    if np.logical_and("_PD" in Pastasfiles[0],
                                      "_NB" in Pastasfiles[1]):

                        # Identifies missing well and index
                        missing = "NL"

                        # If PD
                        if num_model == 0:

                            proxy_name = "Proxy BK"

                        # If NB
                        else:

                            proxy_name = "Proxy NL"

                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30, return_warmup=False)
                        temp = temp.rename(proxy_name)
                        well_data.append(temp)

                        # Loads model, PD
                        temp = temp.rename(curr_well)
                        well_data.append(temp)

                    # If only two files
                    # If first file in sorted list is PD, next and NL
                    # missing BK and NB
                    if np.logical_and("_PD" in Pastasfiles[0],
                                      "_NL" in Pastasfiles[1]):

                        # Identifies missing well and index
                        missing = "NB"

                        # If PD
                        if num_model == 0:

                            proxy_name = "Proxy BK"
                            temp = model.simulate(tmin="1950", tmax=tmax,
                                                  warmup=365*30,
                                                  return_warmup=False)
                            temp = temp.rename(proxy_name)
                            well_data.append(temp)

                            # Loads model, PD
                            temp = temp.rename(curr_well)
                            well_data.append(temp)

                        # If NB
                        else:

                            proxy_name = "Proxy NB"
                            temp = model.simulate(tmin="1950", tmax=tmax,
                                                  warmup=365*30,
                                                  return_warmup=False)
                            temp = temp.rename(curr_well)
                            well_data.append(temp)

                            # Loads model, NL as proxy for NB
                            temp = temp.rename(proxy_name)
                            well_data.append(temp)

                # If only two files
                # If first file in sorted list is NL, next and NB
                # missing BK and PD

                    if np.logical_and("_NL" in Pastasfiles[0],
                                      "_NB" in Pastasfiles[1]):

                        # Identifies missing well and index
                        missing = "PD"

                        # If NL
                        if num_model == 0:

                            temp = model.simulate(tmin="1950", tmax=tmax,
                                                  warmup=365*30,
                                                  return_warmup=False)
                            temp = temp.rename("Proxy BK")
                            well_data.append(temp)

                            # Loads model, NL as proxy for PD, BK
                            temp = temp.rename("Proxy PD")
                            well_data.append(temp)

                            # Loads model, NL
                            temp = temp.rename(curr_well)
                            well_data.append(temp)

                        # NB
                        else:

                            temp = model.simulate(tmin="1950", tmax=tmax,
                                                  warmup=365*30,
                                                  return_warmup=False)
                            temp = temp.rename(curr_well)
                            well_data.append(temp)

            # If only 1 file
            # missing others
            if lenfiles == 1:

                missing = "others"

                # For each well model
                for num_model in range(len(well_names)):

                    # Loads model, PD as proxy for BK
                    model = models[num_model]

                    curr_well = well_names[num_model]

                    # If initial run
                    if init == 1:
                        obs_data.append(model.observations())  # Saves obs data

                        if mode != "my_pump" and mode != "my_sub":
                            # Saves min pastas values
                            min_ = model.parameters.loc[:, "pmin"].values
                            # rid of nan
                            # min_[np.isnan(min_)] = model.parameters[
                            #     "optimal"][np.isnan(min_)] - 10
                            min_pastas.append(min_)
                            # Saves max pastas values
                            max_ = model.parameters.loc[:, "pmax"].values
                            # max_[np.isnan(max_)] = model.parameters[
                            #     "optimal"][np.isnan(max_)] + 10
                            max_pastas.append(max_)

                            # Randomize initial starting points for parameter distribution
                            # startparam = []
                            # for n_min, n_max in zip(min_, max_):
                            #     startparam.append([
                            #         np.random.uniform(n_min, n_max, size=1)][0][0])

                            # Saves parameter ensemble
                            params = get_parameter_sample(ne, model, dist, std_mult)
                            param_ens.append(params)

                    # If changing pumping scenario
                    if pumpflag == 1:

                        # If providing optimal parameters from ESMDA
                        if initoptiparam is None:

                            # If pumping time series given
                            if pump_series is not None:
                                model = pastas_setparam(model,
                                                        well_name=curr_well,
                                                        pump_series=pump_series[
                                                            curr_well])

                            # If path and sheet given instead
                            elif pump_path is not None:
                                # Updating model with new pumping scenario
                                model = pastas_setparam(model, pump_path=pump_path,
                                                        pump_sheet=pump_sheet)

                        else:

                            # If pumping time series given
                            if pump_path is None:
                                if pump_series is not None:

                                    model = pastas_setparam(model,
                                                            initoptiparam=initoptiparam.
                                                            loc[
                                                                curr_well],
                                                            pump_series=pump_series[
                                                                curr_well],
                                                            well_name=curr_well,)

                                else:
                                    model = pastas_setparam(model,
                                                            initoptiparam=initoptiparam.
                                                            loc[
                                                                curr_well],
                                                            well_name=curr_well,)

                            # If path and sheet given instead
                            else:

                                # Updating model with new pumping scenario
                                model = pastas_setparam(model,
                                                        initoptiparam=initoptiparam.
                                                        loc[
                                                            curr_well],
                                                        well_name=curr_well,
                                                        pump_path=pump_path,
                                                        pump_sheet=pump_sheet)
                    # Saving new models if new
                    newmodels.append(model)

                    # Identifies missing well and index

                    if "BK" in curr_well:
                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30,
                                              return_warmup=False)
                        temp = temp.rename(curr_well)
                        well_data.append(temp)

                        # Others proxy
                        temp = temp.rename("Proxy PD")
                        well_data.append(temp)

                        # Others proxy
                        temp = temp.rename("Proxy NL")
                        well_data.append(temp)

                        # Others proxy
                        temp = temp.rename("Proxy NB")
                        well_data.append(temp)

                    elif "PD" in curr_well:

                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30,
                                              return_warmup=False)
                        temp = temp.rename("Proxy BK")
                        well_data.append(temp)

                        # Only head value as proxy for others that are
                        # missing
                        temp = temp.rename(curr_well)
                        well_data.append(temp)

                        # Proxy
                        temp = temp.rename("Proxy NL")
                        well_data.append(temp)

                        # Proxy
                        temp = temp.rename("Proxy NB")
                        well_data.append(temp)

                    elif "NL" in curr_well:

                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30,
                                              return_warmup=False)
                        temp = temp.rename("Proxy BK")
                        well_data.append(temp)

                        # Proxy
                        temp = temp.rename("Proxy PD")
                        well_data.append(temp)

                        # Only head value as proxy for others that are
                        # missing
                        temp = temp.rename(curr_well)
                        well_data.append(temp)

                        # Proxy
                        temp = temp.rename("Proxy NB")
                        well_data.append(temp)

                    elif "NB" in curr_well:

                        temp = model.simulate(tmin="1950", tmax=tmax,
                                              warmup=365*30,
                                              return_warmup=False)
                        temp = temp.rename("Proxy BK")
                        well_data.append(temp)

                        # Proxy
                        temp = temp.rename("Proxy PD")
                        well_data.append(temp)

                        # Proxy
                        temp = temp.rename("Proxy NL")
                        well_data.append(temp)

                        # Only head value as proxy for others that are
                        # missing
                        temp = temp.rename(curr_well)
                        well_data.append(temp)

        # No missing wells
        else:

            missing = None

    # If using only available heads
    else:
        missing = None

        # Needs all four wells if proxyflag is not on
        if lenfiles < 4:

            sys.exit("Needs all four wells if proxyflag is not on")

    # If there is no missing data
    if isinstance(missing, type(None)):

        # For each well model
        for num_model in range(len(well_names)):

            # Loads model, PD as proxy for BK
            model = models[num_model]

            curr_well = well_names[num_model]

            # If initial run
            if init == 1:
                obs_data.append(model.observations())  # Saves obs data
                if mode != "my_pump" or mode != "my_sub":
                    # Saves min pastas values
                    min_ = model.parameters.loc[:, "pmin"].values
                    # rid of nan
                    # min_[np.isnan(min_)] = model.parameters[
                    #     "optimal"][np.isnan(min_)] - 10
                    min_pastas.append(min_)
                    # Saves max pastas values
                    max_ = model.parameters.loc[:, "pmax"].values
                    # max_[np.isnan(max_)] = model.parameters[
                    #     "optimal"][np.isnan(max_)] + 10
                    max_pastas.append(max_)

                    # Randomize initial starting points for parameter distribution
                    # startparam = []
                    # for n_min, n_max in zip(min_, max_):
                    #     startparam.append([
                    #         np.random.uniform(n_min, n_max, size=1)][0][0])

                    # Saves parameter ensemble
                    params = get_parameter_sample(ne, model, dist, std_mult)
                    param_ens.append(params)

            # If changing pumping scenario
            if pumpflag == 1:

                # If providing optimal parameters from ESMDA
                if initoptiparam is None:

                    # If pumping time series given
                    if pump_series is not None:
                        model = pastas_setparam(model,
                                                well_name=curr_well,
                                                pump_series=pump_series[
                                                    curr_well])

                    # If path and sheet given instead
                    elif pump_path is not None:
                        # Updating model with new pumping scenario
                        model = pastas_setparam(model, pump_path=pump_path,
                                                pump_sheet=pump_sheet)

                else:

                    # If pumping time series given
                    if pump_path is None:
                        if pump_series is not None:
                            model = pastas_setparam(model,
                                                    initoptiparam=initoptiparam.
                                                    loc[
                                                        curr_well],
                                                    pump_series=pump_series[
                                                        curr_well],
                                                    well_name=curr_well,)

                        else:
                            model = pastas_setparam(model,
                                                    initoptiparam=initoptiparam.
                                                    loc[
                                                        curr_well],
                                                    well_name=curr_well,)

                    # If path and sheet given instead
                    else:

                        # Updating model with new pumping scenario
                        model = pastas_setparam(model,
                                                initoptiparam=initoptiparam.loc[
                                                    curr_well],
                                                well_name=curr_well,
                                                pump_path=pump_path,
                                                pump_sheet=pump_sheet)

            # Saving new models if new
            newmodels.append(model)

            temp = model.simulate(tmin="1950", tmax=tmax,
                                  warmup=365*30,
                                  return_warmup=False)
            temp = temp.rename(curr_well)
            well_data.append(temp)

    # Well data with matching dates only
    well_data_dates = functools.reduce(lambda left, right:
                                       pd.merge(left, right,
                                                left_index=True,
                                                right_index=True),
                                       well_data)
    well_data_dates = well_data_dates[
        (well_data_dates.index.year >= int(tmin)) &
        (well_data_dates.index.year <= int(tmax))]

    # All well data
    all_well4_data = functools.reduce(lambda left, right:
                                      pd.concat([left, right], axis=1,
                                                sort=True),
                                      well_data)

    # If initial run
    if init == 1:
        # All obs data
        all_obs_data = functools.reduce(lambda left, right:
                                        pd.concat([left, right], axis=1,
                                                  sort=True),
                                        obs_data)

        # well_data_dates - Well data with matching dates only
        # all_well4_data - all well data no matter the date
        # all_obs_data - all well data no matter the date
        # min_pastas - max of pastas parameter values
        # max_pastas - min of pastas parameter values
        return well_data_dates, all_well4_data, all_obs_data, np.array(min_pastas), \
            np.array(max_pastas), param_ens, well_names, newmodels

    else:

        return well_data_dates, all_well4_data


# %%###########################################################################
# Sets initial condition for subsidence run
##############################################################################

# Assuming has data for all four aquifers
# Assuming conceptual model of clay above BK, between BK and PD, PD and NL, NL
# and NB for a total of 4 clay layers.
def set_ic(headb, headt, i, mode, fullheadt, fullheadb,
           tmin, tmax, SS_data, wellnest, aq_namet, aq_nameb,
           Kv_cl, Sskv_cl, Sske_cl, Sske_aqt, Sske_aqb, CC, Nz,
           Thick_cl, nclay, Thick_aqt, Thick_aqb, Nt,
           all_well_data=None):
    """Runs groundwater models for aquifers to set initial conditions for clay
    layers

    headb - head of bottom aquifer
    headt - head of top aquifer
    i - current clay layer (1 - 4) where 1 is the top clay layer
    mode - raw groundwater data or time series from pastas (raw needs to
    to be interpolated). options: raw, pastas
    all_well_data - raw observed groundwater data not within tmin and tmax
    fullheadt - all Pastas simulated groundwater data from top aquifer despite the
    date
    fullheadb - all Pastas simulated groundwater data from bottom aquifer despite
    the date
    tmin, tmax - (str) minimum and maximum year to calculate sub
    SS_data - steady state heads relative to land surface from coastal dem 2.1. SS
    heads taken from MODFLOW model. Only used with no data from 1950 onwards but
    Pastas was used to simulate from 1950 onwards. Shouldn't be used but is an
    option
    wellnest - string of well nest name
    aq_namet - name of top aquifer
    aq_nameb - name of bottom aquifer
    Kv_cl - value of vertical hydraulic conductivity of clay layer
    Sskv_cl - value of inelastic specific storage of clay layer
    Sske_cl - value of elastic specific storage of clay layer
    Sske_aqt - value of elastic specific storage of top aquifer layer
    Sske_aqb - value of elastic specific storage of bottom aquifer layer
    CC - convergence criteria
    Nz - number of nodes in the z direction
    Thick_cl - Thickness of clay layer
    nclay - number of clay model layers
    Thick_aqt - Thickness of top aquifer
    Thick_aqb - thickness of bottom aquifer
    Nt - number of time steps

    Returns
    t_ic - time for spin up run
    h_ic - head for spin up run for clay model layers
    """

    # Create daily time series
    df = pd.DataFrame(index=pd.date_range("1950-01-01",
                                          headb.index[0],
                                          freq="d"))

    # time
    # Creating time time series [0: len of time series]
    timet_ic = np.arange(len(df.index))

    # If not the first clay layer
    if i != 1:

        # First official model aquifer head in top and bottom
        headt1 = headt.iloc[0]
        headb1 = headb.iloc[0]

        # Getting subset of dates that are before tmin to be used
        # in linear interpolation of head
        # Top
        if mode == "raw":
            subsetdate_t = all_well_data.index[np.logical_and(
                ~all_well_data.iloc[:, i-2].isna(),
                all_well_data.index < headt.index[0])]
            interpdata = all_well_data.loc[
                subsetdate_t].iloc[:, i-2]
        elif mode == "Pastas":
            subsetdate_t = fullheadt.index[np.logical_and(
                ~fullheadt.isna(),
                fullheadt.index.year < int(tmin))]
            interpdata = fullheadt.loc[subsetdate_t]

        # Getting subset of index of those dates that are before
        # tmin to be used in linear interpolation of head
        subsetindex_t = []

        for j in range(len(subsetdate_t)):
            subsetindex_t = np.append(subsetindex_t,
                                      np.flatnonzero(
                                          df.index ==
                                          subsetdate_t[j]))

        # If no earlier GW obs before model start
        if len(subsetindex_t) == 0:

            # Two values and will interpolate between
            timet2_ic = [0, timet_ic[-1]]
            headt2_ic = [SS_data.loc[wellnest, aq_namet], headt1]

        # If there are GW obs before model start, uses it for
        # llnear interpolation with SS heads
        else:

            # Converting to int
            subsetindex_t = subsetindex_t.astype(int)

            # Values and will interpolate between; time for
            # interpolation
            timet2_ic = np.insert(subsetindex_t, 0, 0)
            timet2_ic = np.append(timet2_ic, timet_ic[-1])

            # SS, head before model run, and first model head
            # Values and will interpolate between
            # Top aquifer
            headt2_ic_subset = interpdata.values
            headt2_ic = np.insert(headt2_ic_subset, 0,
                                  SS_data.loc[wellnest, aq_namet])
            headt2_ic = np.append(headt2_ic, headt1)

        # Bottom
        if mode == "raw":
            subsetdate_b = all_well_data.index[np.logical_and(
                ~all_well_data.iloc[:, i-1].isna(),
                all_well_data.index < headb.index[0])]
            interpdata = all_well_data.loc[
                subsetdate_b].iloc[:, i-1]
        elif mode == "Pastas":
            subsetdate_b = fullheadb.index[np.logical_and(
                ~fullheadb.isna(),
                fullheadb.index.year < int(tmin))]
            interpdata = fullheadb.loc[subsetdate_b]

        # Getting subset of index of those dates that are before
        # tmin to be used in linear interpolation of head
        subsetindex_b = []

        for j in range(len(subsetdate_b)):
            subsetindex_b = np.append(subsetindex_b,
                                      np.flatnonzero(
                                          df.index ==
                                          subsetdate_b[j]))

        # If no earlier GW obs before model start
        if len(subsetindex_b) == 0:
            # Two values and will interpolate between
            timeb2_ic = [0, timet_ic[-1]]
            headb2_ic = [SS_data.loc[wellnest, aq_nameb], headb1]

        # If there are GW obs before model start, uses it for
        # llnear interpolation with SS heads
        else:

            # Converting to int
            subsetindex_b = subsetindex_b.astype(int)

            # Values and will interpolate between; time for
            # interpolation
            timeb2_ic = np.insert(subsetindex_b, 0, 0)
            timeb2_ic = np.append(timeb2_ic, timet_ic[-1])

            # SS, head before model run, and first model head
            # Values and will interpolate between
            # Bottom aquifer
            headb2_ic_subset = interpdata.values
            headb2_ic = np.insert(headb2_ic_subset, 0,
                                  SS_data.loc[wellnest, aq_nameb])
            headb2_ic = np.append(headb2_ic, headb1)

        # Interpolating
        headb_ic = pd.Series(np.interp(timet_ic, timeb2_ic,
                                       headb2_ic))  # Linear
        headb_ic.set_index = df.index
        headt_ic = pd.Series(np.interp(timet_ic, timet2_ic,
                                       headt2_ic))  # Linear
        headt_ic.set_index = df.index

        # Using Pastas constant d for initial condition
        # Linearly interpolated between top and bottom
        # constant d
        spacing = np.linspace(0, Nz+1, num=Nz+2, endpoint=True)
        constant_d_ic = np.interp(spacing,
                                  [0, Nz+1],
                                  [SS_data.loc[wellnest, aq_namet],
                                   SS_data.loc[wellnest, aq_nameb]])

    # If top clay layer i == 1
    else:
        # Last spin up run is the first value in the first model
        # run
        # First official model aquifer head in top and bottom
        headb1 = headb.iloc[0]

        # Getting subset of dates that are before tmin to be used
        # in linear interpolation of head
        # Bottom
        if mode == "raw":
            subsetdate_b = all_well_data.index[np.logical_and(
                ~all_well_data.iloc[:, i-1].isna(),
                all_well_data.index < headb.index[0])]
            interpdata = all_well_data.loc[subsetdate_b].iloc[:, i]
        elif mode == "Pastas":
            subsetdate_b = fullheadb.index[np.logical_and(
                ~fullheadb.isna(),
                fullheadb.index.year < int(tmin))]
            interpdata = fullheadb.loc[subsetdate_b]

        # Getting subset of index of those dates that are before
        # tmin to be used in linear interpolation of head
        subsetindex_b = []

        for j in range(len(subsetdate_b)):
            subsetindex_b = np.append(subsetindex_b,
                                      np.flatnonzero(
                                          df.index ==
                                          subsetdate_b[j]))

        # If no earlier GW obs before model start
        if len(subsetindex_b) == 0:
            # Two values and will interpolate between
            timeb2_ic = [0, timet_ic[-1]]
            headb2_ic = [SS_data.loc[wellnest, aq_nameb], headb1]

        # If there are GW obs before model start, uses it for
        # llnear interpolation with SS heads
        else:

            # Converting to int
            subsetindex_b = subsetindex_b.astype(int)

            # Values and will interpolate between; time for
            # interpolation
            timeb2_ic = np.insert(subsetindex_b, 0, 0)
            timeb2_ic = np.append(timeb2_ic, timet_ic[-1])

            # SS, head before model run, and first model head
            # Values and will interpolate between
            # Bottom aquifer
            headb2_ic_subset = interpdata.values
            headb2_ic = np.insert(headb2_ic_subset, 0,
                                  SS_data.loc[wellnest, aq_nameb])
            headb2_ic = np.append(headb2_ic, headb1)

        # Interpolating
        headb_ic = pd.Series(np.interp(timet_ic, timeb2_ic,
                                       headb2_ic))  # Linear
        headb_ic.set_index = df.index

        headt_ic = None

        # Using Pastas constant d for initial condition
        # Linearly interpolated between top and bottom
        # constant d
        spacing = np.linspace(0, Nz+1, num=Nz+2, endpoint=True)
        constant_d_ic = np.interp(spacing,
                                  [0, Nz+1],
                                  [0, SS_data.loc[wellnest, aq_nameb]])

    # print(wellnest, " Clay " + str(i) + " Initial Condition\n")
    # Calculates sub
    # Returns interpolated t, cum sub total, interp top head, bot
    # head, cum sub inelastic, head matrix with top and bottom row
    # as top and bottom aquifer (row is node, column is time)
    t_ic, _, _, _, _, h_ic = \
        calc_deformation(timet_ic, headt_ic, headb_ic, Kv_cl,
                         Sskv_cl, Sske_cl, Sske_sandt=Sske_aqt,
                         Sske_sandb=Sske_aqb, claythick=Thick_cl,
                         nclay=nclay, sandthickt=Thick_aqt,
                         sandthickb=Thick_aqb,
                         Nz=Nz, CC=CC, Nt=Nt,
                         ic=constant_d_ic)

    # t_ic - time for spin up run
    # h_ic - head for spin up run for clay model layers
    return t_ic, h_ic


# %%###########################################################################
# Runs the bulk of code of the subsidence model for the four clay layers
##############################################################################

# Assuming has data for all four aquifers
# Assuming conceptual model of clay above BK, between BK and PD, PD and NL, NL
# and NB for a total of 4 clay layers.
def run_sub(num_clay, all_well4_data, well_data_dates, mode,
            tmin, tmax, SS_data, wellnest, K_data, Sskv_data, Sske_data, CC, Nz,
            Thick_data, ic_run, sub_total, subv_total, all_results,
            well_data=None, user_ic=None):
    """Runs code for bulk of subsidence modeling

    num_clay - number of clay layers
    Thick_data - thickness data for each well nest and layer
    all_well4_data - all well data no matter the date
    well_data_dates - well data with only overlapping dates
    mode - raw groundwater data or time series from pastas (raw needs to
    to be interpolated). options: raw, pastas
    all_well_data - raw observed groundwater data not within tmin and tmax
    tmin, tmax - (str) minimum and maximum year to calculate sub
    K_data - vertical hydraulic conductivity data for each well nest and layer
    Sske_data - elastic specific storage data for each well nest and layer
    Sskv_data - inelastic specific storage data for each well nest and layer
    SS_data - steady state heads relative to land surface from coastal dem 2.1. SS
    heads taken from MODFLOW model. Only used with no data from 1950 onwards but
    Pastas was used to simulate from 1950 onwards. Shouldn't be used but is an
    option
    wellnest - string of well nest name
    ic_run - flag if running spin up run
    sub_total - list of lists (stores results for total subsidence)
    subv_total - list of lists (stores results for inelastic sub)
    all_results - list of lists (stores all results)
    user_ic - if user provides t_ic (index 0) and h_ic (index 1) (list of lists i.e. t_ic
    contains all the t_ic for multiple clays)

    Returns
    sub_total - list of lists (stores results for total subsidence)
    subv_total - list of lists (stores results for inelastic sub)
    all_results - list of lists (stores all results)
    """

    # Keeps track of current z (bottom of layer)
    curr_z = 0

    # For each model
    for i in range(1, num_clay+1):

        # print(wellnest, " Clay " + str(i))

        # Specifies index and clay layer names
        # VSC = very soft clay, MSC = medium stiff clay, SC = stiff
        # clay, HC = hard clay

        # If clay layer BK aquifer
        if i == 1:
            clay_name = "VSC"
            aq_namet = "BK"
            aq_nameb = "BK"

        # If clay layer between BK and PD aquifer
        elif i == 2:
            clay_name = "MSC"
            aq_namet = "BK"
            aq_nameb = "PD"

        # If clay layer between PD and NL aquifer
        elif i == 3:
            clay_name = "SC"
            aq_namet = "PD"
            aq_nameb = "NL"

        # If clay layer between NL and NB aquifer
        elif i == 4:
            clay_name = "HC"
            aq_namet = "NL"
            aq_nameb = "NB"

        # Thickness data, thickness for the clay layer, and  top and
        # bottom aquifer
        Thick_cl = Thick_data.loc[wellnest, clay_name]
        Thick_aqb = Thick_data.loc[wellnest, aq_nameb]
        Thick_aqt = Thick_data.loc[wellnest, aq_namet]

        # Time for both aquifers is the same
        # If clay layer above BK, no aquifer above it
        if i == 1:

            if mode == "Pastas":

                # BK head
                # Only bottom aquifer
                fullheadt = None
                fullheadb = all_well4_data.iloc[:, i-1]
                headb = well_data_dates.iloc[:, i-1]

            elif mode == "raw":

                # BK head
                # No top aquifer, only bottom aquifer
                headb = well_data.iloc[:, i-1]

            # No top aquifer
            headt = None

            # Thickness/Specific storage of top aquifer is 0 because
            # it doesn't exist
            Thick_aqt = 0
            Sske_aqt = 0

        # All other clay layers not first or last
        elif i != 4:

            Thick_aqb /= 2  # NB aquifer not halved.
            # Not simulating clay below it

        # If not first aquifer
        if i != 1:

            if mode == "Pastas":

                fullheadb = all_well4_data.iloc[:, i-1]
                headb = well_data_dates.iloc[:, i-1]

                fullheadt = all_well4_data.iloc[:, i-2]
                headt = well_data_dates.iloc[:, i-2]

            elif mode == "raw":

                headb = well_data.iloc[:, i-1]
                headt = well_data.iloc[:, i-2]

            Sske_aqt = Sske_data.loc[wellnest, aq_namet]

        # Creating time time series [0: len of time series]
        timet = np.arange(len(headb.index))

        # Thickness of top aquifer needs to be halved
        # For all clay layers (top will be zero even if halved)
        Thick_aqt /= 2

        # Specific storage for clays, needed for DELAY CALCULATIONS
        # Inelastic (v) and elastic (e)
        Sskv_cl = Sskv_data.loc[wellnest, clay_name]
        Sske_cl = Sske_data.loc[wellnest, clay_name]
        Sske_aqb = Sske_data.loc[wellnest, aq_nameb]

        # Kv for clays (m/day)
        # Assuming Kv = Kh
        # Using Chula value for BK clay for all clay values as starting
        Kv_cl = K_data.loc[wellnest, clay_name]

        # Number of clay layers
        nclay = 1

        # Number of time steps
        Nt = 100

        # z distribution, not used for calculation
        # Only for plotting adnd reference
        # mesh points in space
        # Current z is at the bottom of the top aquifer
        curr_z += Thick_aqt * 2

        # Z distribution from bottom of top aq to
        # bottom of clay
        dz = Thick_cl/Nz
        z = np.arange(curr_z + dz/2,
                      curr_z + Thick_cl+dz/2,
                      dz)
        # Current z updated to now bottom of clay
        z = np.insert(z, 0, curr_z)
        curr_z += Thick_cl
        z = np.append(z, curr_z)

        # If running transient simulation before model run
        # to get clay heads to where they need to be
        if ic_run:

            # user provided t and h from initial condition interpolation
            if user_ic:
                t_ic = user_ic[0][i-1]
                h_ic = np.array(user_ic[1][i-1])

                # Calculates sub
                # Returns interpolated t, cum sub total, interp top head, bot
                # head, cum sub inelastic, head matrix with top and bottom row
                # as top and bottom aquifer (row is node, column is time)
                interp_t, sub, boundaryt, boundaryb, sub_v, h = \
                    calc_deformation(timet, headt, headb, Kv_cl,
                                     Sskv_cl, Sske_cl, Sske_sandt=Sske_aqt,
                                     Sske_sandb=Sske_aqb, claythick=Thick_cl,
                                     nclay=nclay, sandthickt=Thick_aqt,
                                     sandthickb=Thick_aqb,
                                     Nz=Nz, CC=CC, Nt=Nt,
                                     ic=h_ic)

            # if user didn't provide, calculates
            else:

                t_ic, h_ic = set_ic(headb, headt, i, mode, fullheadt,
                                    fullheadb, tmin, tmax, SS_data, wellnest,
                                    aq_namet, aq_nameb, Kv_cl, Sskv_cl, Sske_cl,
                                    Sske_aqt, Sske_aqb, CC, Nz, Thick_cl, nclay,
                                    Thick_aqt, Thick_aqb, Nt)

                # Calculates sub
                # Returns interpolated t, cum sub total, interp top head, bot
                # head, cum sub inelastic, head matrix with top and bottom row
                # as top and bottom aquifer (row is node, column is time)
                interp_t, sub, boundaryt, boundaryb, sub_v, h = \
                    calc_deformation(timet, headt, headb, Kv_cl,
                                     Sskv_cl, Sske_cl, Sske_sandt=Sske_aqt,
                                     Sske_sandb=Sske_aqb, claythick=Thick_cl,
                                     nclay=nclay, sandthickt=Thick_aqt,
                                     sandthickb=Thick_aqb,
                                     Nz=Nz, CC=CC, Nt=Nt,
                                     ic=h_ic[:, -1])

        # If not running to get initial condition
        else:

            # Calculates sub
            # Returns interpolated t, cum sub total, interp top head, bot
            # head, cum sub inelastic, head matrix with top and bottom row
            # as top and bottom aquifer (row is node, column is time)
            interp_t, sub, boundaryt, boundaryb, sub_v, h = \
                calc_deformation(timet, headt, headb, Kv_cl,
                                 Sskv_cl, Sske_cl, Sske_sandt=Sske_aqt,
                                 Sske_sandb=Sske_aqb, claythick=Thick_cl,
                                 nclay=nclay, sandthickt=Thick_aqt,
                                 sandthickb=Thick_aqb,
                                 Nz=Nz, CC=CC, Nt=Nt)
        # Well names
        if mode == "Pastas":

            # Getting well name from Pastas model file name
            well_name = all_well4_data.columns[i-1]

        elif mode == "raw":

            well_name = well_data.columns[i-1]

        # Adds subsidence to total of all clay
        # Stores records as wellnest, well, data in list
        sub_total.append([wellnest, well_name,
                          interp_t, sub])
        subv_total.append([wellnest, well_name,
                           interp_t, sub_v])

        # If running transient simulation before model run
        # to get clay heads to where they need to be
        if ic_run:

            # Saves heads in clay nodes, z distribution
            # time original (original time series (0:len(date))), date
            # Saves initial condition head and initial condition time
            all_results.append([wellnest, well_name,
                                timet, headb.index, h, z, t_ic, h_ic])

        else:

            # Saves heads in clay nodes, z distribution
            # time original (original time series (0:len(date))), date
            all_results.append([wellnest, well_name,
                                timet, headb.index, h, z])

    return sub_total, subv_total, all_results


# %%###########################################################################
# ESMDA implementation
##############################################################################

# Random number generator
rng = np.random.default_rng()


def forward_sub(p, wellnestlist,
                mode, tmin, tmax,
                Thick_data, K_data, Sskv_data,
                Sske_data, CC,
                Nz, num_clay, all_well4_data,
                well_data_dates,
                wellnest,
                SS_data, p_multop,
                ic_run, return_sub, dobs,
                gw_obs_indices,
                well_names, hidstate_t=None):
    """Running pastas and subsidence model for sub only

    p - list of param multiplers (1 for Sskv, 2 for Sske, 3 for K)
    wellnestlist - list of wellnest to calculate subsidence for
    mode - raw groundwater data or time series from pastas (raw needs to
    to be interpolated). options: raw, pastas
    tmin, tmax - (str) minimum and maximum year to calculate sub
    CC - convergence criteria
    Nz - number of nodes in the z direction
    ic_run - True or false to generate initial condition run for clays
    proxyflag - 1 if using available heads as proxy for missing heads
    pumpflag - 1 if changing pumping scenario for Pastas
    model_path - path to python models
    pump_path - path to pumping excel sheet
    pump_sheet - sheet of specific pumping scenario
    califlag - 1 if calibrating subsidence model
    dobs - observations
    hidstate_t - initial hidden state time given by user

    The data sets have specific names for clays and aquifers
    Thick_data - thickness of clay and aquifers
    K_data - vertical hydraulic conductivity of clay and aquifers
    Sskv - inelastic specific storage term
    Sske - elastic specific storage term

    Returns
    y - subsidence (cm/yr) time series
    """

    # Preallocation
    # Head time series for each  node
    all_results = []

    # Subsidence sum for all clay layers
    sub_total = []

    # Inelastic subsidence sum for all clay layers
    subv_total = []

    # Temporary assignments
    temp_Sskv = Sskv_data.copy()
    temp_Sske = Sske_data.copy()
    temp_K = K_data.copy()

    # If working with multipliers
    if p_multop[0]:

        # If only calibrating Sskv and Sske
        if p_multop[1] == "Sskv":

            temp_Sskv.loc[wellnestlist[0]] = np.array(
                Sskv_data.loc[wellnestlist[0]]) * np.array(np.exp(p[0]))

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + 1):((numt+1)*len_nodes+1)]
                    for numt in range(len(hidstate_t))]

        elif p_multop[1] == "Sske":

            temp_Sske.iloc[0, ::2] = Sskv_data.loc[wellnestlist[0]].iloc[::2] * \
                np.array(np.exp(p[0]))
            temp_Sske.iloc[0, 1::2] = temp_Sske.iloc[0, 0::2] / 10

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + 1):((numt+1)*len_nodes+1)]
                    for numt in range(len(hidstate_t))]

        # If calibrating K
        elif p_multop[1] == "K":

            temp_K.loc[wellnestlist[0]] = np.array(
                K_data.loc[wellnestlist[0]]) * np.array(np.exp(p[0]))

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + 1):((numt+1)*len_nodes+1)]
                    for numt in range(len(hidstate_t))]

        # If calibrating Sskv and Sske and K
        elif p_multop[1] == "all":

            # BASIN HOPPING
            temp_Sskv.loc[wellnestlist[0]] = np.array(
                Sskv_data.loc[wellnestlist[0]]) * np.array(np.exp(p[0]))

            temp_Sske.iloc[0, ::2] = Sskv_data.loc[wellnestlist[0]].iloc[::2] * \
                np.array(np.exp(p[1]))
            temp_Sske.iloc[0, 1::2] = temp_Sske.iloc[0, 0::2] / 10

            temp_K.loc[wellnestlist[0]] = np.array(
                K_data.loc[wellnestlist[0]]) * np.array(np.exp(p[2]))

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + 3):((numt+1)*len_nodes+3)]
                    for numt in range(len(hidstate_t))]

        # If calibrating Sskv and Sske
        elif p_multop[1] == "Ss":

            # BASIN HOPPING
            temp_Sskv.loc[wellnestlist[0]] = np.array(
                Sskv_data.loc[wellnestlist[0]]) * np.array(np.exp(p[0]))

            temp_Sske.iloc[0, ::2] = Sskv_data.loc[wellnestlist[0]].iloc[::2] * \
                np.array(np.exp(p[1]))
            temp_Sske.iloc[0, 1::2] = temp_Sske.iloc[0, 0::2] / 10

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + 2):((numt+1)*len_nodes+2)]
                    for numt in range(len(hidstate_t))]

        # If calibrating Sskv and K
        elif p_multop[1] == "SsK":

            # BASIN HOPPING
            temp_Sskv.loc[wellnestlist[0]] = np.array(
                Sskv_data.loc[wellnestlist[0]]) * np.array(np.exp(p[0]))

            temp_K.loc[wellnestlist[0]] = np.array(
                K_data.loc[wellnestlist[0]]) * np.array(np.exp(p[1]))

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + 2):((numt+1)*len_nodes+2)]
                    for numt in range(len(hidstate_t))]

        # Catching mistake
        else:

            sys.exit("Sorry: wrong settings")

    # Working with non multipliers and calculating each layer value
    else:

        temp_Sskv.loc[wellnestlist[0]] = p[0:8]
        temp_Sske.loc[wellnestlist[0]] = p[8:16]
        temp_K.loc[wellnestlist[0]] = p[16:]

    # Calculates subsidence without reloading Pastas models
    sub_total, subv_total, all_results = run_sub(num_clay, all_well4_data,
                                                 well_data_dates, mode,
                                                 tmin, tmax, SS_data,
                                                 wellnest,
                                                 temp_K, temp_Sskv, temp_Sske,
                                                 CC, Nz, Thick_data, ic_run,
                                                 sub_total, subv_total,
                                                 all_results, user_ic=[hidstate_t, h_ic])

    # Post process data
    sub_total, subv_total, ann_sub, \
        avgsub = bkk_postproc(wellnestlist,
                              sub_total,
                              subv_total,
                              all_results)

    # preparation
    daterange = pd.date_range(dt.datetime(1978, 12, 31), periods=43,
                              freq="Y").tolist()
    df = pd.DataFrame(daterange, columns=["date"])

    # annual data in cm
    plot_data = df.merge(ann_sub[0][1]*100, left_on=df.date,
                         right_on=ann_sub[0][1].index,
                         how="left")

    # Renaming for other merge
    plot_data = plot_data.rename(columns={"key_0": "key0"})

    # Filling na with 0
    plot_data = plot_data.fillna(0)

    # Benchamrks already in cm
    plot_data = plot_data.merge(dobs,
                                left_on=plot_data.key0,
                                right_on=pd.to_datetime(dobs.to_frame().index),
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
    d_pred = plot_data.AnnRates[landlevel != 0]

    return d_pred


def forwardmy_gwparamSELECT(p, models,
                            dobs,
                            gw_obs_indices,
                            pump_path=None,
                            pump_sheet=None,
                            params_i=None,
                            other_params=None,
                            others_i=None):
    """Running pastas model for py esmda with only some parameters

    p - list of param multiplers (1 for Sskv, 2 for Sske, 3 for K)
    wellnestlist - list of wellnest to calculate subsidence for
    mode - raw groundwater data or time series from pastas (raw needs to
    to be interpolated). options: raw, pastas
    tmin, tmax - (str) minimum and maximum year to calculate sub
    CC - convergence criteria
    Nz - number of nodes in the z direction
    ic_run - True or false to generate initial condition run for clays
    proxyflag - 1 if using available heads as proxy for missing heads
    pumpflag - 1 if changing pumping scenario for Pastas
    model_path - path to python models
    pump_path - path to pumping excel sheet
    pump_sheet - sheet of specific pumping scenario
    califlag - 1 if calibrating subsidence model
    dobs - observations
    params_i - indices of parameters to calibrate
    other_params - other parameter values to not calibrate
    other_i - other index

    The data sets have specific names for clays and aquifers
    Thick_data - thickness of clay and aquifers
    K_data - vertical hydraulic conductivity of clay and aquifers
    Sskv - inelastic specific storage term
    Sske - elastic specific storage term

    Returns
    y - subsidence (cm/yr) time series
    """

    # Initiate an array of predicted results.
    d_pred = np.zeros([len(dobs)])

    # Length of parameters
    n_param = len(params_i)

    # Initiate an array of predicted results.
    for model_i, model in enumerate(models):

        # parameter names
        param_name = model.parameters.index[params_i].values

        for param_i, param in enumerate(param_name):

            model.set_parameter(name=param,
                                initial=p[model_i*n_param+param_i],
                                optimal=p[model_i*n_param+param_i])

        if n_param != 4:

            # Other names
            other_name = model.parameters.index[others_i].values

            for other_i, other in enumerate(other_name):

                model.set_parameter(name=other,
                                    initial=other_params[model_i][others_i[other_i]],
                                    optimal=other_params[model_i][others_i[other_i]])

        tmin = str(gw_obs_indices[model_i][0].year)
        # Give an extra year since it stops at the start of the year
        tmax = str(gw_obs_indices[model_i][-1].year + 1)
        # Newly simulated head
        d_pred_temp = model.simulate(tmin=tmin, tmax=tmax)

        # Stripping d_pred to match dobs
        # Matching obs to predicted
        df = pd.DataFrame(gw_obs_indices[model_i],
                          columns=["date"], index=gw_obs_indices[model_i])

        len_ = len(gw_obs_indices[model_i])

        # annual data in cm
        plot_data = df.merge(d_pred_temp, left_on=df.index,
                             right_on=d_pred_temp.index,
                             how="left")

        # Renaming for second merge
        plot_data = plot_data.rename(columns={"key_0": "key0"})

        if model_i == 0:

            index0 = 0
            nd_temp_1 = 0 + len_

        else:

            index0 += len(gw_obs_indices[model_i-1])
            nd_temp_1 += len_

        # Renaming for ESMDA code later on
        plot_data = plot_data.rename(columns={"Simulation": "AnnRates"})

        d_pred[index0:nd_temp_1] = plot_data.AnnRates
        # Renaming for second merge
        plot_data = plot_data.rename(columns={"AnnRates ": "Simulation0"})

    return d_pred


def forward_gwparam_subSELECT_pump_my(p, wellnestlist, Pastasfiles, lenfiles,
                                      proxyflag, model_path, pumpflag,
                                      mode, tmin, tmax,
                                      Thick_data, K_data, Sskv_data,
                                      Sske_data, CC,
                                      Nz, num_clay, all_well4_data,
                                      well_data_dates,
                                      wellnest,
                                      SS_data, p_multop,
                                      ic_run, return_sub, dobs,
                                      gw_obs_indices, models,
                                      well_names,
                                      n_pump, annual_date_pump,
                                      daily_date_pump, nd_sub, nd_gw,
                                      initoptiparam=None,
                                      pump_path=None,
                                      pump_sheet=None,
                                      pump_series=None,
                                      hidstate_t=None, hidstate=None,
                                      params_i=None,
                                      other_params=None,
                                      others_i=None):
    """Running pastas and subsidence model for esmda
    Selected parameters from gw and subsidence + pumping + clay heads


    p - list of param multiplers (1 for Sskv, 2 for Sske, 3 for K)
    wellnestlist - list of wellnest to calculate subsidence for
    mode - raw groundwater data or time series from pastas (raw needs to
    to be interpolated). options: raw, pastas
    tmin, tmax - (str) minimum and maximum year to calculate sub
    CC - convergence criteria
    Nz - number of nodes in the z direction
    ic_run - True or false to generate initial condition run for clays
    proxyflag - 1 if using available heads as proxy for missing heads
    pumpflag - 1 if changing pumping scenario for Pastas
    model_path - path to python models
    pump_path - path to pumping excel sheet
    pump_sheet - sheet of specific pumping scenario
    califlag - 1 if calibrating subsidence model
    dobs - observations
    hidstate_t - initial hidden state time given by user
    hidstate - initial hidden state given by user

    The data sets have specific names for clays and aquifers
    Thick_data - thickness of clay and aquifers
    K_data - vertical hydraulic conductivity of clay and aquifers
    Sskv - inelastic specific storage term
    Sske - elastic specific storage term

    Returns
    y - subsidence (cm/yr) time series
    """

    # Initiate an array of predicted results.
    d_pred = np.zeros([len(dobs)])

    # Preallocation
    # Head time series for each  node
    all_results = []

    # Subsidence sum for all clay layers
    sub_total = []

    # Inelastic subsidence sum for all clay layers
    subv_total = []

    # Optimal params for running ESMDA
    # Number of wells
    num_wells = len(models)

    # Length of parameters
    n_param = len(params_i)

    # list to dataframe for pumping
    pumping_data = []

    # Preallocation of new models with new parameters
    newmodels = []

    # Temporary assignments
    temp_Sskv = Sskv_data.copy()
    temp_Sske = Sske_data.copy()
    temp_K = K_data.copy()
    # If working with multipliers
    if p_multop[0]:

        # If only calibrating Sskv and Sske
        if p_multop[1] == "Sskv":

            temp_Sskv.loc[wellnestlist[0]] = np.array(
                Sskv_data.loc[wellnestlist[0]]) * np.array(
                    np.exp(p[num_wells*n_param+n_pump]))

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + (
                num_wells*n_param+n_pump+1)):(
                    (numt+1)*len_nodes+(
                        num_wells*n_param+n_pump+1))] for numt in range(len(hidstate_t))]

        elif p_multop[1] == "Sske":

            temp_Sske.iloc[0, ::2] = Sskv_data.loc[wellnestlist[0]].iloc[::2] * \
                np.array(np.exp(p[num_wells*n_param+n_pump]))
            temp_Sske.iloc[0, 1::2] = temp_Sske.iloc[0, 0::2] / 10

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + 1):((numt+1)*len_nodes+1)]
                    for numt in range(len(hidstate_t))]

        # If calibrating K
        elif p_multop[1] == "K":

            temp_K.loc[wellnestlist[0]] = np.array(
                K_data.loc[wellnestlist[0]]) * np.array(
                    np.exp(p[num_wells*n_param+n_pump]))

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + 1):((numt+1)*len_nodes+1)]
                    for numt in range(len(hidstate_t))]

        # If calibrating Sskv and Sske and K
        elif p_multop[1] == "all":

            # BASIN HOPPING
            temp_Sskv.loc[wellnestlist[0]] = np.array(
                Sskv_data.loc[wellnestlist[0]]) * np.array(
                    np.exp(p[num_wells*n_param+n_pump]))

            temp_Sske.iloc[0, ::2] = Sskv_data.loc[wellnestlist[0]].iloc[::2] * \
                np.array(np.exp(p[num_wells*n_param+n_pump+1]))
            temp_Sske.iloc[0, 1::2] = temp_Sske.iloc[0, 0::2] / 10

            temp_K.loc[wellnestlist[0]] = np.array(
                K_data.loc[wellnestlist[0]]) * np.array(
                    np.exp(p[num_wells*n_param+n_pump+2]))

            # Number of clay nodes
            len_nodes = Nz + 2
            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + (
                num_wells*n_param+n_pump+3)):((
                    numt+1)*len_nodes+(
                        num_wells*n_param+n_pump+3))] for numt in range(len(hidstate_t))]

        # If calibrating Sskv and Sske
        elif p_multop[1] == "Ss":

            # BASIN HOPPING
            temp_Sskv.loc[wellnestlist[0]] = np.array(
                Sskv_data.loc[wellnestlist[0]]) * np.array(np.exp(p[0]))

            temp_Sske.iloc[0, ::2] = Sskv_data.loc[wellnestlist[0]].iloc[::2] * \
                np.array(np.exp(p[1]))
            temp_Sske.iloc[0, 1::2] = temp_Sske.iloc[0, 0::2] / 10

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + 2):((numt+1)*len_nodes+2)]
                    for numt in range(len(hidstate_t))]

        # If calibrating Sskv and K
        elif p_multop[1] == "SsK":

            # BASIN HOPPING
            temp_Sskv.loc[wellnestlist[0]] = np.array(
                Sskv_data.loc[wellnestlist[0]]) * np.array(
                    np.exp(p[num_wells*n_param+n_pump]))

            temp_Sske.loc[wellnest][::2] = temp_Sskv.loc[wellnest][::2] * .15
            temp_Sske.loc[wellnest][1::2] = temp_Sske.loc[wellnest][0::2] / 10

            temp_K.loc[wellnestlist[0]] = np.array(
                K_data.loc[wellnestlist[0]]) * np.array(
                    np.exp(p[num_wells*n_param+n_pump+1]))

            # Each hidden state value METHOD
            # Number of clay nodes
            # len_nodes = Nz + 2

            # # Gets initial clay heads
            # h_ic = [p[(numt*len_nodes + (
            #     num_wells*n_param+n_pump+2)):((
            #         numt+1)*len_nodes+(
            #             num_wells*n_param+n_pump+2))] for numt in range(len(hidstate_t))]

            # Multiplier of hidden state METHOD
            # len_clay = len(hidstate)
            # # Alternative method for clay heads
            # # Usign multiplier instead of individual values to
            # # maintain pattern from diffusion model
            # h_ic = []
            # [h_ic.append(pd.Series(hidstate[numhid] * p[(numhid + (
            #     num_wells*n_param+n_pump+2))])) for numhid in range(len_clay)]

    # Initiate an array of predicted results.
    for model_i, model in enumerate(models):

        # parameter names
        param_name = model.parameters.index[params_i].values

        for param_i, param in enumerate(param_name):

            model.set_parameter(name=param,
                                initial=p[model_i*n_param+param_i],
                                optimal=p[model_i*n_param+param_i])

        if n_param != 4:
            # Other names
            other_name = model.parameters.index[others_i].values

            for other_i, other in enumerate(other_name):
                # Set to optimal parameter
                # model.set_parameter(name=other,
                #                     initial=other_params[model_i][others_i[other_i]],
                #                     optimal=other_params[model_i][others_i[other_i]])
                # Set to the truth
                if other == "well_A":

                    model.set_parameter(name=other,
                                        initial=-.1,
                                        optimal=-.1)
                elif other == "constant_d":

                    model.set_parameter(name=other,
                                        initial=2,
                                        optimal=2)

        # Saving new models with new parameters
        newmodels.append(model)

    # Saving pumping data
    pumping_data.append(p[num_wells*n_param:num_wells*n_param+n_pump])

    # Pumping ensemble
    pumping_series = pd.DataFrame(pumping_data)
    pumping_series = pumping_series.transpose()
    pumping_series.index = annual_date_pump.index
    df = pd.DataFrame(index=daily_date_pump.index)
    df = pd.concat([df, pumping_series], join="outer",
                   keys=["Date", "Date"], axis=1)
    df.columns = df.columns.droplevel()

    # Interpolating pumping data
    pump_interp = df.interpolate(method="cubic")
    pump_interp = pump_interp.dropna()

    # Repeat pumping if num_wells > 1
    if num_wells > 1:

        pump_interp = pd.concat(
            [pump_interp] * (num_wells), axis=1, ignore_index=True)

    pump_interp.columns = well_names
    well_data_dates, \
        all_well4_data = load_Pastas_ESMDA(Pastasfiles,
                                           lenfiles,
                                           proxyflag,
                                           newmodels, well_names,
                                           model_path,
                                           pumpflag,
                                           tmin, tmax, SS_data,
                                           init=0, ne=None,
                                           initoptiparam=None,
                                           pump_path=None,
                                           pump_sheet=None,
                                           pump_series=pump_interp
                                           )

    # Calculates subsidence without reloading Pastas models
    sub_total, subv_total, all_results = run_sub(num_clay, all_well4_data,
                                                 well_data_dates, mode,
                                                 tmin, tmax, SS_data,
                                                 wellnest,
                                                 temp_K, temp_Sskv, temp_Sske,
                                                 CC, Nz, Thick_data, ic_run,
                                                 sub_total, subv_total,
                                                 all_results, user_ic=None)
    #                                             [hidstate_t, h_ic])

    # Post process data
    sub_total, subv_total, ann_sub, \
        avgsub = bkk_postproc(wellnestlist,
                              sub_total,
                              subv_total,
                              all_results)

    # preparation
    daterange = pd.date_range(dt.datetime(1978, 12, 31), periods=43,
                              freq="Y").tolist()
    df = pd.DataFrame(daterange, columns=["date"])

    # annual data in cm
    plot_data = df.merge(ann_sub[0][1]*100, left_on=df.date,
                         right_on=ann_sub[0][1].index,
                         how="left")

    # Renaming for other merge
    plot_data = plot_data.rename(columns={"key_0": "key0"})

    # Filling na with 0
    plot_data = plot_data.fillna(0)

    # Benchamrks already in cm
    plot_data = plot_data.merge(dobs.iloc[0:nd_sub].rename("dobs_sub"),
                                left_on=plot_data.key0,
                                right_on=pd.to_datetime(
                                    dobs.iloc[0:nd_sub].index),
                                how="left")
    # Renaming for other merge
    plot_data = plot_data.rename(columns={"key_0": "key1"})

    # Filling na with 0
    plot_data = plot_data.fillna(0)

    plot_data = plot_data.dropna()
    landlevel = plot_data[
        plot_data.columns[
            plot_data.columns.str.contains(
                "dobs_sub")].item()]

    d_pred[0:nd_sub] = plot_data.AnnRates[landlevel != 0]

    # No proxy head
    well_data_temp = all_well4_data[well_names]

    # Newly simulated head
    for well_num in range(len(well_names)):

        d_pred_temp = well_data_temp.iloc[:, well_num]
        d_pred_temp.name = "Simulation"

        # Stripping d_pred to match dobs
        # Matching obs to predicted
        df = pd.DataFrame(gw_obs_indices[well_num+1], columns=["date"],
                          index=gw_obs_indices[well_num+1])
        len_ = len(gw_obs_indices[well_num+1])  # +1 because first is sub
        # annual data in cm

        temp = df.merge(d_pred_temp, left_on=pd.to_datetime(
            df.index), right_on=d_pred_temp.index, how="left")

        # Renaming for second merge
        temp = temp.rename(columns={"key_0": "key0"})

        if well_num == 0:

            index0 = nd_sub
            nd_temp_1 = nd_sub + len_

        else:

            index0 += len(gw_obs_indices[well_num+1-1])
            nd_temp_1 += len_

        d_pred[index0:nd_temp_1] = temp.Simulation

        # Renaming for second merge
        temp = temp.rename(columns={"Simulation": "Simulation0"})

    return d_pred


def forward_all_my(p, wellnestlist, Pastasfiles, lenfiles,
                   proxyflag, model_path, pumpflag,
                   mode, tmin, tmax,
                   Thick_data, K_data, Sskv_data,
                   Sske_data, CC,
                   Nz, num_clay, all_well4_data,
                   well_data_dates,
                   wellnest,
                   SS_data, p_multop,
                   ic_run, return_sub, dobs,
                   gw_obs_indices, models,
                   well_names,
                   n_pump, annual_date_pump,
                   daily_date_pump, nd_sub, nd_gw,
                   initoptiparam=None,
                   pump_path=None,
                   pump_sheet=None,
                   pump_series=None):
    """Running pastas and subsidence model for py esmda

    p - list of param multiplers (1 for Sskv, 2 for Sske, 3 for K)
    wellnestlist - list of wellnest to calculate subsidence for
    mode - raw groundwater data or time series from pastas (raw needs to
    to be interpolated). options: raw, pastas
    tmin, tmax - (str) minimum and maximum year to calculate sub
    CC - convergence criteria
    Nz - number of nodes in the z direction
    ic_run - True or false to generate initial condition run for clays
    proxyflag - 1 if using available heads as proxy for missing heads
    pumpflag - 1 if changing pumping scenario for Pastas
    model_path - path to python models
    pump_path - path to pumping excel sheet
    pump_sheet - sheet of specific pumping scenario
    califlag - 1 if calibrating subsidence model
    dobs - observations

    The data sets have specific names for clays and aquifers
    Thick_data - thickness of clay and aquifers
    K_data - vertical hydraulic conductivity of clay and aquifers
    Sskv - inelastic specific storage term
    Sske - elastic specific storage term

    Returns
    y - subsidence (cm/yr) time series
    """

    # Initiate an array of predicted results.
    d_pred = np.zeros([len(dobs)])

    # Preallocation
    # Head time series for each  node
    all_results = []

    # Subsidence sum for all clay layers
    sub_total = []

    # Inelastic subsidence sum for all clay layers
    subv_total = []

    # Temporary assignments
    temp_Sskv = Sskv_data.copy()
    temp_Sske = Sske_data.copy()
    temp_K = K_data.copy()
    # BASIN HOPPING
    temp_Sskv.loc[wellnestlist[0]] = np.array(
        Sskv_data.loc[wellnestlist[0]]) * np.array(p.iloc[-3])

    temp_Sske.iloc[0, ::2] = Sskv_data.loc[
        wellnestlist[0]].iloc[::2] * \
        np.array(p.iloc[-2])
    temp_Sske.iloc[0, 1::2] = temp_Sske.iloc[0, 0::2] / 10

    temp_K.loc[wellnestlist[0]] = np.array(
        K_data.loc[wellnestlist[0]]) * np.array(p.iloc[-1])

    # Optimal params for running ESMDA
    # Number of wells
    num_well = len(well_names)

    # Optimal params for running ESMDA
    initopti_data = []
    pumping_data = []
    for well_ in range(num_well):
        initopti_data.append(p[well_*4:(well_+1)*4])
        pumping_data.append(p[num_well*4:num_well*4+n_pump])

    initoptiparam = pd.DataFrame(list(map(np.ravel, initopti_data)))
    initoptiparam.index = well_names

    # Pumping ensemble
    pumping_series = pd.DataFrame(pumping_data)
    pumping_series = pumping_series.transpose()
    pumping_series.index = annual_date_pump.index
    df = pd.DataFrame(index=daily_date_pump.index)
    df = pd.concat([df, pumping_series], join="outer",
                   keys=["Date", "Date"], axis=1)
    df.columns = df.columns.droplevel()

    # Interpolating pumping data
    pump_interp = df.interpolate(method="cubic")
    pump_interp = pump_interp.dropna()
    pump_interp.columns = well_names
    well_data_dates, \
        all_well4_data = load_Pastas_ESMDA(Pastasfiles,
                                           lenfiles,
                                           proxyflag,
                                           models, well_names,
                                           model_path,
                                           pumpflag,
                                           tmin, tmax, SS_data,
                                           init=0, ne=None,
                                           initoptiparam=initoptiparam,
                                           pump_path=None,
                                           pump_sheet=None,
                                           pump_series=pump_interp
                                           )

    # Calculates subsidence without reloading Pastas models
    sub_total, subv_total, all_results = run_sub(num_clay, all_well4_data,
                                                 well_data_dates, mode,
                                                 tmin, tmax, SS_data,
                                                 wellnest,
                                                 temp_K, temp_Sskv, temp_Sske,
                                                 CC, Nz, Thick_data, ic_run,
                                                 sub_total, subv_total,
                                                 all_results)

    # Post process data
    sub_total, subv_total, ann_sub, \
        avgsub = bkk_postproc(wellnestlist,
                              sub_total,
                              subv_total,
                              all_results)

    # preparation
    daterange = pd.date_range(dt.datetime(1978, 12, 31), periods=43,
                              freq="Y").tolist()
    df = pd.DataFrame(daterange, columns=["date"])

    # annual data in cm
    plot_data = df.merge(ann_sub[0][1]*100, left_on=df.date,
                         right_on=ann_sub[0][1].index,
                         how="left")

    # Renaming for other merge
    plot_data = plot_data.rename(columns={"key_0": "key0"})

    # Filling na with 0
    plot_data = plot_data.fillna(0)

    # Benchamrks already in cm
    plot_data = plot_data.merge(dobs.iloc[0:nd_sub].rename("dobs_sub"),
                                left_on=plot_data.key0,
                                right_on=pd.to_datetime(
                                    dobs.iloc[0:nd_sub].index),
                                how="left")
    # Renaming for other merge
    plot_data = plot_data.rename(columns={"key_0": "key1"})

    # Filling na with 0
    plot_data = plot_data.fillna(0)

    plot_data = plot_data.dropna()
    landlevel = plot_data[
        plot_data.columns[
            plot_data.columns.str.contains(
                "dobs_sub")].item()]

    d_pred[0:nd_sub] = plot_data.AnnRates[landlevel != 0]

    # No proxy head
    well_data_temp = all_well4_data[well_names]

    # Newly simulated head
    for well_num in range(len(well_names)):

        d_pred_temp = well_data_temp.iloc[:, well_num]
        d_pred_temp.name = "Simulation"

        # Stripping d_pred to match dobs
        # Matching obs to predicted
        df = pd.DataFrame(gw_obs_indices[well_num+1], columns=["date"],
                          index=gw_obs_indices[well_num+1])
        len_ = len(gw_obs_indices[well_num+1])  # +1 because first is sub
        # annual data in cm

        temp = df.merge(d_pred_temp, left_on=pd.to_datetime(
            df.index), right_on=d_pred_temp.index, how="left")

        # Renaming for second merge
        temp = temp.rename(columns={"key_0": "key0"})

        if well_num == 0:

            index0 = nd_sub
            nd_temp_1 = nd_sub + len_

        else:

            index0 += len(gw_obs_indices[well_num+1-1])
            nd_temp_1 += len_

        d_pred[index0:nd_temp_1] = temp.Simulation

        # Renaming for second merge
        temp = temp.rename(columns={"Simulation": "Simulation0"})

    return d_pred


def forwardmy_gwparamSELECT_pump(p, models,
                                 dobs,
                                 gw_obs_indices,
                                 Pastasfiles,
                                 lenfiles, model_path, proxyflag,
                                 well_names, pumpflag, m_ensemble,
                                 n_pump, annual_date_pump,
                                 daily_date_pump, SS_data,
                                 pump_path=None,
                                 pump_sheet=None,
                                 params_i=None,
                                 other_params=None,
                                 others_i=None):
    """Running pastas model for py esmda with only some parameters
    and pumping

    p - list of param multiplers (1 for Sskv, 2 for Sske, 3 for K)
    wellnestlist - list of wellnest to calculate subsidence for
    models - list of pastas model instances
    dobs - observations
    gw_obs_indices - list of groundwater obs index
    well_names - list of well names
    n_pump - number of pumping rates
    annual_date_pump - annual dates for pumping
    daily_date_pump - daily dates for pumping
    params_i - indices of parameters to calibrate
    other_params - other parameter values to not calibrate
    other_i - other index

    The data sets have specific names for clays and aquifers
    Thick_data - thickness of clay and aquifers
    K_data - vertical hydraulic conductivity of clay and aquifers
    Sskv - inelastic specific storage term
    Sske - elastic specific storage term

    Returns
    y - subsidence (cm/yr) time series
    """

    # Initiate an array of predicted results.
    d_pred = np.zeros([len(dobs)])

    # Length of parameters
    n_param = len(params_i)

    # list to dataframe for function
    initopti_data = []
    pumping_data = []

    # Number of wells/models
    num_wells = len(models)

    # Initiate an array of predicted results.
    for model_i, model in enumerate(models):

        pumping_data.append(p[num_wells*n_param:num_wells*n_param+n_pump])

        # parameter names
        param_name = model.parameters.index[params_i].values

        tmin = str(gw_obs_indices[model_i][0].year)
        # Give an extra year since it stops at the start of the year
        tmax = str(gw_obs_indices[model_i][-1].year + 1)

        for param_i, param in enumerate(param_name):

            model.set_parameter(name=param,
                                initial=p[model_i*n_param+param_i],
                                optimal=p[model_i*n_param+param_i])

        if n_param != 4:

            # Other names
            other_name = model.parameters.index[others_i].values

            for other_i, other in enumerate(other_name):

                model.set_parameter(name=other,
                                    initial=other_params[model_i][others_i[other_i]],
                                    optimal=other_params[model_i][others_i[other_i]])
        models[model_i] = model

    initoptiparam = pd.DataFrame(list(map(np.ravel, initopti_data)))
    initoptiparam.index = well_names

    # Pumping ensemble
    pumping_series = pd.DataFrame(pumping_data)
    pumping_series = pumping_series.transpose()
    pumping_series.index = annual_date_pump.index
    df = pd.DataFrame(index=daily_date_pump.index)
    df = pd.concat([df, pumping_series], join="outer",
                   keys=["Date", "Date"], axis=1)
    df.columns = df.columns.droplevel()

    # Interpolating pumping data
    pump_interp = df.interpolate(method="cubic")
    pump_interp = pump_interp.dropna()
    pump_interp.columns = well_names

    well_data_dates, \
        all_well4_data = load_Pastas_ESMDA(Pastasfiles,
                                           lenfiles,
                                           proxyflag,
                                           models, well_names,
                                           model_path,
                                           pumpflag,
                                           tmin, tmax, SS_data,
                                           init=0, ne=m_ensemble,
                                           initoptiparam=None,
                                           pump_path=None,
                                           pump_sheet=None,
                                           pump_series=pump_interp
                                           )
    # No proxy head
    well_data_temp = all_well4_data[well_names]

    # Initiate an array of predicted results.
    for model_i, model in enumerate(models):

        # Newly simulated head
        d_pred_temp = well_data_temp.iloc[:, model_i]
        d_pred_temp.name = "Simulation"

        # Stripping d_pred to match dobs
        # Matching obs to predicted
        df = pd.DataFrame(gw_obs_indices[model_i],
                          columns=["date"], index=gw_obs_indices[model_i])

        len_ = len(gw_obs_indices[model_i])

        # annual data in cm
        plot_data = df.merge(d_pred_temp, left_on=df.index,
                             right_on=d_pred_temp.index,
                             how="left")

        # Renaming for second merge
        plot_data = plot_data.rename(columns={"key_0": "key0"})

        if model_i == 0:

            index0 = 0
            nd_temp_1 = 0 + len_

        else:

            index0 += len(gw_obs_indices[model_i-1])
            nd_temp_1 += len_

        # Renaming for ESMDA code later on
        plot_data = plot_data.rename(columns={"Simulation": "AnnRates"})

        d_pred[index0:nd_temp_1] = plot_data.AnnRates
        # Renaming for second merge
        plot_data = plot_data.rename(columns={"AnnRates ": "Simulation0"})

    return d_pred


def forwardls_gwparamSELECT_pump(p, models,
                                 dobs,
                                 gw_obs_indices,
                                 Pastasfiles,
                                 lenfiles, model_path, proxyflag,
                                 well_names, pumpflag, m_ensemble,
                                 n_pump, annual_date_pump,
                                 daily_date_pump, SS_data, obs_var,
                                 pump_path=None,
                                 pump_sheet=None,
                                 params_i=None,
                                 other_params=None,
                                 others_i=None):
    """Running pastas model for least squares with only some parameters
    and pumping

    p - list of param multiplers (1 for Sskv, 2 for Sske, 3 for K)
    wellnestlist - list of wellnest to calculate subsidence for
    models - list of pastas model instances
    dobs - observations
    gw_obs_indices - list of groundwater obs index
    well_names - list of well names
    n_pump - number of pumping rates
    annual_date_pump - annual dates for pumping
    daily_date_pump - daily dates for pumping
    params_i - indices of parameters to calibrate
    other_params - other parameter values to not calibrate
    other_i - other index

    The data sets have specific names for clays and aquifers
    Thick_data - thickness of clay and aquifers
    K_data - vertical hydraulic conductivity of clay and aquifers
    Sskv - inelastic specific storage term
    Sske - elastic specific storage term

    Returns
    y - subsidence (cm/yr) time series
    """

    v = p.valuesdict()

    # Initiate an array of predicted results.
    d_pred = np.zeros([len(dobs)])

    # Length of parameters
    n_param = len(params_i)

    # list to dataframe for function
    initopti_data = []
    pumping_data = []

    # Initiate an array of predicted results.
    for model_i, model in enumerate(models):

        # parameter names
        param_name = model.parameters.index[params_i].values

        tmin = str(gw_obs_indices[model_i][0].year)
        # Give an extra year since it stops at the start of the year
        tmax = str(gw_obs_indices[model_i][-1].year + 1)

        for param_i, param in enumerate(param_name):

            model.set_parameter(name=param,
                                initial=v[param],
                                optimal=v[param])

        if n_param != 4:

            # Other names
            other_name = model.parameters.index[others_i].values

            for other_i, other in enumerate(other_name):

                model.set_parameter(name=other,
                                    initial=other_params[model_i][others_i[other_i]],
                                    optimal=other_params[model_i][others_i[other_i]])
        models[model_i] = model

    initoptiparam = pd.DataFrame(list(map(np.ravel, initopti_data)))
    initoptiparam.index = well_names

    # Pumping ensemble
    n_pump = len(annual_date_pump)

    for pump_i in range(n_pump):

        pumping_data.append(v["pump"+str(pump_i)])

    pumping_series = pd.DataFrame(pumping_data)
    pumping_series.index = annual_date_pump.index
    df = pd.DataFrame(index=daily_date_pump.index)
    df = pd.concat([df, pumping_series], join="outer",
                   keys=["Date", "Date"], axis=1)
    df.columns = df.columns.droplevel()

    # Interpolating pumping data
    pump_interp = df.interpolate(method="cubic")
    pump_interp = pump_interp.dropna()
    pump_interp.columns = well_names

    well_data_dates, \
        all_well4_data = load_Pastas_ESMDA(Pastasfiles,
                                           lenfiles,
                                           proxyflag,
                                           models, well_names,
                                           model_path,
                                           pumpflag,
                                           tmin, tmax, SS_data,
                                           init=0, ne=m_ensemble,
                                           initoptiparam=None,
                                           pump_path=None,
                                           pump_sheet=None,
                                           pump_series=pump_interp
                                           )
    # No proxy head
    well_data_temp = all_well4_data[well_names]

    # Initiate an array of predicted results.
    for model_i, model in enumerate(models):

        # Newly simulated head
        d_pred_temp = well_data_temp.iloc[:, model_i]
        d_pred_temp.name = "Simulation"

        # Stripping d_pred to match dobs
        # Matching obs to predicted
        df = pd.DataFrame(gw_obs_indices[model_i],
                          columns=["date"], index=gw_obs_indices[model_i])

        len_ = len(gw_obs_indices[model_i])

        # annual data in cm
        plot_data = df.merge(d_pred_temp, left_on=df.index,
                             right_on=d_pred_temp.index,
                             how="left")

        # Renaming for second merge
        plot_data = plot_data.rename(columns={"key_0": "key0"})

        if model_i == 0:

            index0 = 0
            nd_temp_1 = 0 + len_

        else:

            index0 += len(gw_obs_indices[model_i-1])
            nd_temp_1 += len_

        # Renaming for ESMDA code later on
        plot_data = plot_data.rename(columns={"Simulation": "AnnRates"})

        d_pred[index0:nd_temp_1] = plot_data.AnnRates
        # Renaming for second merge
        plot_data = plot_data.rename(columns={"AnnRates ": "Simulation0"})

    # Calculating residuals
    resid = d_pred - dobs

    # Cost function
    # Jx = .5 * resid.T @ np.linalg.pinv(ce) @ resid
    Jx = resid/obs_var

    return Jx


def forward_gwparam_subSELECT_pump_ls(p, wellnestlist, Pastasfiles, lenfiles,
                                      proxyflag, model_path, pumpflag,
                                      mode, tmin, tmax,
                                      Thick_data, K_data, Sskv_data,
                                      Sske_data, CC,
                                      Nz, num_clay, all_well4_data,
                                      well_data_dates,
                                      wellnest,
                                      SS_data, p_multop,
                                      ic_run, return_sub, dobs,
                                      gw_obs_indices, models,
                                      well_names,
                                      n_pump, annual_date_pump,
                                      daily_date_pump, nd_sub, nd_gw, obs_var,
                                      initoptiparam=None,
                                      pump_path=None,
                                      pump_sheet=None,
                                      pump_series=None,
                                      hidstate_t=None, hidstate=None,
                                      params_i=None,
                                      other_params=None,
                                      others_i=None):
    """Running pastas and subsidence model for esmda
    Selected parameters from gw and subsidence + pumping + clay heads


    p - list of param multiplers (1 for Sskv, 2 for Sske, 3 for K)
    wellnestlist - list of wellnest to calculate subsidence for
    mode - raw groundwater data or time series from pastas (raw needs to
    to be interpolated). options: raw, pastas
    tmin, tmax - (str) minimum and maximum year to calculate sub
    CC - convergence criteria
    Nz - number of nodes in the z direction
    ic_run - True or false to generate initial condition run for clays
    proxyflag - 1 if using available heads as proxy for missing heads
    pumpflag - 1 if changing pumping scenario for Pastas
    model_path - path to python models
    pump_path - path to pumping excel sheet
    pump_sheet - sheet of specific pumping scenario
    califlag - 1 if calibrating subsidence model
    dobs - observations
    hidstate_t - initial hidden state time given by user
    hidstate - initial hidden state given by user

    The data sets have specific names for clays and aquifers
    Thick_data - thickness of clay and aquifers
    K_data - vertical hydraulic conductivity of clay and aquifers
    Sskv - inelastic specific storage term
    Sske - elastic specific storage term

    Returns
    y - subsidence (cm/yr) time series
    """

    v = p.valuesdict()

    # Obs covariance
    # [nd x nd]
    # ce = np.diag(obs_var)

    # Initiate an array of predicted results.
    d_pred = np.zeros([len(dobs)])

    # Preallocation
    # Head time series for each  node
    all_results = []

    # Subsidence sum for all clay layers
    sub_total = []

    # Inelastic subsidence sum for all clay layers
    subv_total = []

    # Optimal params for running ESMDA
    # Number of wells
    num_wells = len(models)

    # Length of parameters
    n_param = len(params_i)

    # list to dataframe for pumping
    pumping_data = []

    # Preallocation of new models with new parameters
    newmodels = []

    # Temporary assignments
    temp_Sskv = Sskv_data.copy()
    temp_Sske = Sske_data.copy()
    temp_K = K_data.copy()
    # If working with multipliers
    if p_multop[0]:

        # If only calibrating Sskv and Sske
        if p_multop[1] == "Sskv":

            temp_Sskv.loc[wellnestlist[0]] = np.array(
                Sskv_data.loc[wellnestlist[0]]) * np.array(
                    np.exp(p[num_wells*n_param+n_pump]))

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + (
                num_wells*n_param+n_pump+1)):(
                    (numt+1)*len_nodes+(
                        num_wells*n_param+n_pump+1))] for numt in range(len(hidstate_t))]

        elif p_multop[1] == "Sske":

            temp_Sske.iloc[0, ::2] = Sskv_data.loc[wellnestlist[0]].iloc[::2] * \
                np.array(np.exp(p[num_wells*n_param+n_pump]))
            temp_Sske.iloc[0, 1::2] = temp_Sske.iloc[0, 0::2] / 10

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + 1):((numt+1)*len_nodes+1)]
                    for numt in range(len(hidstate_t))]

        # If calibrating K
        elif p_multop[1] == "K":

            temp_K.loc[wellnestlist[0]] = np.array(
                K_data.loc[wellnestlist[0]]) * np.array(
                    np.exp(p[num_wells*n_param+n_pump]))

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + 1):((numt+1)*len_nodes+1)]
                    for numt in range(len(hidstate_t))]

        # If calibrating Sskv and Sske and K
        elif p_multop[1] == "all":

            # BASIN HOPPING
            temp_Sskv.loc[wellnestlist[0]] = np.array(
                Sskv_data.loc[wellnestlist[0]]) * np.array(
                    np.exp(p[num_wells*n_param+n_pump]))

            temp_Sske.iloc[0, ::2] = Sskv_data.loc[wellnestlist[0]].iloc[::2] * \
                np.array(np.exp(p[num_wells*n_param+n_pump+1]))
            temp_Sske.iloc[0, 1::2] = temp_Sske.iloc[0, 0::2] / 10

            temp_K.loc[wellnestlist[0]] = np.array(
                K_data.loc[wellnestlist[0]]) * np.array(
                    np.exp(p[num_wells*n_param+n_pump+2]))

            # Number of clay nodes
            len_nodes = Nz + 2
            # Gets initial clay heads
            h_ic = [p[(numt*len_nodes + (
                num_wells*n_param+n_pump+3)):((
                    numt+1)*len_nodes+(
                        num_wells*n_param+n_pump+3))] for numt in range(len(hidstate_t))]

        # If calibrating Sskv and Sske
        elif p_multop[1] == "Ss":

            # BASIN HOPPING
            temp_Sskv.loc[wellnestlist[0]] = np.array(
                Sskv_data.loc[wellnestlist[0]]) * np.array(np.exp(p[0]))

            temp_Sske.iloc[0, ::2] = Sskv_data.loc[wellnestlist[0]].iloc[::2] * \
                np.array(np.exp(p[1]))
            temp_Sske.iloc[0, 1::2] = temp_Sske.iloc[0, 0::2] / 10

            # Number of clay nodes
            len_nodes = Nz + 2

            # Gets initial clay heads
            # h_ic = [p[(numt*len_nodes + 2):((numt+1)*len_nodes+2)]
            #         for numt in range(len(hidstate_t))]

        # If calibrating Sskv and K
        elif p_multop[1] == "SsK":

            # BASIN HOPPING
            temp_Sskv.loc[wellnestlist[0]] = np.array(
                Sskv_data.loc[wellnestlist[0]]) * np.array(
                    np.exp(v[p_multop[1]+str(0)]))

            temp_Sske.loc[wellnest][::2] = temp_Sskv.loc[wellnest][::2] * .15
            temp_Sske.loc[wellnest][1::2] = temp_Sske.loc[wellnest][0::2] / 10

            temp_K.loc[wellnestlist[0]] = np.array(
                K_data.loc[wellnestlist[0]]) * np.array(
                    np.exp(v[p_multop[1]+str(1)]))

            # Each hidden state value METHOD
            # Number of clay nodes
            # len_nodes = Nz + 2

            # # Gets initial clay heads
            # h_ic = [p[(numt*len_nodes + (
            #     num_wells*n_param+n_pump+2)):((
            #         numt+1)*len_nodes+(
            #             num_wells*n_param+n_pump+2))] for numt in range(len(hidstate_t))]

            # Multiplier of hidden state METHOD
            # len_clay = len(hidstate)
            # # Alternative method for clay heads
            # # Usign multiplier instead of individual values to
            # # maintain pattern from diffusion model
            # h_ic = []
            # [h_ic.append(pd.Series(hidstate[numhid] * p[(numhid + (
            #     num_wells*n_param+n_pump+2))])) for numhid in range(len_clay)]

    # Initiate an array of predicted results.
    for model_i, model in enumerate(models):

        # parameter names
        param_name = model.parameters.index[params_i].values

        for param_i, param in enumerate(param_name):

            model.set_parameter(name=param,
                                initial=v[param+str(model_i)],
                                optimal=v[param+str(model_i)])

        if n_param != 4:

            # Other names
            other_name = model.parameters.index[others_i].values

            for other_i, other in enumerate(other_name):
                # Sets to optimal param
                # model.set_parameter(name=other,
                #                     initial=other_params[model_i][others_i[other_i]],
                #                     optimal=other_params[model_i][others_i[other_i]])
                # Sets to true value
                # Set to the truth
                if other == "well_A":

                    model.set_parameter(name=other,
                                        initial=-.1,
                                        optimal=-.1)
                elif other == "constant_d":

                    model.set_parameter(name=other,
                                        initial=2,
                                        optimal=2)

        # Saving new models with new parameters
        newmodels.append(model)

    # Pumping ensemble
    n_pump = len(annual_date_pump)

    for pump_i in range(n_pump):

        pumping_data.append(v["pump"+str(pump_i)])

    # Pumping ensemble
    pumping_series = pd.DataFrame(pumping_data)
    pumping_series.index = annual_date_pump.index
    df = pd.DataFrame(index=daily_date_pump.index)
    df = pd.concat([df, pumping_series], join="outer",
                   keys=["Date", "Date"], axis=1)
    df.columns = df.columns.droplevel()

    # Interpolating pumping data
    pump_interp = df.interpolate(method="cubic")
    pump_interp = pump_interp.dropna()

    # Repeat pumping if num_wells > 1
    if num_wells > 1:

        pump_interp = pd.concat(
            [pump_interp] * (num_wells), axis=1, ignore_index=True)

    pump_interp.columns = well_names
    well_data_dates, \
        all_well4_data = load_Pastas_ESMDA(Pastasfiles,
                                           lenfiles,
                                           proxyflag,
                                           newmodels, well_names,
                                           model_path,
                                           pumpflag,
                                           tmin, tmax, SS_data,
                                           init=0, ne=None,
                                           initoptiparam=None,
                                           pump_path=None,
                                           pump_sheet=None,
                                           pump_series=pump_interp
                                           )

    # Calculates subsidence without reloading Pastas models
    sub_total, subv_total, all_results = run_sub(num_clay, all_well4_data,
                                                 well_data_dates, mode,
                                                 tmin, tmax, SS_data,
                                                 wellnest,
                                                 temp_K, temp_Sskv, temp_Sske,
                                                 CC, Nz, Thick_data, ic_run,
                                                 sub_total, subv_total,
                                                 all_results, user_ic=None)
    #                                            [hidstate_t, h_ic])

    # Post process data
    sub_total, subv_total, ann_sub, \
        avgsub = bkk_postproc(wellnestlist,
                              sub_total,
                              subv_total,
                              all_results)

    # preparation
    daterange = pd.date_range(dt.datetime(1978, 12, 31), periods=43,
                              freq="Y").tolist()
    df = pd.DataFrame(daterange, columns=["date"])

    # annual data in cm
    plot_data = df.merge(ann_sub[0][1]*100, left_on=df.date,
                         right_on=ann_sub[0][1].index,
                         how="left")

    # Renaming for other merge
    plot_data = plot_data.rename(columns={"key_0": "key0"})

    # Filling na with 0
    plot_data = plot_data.fillna(0)

    # Benchamrks already in cm
    plot_data = plot_data.merge(dobs.iloc[0:nd_sub].rename("dobs_sub"),
                                left_on=plot_data.key0,
                                right_on=pd.to_datetime(
                                    dobs.iloc[0:nd_sub].index),
                                how="left")
    # Renaming for other merge
    plot_data = plot_data.rename(columns={"key_0": "key1"})

    # Filling na with 0
    plot_data = plot_data.fillna(0)

    plot_data = plot_data.dropna()
    landlevel = plot_data[
        plot_data.columns[
            plot_data.columns.str.contains(
                "dobs_sub")].item()]

    d_pred[0:nd_sub] = plot_data.AnnRates[landlevel != 0]

    # No proxy head
    well_data_temp = all_well4_data[well_names]

    # Newly simulated head
    for well_num in range(len(well_names)):

        d_pred_temp = well_data_temp.iloc[:, well_num]
        d_pred_temp.name = "Simulation"

        # Stripping d_pred to match dobs
        # Matching obs to predicted
        df = pd.DataFrame(gw_obs_indices[well_num+1], columns=["date"],
                          index=gw_obs_indices[well_num+1])
        len_ = len(gw_obs_indices[well_num+1])  # +1 because first is sub
        # annual data in cm

        temp = df.merge(d_pred_temp, left_on=pd.to_datetime(
            df.index), right_on=d_pred_temp.index, how="left")

        # Renaming for second merge
        temp = temp.rename(columns={"key_0": "key0"})

        if well_num == 0:

            index0 = nd_sub
            nd_temp_1 = nd_sub + len_

        else:

            index0 += len(gw_obs_indices[well_num+1-1])
            nd_temp_1 += len_

        d_pred[index0:nd_temp_1] = temp.Simulation

        # Renaming for second merge
        temp = temp.rename(columns={"Simulation": "Simulation0"})

    # Calculating residuals
    # All
    resid = d_pred - dobs

    # # Cost function
    # # Jx = .5 * resid.T @ np.linalg.pinv(ce) @ resid
    Jx = resid/obs_var

    # GW/ SUb separate
    # sub_resid = (d_pred[0:nd_sub] - dobs.iloc[0:nd_sub])/(
    #     obs_var[0:nd_sub] * nd_sub)
    # gw_resid = (d_pred[nd_sub:] - dobs.iloc[nd_sub:])/(
    #     obs_var[nd_sub:] * nd_gw)
    # Jx = pd.concat([sub_resid, gw_resid])

    return Jx


def forwardmy_pump(p, models,
                   dobs,
                   gw_obs_indices,
                   Pastasfiles,
                   lenfiles, model_path, proxyflag,
                   well_names, pumpflag, m_ensemble,
                   n_pump, annual_date_pump,
                   daily_date_pump, SS_data,
                   pump_path=None,
                   pump_sheet=None,
                   ):
    """Running pastas model for py esmda with pumping

    p - list of param multiplers (1 for Sskv, 2 for Sske, 3 for K)
    wellnestlist - list of wellnest to calculate subsidence for
    models - list of pastas model instances
    dobs - observations
    gw_obs_indices - list of groundwater obs index
    well_names - list of well names
    n_pump - number of pumping rates
    annual_date_pump - annual dates for pumping
    daily_date_pump - daily dates for pumping
    params_i - indices of parameters to calibrate
    other_params - other parameter values to not calibrate
    other_i - other index

    The data sets have specific names for clays and aquifers
    Thick_data - thickness of clay and aquifers
    K_data - vertical hydraulic conductivity of clay and aquifers
    Sskv - inelastic specific storage term
    Sske - elastic specific storage term

    Returns
    y - subsidence (cm/yr) time series
    """

    # Initiate an array of predicted results.
    d_pred = np.zeros([len(dobs)])

    # list to dataframe for function
    pumping_data = []

    # Initiate an array of predicted results.
    for model_i, model in enumerate(models):

        pumping_data.append(p)

        tmin = str(gw_obs_indices[model_i][0].year)
        # Give an extra year since it stops at the start of the year
        tmax = str(gw_obs_indices[model_i][-1].year + 1)

    # Pumping ensemble
    pumping_series = pd.DataFrame(pumping_data)
    pumping_series = pumping_series.transpose()
    pumping_series.index = annual_date_pump.index
    df = pd.DataFrame(index=daily_date_pump.index)
    df = pd.concat([df, pumping_series], join="outer",
                   keys=["Date", "Date"], axis=1)
    df.columns = df.columns.droplevel()

    # Interpolating pumping data
    pump_interp = df.interpolate(method="cubic")
    pump_interp = pump_interp.dropna()
    pump_interp.columns = well_names

    well_data_dates, \
        all_well4_data = load_Pastas_ESMDA(Pastasfiles,
                                           lenfiles,
                                           proxyflag,
                                           models, well_names,
                                           model_path,
                                           pumpflag,
                                           tmin, tmax, SS_data,
                                           init=0, ne=m_ensemble,
                                           initoptiparam=None,
                                           pump_path=None,
                                           pump_sheet=None,
                                           pump_series=pump_interp
                                           )
    # No proxy head
    well_data_temp = all_well4_data[well_names]

    # Initiate an array of predicted results.
    for model_i, model in enumerate(models):

        # Newly simulated head
        d_pred_temp = well_data_temp.iloc[:, model_i]
        d_pred_temp.name = "Simulation"

        # Stripping d_pred to match dobs
        # Matching obs to predicted
        df = pd.DataFrame(gw_obs_indices[model_i],
                          columns=["date"], index=gw_obs_indices[model_i])

        len_ = len(gw_obs_indices[model_i])

        # annual data in cm
        plot_data = df.merge(d_pred_temp, left_on=df.index,
                             right_on=d_pred_temp.index,
                             how="left")

        # Renaming for second merge
        plot_data = plot_data.rename(columns={"key_0": "key0"})

        if model_i == 0:

            index0 = 0
            nd_temp_1 = 0 + len_

        else:

            index0 += len(gw_obs_indices[model_i-1])
            nd_temp_1 += len_

        # Renaming for ESMDA code later on
        plot_data = plot_data.rename(columns={"Simulation": "AnnRates"})

        d_pred[index0:nd_temp_1] = plot_data.AnnRates
        # Renaming for second merge
        plot_data = plot_data.rename(columns={"AnnRates ": "Simulation0"})

    return d_pred


def esmda_gw_my(fwd, mprior, m_bounds, alphas, dobs, obs_var, par_var,
                args=(), kwargs={},
                return_steps=False,
                return_cost=False):
    """Simple ES-MDA function for pastas.

    fwd - forward model
    plot_data - observations matched with predicted values
    mprior - prior parameters from initial or previous ensemble member
    m_bounds - boundary of parameters 2 x number of param
    alphas - tweaking parameters in esmda implementation
    dobs - observations
    obs_var - observation variance
    return_cost - returns cost function

    """

    # Check alphas
    if np.abs(np.sum(1/alphas) - 1.0) > 0.001:
        raise ValueError(f"sum(1/alpha) != 1:  {np.sum(1/alphas):.2f}")

    # Get shapes
    nm, ne = mprior.shape
    nd = dobs.size

    # Obs covariance
    # [nd x nd]
    ce = np.diag(obs_var)

    # Par covariance
    if len(par_var) == 1:
        cm = np.diag((np.ones(nm) * par_var) ** 2)

    else:
        cm = np.diag(np.array(par_var) ** 2)

    # Pre-allocate output if return_steps
    if return_steps:
        mprior_out = np.zeros((alphas.size+1, nm, ne))
        dpred_out = np.zeros((alphas.size+1, nd, ne))

    # Preallocation
    # Saving annual subsidence time series from 1978-2020
    # For each ensemble
    dpred_ne = np.empty((nd, ne))
    # Kalman = np.empty((alphas.size, nm, nd))
    # like = np.empty((alphas.size, nd, ne))
    cost = np.empty((alphas.size, ne, ne))

    # Prior Run
    # Loop over ensemble members
    for j in range(ne):
        print("\nPrior Run" +
              "\nEnsemble Member #: " + str(j+1))

        # Run the ensemble from time zero8
        dpred = fwd(mprior.iloc[:, j], *args, **kwargs)

        dpred_ne[:, j] = dpred.copy()

        # Store if required
        if return_steps:
            mprior_out[0, :, :] = mprior.copy()
            dpred_out[0, :, :] = dpred_ne.copy()

    # Loop over alphas
    for i, alpha in enumerate(alphas):

        # Perturb the observation vector
        # Element by element multiplication
        # [nd] + cons * [nd] * [ne x nd] => transposed = [nd x ne]
        duc = np.transpose(dobs.values + np.sqrt(alpha) * np.sqrt(obs_var) *
                           np.random.randn(ne, nd))

        # Cov: using mean of the ensemble member parameters
        # [nd x ne]
        repeat = dpred_ne.mean(axis=1)
        repeat = np.array([repeat] * ne).T
        deltaD = dpred_ne - repeat

        # Cov: using mean of the ensemble member data
        # [nm x ne]
        repeat = mprior.mean(axis=1)
        repeat = np.array([repeat] * ne).T
        deltaM = mprior - repeat

        # Kalman Gain
        # [nm x ne][ne x nd][[nd x ne][ne x nd] + [nd x nd]] => [nm x nd]
        K = (deltaM @ deltaD.T) @ np.linalg.pinv((deltaD @ deltaD.T +
                                                 alpha * (ne-1) * ce))

        # Saving Kalman Gain
        # Kalman[i, :, :] = K

        like_temp = duc-dpred_ne

        # Update ensemble paramters
        # Update mprior; (nm x ne) + (nm x nd) x (nd x ne) = nm x ne (Equation 1)
        mnew = mprior + K @ (like_temp)

        # Saving likelihood
        # like[i, :, :] = duc

        # Clipping according to parameter number
        # Clip the ensemble members to bounds
        for m in range(nm):

            mnew.iloc[m, :] = np.clip(mnew.iloc[m, :],
                                      m_bounds[0][0][m],
                                      m_bounds[1][0][m])

        # Compute and store cost
        if return_cost:

            # Saving cost function
            # costmat = .5 * (mnew - mprior).T @ (
            #     (1/deltaM) * (mnew - mprior)) + .5 * (
            #         (dpred_ne - duc).T @ (alpha*(1/deltaD) * (dpred_ne - duc)))
            # costmat = .5 * (mnew - mprior).T @ np.linalg.pinv(
            #     deltaM @ deltaM.T) @ (mnew - mprior) + \
            #     .5 * (dpred_ne - np.tile(
            #         dobs.T, (ne, 1)).T).T @ alpha*np.linalg.pinv(deltaD @ deltaD.T) @ (
            #             dpred_ne - np.tile(dobs.T, (ne, 1)).T)
            A_term = .5 * (mnew - mprior).T @ np.linalg.pinv(
                cm) @ (mnew - mprior)
            B_term = .5 * (dpred_ne - duc).T @ np.linalg.pinv(alpha * ce) @ (
                dpred_ne - duc)
            costmat = A_term + B_term

            cost[i, :] = costmat

        mprior = mnew.copy()

        # Preallocation
        # Saving annual subsidence time series from 1978-2020
        # For each ensemble
        dpred_ne = np.empty((nd, ne))

        # Loop over ensemble members
        for j in range(ne):

            print("\nESMDA Run #: " + str(i+1) +
                  "\nEnsemble Member #: " + str(j+1))
            # Run the ensemble from time zero
            dpred = fwd(mprior.iloc[:, j], *args, **kwargs)

            dpred_ne[:, j] = dpred

        # Store if required
        if return_steps:
            mprior_out[i+1, :, :] = mprior.copy()
            dpred_out[i+1, :, :] = dpred_ne.copy()

    # Output options
    # return {'dpred': dpred_out, 'mprior': mprior_out,
    #         'K': Kalman, 'like': like}
    # Store if required
    if return_steps and not return_cost:
        return {'dpred': dpred_out, 'mprior': mprior_out}

    # Compute and store chi
    if return_cost:
        return {'dpred': dpred_out, 'mprior': mprior_out,
                'cost': cost}


def esmda_sub_my(fwd, mprior, m_bounds, alphas, dobs, obs_var,
                 par_var,
                 args=(), kwargs={},
                 na_win=1,
                 return_steps=False,
                 return_cost=False):
    """Simple ES-MDA function for sub.

    fwd - forward model
    plot_data - observations matched with predicted values
    mprior - prior parameters from initial or previous ensemble member
    m_bounds - boundary of parameters 2 x number of param
    alphas - tweaking parameters in esmda implementation
    na_win - number of assimilation windows
    dobs - observations
    obs_var - observation variance

    """

    # Check alphas
    if np.abs(np.sum(1/alphas) - 1.0) > 0.001:
        raise ValueError(f"sum(1/alpha) != 1:  {np.sum(1/alphas):.2f}")

    # Get shapes
    nm, ne = mprior.shape
    nd = dobs.size

    # Obs covariance
    # [nd x nd]
    ce = np.diag(obs_var)

    # Par covariance
    cm = np.diag((np.ones(nm) * par_var) ** 2)

    # Pre-allocate output if return_steps
    if return_steps:

        mprior_out = np.zeros((alphas.size+1, nm, ne))
        dpred_out = np.zeros((alphas.size+1, nd, ne))

    # Preallocation
    # Saving annual subsidence time series from 1978-2020
    # For each ensemble
    dpred_ne = np.empty((nd, ne))
    # Kalman = np.empty((alphas.size, nm, nd))
    # like = np.empty((alphas.size, nd, ne))
    cost = np.empty((alphas.size, ne, ne))

    # Window obs
    # Prior Run
    # Loop over ensemble members
    for j in range(ne):
        print("\nPrior Run" +
              "\nEnsemble Member #: " + str(j+1))

        # Run the ensemble from time zero8
        dpred = fwd(mprior.iloc[:, j], *args, **kwargs)

        dpred_ne[:, j] = dpred.copy()

        # Store if required
        if return_steps:

            # If the number of assimilation windows = 1
            if na_win == 1:
                mprior_out[0, :, :] = mprior.copy()
                dpred_out[0, :, :] = dpred_ne.copy()

            else:
                mprior_out[0, :, :] = mprior.copy()
                dpred_out[0, :, :] = dpred_ne.copy()

    # Loop over alphas
    for i, alpha in enumerate(alphas):

        # Perturb the observation vector
        # Element by element multiplication
        # [nd] + cons * [nd] * [ne x nd] => transposed = [nd x ne]
        # duc = np.transpose(dobs.values + np.sqrt(alpha) * np.sqrt(obs_var) *
        #                    np.random.randn(ne, nd))
        duc = np.transpose(dobs.values + np.zeros([ne, nd]))
        # Cov: using mean of the ensemble member parameters
        # [nd x ne]
        repeat = dpred_ne.mean(axis=1)
        repeat = np.array([repeat] * ne).T
        deltaD = dpred_ne - repeat

        # Cov: using mean of the ensemble member data
        # [nm x ne]
        repeat = mprior.mean(axis=1)
        repeat = np.array([repeat] * ne).T
        deltaM = mprior - repeat

        # Kalman Gain
        # [nm x ne][ne x nd][[nd x ne][ne x nd] + [nd x nd]] => [nm x nd]
        K = (deltaM @ deltaD.T) @ np.linalg.pinv((deltaD @ deltaD.T +
                                                 alpha * (ne-1) * ce))

        # Saving Kalman Gain
        # Kalman[i, :, :] = K

        like_temp = duc-dpred_ne

        # Update ensemble paramters
        # Update mprior; (nm x ne) + (nm x nd) x (nd x ne) = nm x ne (Equation 1)
        mnew = mprior + K @ (like_temp)

        # Saving likelihood
        # like[i, :, :] = duc

        # Clipping according to parameter number
        # Clip the ensemble members to bounds
        for m in range(nm):

            mnew.iloc[m, :] = np.clip(mnew.iloc[m, :],
                                      m_bounds[m][0],
                                      m_bounds[m][1])

        # Compute and store cost
        if return_cost:

            # Saving cost function
            # costmat = .5 * (mnew - mprior).T @ (
            #     (1/deltaM) * (mnew - mprior)) + .5 * (
            #         (dpred_ne - duc).T @ (alpha*(1/deltaD) * (dpred_ne - duc)))
            # costmat = .5 * (mnew - mprior).T @ np.linalg.pinv(
            #     deltaM @ deltaM.T) @ (mnew - mprior) + \
            #     .5 * (dpred_ne - np.tile(
            #         dobs.T, (ne, 1)).T).T @ alpha*np.linalg.pinv(deltaD @ deltaD.T) @ (
            #             dpred_ne - np.tile(dobs.T, (ne, 1)).T)
            A_term = .5 * (mnew - mprior).T @ np.linalg.pinv(
                cm) @ (mnew - mprior)
            B_term = .5 * (dpred_ne - duc).T @ np.linalg.pinv(alpha * ce) @ (
                dpred_ne - duc)
            costmat = A_term + B_term

            # cost[i, :] = np.diagonal(costmat)
            cost[i, :] = costmat

        mprior = mnew.copy()

        # Preallocation
        # Saving annual subsidence time series from 1978-2020
        # For each ensemble
        dpred_ne = np.empty((nd, ne))

        # Loop over ensemble members
        for j in range(ne):

            print("\nESMDA Run #: " + str(i+1) +
                  "\nEnsemble Member #: " + str(j+1))
            # Run the ensemble from time zero
            dpred = fwd(mprior.iloc[:, j], *args, **kwargs)

            dpred_ne[:, j] = dpred.copy()

        # Compute and store chi
        if return_cost:

            # Cov: using mean of the ensemble member parameters
            # [nd x ne]
            repeat = dpred_ne.mean(axis=1)
            repeat = np.array([repeat] * ne).T
            # deltagz = dpred_ne - repeat

            # Predicted - obs
            like_temp = duc-dpred_ne

            # chi_sq = (like_temp).T @ np.linalg.pinv((ce + deltagz @ deltagz.T)) @ \
            #     (like_temp)

        # Store if required
        if return_steps:

            mprior_out[i+1, :, :] = mprior.copy()
            dpred_out[i+1, :, :] = dpred_ne.copy()

    # Output options
    # return {'dpred': dpred_out, 'mprior': mprior_out,
    #         'K': Kalman, 'like': like}
    # Store if required
    if return_steps and not return_cost:
        return {'dpred': dpred_out, 'mprior': mprior_out}

    # Compute and store cost function
    if return_cost:
        return {'dpred': dpred_out, 'mprior': mprior_out,
                'cost': cost}


# Running particle filter
def PF(fwd, mprior, m_bounds, dobs, like_sigma,
       args=(), kwargs={},
       return_steps=False,
       return_cost=False):
    """

    Parameters
    ----------
    fwd : model
        forward model generating predictions
    mprior : dataframe
        particle and parameters to estimate
    m_bounds : list of zipped min and max
        boundaries of parameters
    alphas : TYPE
        DESCRIPTION.
    dobs : dataframe
        observations
    like_sigma : int
        standard deviation of likelihood
    return_steps : true or false

    Returns
    -------
    weights
        weights of each particle

    prior_weighted
        weighted groundwater head

    prior
        prior parameters used

    """
    # Get shapes
    nm, ne = mprior.shape
    nd = dobs.size

    # Pre-allocate output if return_steps
    if return_steps:
        weights = np.zeros((1, ne))
        prior_weighted = np.zeros((nd, ne))

    # Preallocation
    # Saving annual subsidence time series from 1978-2020
    # For each ensemble
    dpred_ne = np.empty((nd, ne))

    # Prior Run
    # Loop over ensemble members
    for j in range(ne):
        print("\nPF Run" +
              "\nEnsemble Member #: " + str(j+1))

        # Run the ensemble from time zero8
        dpred = fwd(mprior.iloc[:, j], *args, **kwargs)

        dpred_ne[:, j] = dpred.copy()

    # Calculating the innovation/likelihood
    # Squared differences between data and prior for each particle
    # Repeat dobs to be the size of nd x ne
    dobs = np.array([dobs] * ne).T
    innov = (dpred_ne - dobs) ** 2
    innov = np.mean(innov, axis=1)

    # Calculation of likelihood
    likelihood = 1 / (1 + (innov/like_sigma**2))

    # Calculation of the weights
    # The weights are the product of each `innov` per particle
    # - Shape innov  : (n_obs, n_part)
    # - Shape weights : (n_part, )
    weights = np.prod(likelihood, axis=0)

    # Normalize weights by total weight
    weights /= np.sum(weights)

    # Calculation of the weighted deformation
    # The weighted deformation is the sum of prior * weights,
    # which is equivalent to the dot product.
    prior_weighted = np.dot(dpred_ne, weights)

    # Output options
    # return {'dpred': dpred_out, 'mprior': mprior_out,
    #         'K': Kalman, 'like': like}
    # Store if required
    if return_steps:
        return {'weights': weights,
                'prior_weighted': prior_weighted,
                'mprior': mprior,
                'likelihood': likelihood}


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

        if option[0] == "normal":

            # For each t
            for t in range(n_pump):

                temp = np.random.normal(ann_pump.iloc[t, 1],
                                        ann_pump.iloc[t, 2])

                # Make sure temp is > 0
                if temp < 0:

                    temp_list.append(0)

                else:

                    temp_list.append(temp)

        mat = pd.DataFrame(temp_list, index=ann_pump.index,
                           columns=[i])
        df = pd.concat([df, mat], join="outer",
                       keys=["Date", "Date"], axis=1)
        df.columns = df.columns.droplevel()

    return df


# %%###########################################################################
# Runs Pastas and subsidence models and saves data
##############################################################################

# Assuming has data for all four aquifers
# Assuming conceptual model of clay above BK, between BK and PD, PD and NL, NL
# and NB for a total of 4 clay layers.
def bkk_subsidence(wellnestlist, mode, tmin, tmax,
                   Thick_data, K_data, Sskv_data, Sske_data, CC, Nz, ic_run,
                   proxyflag, pumpflag, califlag, esmdaflag, esmdaindex=None,
                   p_multop=None, return_sub=None, na=None, na_win=None, ne=None,
                   obs_error=None, par_error=None, like_sigma=None,
                   model_path=None, pump_path=None, pump_series=None,
                   pump_sheet=None, initoptiparam=None, user_obs=None,
                   user_obs_indices=None,
                   dist="norm", pump_ens=None,
                   annual_pump=None, listdaily_pump=None,
                   user_models=None, init_state=None):
    """Calculate sub for four clay layers and four confined aquifers.

    wellnestlist - list of wellnest to calculate subsidence for
    mode - raw groundwater data or time series from pastas (raw needs to
    to be interpolated). options: raw, pastas
    tmin, tmax - (str) minimum and maximum year to calculate sub
    CC - convergence criteria
    Nz - number of nodes in the z direction
    ic_run - True or false to generate initial condition run for clays
    proxyflag - 1 if using available heads as proxy for missing heads
    pumpflag - 1 if changing pumping scenario for Pastas
    model_path - path to python models
    pump_path - path to pumping excel sheet
    pumping_series - df of pumping rates
    pump_sheet - sheet of specific pumping scenario
    califlag - 1 if calibrating subsidence model
    esmdaflag - 1 if calibrating subsidence model with esmda
    na - number of data assimilation steps
    ne - number of ensemble members
    obs_error - observation stand deviation
    par_error - parameter stand deviation
    p_multop - [paramaeter multiplier option (True, False), and mode
                either "Sskv" or "K"]
    return_sub - return subsidence True or False
    initoptimparam - optimal parameter of Pastas model provided
    user_dobs - provided obs (may be synthetic)
    user_obs_indicies - indicies of obs
    dist - normal distriubtion intial
    user_models - user providing pastas models
    listdaily_pump - daily dates for pumping
    annual_pump - annual dates for pumping
    pumping_ens - df with ens of pumping time series
    like_sigma - SD of likelihood dist
    init_state - user provided initial state for state estimation (DA)

    The data sets have specific names for clays and aquifers
    Thick_data - thickness of clay and aquifers
    K_data - vertical hydraulic conductivity of clay and aquifers
    Sskv - inelastic specific storage term
    Sske - elastic specific storage term
    Returns
    all_total - list of lists: all subsidence data (total and inelastic) for
    # each clay layer
    sub_total - list of lists: sub total for all four clay layers (m)
    subv_total - list of lists: inelastic sub total for all four clay layers
    # (m)
    """
    # Preallocation
    # Head time series for each  node
    all_results = []

    # Subsidence sum for all clay layers
    sub_total = []

    # Inelastic subsidence sum for all clay layers
    subv_total = []

    # CORRECTING GW HEAD DATA TO LAND SURFACE (COASTAL DEM 2.1)
    landsurf_path = os.path.join(os.path.abspath("inputs"),
                                 "LandSurfElev_GWWellLocs.xlsx")

    # Each well nest has its own Ss and K sheet
    landsurf_data = pd.read_excel(landsurf_path,
                                  sheet_name="2.1",
                                  usecols="C:F",
                                  index_col=0)

    # If running transient simulation before model run
    # to get clay heads to where they need to be
    if ic_run:

        SS_path = os.path.join(os.path.abspath("inputs"),
                               "SS_Head_GWWellLocs.xlsx")

        # Each well nest has its own Ss and K sheet
        SS_data = pd.read_excel(SS_path,
                                sheet_name="SS_Py",
                                index_col=0)

    # For each well nest in the list
    for wellnest in wellnestlist:

        # If calculating subsidence from raw groundwater data
        if mode == "raw":

            # Preprocesses wellnest groundwater data and returns dataframe
            # with data of matching dates after interpolation
            # Returns data within tmin and tmax in well_data, and data not
            # within tmin and tmax in all_well_data
            well_data, all_well_data = bkk_wellnest_preproc(wellnest,
                                                            tmin, tmax,
                                                            proxyflag)

            # Correcting obs GW to land surface
            well_data += (landsurf_data.RASTERVALU.loc[wellnest])
            all_well_data += (landsurf_data.RASTERVALU.loc[wellnest])

            # Number clay layers
            num_clay = len(well_data.columns)

        elif mode == "Pastas":

            # Get Pastas model file names for each wellnest (Should have four
            # files for each aquifer)
            Pastasfiles = [filename for filename in os.listdir(model_path)
                           if filename.startswith(wellnest) &
                           filename.endswith(".pas")]

            # Reordering from shallowest to deepest aquifer
            # Reorder well list to shallow to deep aquifers
            # BK, PD, NL, NB
            Pastasfiles = [x for y in ["_BK", "_PD", "_NL", "_NB"]
                           for x in Pastasfiles if y in x]
            lenfiles = len(Pastasfiles)

            # If no models provided
            if not user_models:
                # Loading models for good
                # models, well names, optimized pastas parameters
                models, well_names, pastas_optparam = load_Pastas_models(
                    Pastasfiles, model_path, SS_data)

            else:

                # Loading models for good
                # models, well names, optimized pastas parameters
                _, well_names, pastas_optparam = load_Pastas_models(
                    Pastasfiles, model_path, SS_data)
                models = user_models

            # If not calibrating
            if califlag == 0:
                well_data_dates, \
                    all_well4_data = load_Pastas(Pastasfiles,
                                                 lenfiles,
                                                 proxyflag, models,
                                                 well_names,
                                                 model_path,
                                                 pumpflag,
                                                 tmin, tmax,
                                                 initoptiparam=initoptiparam,
                                                 pump_path=pump_path,
                                                 pump_sheet=pump_sheet,
                                                 pump_series=pump_series
                                                 )

            # If calibrating
            else:

                if len(par_error) > 4:
                    # Taking parameter confidence of gw
                    par_conf_gw = par_error[3:]
                    std_mult = par_conf_gw

                else:

                    std_mult = par_error
                # Number of ensemble members
                ne_pump = ne
                well_data_dates, all_well4_data, all_obs_data, min_pastas, \
                    max_pastas, param_ens, well_names, models = load_Pastas_ESMDA(
                        Pastasfiles,
                        lenfiles,
                        proxyflag,
                        models,
                        well_names,
                        model_path,
                        pumpflag,
                        tmin, tmax,
                        SS_data, init=1, std_mult=std_mult,
                        ne=ne_pump,
                        initoptiparam=initoptiparam,
                        pump_path=pump_path,
                        pump_sheet=pump_sheet,
                        pump_series=pump_series,
                        dist=dist, mode=esmdaflag
                    )

            num_clay = 4

        # If calibrating
        if califlag == 1:

            # If calibrating mroe than one well nest, won't work
            if len(wellnestlist) > 1:

                sys.exit("Calibration must be done with one well nest at a time.")

            # LOADING OBS AND BOUNDS DATA
            # For each wellnest in list
            # num_well is the index, wellnest = name
            wellnest = wellnestlist[0]

            # If user doesn't provide obs
            if user_obs is None:
                dobs = pd.Series(np.empty(1, dtype=object))
                # Saving gw obs index
                gw_obs_indices = []
                # Adding subsidence observations to groundwater observations
                for welli in range(len(all_obs_data.columns)):
                    data = all_obs_data.iloc[:, welli].dropna()
                    gw_obs_indices.append(data.index)
                    dobs = pd.concat([dobs, data])
                dobs = dobs[1:]

            # If user provides obs
            else:

                dobs = user_obs
                gw_obs_indices = user_obs_indices

            # Generate alphas
            # Alphas
            cov_obs_inflation_geo = 1.2
            cov_obs_inflation_factors = [1.1]  # list[float]
            for i in range(1, na):

                cov_obs_inflation_factors.append(
                    cov_obs_inflation_factors[i - 1] / cov_obs_inflation_geo
                )

            scaling_factor = np.sum(1 /
                                    np.array(cov_obs_inflation_factors))
            # : float
            cov_obs_inflation_factors = [
                alpha * scaling_factor
                for alpha in cov_obs_inflation_factors
            ]
            alphas = np.array(cov_obs_inflation_factors)

            # If pastas parameter calibration with my esmda func
            if "my_gwparamSELECT" == esmdaflag:

                nd_sub = 0

                # Mprior
                mprior_init = pd.DataFrame()

                # Parameter indexing
                param_index = esmdaindex

                # Initial samples of parameters
                # should be [nd x ne]
                for param_i in param_ens:

                    mprior = param_i.iloc[:, param_index].T
                    mprior = pd.DataFrame(mprior)
                    mprior_init = pd.concat([mprior_init.reset_index(drop=True),
                                             mprior.reset_index(drop=True)],
                                            axis=0,
                                            ignore_index=True)

                # Other index
                range_ = set(range(0, len(pastas_optparam[0])))
                other_i = np.array(list(
                    set(param_index).symmetric_difference(range_)))

                # Number of observations
                nd = dobs.size
                nd_gw = nd - nd_sub
                # Associated standard deviation: ones (for this scenario)
                obs_error_sub = 1.5
                obs_error_gw = obs_error
                obs_std = np.append(np.ones(nd_sub) * obs_error_sub,
                                    np.ones(nd_gw) * obs_error_gw)
                obs_var = obs_std**2

                # Boundaries of parameters
                mins = np.append(min_pastas[:, param_index], np.zeros(1))
                maxs = np.append(max_pastas[:, param_index], np.zeros(1))
                mins = mins[:-1]
                maxs = maxs[:-1]
                m_bounds = list(zip([mins, maxs]))

                # Running ESMDA
                out = esmda_gw_my(forwardmy_gwparamSELECT, mprior_init, m_bounds,
                                  alphas, dobs, obs_var,
                                  args=(models, dobs, gw_obs_indices, None, None,
                                        param_index, pastas_optparam,
                                        other_i),
                                  return_steps=True,
                                  return_cost=True)

                # Returns ESMDA solver and pastas models
                return out, models, well_names

            # If parameter calibration and pumping estimation
            elif "my_gwparamSELECT_pump" == esmdaflag:

                nd_sub = 0

                # Mprior
                mprior_init = pd.DataFrame()

                # Parameter indexing
                param_index = esmdaindex

                # Renaming pumping ensemble
                pumping_ens = pump_ens

                # Initial samples of parameters
                # should be [nd x ne]
                for param_i in param_ens:

                    mprior = param_i.iloc[:, param_index].T
                    mprior = pd.DataFrame(mprior)
                    mprior_init = pd.concat([mprior_init.reset_index(drop=True),
                                             mprior.reset_index(drop=True)],
                                            axis=0,
                                            ignore_index=True)

                # Other index
                range_ = set(range(0, len(pastas_optparam[0])))
                other_i = np.array(list(
                    set(param_index).symmetric_difference(range_)))

                mprior_init = pd.concat([mprior_init.reset_index(drop=True),
                                         pumping_ens.reset_index(drop=True)],
                                        axis=0,
                                        ignore_index=True)

                # Bounds; Pastas, pumping; sub
                n_pump = len(annual_pump)

                # Saving confidence of pumping
                par_pump = np.ones(n_pump) * 30

                # Confidence of gw
                # Taking parameter confidence of gw
                par_conf_gw = np.array(par_error)

                # Adding pastas param and pumping
                if len(models) == 1:
                    # Saving all confidence
                    par_error = par_conf_gw.tolist() + par_pump.tolist()
                    mins = np.append(min_pastas[0, param_index], np.zeros(n_pump))

                else:
                    par_error = par_conf_gw.tolist()
                    # Saving all confidence
                    for n_models in range(len(models)-1):
                        par_error += par_conf_gw.tolist()

                    mins = np.append(min_pastas[0, param_index],
                                     min_pastas[0, param_index])
                    for n_models in range(len(models)-2):

                        mins = np.append(mins, min_pastas[0, param_index])
                    mins = np.append(mins, np.zeros(n_pump))

                mins[np.isnan(mins)] = -30000
                # Adding hidden state
                # mins = np.append(mins, state_mins)
                # Adding pastas param and pumping
                if len(models) == 1:
                    maxs = np.append(max_pastas[0, param_index], np.repeat(5E2, n_pump))
                else:
                    maxs = np.append(max_pastas[0, param_index],
                                     max_pastas[0, param_index])
                    for n_models in range(len(models)-2):

                        maxs = np.append(maxs, max_pastas[0, param_index])
                    maxs = np.append(maxs, np.repeat(5E2, n_pump))

                maxs[np.isnan(maxs)] = 10000

                # Adding hidden state
                # maxs = np.append(maxs, state_maxs)
                m_bounds = list(zip((mins, maxs)))

                # Number of observations
                nd_gw = dobs.size
                nd = nd_gw
                # Associated standard deviation: ones (for this scenario)
                obs_error_gw = obs_error
                obs_std = np.ones(nd_gw) * obs_error_gw
                obs_var = obs_std**2

                # Running ESMDA
                out = esmda_gw_my(forwardmy_gwparamSELECT_pump, mprior_init, m_bounds,
                                  alphas, dobs, obs_var, par_error,
                                  args=(models, dobs, gw_obs_indices, Pastasfiles,
                                        lenfiles, model_path, proxyflag,
                                        well_names, pumpflag, ne, n_pump,
                                        annual_pump, listdaily_pump, SS_data,
                                        None, None,
                                        param_index, pastas_optparam,
                                        other_i),
                                  return_steps=True,
                                  return_cost=False)

                # Returns ESMDA solver and pastas models
                return out, models, well_names

            # If parameter calibration and pumping estimation Least squares
            elif "ls_gwparamSELECT_pump" == esmdaflag:

                # Mprior
                mprior_init = pd.DataFrame()

                # Parameter indexing
                param_index = esmdaindex

                # Renaming pumping ensemble
                pumping_ens = pump_ens

                params = lmfit.Parameters()

                # Initial samples of Pastas parameters
                # should be [nd x ne]
                for param_i in param_ens:

                    for i, pasta_param in enumerate(param_i.columns[param_index]):

                        # If only one well
                        if len(models) == 1:
                            params.add(
                                name=pasta_param,
                                value=param_i[pasta_param].values[0],
                                min=min_pastas[0, param_index[i]],
                                max=max_pastas[0, param_index[i]])

                        # If more models
                        else:

                            # Adding parameters for each model
                            for n_models in range(len(models)-1):
                                params.add(
                                    name=pasta_param,
                                    value=param_i[pasta_param].values[0],
                                    min=min_pastas[0, param_index[i]],
                                    max=max_pastas[0, param_index[i]])

                # Other index
                range_ = set(range(0, len(pastas_optparam[0])))
                other_i = np.array(list(
                    set(param_index).symmetric_difference(range_)))

                # Number of pumping
                n_pump = len(annual_pump)
                # For each pumping ens
                for pump_i in range(n_pump):
                    params.add(
                        name="pump"+str(pump_i),
                        value=pumping_ens.iloc[pump_i, 0],
                        min=0,
                        max=500)

                # Number of observations
                nd = dobs.size
                nd_gw = nd

                # Associated standard deviation: ones (for this scenario)
                obs_error_gw = obs_error
                obs_std = np.ones(nd_gw) * obs_error_gw
                obs_var = obs_std**2

                # Running ESMDA
                out = lmfit.minimize(forwardls_gwparamSELECT_pump, params,
                                     method='least_squares',
                                     args=(models, dobs, gw_obs_indices, Pastasfiles,
                                           lenfiles, model_path, proxyflag,
                                           well_names, pumpflag, ne, n_pump,
                                           annual_pump, listdaily_pump, SS_data, obs_var,
                                           None, None,
                                           param_index, pastas_optparam,
                                           other_i))

                # Returns ESMDA solver and pastas models
                return out, models, well_names

            # If parameter calibration and pumping estimation
            elif "my_gwparam_subSELECT_pump" == esmdaflag:

                # Saves parameter std
                # par_conf = par_error.copy()

                # Mprior
                mprior_init = pd.DataFrame()

                # Parameter indexing
                param_index = esmdaindex

                # Renaming pumping ensemble
                pumping_ens = pump_ens

                # Initial samples of parameters
                # should be [nm x ne]
                for param_i in param_ens:

                    mprior = param_i.iloc[:, param_index].T
                    mprior = pd.DataFrame(mprior)
                    mprior_init = pd.concat([mprior_init.reset_index(drop=True),
                                             mprior.reset_index(drop=True)],
                                            axis=0,
                                            ignore_index=True)

                # Other index
                range_ = set(range(0, len(pastas_optparam[0])))
                other_i = np.array(list(
                    set(param_index).symmetric_difference(range_)))

                mprior_init = pd.concat([mprior_init.reset_index(drop=True),
                                         pumping_ens.reset_index(drop=True)],
                                        axis=0,
                                        ignore_index=True)

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

                # Only subsidence parameters
                if p_multop[1] == "all":
                    nm = 3
                    p_mult = []
                    [p_mult.append(random.uniform(m_bounds[x][0], m_bounds[x][1]))
                     for x in range(nm)]
                    par_error_sub = par_error[0]
                # Only subsidence parameters
                elif p_multop[1] == "Ss":
                    nm = 2
                    p_mult = []
                    [p_mult.append(random.uniform(m_bounds[x][0], m_bounds[x][1]))
                     for x in range(nm)]

                    par_error_sub = par_error[:2]

                elif p_multop[1] == "Sskv":

                    p_mult = [random.uniform(
                        m_bounds[0][0], m_bounds[0][1])]

                    par_error_sub = par_error[0]/10

                elif p_multop[1] == "Sske":

                    p_mult = [random.uniform(
                        m_bounds[1][0], m_bounds[1][1])]

                    par_error_sub = par_error[1]

                elif p_multop[1] == "K":

                    p_mult = [random.uniform(
                        m_bounds[2][0], m_bounds[2][1])]

                    par_error_sub = par_error[2]

                # Only Sskv and K
                elif p_multop[1] == "SsK":
                    nm = 2
                    p_mult = []
                    [p_mult.append(random.uniform(m_bounds[x][0], m_bounds[x][1]))
                     for x in np.arange(0, nm+1, nm)]

                    par_error_sub = par_error[0:nm+1:nm]

                # Number of model parameters
                nm = len(p_mult)

                # Prior: Let's start average of the bounds
                mprior_mean = p_mult

                mprior_mean = np.array(mprior_mean)

                # mprior_mean = np.ones(nm) * random.uniform(m_bounds[x][0], m_bounds[x][1])

                # Parameter error
                par_std = par_error_sub

                # Prior distribution of parameters
                param_init = rng.normal(loc=mprior_mean, scale=par_std,
                                        size=(ne, nm)).T

                param_init = pd.DataFrame(param_init)

                # Adding sub parameters to other initial parameters
                mprior_init = pd.concat([mprior_init.reset_index(drop=True),
                                         param_init],
                                        axis=0,
                                        ignore_index=True)

                # Preallocation for hidden state mins maxs
                # state_mins = []
                # state_maxs = []

                # # Prior distribution of hidden state
                # par_hid = []
                # for hid in init_state[1]:

                #     # Using each clay node individually METHOD
                #     # For each clay node, a new value
                #     state_init = rng.normal(loc=hid, scale=0.5,
                #                             size=(ne, len(hid))).T

                #     # Adding to prior
                #     state_init = pd.DataFrame(state_init)
                #     mprior_init = pd.concat([mprior_init, state_init], axis=0,
                #                             ignore_index=True)
                #     # State boundaries
                #     state_min = np.empty(len(hid))
                #     state_min[:] = np.nan
                #     state_mins.append(state_min)
                #     state_max = np.empty(len(hid))
                #     state_max[:] = np.nan
                #     state_maxs.append(state_max)

                #     # Saving confidence
                #     par_hid.extend(np.ones(len(hid)) * 30)

                #     # Using ERROR vs entire clay METHOD
                #     state_init = rng.normal(loc=1, scale=0.25,
                #                             size=(ne, 1)).T
                #     # Adding to prior
                #     state_init = pd.DataFrame(state_init)
                #     mprior_init = pd.concat([mprior_init, state_init], axis=0,
                #                             ignore_index=True)

                #     # State boundaries
                #     state_min = np.empty(1)
                #     state_min[:] = 0
                #     state_mins.append(state_min)
                #     state_max = np.empty(1)
                #     state_max[:] = 10
                #     state_maxs.append(state_max)

                #     # Saving confidence
                #     par_hid.extend([30])

                # Bounds; Pastas, pumping; sub
                n_pump = len(annual_pump)

                # Saving confidence of pumping
                par_pump = np.ones(n_pump) * 30

                # Confidence of gw
                # Taking parameter confidence of gw
                par_conf_gw = np.array(par_error)[np.array(param_index+3)]

                # Adding pastas param and pumping
                if len(models) == 1:
                    # Saving all confidence
                    par_error = par_conf_gw.tolist() + par_pump.tolist() + \
                        par_error_sub
                    mins = np.append(min_pastas[0, param_index], np.zeros(n_pump))

                else:
                    par_error = par_conf_gw.tolist()
                    # Saving all confidence
                    for n_models in range(len(models)-1):
                        par_error += par_conf_gw.tolist()

                    par_error += par_pump.tolist() + par_error_sub

                    mins = np.append(min_pastas[0, param_index],
                                     min_pastas[0, param_index])
                    for n_models in range(len(models)-2):

                        mins = np.append(mins, min_pastas[0, param_index])
                    mins = np.append(mins, np.zeros(n_pump))

                # Adding sub param
                mins = np.append(mins, parambound.loc[
                    parambound.Wellnest == wellnestlist[0]].iloc[0:nm, 0])
                mins[np.isnan(mins)] = -30000
                # Adding hidden state
                # mins = np.append(mins, state_mins)
                # Adding pastas param and pumping
                if len(models) == 1:
                    maxs = np.append(max_pastas[0, param_index], np.repeat(5E2, n_pump))
                else:
                    maxs = np.append(max_pastas[0, param_index],
                                     max_pastas[0, param_index])
                    for n_models in range(len(models)-2):

                        maxs = np.append(maxs, max_pastas[0, param_index])
                    maxs = np.append(maxs, np.repeat(5E2, n_pump))

                # Adding sub param
                maxs = np.append(maxs, parambound.loc[
                    parambound.Wellnest == wellnestlist[0]].iloc[0:nm, 1])
                maxs[np.isnan(maxs)] = 10000

                # Adding hidden state
                # maxs = np.append(maxs, state_maxs)
                m_bounds = list(zip((mins, maxs)))

                # Number of observations
                nd = dobs.size
                nd_sub = len(gw_obs_indices[0])
                nd_gw = nd - nd_sub
                # Associated standard deviation: ones (for this scenario)
                if len(models) > 1:
                    # multiplier = .5  # Synthetic
                    multiplier = 1  # Real case
                    obs_error_sub = 1.5 * multiplier

                else:
                    obs_error_sub = 1.5
                obs_error_gw = obs_error
                obs_std = np.append(np.ones(nd_sub) * obs_error_sub,
                                    np.ones(nd_gw) * obs_error_gw)
                obs_var = obs_std**2

                # Running ESMDA
                out = esmda_gw_my(forward_gwparam_subSELECT_pump_my,
                                  mprior_init, m_bounds,
                                  alphas, dobs, obs_var, par_error,
                                  args=(wellnestlist, Pastasfiles,
                                        lenfiles, proxyflag,
                                        model_path, pumpflag,
                                        mode, tmin, tmax,
                                        Thick_data, K_data, Sskv_data,
                                        Sske_data, CC,
                                        Nz, num_clay, all_well4_data,
                                        well_data_dates,
                                        wellnest,
                                        SS_data, p_multop,
                                        ic_run, return_sub, dobs,
                                        gw_obs_indices, models,
                                        well_names, n_pump,
                                        annual_pump,
                                        listdaily_pump,
                                        nd_sub, nd_gw, None, None, None, None,
                                        init_state[0], init_state[1],
                                        param_index,
                                        pastas_optparam,
                                        other_i),
                                  return_steps=True,
                                  return_cost=True)

                # Returns ESMDA solver and pastas models
                return out, models, well_names

            # If parameter calibration and pumping estimation Least squares
            elif "ls_gwparam_subSELECT_pump" == esmdaflag:

                # Mprior
                mprior_init = pd.DataFrame()

                # Parameter indexing
                param_index = esmdaindex

                # Renaming pumping ensemble
                pumping_ens = pump_ens

                params = lmfit.Parameters()

                # Initial samples of Pastas parameters
                # should be [nd x ne]
                for param_i in param_ens:

                    for i, pasta_param in enumerate(param_i.columns[param_index]):

                        # Adding parameters for each model
                        for n_models in range(len(models)):
                            params.add(
                                name=pasta_param+str(n_models),
                                value=param_i[pasta_param].values[0],
                                min=min_pastas[0, param_index[i]],
                                max=max_pastas[0, param_index[i]])

                # Other index
                range_ = set(range(0, len(pastas_optparam[0])))
                other_i = np.array(list(
                    set(param_index).symmetric_difference(range_)))

                # Number of pumping
                n_pump = len(annual_pump)
                # For each pumping ens
                for pump_i in range(n_pump):
                    params.add(
                        name="pump"+str(pump_i),
                        value=pumping_ens.iloc[pump_i, 0],
                        min=0,
                        max=500)

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

                # Only subsidence parameters
                if p_multop[1] == "all":
                    nm = 3
                    p_mult = []
                    [p_mult.append(random.uniform(m_bounds[x][0], m_bounds[x][1]))
                     for x in range(nm)]
                    par_error_sub = par_error[0]
                # Only subsidence parameters
                elif p_multop[1] == "Ss":
                    nm = 2
                    p_mult = []
                    [p_mult.append(random.uniform(m_bounds[x][0], m_bounds[x][1]))
                     for x in range(nm)]

                    par_error_sub = par_error[:2]

                elif p_multop[1] == "Sskv":

                    p_mult = [random.uniform(
                        m_bounds[0][0], m_bounds[0][1])]

                    par_error_sub = par_error[0]/10

                elif p_multop[1] == "Sske":

                    p_mult = [random.uniform(
                        m_bounds[1][0], m_bounds[1][1])]

                    par_error_sub = par_error[1]

                elif p_multop[1] == "K":

                    p_mult = [random.uniform(
                        m_bounds[2][0], m_bounds[2][1])]

                    par_error_sub = par_error[2]

                # Only Sskv and K
                elif p_multop[1] == "SsK":
                    nm = 2
                    p_mult = []
                    [p_mult.append(random.uniform(m_bounds[x][0], m_bounds[x][1]))
                     for x in np.arange(0, nm+1, nm)]

                    par_error_sub = par_error[0:nm+1:nm]

                # Number of model parameters
                nm = len(p_mult)

                # Prior: Let's start average of the bounds
                mprior_mean = p_mult

                mprior_mean = np.array(mprior_mean)

                # Parameter error
                par_std = par_error_sub

                # Prior distribution of parameters
                param_init = rng.normal(loc=mprior_mean, scale=par_std,
                                        size=(ne, nm)).T

                # For each pumping ens
                for sub_param_i in range(nm):
                    params.add(
                        name=p_multop[1]+str(sub_param_i),
                        value=param_init[sub_param_i, 0],
                        min=parambound.loc[
                            parambound.Wellnest == wellnestlist[0]].iloc[sub_param_i, 0],
                        max=parambound.loc[
                            parambound.Wellnest == wellnestlist[0]].iloc[sub_param_i, 1])

                # Number of observations
                nd = dobs.size
                nd_sub = len(gw_obs_indices[0])
                nd_gw = nd - nd_sub
                # Associated standard deviation: ones (for this scenario)
                obs_error_sub = 1.5
                obs_error_gw = obs_error
                obs_std = np.append(np.ones(nd_sub) * obs_error_sub,
                                    np.ones(nd_gw) * obs_error_gw)
                obs_var = obs_std**2

                # Running ESMDA
                out = lmfit.minimize(forward_gwparam_subSELECT_pump_ls, params,
                                     method='least_squares',
                                     args=(wellnestlist, Pastasfiles,
                                           lenfiles, proxyflag,
                                           model_path, pumpflag,
                                           mode, tmin, tmax,
                                           Thick_data, K_data, Sskv_data,
                                           Sske_data, CC,
                                           Nz, num_clay, all_well4_data,
                                           well_data_dates,
                                           wellnest,
                                           SS_data, p_multop,
                                           ic_run, return_sub, dobs,
                                           gw_obs_indices, models,
                                           well_names, n_pump,
                                           annual_pump,
                                           listdaily_pump,
                                           nd_sub, nd_gw, obs_var,
                                           None, None, None, None,
                                           None, None,
                                           param_index,
                                           pastas_optparam,
                                           other_i))

                # Returns ESMDA solver and pastas models
                return out, models, well_names

            # If parameter calibration and pumping estimation
            elif "my_pump" == esmdaflag:

                nd_sub = 0

                # Mprior
                mprior_init = pd.DataFrame()

                # Parameter indexing
                param_index = esmdaindex

                # Renaming pumping ensemble
                pumping_ens = pump_ens

                mprior_init = pd.concat([mprior_init.reset_index(drop=True),
                                         pumping_ens.reset_index(drop=True)],
                                        axis=0,
                                        ignore_index=True)

                # Bounds
                n_pump = len(pumping_ens)
                mins = np.zeros(n_pump)
                maxs = np.repeat(4E2, n_pump)
                m_bounds = list(zip([mins, maxs]))

                # Number of observations
                nd = dobs.size
                nd_gw = nd - nd_sub
                # Associated standard deviation: ones (for this scenario)
                obs_error_sub = 1.5
                obs_error_gw = obs_error
                obs_std = np.append(np.ones(nd_sub) * obs_error_sub,
                                    np.ones(nd_gw) * obs_error_gw)
                obs_var = obs_std**2

                # Running ESMDA
                out = esmda_gw_my(forwardmy_pump, mprior_init, m_bounds,
                                  alphas, dobs, obs_var, par_error,
                                  args=(models, dobs, gw_obs_indices, Pastasfiles,
                                        lenfiles, model_path, proxyflag,
                                        well_names, pumpflag, ne, n_pump,
                                        annual_pump, listdaily_pump, SS_data,
                                        None, None),
                                  return_steps=True,
                                  return_cost=True)

                # Returns ESMDA solver and pastas models
                return out, models, well_names

            # If sub parameter cali
            elif "my_sub" == esmdaflag:

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

                # Only subsidence parameters
                if p_multop[1] == "all":
                    nm = 3
                    p_mult = []
                    [p_mult.append(random.uniform(m_bounds[x][0], m_bounds[x][1]))
                     for x in range(nm)]

                # Only subsidence parameters
                elif p_multop[1] == "Ss":
                    nm = 2
                    p_mult = []
                    [p_mult.append(random.uniform(m_bounds[x][0], m_bounds[x][1]))
                     for x in range(nm)]

                    par_error = par_error[:2]

                # Only Sskv and K
                elif p_multop[1] == "SsK":
                    nm = 2
                    p_mult = []
                    [p_mult.append(random.uniform(m_bounds[x][0], m_bounds[x][1]))
                     for x in np.arange(0, nm+1, nm)]

                    par_error = par_error[0:nm+1:nm]

                    par_error = par_error[0:nm+1:nm]

                elif p_multop[1] == "Sskv":

                    p_mult = [random.uniform(
                        m_bounds[0][0], m_bounds[0][1])]

                    par_error = par_error[0]

                elif p_multop[1] == "Sske":

                    p_mult = [random.uniform(
                        m_bounds[1][0], m_bounds[1][1])]

                    par_error = par_error[1]

                elif p_multop[1] == "K":

                    p_mult = [random.uniform(
                        m_bounds[2][0], m_bounds[2][1])]

                    par_error = par_error[2]

                # Number of model parameters
                nm = len(p_mult)

                # Prior: Let's start average of the bounds
                mprior_mean = p_mult

                mprior_mean = np.array(mprior_mean)

                # mprior_mean = np.ones(nm) * random.uniform(m_bounds[x][0], m_bounds[x][1])

                # Parameter error
                par_std = par_error

                # Prior distribution of parameters
                param_init = rng.normal(loc=mprior_mean, scale=par_std,
                                        size=(ne, nm)).T

                # Preallocation of empty df
                mprior_init = pd.DataFrame(param_init)

                # Prior distribution of hidden state
                for hid in init_state[1]:

                    state_init = rng.normal(loc=hid, scale=5,
                                            size=(ne, len(hid))).T

                    state_init = pd.DataFrame(state_init)
                    mprior_init = pd.concat([mprior_init, state_init], axis=0,
                                            ignore_index=True)
                    # State boundaries
                    state_min = np.empty(len(hid))
                    state_min[:] = np.nan
                    state_max = np.empty(len(hid))
                    state_max[:] = np.nan

                    m_bounds += list(zip(state_min, state_max))

                # # Uniform
                # mprior_init = np.empty([nm, ne])

                # for x in range(nm):
                #     if x == 0:
                #         mprior_init[x] = rng.uniform(m_bounds[x][0],
                #                                      m_bounds[x][1]/10,
                #                                      size=ne)
                #     else:
                #         mprior_init[x] = rng.uniform(m_bounds[x][0],
                #                                      m_bounds[x][1],
                #                                      size=ne)

                mprior_init = pd.DataFrame(mprior_init)

                # If number of assimilation windows > 1
                # Preallocation and separating obs into different windows
                if na_win > 1:

                    # Where subsets of obs kept
                    obs_df = np.array_split(dobs, na_win)

                    outs = []

                    # Number of observations
                    for win in range(na_win):

                        nd_sub = len(obs_df[win])

                        # Associated standard deviation: ones (for this scenario)
                        obs_error_sub = obs_error
                        obs_std = np.ones(nd_sub) * obs_error_sub
                        obs_var = obs_std**2

                        # Running ESMDA
                        out = esmda_sub_my(forward_sub, mprior_init, m_bounds,
                                           alphas, obs_df[win], obs_var, par_error,
                                           args=(wellnestlist,
                                                 mode, tmin, tmax,
                                                 Thick_data, K_data, Sskv_data,
                                                 Sske_data, CC,
                                                 Nz, num_clay, all_well4_data,
                                                 well_data_dates,
                                                 wellnest,
                                                 SS_data, p_multop,
                                                 ic_run, return_sub, obs_df[win],
                                                 obs_df[win].index,
                                                 well_names, init_state[0]),
                                           na_win=na_win,
                                           return_steps=True,
                                           return_cost=True)
                        outs.append(out)

                        # New mprior_init from old assimilation window
                        mprior_init = pd.DataFrame(out["mprior"][-1, :, :])

                    # Returns ESMDA solver and pastas models
                    return outs, models

                else:
                    nd_sub = len(dobs)

                    # Associated standard deviation: ones (for this scenario)
                    obs_error_sub = obs_error
                    obs_std = np.ones(nd_sub) * obs_error_sub
                    obs_var = obs_std**2

                    # Running ESMDA
                    out = esmda_sub_my(forward_sub, mprior_init, m_bounds,
                                       alphas, dobs, obs_var, par_error,
                                       args=(wellnestlist,
                                             mode, tmin, tmax,
                                             Thick_data, K_data, Sskv_data,
                                             Sske_data, CC,
                                             Nz, num_clay, all_well4_data,
                                             well_data_dates,
                                             wellnest,
                                             SS_data, p_multop,
                                             ic_run, return_sub, dobs,
                                             gw_obs_indices,
                                             well_names, init_state[0]),
                                       na_win=na_win,
                                       return_steps=True,
                                       return_cost=True)

                    # Returns ESMDA solver and pastas models
                    return out, models

            # If sub parameter cali
            elif "my_all" == esmdaflag:

                # Number of sub obs
                nd_sub = len(gw_obs_indices[0])

                # Annual pumping data (mean), std
                pumppath = os.path.join(os.path.abspath("inputs"),
                                        "BasinPumping_Annual_ESMDA.xlsx")
                pumpsheet = "EstTotalPump_54-60"
                annual_pump = pd.read_excel(pumppath, sheet_name=pumpsheet,
                                            index_col=0, parse_dates=["Date"])

                # List of daily dates to be interpolated
                pumpsheet = "EstTotalPump_54-60_Int50"
                listdaily_pump = pd.read_excel(pumppath, sheet_name=pumpsheet,
                                               index_col=0,
                                               parse_dates=["Date"])

                # Generating ensemble of time series
                option = ["normal", .99]
                pumping_ens = generate_pumping_ens(annual_pump, ne_pump, option)

                # Mprior
                mprior_init = pd.DataFrame()
                for parami in param_ens:

                    mprior = parami.T
                    mprior = pd.DataFrame(mprior)
                    mprior_init = pd.concat([mprior_init.reset_index(drop=True),
                                             mprior.reset_index(drop=True)],
                                            axis=0,
                                            ignore_index=True)

                mprior_init = pd.concat([mprior_init.reset_index(drop=True),
                                         pumping_ens.reset_index(drop=True)],
                                        axis=0,
                                        ignore_index=True)

                # Number of subsidence parameters
                # PARAMETER BOUNDARIES!
                # parambound_path = os.path.join(os.path.abspath("inputs"),
                #                                "SUBParametersCali.xlsx")
                parambound_path = os.path.join(os.path.abspath("inputs"),
                                               "SUBParametersPriortoManual.xlsx")

                parambound = pd.read_excel(parambound_path,
                                           sheet_name="bounds_mult",
                                           index_col=0)
                parambound = pd.DataFrame(parambound)

                sub_m_bounds = list(zip(parambound.loc[
                    parambound.Wellnest == wellnestlist[0]].iloc[0:3, 0],
                    parambound.loc[
                        parambound.Wellnest == wellnestlist[0]].iloc[0:3, 1]))

                # Number of sub param
                sub_nm = len(parambound.loc[
                    parambound.Wellnest == wellnestlist[0]].iloc[0:3, 0])
                for x in range(sub_nm):

                    test = rng.uniform(sub_m_bounds[x][0],
                                       sub_m_bounds[x][1],
                                       size=ne_pump)
                    test = test[None, :]
                    mprior_init = np.append(mprior_init,
                                            test, axis=0)

                mprior_init = pd.DataFrame(mprior_init)

                # Bounds
                n_pump = len(annual_pump)
                mins = np.append(min_pastas, np.zeros(n_pump))
                mins = np.append(mins, parambound.loc[
                    parambound.Wellnest == wellnestlist[0]].iloc[0:3, 0])
                mins[np.isnan(mins)] = -30000
                maxs = np.append(max_pastas, np.repeat(5E2, n_pump))
                maxs = np.append(maxs, parambound.loc[
                    parambound.Wellnest == wellnestlist[0]].iloc[0:3, 1])
                maxs[np.isnan(maxs)] = 10000
                m_bounds = list(zip((mins, maxs)))

                # Number of observations
                nd = dobs.size
                nd_gw = nd - nd_sub
                # Associated standard deviation: ones (for this scenario)
                obs_error_sub = 3
                obs_error_gw = 10
                obs_std = np.append(np.ones(nd_sub) * obs_error_sub,
                                    np.ones(nd_gw) * obs_error_gw)
                obs_var = obs_std**2

                # Alphas
                cov_obs_inflation_geo = 1.2
                cov_obs_inflation_factors = [1.1]  # list[float]
                for i in range(1, na):

                    cov_obs_inflation_factors.append(
                        cov_obs_inflation_factors[i - 1] / cov_obs_inflation_geo
                    )

                scaling_factor = np.sum(1 /
                                        np.array(cov_obs_inflation_factors))
                # : float
                cov_obs_inflation_factors = [
                    alpha * scaling_factor
                    for alpha in cov_obs_inflation_factors
                ]
                alphas = np.array(cov_obs_inflation_factors)

                out = esmda_gw_my(forward_all_my, mprior_init, m_bounds,
                                  alphas, dobs,
                                  obs_var,
                                  args=(wellnestlist, Pastasfiles,
                                        lenfiles, proxyflag,
                                        model_path, pumpflag,
                                        mode, tmin, tmax,
                                        Thick_data, K_data, Sskv_data,
                                        Sske_data, CC,
                                        Nz, num_clay, all_well4_data,
                                        well_data_dates,
                                        wellnest,
                                        SS_data, p_multop,
                                        ic_run, return_sub, dobs,
                                        gw_obs_indices, models,
                                        well_names, n_pump,
                                        annual_pump,
                                        listdaily_pump,
                                        nd_sub, nd_gw),
                                  return_steps=True)
                return out, models

            elif "PF" == esmdaflag:

                # Mprior
                mprior_init = pd.DataFrame()

                # Parameter indexing
                param_index = esmdaindex

                # Initial samples of parameters
                # should be [nd x ne]
                for param_i in param_ens:

                    mprior = param_i.iloc[:, param_index].T
                    mprior = pd.DataFrame(mprior)
                    mprior_init = pd.concat([mprior_init.reset_index(drop=True),
                                             mprior.reset_index(drop=True)],
                                            axis=0,
                                            ignore_index=True)

                # Other index
                range_ = set(range(0, len(pastas_optparam[0])))
                other_i = np.array(list(
                    set(param_index).symmetric_difference(range_)))

                # Boundaries of parameters
                mins = np.append(min_pastas[:, param_index], np.zeros(1))
                maxs = np.append(max_pastas[:, param_index], np.zeros(1))
                mins = mins[:-1]
                maxs = maxs[:-1]
                m_bounds = list(zip([mins, maxs]))

                # Running ESMDA
                out = PF(forwardmy_gwparamSELECT, mprior_init, m_bounds,
                         dobs, like_sigma,
                         args=(models, dobs, gw_obs_indices, None, None,
                               param_index, pastas_optparam,
                               other_i),
                         return_steps=True,
                         return_cost=True)

                # Returns ESMDA solver and pastas models
                return out, models, well_names

        # A regular subsidence model run
        else:

            sub_total, subv_total, all_results = run_sub(num_clay,
                                                         all_well4_data,
                                                         well_data_dates, mode,
                                                         tmin, tmax, SS_data,
                                                         wellnest,
                                                         K_data, Sskv_data,
                                                         Sske_data,
                                                         CC, Nz, Thick_data,
                                                         ic_run,
                                                         sub_total, subv_total,
                                                         all_results)

    # Returns heads in clay nodes, z dist, cum sub time series for each
    # well, cum inelastic sub time series for each well, original time
    # step
    return all_results, sub_total, subv_total


# %%###########################################################################
# Post processes data
##############################################################################

# Need to downsample sub data into daily
def bkk_postproc(wellnestlist, sub_total, subv_total, all_results):
    """Take results of calcs, cleans it, reinterpolates to original date.

    wellnestlist - list of wellnest to calculate subsidence for
    sub_total - list of lists: wellnest, well, interp t, cum sub results (m)
    subv_total - list of lists: wellnest, well, interp t, cum sub inelastic
    results (m)
    all_results - lists of lists: wellnest, well, time original (original time
    series (0:len(date))), date, heads in clay nodes, z distribution, time in
    spin up period, heads in clay nodes in spin up (last two if run
    spin up period)

    Returns: sub_total - list of lists: reinterpolated sub_total cum sub
    results (m)
    [4] index
    subvtotal - list of lists:  reinterpolated subv_total cum sub inelastic
    results (m): [4] index
    annual_data_all - lists of lists of annual total sub for all 4 clay for
    each well nest
    avg_sub_perc - lists of lists of average subsidence percentage for each well
    from total subsidence across all time steps
    """
    # Preallocation
    # Saves annual total sub for all 4 clay for all wellnests
    annual_data_all = []

    # Saves sub percentage from each well
    avg_sub_perc = []

    # For each well nest in the list
    # num_well is the index, wellnest = name
    for num_well, wellnest in enumerate(wellnestlist):

        # Total sub for all 4 clay for one wellnest
        # Resets after one well nest completed
        cumsum_4cl = np.zeros([4, len(all_results[num_well*4][2])])

        # Assumes four wells for each wellnest
        # For each well in each wellnest
        # num_well*4+i guarantees well indexing within each well nest
        for i in range(4):

            # time original (original time series (0:len(date)))
            t_og = all_results[num_well*4+i][2]

            # Reinterpolated to original time series
            # [2] has the model time series, [3] has the cum sub results
            sub_total[num_well*4+i].append(np.interp(t_og,
                                                     sub_total[num_well*4+i][2],
                                                     sub_total[num_well*4+i][3]))

            # Adding this new interpolation to cum_sum4cl
            # Each i for each well in one well nest
            cumsum_4cl[i] += sub_total[num_well*4+i][4]

            # Reinterpolated to original time series
            # [2] has the model time series, [3] has the cum sub results
            subv_total[num_well*4+i].append(np.interp(t_og,
                                                      subv_total[num_well*4+i][2],
                                                      subv_total[num_well*4+i][3]))

        #  original date series
        date = all_results[num_well*4][3]
        # For each well nest, creating new data frame for the cum total sub sum
        df_data = {"Date": date, "CumTotSum": np.sum(cumsum_4cl, axis=0)}
        df = pd.DataFrame(df_data,
                          columns=["Date", "CumTotSum"],
                          index=date)
        df["month"] = df.index.month
        df["day"] = df.index.day

        # Resampling to each year
        annual_data = df.CumTotSum[(df["month"] == 12) &
                                   (df["day"] == 31)].to_frame()
        annual_data["year"] = annual_data.index.year

        # # IMPORTANT INFO
        # # For benchmark measurements, the first year is 0, the second year is
        # # the compaction rate over that first year.
        # # For implicit Calc, the first year has a compaction rate over that
        # # year, so need to move Implicit values down one to match benchmark
        # # measurements.
        # # Index has the right years
        # First data value is the previous year at 0 compaction
        # Adds an extra year to the end
        annual_data.loc[annual_data.index[-1] + pd.offsets.DateOffset(years=1)] = 0
        annual_data = annual_data.shift(1)  # Shifts all values down one year
        annual_data.iloc[0] = 0  # Sets first value as 0
        annual_data.index = annual_data.index.shift(-12, freq="M")

        # Adding annual rates
        annual_data["AnnRates"] = annual_data.CumTotSum.diff()

        # Saving annual data for all well nests
        annual_data_all.append([wellnest, annual_data])

        # Looking at sub percentages for each well
        all_clay = np.sum(cumsum_4cl, axis=0)
        for j in range(4):
            avg_sub_perc.append([wellnest, j,
                                 (np.average(np.divide(cumsum_4cl[j, 1:],
                                                       all_clay[1:])))])
    # Returning
    return sub_total, subv_total, annual_data_all, avg_sub_perc
