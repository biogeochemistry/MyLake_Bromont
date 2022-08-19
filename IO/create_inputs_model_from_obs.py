#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : Marianne Cote
# Created Date: 2022-03-22
# version ='1.0'
# ---------------------------------------------------------------------------
""" Script formatting the Raw data into correct format for MyLake"""
# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
import numpy as np
import pandas as pd
import h5py
import datetime
import os

# ---------------------------------------------------------------------------
# Global Variables
# ---------------------------------------------------------------------------
variables = ['clt', 'hurs', 'tas', 'rsds', 'ps', 'pr', 'sfcWind']
lakes_dict = {"Bromont": "Bromont"}

rawdata_directory = r"../Raw_data/Data"
observation_directory = r"../obs"
input_directory = r"../IO"

depth_resolution: float = 1  # MC 2022-03-04 default: 1.0 NOTE: don't change this unless you know what you are doing.


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

def get_bathymetry_formatted(raw_bathymetry_depths: list, raw_bathymetry_areas: list,
                             resolution: float = 1.0):
    """
    Used the bathymetry measured and generate a bathymetry in function of the depth resolution specify. 
    For the range of depth created in function of the depth resollution, the script interpolates the area for the depth 
    without measured.

    :param raw_bathymetry_depths:   List of depths in the water column with area value
    :type raw_bathymetry_depths:    list[float]

    :param raw_bathymetry_areas:    List of the area for each depth given
    :type raw_bathymetry_areas:     list[float]

    :param resolution:              Resolution for the bathymetry wanted. The value represente the height between depth.
                                    The resolution is by default set to 1 meter.
    :type resolution:               float

    :return:                        depth_levels:   List of the depth for the new range using the given depth resolution
                                    area_levels:    List of the area for each depth given in the variable depth_levels
    """
    # Create depth_levels in function of maximum depth and depth resolution
    if raw_bathymetry_depths[-1] == int(raw_bathymetry_depths[-1]):
        max_depth = raw_bathymetry_depths[-1]
    elif raw_bathymetry_depths[-1] == int(raw_bathymetry_depths[-1]) + 0.5:
        max_depth = raw_bathymetry_depths[-1]
    elif raw_bathymetry_depths[-1] > int(raw_bathymetry_depths[-1]) + 0.5:
        max_depth = int(raw_bathymetry_depths[-1]) + 0.5
    else:
        max_depth = int(raw_bathymetry_depths[-1])

    depth_levels = np.arange(0, max_depth + resolution, resolution)

    # Get Area for each depth levels in function of depth_resolution
    area_levels = [None] * len(depth_levels)
    for i in range(0, len(depth_levels)):
        i: int = int(i)
        depth = float(depth_levels[i])
        if depth in raw_bathymetry_depths:
            depth_index = [i for i, val in enumerate(raw_bathymetry_depths) if val == depth]
            if len(depth_index) == 1:
                area_levels[i] = raw_bathymetry_areas[depth_index[0]]
            else:
                n = 0
                sum_area = 0
                for pos in depth_index:
                    sum_area += raw_bathymetry_areas[pos]
                    n += 1
                area_levels[i] = sum_area / n
        else:
            xa, xb, ya, yb = False, False, False, False
            for x_depth in range(0, len(raw_bathymetry_depths)):
                if raw_bathymetry_depths[x_depth] < depth:
                    xa = raw_bathymetry_depths[x_depth]
                    ya = raw_bathymetry_areas[x_depth]
                else:
                    xb = raw_bathymetry_depths[x_depth]
                    yb = raw_bathymetry_areas[x_depth]
                    break
            if xa is False:
                area_levels[i] = raw_bathymetry_areas[0]
            elif xb is False:
                area_levels[i] = raw_bathymetry_areas[-1]
            else:
                area_levels[i] = findYPoint(xa, xb, ya, yb, depth)

    return depth_levels, area_levels


def get_initial_concentration(depth_levels: list, obs_concentration: pd.DataFrame,
                              date_with_initial_concentration: str = "2021-01-01"):
    """
    get initinal concentration into the variable needed for the modelisation.
    It included conversion into the units process by the model, creation of variation from data observed and generate
    variable with default values if not measured. It also interpolate the data from depth without data.

    :param depth_levels:                    List of the depth in the water column that will be simulated by the model
    :type depth_levels:                     list

    :param obs_concentration:               Table for the measured variable value by time.
                                            First row is the column names:
                                                [Date, Depth, Temp, Turbidity, TP, TDP, Chl_a, DOC, DO, pH]
                                            Data format by column:
                                                Date:       Sample date, string in format YYYY-MM-DD
                                                Depth:      Sample Depth (by meters), float
                                                Temp:       Temperature measure (°C), float
                                                Turbidity:  Suspended sediment in the water (mg/l), float
                                                TP:         Total Phosphorus (µg/l), float
                                                TDP:        Total Dissolved Phosphorus (µg/l), float
                                                PO4_P:      Phosphate (µg/l)
                                                Chl_a:      Chlorophyll a (mg/m3), float
                                                DOC:        Dissolved Organic Carbon (mg/l), float
                                                DO:         Dissolved Oxygen (mg/l), float
                                                pH:         pH (-), float
    :type obs_concentration:                Pandas DataFrame

    :param date_with_initial_concentration: date (YYYY-MM-DD) with observation that will be used to set the first day
                                            for simulation.
    :type date_with_initial_concentration:  str

    :return:                                Lists of value for the variables (float):
                                                Tz:         Temperature measure (°C)
                                                POCz:       Suspended inorganic matter (turbidity) profile (mg/m3)
                                                TPz:        Total P profile (mg/m3)
                                                DOPz:       Dissolved organic P profile(DOP = TDP - PO4_P)(mg/m3)
                                                Chlz:       Chlorophyll a profile (mg/m3)
                                                DOCz:       DOC profile (mg/m3)
                                                O2z:        Oxygen profile (mg/m3)
                                                pHz:         pH (-)
                                                Hice:       ***ESTIMATED*** ice thicknesses (m)
                                                Hsnow:      ***ESTIMATED*** snow thicknesses (m)
                                                POPz:       POP = TP - TDP(mg/m3)
    """

    # Set variable list
    dict_variable_names = {"Temp": [None] * len(depth_levels),
                           "DO": [None] * len(depth_levels)}

    #corrects depth if negative
    obs_concentration.loc[obs_concentration["Depth"]< 0, ["Depth"]] += depth_levels[-1]
    # Get values for the initial date
    inital_concentrations = obs_concentration.loc[
        obs_concentration['Date'] == date_with_initial_concentration].sort_values(by=['Depth'])

    for variable in dict_variable_names.keys():
        not_empty_cell = inital_concentrations.loc[inital_concentrations[variable].notna()]
        not_empty_cell.loc[(not_empty_cell[variable] == '<1'), variable] = 1
        concentration_data = [float(x) for x in list(not_empty_cell[variable])]
        observations_depths = list(not_empty_cell['Depth'])
        for i in range(0, len(depth_levels)):
            depth = depth_levels[i]
            if depth in observations_depths:
                depth_index = [i for i, val in enumerate(observations_depths) if val == depth]
                if len(depth_index) == 1:
                    dict_variable_names[variable][i] = concentration_data[depth_index[0]]
                else:
                    n = 0
                    sum_area = 0
                    for pos in depth_index:
                        sum_area += concentration_data[pos]
                        n += 1
                    dict_variable_names[variable][i] = sum_area / n
            else:
                xa, xb, ya, yb = False, False, False, False
                for x_depth in range(0, len(observations_depths)):
                    if observations_depths[x_depth] < depth:
                        xa = observations_depths[x_depth]
                        ya = concentration_data[x_depth]
                    else:
                        xb = observations_depths[x_depth]
                        yb = concentration_data[x_depth]
                        break
                if xa is False:
                    dict_variable_names[variable][i] = concentration_data[0]
                elif xb is False:
                    dict_variable_names[variable][i] = concentration_data[-1]
                else:
                    dict_variable_names[variable][i] = findYPoint(xa, xb, ya, yb, depth)

    Tz, O2z = dict_variable_names["Temp"], [element * 1000 for element in dict_variable_names["DO"]]# *1000 to convert mg/l to mg/m3


    return Tz, O2z


def findYPoint(xa: float, xb: float, ya: float, yb: float, xc: float):
    """
    Function used to calculate the temperature at a depth non simulated by the model. MyLake simulates the temperature
    at each meter (starting at 0.5) and this function permit to compare the temperature at the same depth that it has
    been measured.

    :param xa:  Closest depth (m) simulated below the wanted depth.
    :type xa:   float

    :param xb:  Closest depth (m) simulated over the wanted depth.
    :type xb:   float

    :param ya:  Temperature (°C) at the depth xa.
    :type ya:   float

    :param yb:  Temperature (°C) at the depth xb.
    :type yb:   float

    :param xc:  Depth (m) at which the temperature is wanted.
    :type xc:   float

    :return:    yc: Temperature (°C) at the depth xc.
    """
    m = (float(ya) - float(yb)) / (float(xa) - float(xb))
    yc = (float(xc) - (float(xb))) * m + float(yb)
    return yc



def nbrleapyears(start, end):  # MC 2018-07-10
    """
    determine the number of leap years in the date range
    :param start: start year
    :param end: end year
    :return: number of leap year between the start and end years
    """
    nbryears = 0
    while start <= end:
        if (start % 4 == 0 and start % 100 != 0) or start % 400 == 0:
            nbryears += 1
        start += 1
    return nbryears


class Simulation_informations:
    """
    Object containing information for the simulation
    """

    def __init__(self,  depth_resolution=1.0, input_folder=r"../IO",
                 observation_folder=r"../obs",
                 output_folder=r"../Postproc_code", start_year=2018, end_year=2022, enable_river_inflow=False):

        # relative path to get the information and where the format will be save
        self.input_folder = input_folder
        if not os.path.exists(self.input_folder):
            os.makedirs(self.input_folder)
        self.observation_folder = observation_folder
        if not os.path.exists(self.observation_folder):
            os.makedirs(self.observation_folder)
        self.output_folder = output_folder
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Simulation information
        self.start_year = start_year
        self.end_year = end_year
        self.depth_resolution = depth_resolution
        self.climate_file = os.path.join(observation_folder, "climate_data.csv")
        self.enable_river_inflow = enable_river_inflow

    def mylakeinput(self, outpath: str, inflows_filename: str):
        """
        create a file containing the informations relatively to Mylake
        :param pA: dictionary of paths to HDF5 files
        :param pB: dictionary of paths to HDF5 files
        :param inflowfile: filename of the inflowfile
        :param outpath: filename where a file of Mylake input will be written
        :type pA: dict
        :type pB: dict
        :type inflowfile: str
        :type outpath: str
        :return: string to be written to a file
        """

        all_days = pd.date_range(pd.datetime(self.start_year, 1, 1), pd.datetime(self.end_year, 12, 31), freq='D')
        meteo = pd.read_csv(self.climate_file)


        meteo["Date"] = pd.to_datetime(meteo['Date'])
        # meteo = meteo.sort_values(by = "Date")

        # mask = (meteo['Date'] >= pd.datetime(self.start_year, 1, 1)) & (
        #         meteo['Date'] <= pd.datetime(self.end_year, 12, 31))
        # meteo = meteo.loc[mask].reset_index(drop=True)
        # meteo = meteo.set_index("Date")
        meteo = meteo.groupby(["Date"]).mean()
        meteo["Date"] = meteo.index
        meteo.index = pd.DatetimeIndex(meteo['Date']).floor('D')
        meteo = meteo.reindex(all_days, fill_value=np.nan)
        # meteo = meteo.reindex(all_days)
        meteo["Date"] = meteo.index
        meteo = meteo.interpolate(method='pad', limit=4, limit_direction= 'forward')


        meteo["month_day"] = meteo['Date'].dt.strftime('%m-%d')
        meteoaverage = meteo.groupby(["month_day"]).mean()
        meteoaverage.reset_index(inplace=True)

        meteo_to_estimate = meteo.dropna(axis=1, how='all')
        # emptyrows, emptycolumns = np.where(pd.isnull(meteo.dropna(axis=1, how='all')))
        for column in meteo_to_estimate.columns[:-2]:
            print(column)
            estimatedvalues = meteoaverage.dropna(axis=0,subset=[column])
            for day in estimatedvalues["month_day"]:
                print(day)
                meteo.loc[(meteo[column].isnull())& (meteo["month_day"] == day),column] =  float(estimatedvalues.loc[estimatedvalues["month_day"] == day,column])


        # for i in range(0,len(emptyrows)):
        #     date = meteo["month_day"][emptyrows[i]]
        #     meteoaverageline = meteoaverage.loc[meteoaverage["month_day"] == date]
        #     print(emptyrows[i], emptycolumns[i], meteo.iloc[int(emptyrows[i]),int(emptycolumns[i])])
        #     meteo.iloc[int(emptyrows[i]),int(emptycolumns[i])] = meteoaverageline.iloc[0,int(emptycolumns[i])]


        meteo = meteo.interpolate('pad').ffill()
        meteo = meteo.interpolate('pad').bfill()

        meteo.to_csv(self.climate_file.replace('.csv','_filled_interpolate.csv'), index=False)

        'Year	Month	Day	GlobalRadiation	CloudCover	AirTemperature	RelativeHumidity	AirPressure	WindSpeed	Precipitation	InflowQ	InflowT	InflowC	POC	InflowTP'
        '	InflowDOP	InflowChla	DOC	DIC	O	NO3	NH4	SO4	Fe2	Ca2	pH	CH4	Fe3	Al3	SiO4	SiO2	diatom	POP'

        if self.enable_river_inflow:
            inflows_data = pd.read_csv(inflows_filename)
            inflows_data = inflows_data.interpolate()

            inflows_data["Date"] = pd.to_datetime(inflows_data['Date'])

            mask = (inflows_data['Date'] >= pd.datetime(self.start_year, 1, 1)) & (
                    inflows_data['Date'] <= pd.datetime(self.end_year, 12, 31))
            inflows_data = inflows_data.loc[mask].reset_index(drop=True)
            inflows_data = inflows_data.loc[inflows_data["flow"] == 'inlet']

            # inflows_data["InflowDOP"] = inflows_data["InflowTDP"]- inflows_data["InflowPO4"] negative values?

            # POP = TP - TDP(mg/m3)
            inflows_data["InflowPOP"] = inflows_data["InflowTP"] - inflows_data["InflowTDP"]

            if inflows_data['Date'][0] >= pd.datetime(self.start_year + 1, 1, 1):
                mask = inflows_data['Date'].dt.year == int(self.start_year + 1)
                copy_year = inflows_data[mask].reset_index(drop=True)
                copy_year['Date'] = copy_year['Date'] - pd.offsets.DateOffset(years=1)
                inflows_data = pd.concat([copy_year, inflows_data], ignore_index=True)
            if inflows_data['Date'][len(inflows_data) - 1] <= pd.datetime(self.end_year - 1, 12, 31):
                mask = inflows_data['Date'].dt.year == int(self.end_year - 1)
                copy_year = inflows_data[mask].reset_index()
                copy_year['Date'] = copy_year['Date'] + pd.offsets.DateOffset(years=1)
                inflows_data = pd.concat([inflows_data, copy_year], ignore_index=True)

            inflows_data.index = pd.DatetimeIndex(inflows_data['Date']).floor('D')

            inflows_data = inflows_data.reindex(all_days)
            inflows_data["Date"] = inflows_data.index
            # inflows_data[:][0] = inflows_data[~inflows_data.isnull().any(axis=1)][0]

            inflows_data = inflows_data.interpolate('index').ffill()
            inflows_data = inflows_data.interpolate('index').bfill()

            allyears = list(inflows_data['Date'].dt.year.unique())
            for year in allyears:
                if year == allyears[-1]:
                    inflows_data.loc[
                        (inflows_data['Date'] >= '%s-12-21' % year) & (inflows_data['Date'] <= '%s-12-31' % year),
                        "InflowQ"] = inflows_data["InflowQ"] / 2
                elif year == allyears[0]:
                    test = inflows_data.loc[
                        (inflows_data['Date'] >= '%s-1-1' % year) & (inflows_data['Date'] <= '%s-03-20' % (year)),
                        "InflowQ"]
                    inflows_data.loc[
                        (inflows_data['Date'] >= '%s-1-1' % year) & (inflows_data['Date'] <= '%s-03-20' % (year)),
                        "InflowQ"] = inflows_data["InflowQ"] / 2
                    inflows_data.loc[
                        (inflows_data['Date'] >= '%s-12-21' % year) & (inflows_data['Date'] <= '%s-12-31' % (year)),
                        "InflowQ"] = inflows_data["InflowQ"] / 2
                else:
                    inflows_data.loc[(inflows_data['Date'] >= '%s-1-1' % year) & (
                                inflows_data['Date'] <= '%s-03-20' % (year + 1)), "InflowQ"] = inflows_data["InflowQ"] / 2

            # change DOC mg/L to DOC mg/m3. 1 mg/L = 0.001 mg/m3
            inflows_data["InflowDOC"] = inflows_data["InflowDOC"] * 0.001

            ndays = (datetime.date(self.end_year + 1, 1, 1) - datetime.date(self.start_year, 1, 1)).days

            repd = [datetime.date(self.start_year, 1, 1) + datetime.timedelta(d) for d in range(0,
                                                                                                ndays)]
            print(repd[-1])
            mlyear = np.array([d.year for d in repd])
            mlmonth = np.array([d.month for d in repd])
            mlday = np.array([d.day for d in repd])
            mlndays = ndays  # 365 + 365 + ndays

            repeati = list(range(ndays))  # list(range(365)) + list(range(365)) + list(range(ndays))
            # print(len(mlyear), len(meteo['rsds'][repeati]), mlndays)
            spacer = np.repeat([0], repeats=ndays)[repeati].reshape((mlndays, 1))
            # stream_Q = np.repeat([2000], repeats = ndays)[repeati].reshape((mlndays, 1))
            # stream_T = np.repeat([10], repeats = ndays)[repeati].reshape((mlndays, 1))
            stream_O = np.repeat([12000], repeats=ndays)[repeati].reshape(
                (mlndays, 1))  # MC 06-01-2018 initial parameters stream_O:8000
            stream_C = np.repeat([0.5], repeats=ndays)[repeati].reshape((mlndays, 1))
            # stream_TP = np.repeat([5], repeats = ndays)[repeati].reshape((mlndays, 1))
            stream_DOP = np.repeat([1], repeats=ndays)[repeati].reshape((mlndays, 1))
            stream_SS = np.repeat([0.01], repeats=ndays)[repeati].reshape((mlndays, 1))
            stream_Chl = np.repeat([0.1], repeats=ndays)[repeati].reshape((mlndays, 1))
            stream_DOC = np.repeat([2000], repeats=ndays)[repeati].reshape(
                (mlndays, 1))  # MC 06-01-2018 initial parameters 8000
            stream_DIC = np.repeat([20000], repeats=ndays)[repeati].reshape((mlndays, 1))
            temporarypath = '%s.temp' % outpath

            np.savetxt(temporarypath,
                       np.concatenate((mlyear.reshape((mlndays, 1)),
                                       mlmonth.reshape((mlndays, 1)),
                                       mlday.reshape((mlndays, 1)),
                                       meteo['Global radiation'][repeati].values.reshape((mlndays, 1)),
                                       meteo['Cloud cover'][repeati].values.reshape((mlndays, 1)),
                                       meteo['Air temperature'][repeati].values.reshape((mlndays, 1)),
                                       meteo['Relative humidity'][repeati].values.reshape((mlndays, 1)),
                                       meteo['Air pressure'][repeati].values.reshape((mlndays, 1)),
                                       # np.repeat([0], repeats = ndays)[repeati].reshape((mlndays, 1)),
                                       meteo['Wind speed'][repeati].values.reshape((mlndays, 1)),
                                       meteo['Precipitation'][repeati].values.reshape((mlndays, 1)),
                                       inflows_data['InflowQ'][repeati].values.reshape((mlndays, 1)),
                                       inflows_data['InflowTemp'][repeati].values.reshape((mlndays, 1)),
                                       stream_C, stream_SS,  # C, SS
                                       inflows_data['InflowTP'][repeati].values.reshape((mlndays, 1)), stream_DOP,
                                       # InflowTP InflowDOP
                                       stream_Chl, inflows_data['InflowDOC'][repeati].values.reshape((mlndays, 1)),
                                       # Chl, DOC
                                       stream_DIC, stream_O, spacer, spacer,  # DIC, O, NO3, NH4
                                       spacer, spacer, spacer,
                                       inflows_data['InflowpH'][repeati].values.reshape((mlndays, 1)),  # SO4, Fe, Ca, PH
                                       spacer, spacer, spacer, spacer,  # CH4, Fe3, Al3, SiO4
                                       spacer, spacer, inflows_data['InflowPOP'][repeati].values.reshape((mlndays, 1))),
                                      axis=1),  # SiO2, diatom
                       fmt=['%i', '%i', '%i',  # yy mm dd
                            '%.4g', '%.2f', '%.2f', '%i', '%i', '%.2f', '%.3f',  # rad, cloud, temp, hum, pres, wind, precip
                            '%.3f', '%.3f', '%.3f', '%.3f',  # InflowQ	InflowT	InflowC	POC
                            '%.3f', '%.3f', '%.3f', '%.3f',  # InflowTP InflowDOP	InflowChla	DOC
                            '%.3f', '%.3f', '%i', '%i',  # DIC	O	NO3	NH4
                            '%i', '%i', '%i', '%i',  # SO4	Fe2	Ca2	pH
                            '%i', '%i', '%i', '%i',  # CH4	Fe3	Al3	SiO4
                            '%i', '%i', '%i'],  # SiO2	diatom	POP
                       delimiter='\t',
                       header='mysterious useless line																																POC = TP - DOP\n'
                              'Year	Month	Day	GlobalRadiation	CloudCover	AirTemperature	RelativeHumidity	AirPressure	WindSpeed	Precipitation	'
                              'InflowQ	InflowT	InflowC	POC	InflowTP'
                              '	InflowDOP	InflowChla	DOC	DIC	O	NO3	NH4	SO4	Fe2	Ca2	pH	CH4	Fe3	Al3	SiO4	SiO2	diatom	POP')

        else:
            ndays = (datetime.date(self.end_year + 1, 1, 1) - datetime.date(self.start_year, 1, 1)).days

            repd = [datetime.date(self.start_year, 1, 1) + datetime.timedelta(d) for d in range(0,
                                                                                                ndays)]
            print(repd[-1])
            mlyear = np.array([d.year for d in repd])
            mlmonth = np.array([d.month for d in repd])
            mlday = np.array([d.day for d in repd])
            mlndays = ndays  # 365 + 365 + ndays

            repeati = list(range(ndays))  # list(range(365)) + list(range(365)) + list(range(ndays))
            # print(len(mlyear), len(meteo['rsds'][repeati]), mlndays)
            spacer = np.repeat([0], repeats=ndays)[repeati].reshape((mlndays, 1))
            stream_Q = np.repeat([2000], repeats = ndays)[repeati].reshape((mlndays, 1))
            stream_T = np.repeat([10], repeats = ndays)[repeati].reshape((mlndays, 1))
            stream_O = np.repeat([12000], repeats=ndays)[repeati].reshape(
                (mlndays, 1))  # MC 06-01-2018 initial parameters stream_O:8000
            stream_C = np.repeat([0.5], repeats=ndays)[repeati].reshape((mlndays, 1))
            stream_TP = np.repeat([5], repeats = ndays)[repeati].reshape((mlndays, 1))
            stream_DOP = np.repeat([1], repeats=ndays)[repeati].reshape((mlndays, 1))
            stream_SS = np.repeat([0.01], repeats=ndays)[repeati].reshape((mlndays, 1))
            stream_Chl = np.repeat([1], repeats=ndays)[repeati].reshape((mlndays, 1))
            stream_DOC = np.repeat([2000], repeats=ndays)[repeati].reshape(
                (mlndays, 1))  # MC 06-01-2018 initial parameters 8000
            stream_DIC = np.repeat([20000], repeats=ndays)[repeati].reshape((mlndays, 1))
            temporarypath = '%s.temp' % outpath

            np.savetxt(temporarypath,
                       np.concatenate((mlyear.reshape((mlndays, 1)),
                                       mlmonth.reshape((mlndays, 1)),
                                       mlday.reshape((mlndays, 1)),
                                       meteo['Global radiation'][repeati].values.reshape((mlndays, 1)),
                                       meteo['Cloud cover'][repeati].values.reshape((mlndays, 1)),
                                       meteo['Air temperature'][repeati].values.reshape((mlndays, 1)),
                                       meteo['Relative humidity'][repeati].values.reshape((mlndays, 1)),
                                       meteo['Air pressure'][repeati].values.reshape((mlndays, 1)),
                                       # np.repeat([0], repeats = ndays)[repeati].reshape((mlndays, 1)),
                                       meteo['Wind speed'][repeati].values.reshape((mlndays, 1)),
                                       meteo['Precipitation'][repeati].values.reshape((mlndays, 1)),
                                       stream_Q,stream_T,
                                       stream_C, stream_SS,  # C, SS
                                       stream_TP, stream_DOP,
                                       # InflowTP InflowDOP
                                       stream_Chl, stream_DOC,
                                       # Chl, DOC
                                       stream_DIC, stream_O, spacer, spacer,  # DIC, O, NO3, NH4
                                       spacer, spacer, spacer,spacer,
                                       # SO4, Fe, Ca, PH
                                       spacer, spacer, spacer, spacer,  # CH4, Fe3, Al3, SiO4
                                       spacer, spacer, spacer),
                                      axis=1),  # SiO2, diatom
                       fmt=['%i', '%i', '%i',  # yy mm dd
                            '%.4g', '%.2f', '%.2f', '%i', '%i', '%.2f', '%.3f',
                            # rad, cloud, temp, hum, pres, wind, precip
                            '%i', '%i', '%.3f', '%.3f',  # InflowQ	InflowT	InflowC	POC
                            '%i', '%.3f', '%.3f', '%i',  # InflowTP InflowDOP	InflowChla	DOC
                            '%.3f', '%.3f', '%i', '%i',  # DIC	O	NO3	NH4
                            '%i', '%i', '%i', '%i',  # SO4	Fe2	Ca2	pH
                            '%i', '%i', '%i', '%i',  # CH4	Fe3	Al3	SiO4
                            '%i', '%i', '%i'],  # SiO2	diatom	POP
                       delimiter='\t',
                       header='mysterious useless line																																POC = TP - DOP\n'
                              'Year	Month	Day	GlobalRadiation	CloudCover	AirTemperature	RelativeHumidity	AirPressure	WindSpeed	Precipitation	'
                              'InflowQ	InflowT	InflowC	POC	InflowTP'
                              '	InflowDOP	InflowChla	DOC	DIC	O	NO3	NH4	SO4	Fe2	Ca2	pH	CH4	Fe3	Al3	SiO4	SiO2	diatom	POP')

        with open(temporarypath) as f:
            with open("%s.txt" % outpath, 'w') as g:
                g.write(f.read().replace('-99999999', 'NaN'))
        os.unlink(temporarypath)

        return True


class Lake:
    """
    Create object Lake with all information relative to the lake analysed.
    """

    def __init__(self, lake_name, simulation_information):
        self.name = lake_name
        self.sim_info = simulation_information
        self.input_folder = os.path.join(simulation_information.input_folder, self.name)
        if not os.path.exists(self.input_folder):
            os.makedirs(self.input_folder)
        self.observation_folder = os.path.join(simulation_information.observation_folder, self.name)
        if not os.path.exists(self.observation_folder):
            os.makedirs(self.observation_folder)
        self.output_folder = os.path.join(simulation_information.output_folder, self.name)
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

    def mylakeinit(self, outpath: str):
        """
            create a file of a lake initiated with a max_depth and area.
            Assumes to have a cone shaped bathymetry curve

            :return: string to be written to an init file of MyLake
        """
        # 5-7-2018 MC

        # Get depth and area
        bathy = pd.read_csv(os.path.join(self.observation_folder, "%s_bathymetry.csv" % self.name))
        observation_concentration = pd.read_csv(
            os.path.join(self.observation_folder, "%s_observation_data.csv" % self.name))
        depth_levels, area_levels = get_bathymetry_formatted(list(bathy['Depth']), list(bathy['Area']),
                                                             depth_resolution)

        # get initial concentration and temperature
        Tz,  O2z = get_initial_concentration(depth_levels,
                                                                                    observation_concentration)

        lines = ['\t'.join(
            [('%.2f' % d), ('%.0f' % a)] + [('%.2f' % tz)] + ['0'] + ['0'] + ['200']+ ['4']+ ['11'] + ['0'] * 6 + [('%.2f' % o2z)] + ['0'] * 8 +
            ['0'] + ['0'] * 6 + ['0'])
            for d, a, tz,  o2z in
            zip(depth_levels, area_levels, Tz, O2z)] # phophore and chlorophyl from "Macrozooplankton and the persistence of the deepchlorophyll maximum in a stratiﬁed lake"(pannard, 2015)

        # lines[0] = lines[0] + '\t0\t0'  # snow and ice, plus 16 dummies
        firstlines = '''skip 
     Z (m)	  Az (m2)	 Tz (deg C)	  Cz(mg/m3)	  POCz (mg/m3)	  TPz (mg/m3)	 DOPz (mg/m3)	    Chlaz (mg/m3)	   DOCz (mg/m3)	    TPz_sed (mg/m3)	 Chlaz_sed (mg/m3)	   Fvol_IM (m3/m3, dry w.)	 Hice (m)	 Hsnow (m)	 O2z (mg/m3)	 DICz (mg/m3)	 NO3z (mg/m3)	 NH4z (mg/m3)	 SO4z (mg/m3)	 HSz (mg/m3)	 H2Sz (mg/m3)	 Fe2z (mg/m3)	 Ca2z (mg/m3)	 pHz (mg/m3)	 CH4aqz (mg/m3)	 Fe3z (mg/m3)	 Al3z (mg/m3)	 FeSz (mg/m3)	 CaCO3z (mg/m3)	 CH4gz (mg/m3)	 POPz (mg/m3)
    '''
        lines = [firstlines] + lines
        with open(outpath, 'w') as f:
            f.write('\n'.join(lines))

        return True

    def create_files(self):
        """

        :param modelid: model used
        :param scenarioid: scenario used
        :param depth: depth used for initiate Mylake (see mylakeinit())
        :param area: area used for initiate Mylake (see mylakeinit())
        :param longitude: longitude coordinate for Mylake (see mylakepar())
        :param latitude: latitude coordinate for Mylake (see mylakepar())
        :return:

        """
        # 5-7-2018 MC

        initp = os.path.join(self.input_folder, 'mylake_initial_concentrations.txt')
        inputp = os.path.join(self.input_folder, 'input_%s' % self.name)

        print("Create init file")
        if self.mylakeinit(outpath=initp):
            print("%s created!" % initp)
        else:
            print("%s was not created" % initp)

        inflows = os.path.join(self.observation_folder, "%s_inflow_data.csv" % self.name)
        if self.sim_info.mylakeinput(inputp, inflows):
            print("%s created!" % inputp)
        else:
            print("%s was not created" % inputp)


if __name__ == "__main__":
    import sys

    # default sys.argv : "Bromont"
    lakes = sys.argv[1:]
    print("****Default information for simulation: ***** \n"
          "bathymetry resolution: 1.0\n"
          "relative input path: \IO \n"
          "relative observation path: \obs \n"
          "relative output path: \Postproc_code \n"
          "Simulation period (including spin-off): from 2018-01-01 to 2022-12-31\n")
    print("***Generate Simulation_information Class****\n")
    sim_info = Simulation_informations()
    print("Class generated\n")

    print("***List of the lakes selected****\n")
    for lake in lakes:
        print("%s \n" % lake)

    for lake in lakes:
        print("***Generate Lake_information Class for lake %s****\n" % lake)
        lake_info = Lake(lake, sim_info)
        print("Class generated\n")
        print("**** input files creation ****\n")
        lake_info.create_files()

    # print("hello")
    # lakes = []
    #
    # create_files()

