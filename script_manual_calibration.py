#!/usr/bin/env python

""" Script for manual calibration - Milsbo lakes project
Script Launches a run of the Mylake model using the MyLake_v2_Vansjo version with the parameter's value giving in input,
shows figures of the comparison of the simulated values (for temperature, oxygen, and chl_a) with the observation for
the lake specifies in input (Nedre or Ovre) and loops until the value "End" is given in input.
For all looping period, a report (as table in CSV, naming as "manual_calibration_{lake name}_{Variable}_{start time}_{end time}.csv")
is generated in the output folder, containing for each round of calibration, the lake name, the value given to all
parameters, the calibration performances using a different index, an overall score from visual comparison
(giving in input by the user, between 1 and 10) and optional note for each run, if giving in input by the user.
The main run asks for the lake wanted (Nedre or Ovre) and the
"""

__author__ = "Marianne Cote"

from pandas.plotting import register_matplotlib_converters
import warnings

warnings.filterwarnings("ignore")
register_matplotlib_converters()

import matplotlib.pyplot as plt
from scipy.stats import linregress
import time
import statsmodels.api as smodels
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import scipy.stats as stats
import numpy as np
from math import sqrt, floor, log10, log
from sklearn.metrics import r2_score, mean_squared_error
import statistics
import csv
import os
import shutil
from datetime import datetime, timedelta
import pandas as pd
import seaborn as sns
from numpy import arange, nan, reshape, sqrt
from IPython import get_ipython
import scipy.io as sio

variable_pos = {'T': 5, 'O2': 6, 'Chl': 13}
num2month = {1: "Jan", 2: "Feb", 3: "Mar", 4: "Apr", 5: "May", 6: "June", 7: "July", 8: "Aug", 9: "Sept", 10: "Oct",
             11: "Nov", 12: "Dec"}
variables_dict = {'T': "Temperature", 'O2': "DO concentration", 'Chl': 'Chlorophyll a'}
matlab_folder = r"C:\Program Files\MATLAB\R2019b\bin\matlab"  # Value by default. need to be ajust to where matlab is install and the matlab version

test_cases = 10
parameters = {"Swa_b0": [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5],
              "Swa_b1": [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5],
              "C_shelter": [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
              "I_ScV": [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8],
              "I_ScT": [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8],
              "Alb_melt_ice": [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
              "Alb_melt_snow": [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]}

# full_lake_list = "Allequash"
regions = {"US": ["Allequash_Lake", "Annie", "Big_Muskellunge_Lake", "Black_Oak_Lake", "Crystal_Lake",
                  "Crystal_Bog", "Delavan", "Falling_Creek_Reservoir", "Fish_Lake", "Great_Pond",
                  "Green_Lake", "Laramie_Lake", "Mendota", "Monona", "Okauchee_Lake",
                  "Sammamish", "Sparkling_Lake", "Sunapee", "Tahoe", "Toolik_Lake",
                  "Trout_Lake", "Trout_Bog", "Two_Sisters_Lake", "Washington", "Wingra"],
           "CH": ["Biel", "Lower_Zurich", "Neuchatel"],
           "PT": ["Alqueva"],
           "FR": ["Annecy", "Bourget", "Geneva"],
           "AU": ["Argyle", "Burley_Griffin", "Mt_Bold"],
           "CA": ["Dickie_Lake", "Eagle_Lake", "Harp_Lake"],
           "SE": ["Ekoln_basin_of_Malaren", "Erken"],
           "UK": ["Esthwaite_Water", "Windermere"],
           "IE": ["Feeagh"],
           "FI": ["Kilpisjarvi", "Kuivajarvi", "Paajarvi"],
           "IL": ["Kinneret"],
           "RW": ["Kivu"],
           "CZ": ["Klicava", "Rimov", "Zlutice"],
           "NO": ["Langtjern"],
           "RU": ["Mozaisk", "Vendyurskoe"],
           "DE": ["Mueggelsee", "Rappbode_Reservoir", "Stechlin"],
           "CN": ["Ngoring"],
           "EE": ["Nohipalo_Mustjarv", "Nohipalo_Valgejarv", "Vortsjarv"],
           "ES": ["Sau_Reservoir"],
           "NZ": ["Rotorua", "Tarawera", "Taupo", "Waahi"]}

full_lake_list = ['Annie', 'Allequash_Lake', 'Alqueva', 'Annecy', 'Argyle',
                  'Biel', 'Big_Muskellunge_Lake', 'Black_Oak_Lake', 'Bourget', 'Burley_Griffin',
                  'Crystal_Bog', 'Crystal_Lake', 'Delavan', 'Dickie_Lake', 'Eagle_Lake',
                  'Ekoln_basin_of_Malaren', 'Erken', 'Esthwaite_Water', 'Falling_Creek_Reservoir', 'Feeagh',
                  'Fish_Lake', 'Geneva', 'Great_Pond', 'Green_Lake', 'Harp_Lake',
                  'Kilpisjarvi', 'Kinneret', 'Kivu', 'Kuivajarvi', 'Langtjern',
                  'Laramie_Lake', 'Lower_Zurich', 'Mendota', 'Monona', 'Mozaisk',
                  'Mt_Bold', 'Mueggelsee', 'Neuchatel', 'Ngoring', 'Nohipalo_Mustjarv', 'Nohipalo_Valgejarv',
                  'Okauchee_Lake', 'Paajarvi', 'Rappbode_Reservoir', 'Rimov', 'Rotorua',
                  'Sammamish', 'Sau_Reservoir', 'Sparkling_Lake', 'Stechlin', 'Sunapee',
                  'Tahoe', 'Tarawera', 'Toolik_Lake', 'Trout_Bog', 'Trout_Lake',
                  'Two_Sisters_Lake', 'Vendyurskoe', 'Vortsjarv', 'Washington', 'Windermere',
                  'Wingra']

temp_by_lake = {'Annie': [1, 17], 'Allequash_Lake': [1, 6], 'Alqueva': [1, 60], 'Annecy': [10, 61],
                "Argyle": [0.7, 26.7],
                'Biel': [5, 70], 'Big_Muskellunge_Lake': [1, 18], 'Black_Oak_Lake': [0.9144, 24.384],
                'Bourget': [9.5, 106], 'Burley_Griffin': [3, 12],
                'Crystal_Bog': [0.5, 1.5], 'Crystal_Lake': [1, 17], 'Delavan': [1, 15], 'Dickie_Lake': [1, 9],
                'Eagle_Lake': [1, 10],
                'Ekoln_basin_of_Malaren': [1, 27.25], 'Erken': [0.5, 20], 'Esthwaite_Water': [0.5, 11.5],
                'Falling_Creek_Reservoir': [1, 8], 'Feeagh': [0.9, 42],
                'Fish_Lake': [1, 17], 'Geneva': [3.6, 263.6], 'Great_Pond': [1, 19], 'Green_Lake': [1, 66],
                'Harp_Lake': [1, 27],
                'Kilpisjarvi': [0.75, 31], 'Kinneret': [2, 34], 'Kivu': [0, 90], 'Kuivajarvi': [0.5, 8],
                'Langtjern': [1, 8],
                'Laramie_Lake': [0.8, 6], 'Lower_Zurich': [1, 100], 'Mendota': [0, 20], 'Monona': [1, 18],
                'Mozaisk': [1, 7],
                'Mt_Bold': [1, 34], 'Mueggelsee': [1, 5], 'Neuchatel': [0, 100], 'Ngoring': [4.5, 19.5],
                'Nohipalo_Mustjarv': [0.5, 5],
                'Nohipalo_Valgejarv': [0.5, 7], 'Okauchee_Lake': [1, 25], 'Paajarvi': [1, 60],
                'Rappbode_Reservoir': [2.82, 52.62], 'Rimov': [0.2, 14.2],
                'Rotorua': [2, 16], 'Sammamish': [1, 20], 'Sau_Reservoir': [2, 40], 'Sparkling_Lake': [0, 17],
                'Stechlin': [0, 55],
                'Sunapee': [1, 9], 'Tahoe': [1.5, 168.1], 'Tarawera': [2, 80], 'Toolik_Lake': [1, 16],
                'Trout_Bog': [1, 7],
                'Trout_Lake': [1, 29], 'Two_Sisters_Lake': [4.572, 15.24], 'Vendyurskoe': [4.58, 11.3],
                'Vortsjarv': [0.5, 3], 'Washington': [1, 55],
                'Windermere': [1, 35], 'Wingra': [1, 3]}

variable_pos = {'T': 5, 'O2': 6, 'Chl': 13}


### Matlab files reading function ###
def load_data(f, sediment=0):
    try:
        mat_contents = sio.loadmat(f)
    except NotImplementedError:
        import hdf5storage
        mat_contents = hdf5storage.loadmat(f)
    MyLake_results = mat_contents['MyLake_results']

    if sediment == 1:
        Sediment_results = mat_contents['Sediment_results']
        return MyLake_results, Sediment_results
    return MyLake_results, []


### Function calculating Performance Index ###
def percent_bias_pbias(obs_list, sims_list):
    """
    Finds the sums of squares for all temperatures listed in the comparison file.
    :param obs_list: A list of observed temperatures.
    :param sims_list: A list of simulated temperatures.
    :return: The result of the sums of squares as a float.
    """
    sums = 0
    for x in range(len(obs_list)):
        if obs_list[x] == 'None':
            continue
        sums += (float(obs_list[x]) - float(sims_list[x]))
    if sum(obs_list) != 0:
        percent = sums * 100 / sum(obs_list)
    else:
        percent = 0
    return percent


def nash_sutcliffe_efficiency_nse(obs_list, sims_list):
    """
    Finds the sums of squares for all temperatures listed in the comparison file.
    :param obs_list: A list of observed temperatures.
    :param sims_list: A list of simulated temperatures.
    :return: The result of the sums of squares as a float.
    """
    d = pd.DataFrame(list(zip(obs_list, sims_list)),
                     columns=['obs', 'sim'])

    d["obs"] = d["obs"].astype(float)
    d["sim"] = d["sim"].astype(float)

    if 1 == 1:  # try:
        sums = 0
        st = 0
        mean = d["obs"].mean()
        for x in range(len(obs_list)):
            if obs_list[x] == 'None':
                continue
            sums += (float(obs_list[x]) - float(sims_list[x])) ** 2
            st += (float(obs_list[x]) - mean) ** 2

        NSE = 1 - (sums / st)
    else:  # except:
        NSE = np.nan

    return NSE


def root_mean_square(obs_list, sims_list):
    """
    Finds the root_mean_square for the temperatures listed in the comparison file.
    :param obs_list: A list of observed temperatures.
    :param sims_list: A list of simulated temperatures.
    :return: The result of the root mean square as a float.
    """

    d = pd.DataFrame(list(zip(obs_list, sims_list)),
                     columns=['obs', 'sim'])

    d["obs"] = d["obs"].astype(float)
    d["sim"] = d["sim"].astype(float)

    if 1 == 1:  # try:
        # results = mean_squared_error(d["obs"], d["sim"])
        # results = mean_squared_error(sqare)
        # results = sqrt(results)
        # resultsnormalise = sqrt(mean_squared_error(d["obs"], d["sim"])) / (max(d["obs"]) - min(d["obs"]))
        results = sqrt(mean_squared_error(d["obs"], d["sim"]))
        resultsnormalise = results / standard_deviation(obs_list)
        mino = np.min(d["obs"])
        maxo = np.max(d["obs"])
        difference = maxo - mino
        normalisermse = results / difference

    else:  # except:
        results = np.nan
        resultsnormalise = np.nan
        normalisermse = np.nan
    return results, resultsnormalise, normalisermse


def r_squared(obs_list, sims_list):
    """
    Find the R squared for the simulations compared to the expected observations
    :param obs_list: A list of observed temperatures.
    :param sims_list: A list of simulated temperatures.
    :return: results of R squared, as float
    """

    x = []
    y = []

    for i in obs_list:
        try:

            x.append(float(i))
            y.append(float(sims_list[obs_list.index(i)]))

        except ValueError:
            continue
        except IndexError:
            break
    try:
        rsquare = r2_score(x, y)
    except:
        rsquare = np.nan
    return rsquare


def sums_of_squares(obs_list, sims_list):
    """
    Finds the sums of squares for all temperatures listed in the comparison file.
    :param obs_list: A list of observed temperatures.
    :param sims_list: A list of simulated temperatures.
    :return: The result of the sums of squares as a float.
    """
    sums = 0
    for x in range(len(obs_list)):
        if obs_list[x] == 'None':
            continue
        sums += (float(obs_list[x]) - float(sims_list[x])) ** 2

    return sums


def standard_deviation(obs_list):
    """
    Find the standard deviation of the observations
    :param obs_list: Type list. The list of observed temperatures.
    :return: The standard deviation of obs_list
    """
    observations = []
    for obs in obs_list:
        try:
            observations.append(float(obs))
        except ValueError:
            continue

    return statistics.stdev(observations)


def rmse_by_sd(obs_list, rmse):
    """
    Divides RMSE of the simulations by the SD of the observations
    :param obs_list: A list of observed temperatures.
    :param rmse: Float
    :return: A float, RMSE / SD
    """
    try:
        results = rmse / standard_deviation(obs_list)
    except ZeroDivisionError:
        results = "Error_Zero_Division"
    return results


### Fuction related to User input ###
def ask_for_matlab_directory():
    """
        function to get the exact matlab directory needed to run the model.
        function is case sensible.
        :return: the matlab directory
        """
    if os.path.exists(r"C:\Program Files\MATLAB\R2019b\bin"):
        print(r"matlab.exe used: C:\Program Files\MATLAB\R2019b\bin\matlab")
        return r"C:\Program Files\MATLAB\R2019b\bin\matlab"
    elif os.path.exists(r"C:\Program Files\MATLAB\R2019a\bin"):
        print(r"matlab.exe used: C:\Program Files\MATLAB\R2019a\bin\matlab")
        return r"C:\Program Files\MATLAB\R2019a\bin\matlab"
    else:
        while True:
            directory = input(r"Enter path to matlab.exe (ex: C:\Program Files\MATLAB\R2019b\bin\matlab) : ")
            if os.path.exists(r"%s" % (directory)):
                print("Path enter does not exist or does not include the '\matlab'. Enter valid path.\n")
                continue
            else:
                # user enter correct lake name
                break

        return r"%s" % (directory)


def aks_for_which_lake_wanted_to_calibrated():
    """
    function to get in input the lake that will be calibrated. For now, options are "Nedre" or "Ovre".
    function is not case sensible.
    :return: the lake name ("Nedre" or "Ovre")
    """
    while True:
        lake_name = input("Enter which lake will be calibrated 'Nedre' or 'Ovre' (not case sensitive): ")
        if lake_name.upper() not in ("NEDRE", "OVRE"):
            print("Lake name giving is not an option, choose between 'Nedre' and 'Ovre'.\n")
            continue
        else:
            # user enter correct lake name
            break

    return "%s%s" % (lake_name[0].upper(), lake_name[1:].lower())


def aks_for_what_is_modeled():
    """

    :return: what_is_calibrated, what_variable_is_calibrated
    """
    while True:
        what_is_calibrated = input(
            "\nType the number of what you want to calibrated:\nThe water column (1), The sediment module (2), or Both (3): ")
        if what_is_calibrated not in ("1", "2", "3"):
            print("%s is not a option, please select '1', '2' or '3'.\n" % what_is_calibrated)
            continue
        else:
            # user enter correct lake name
            what_is_calibrated = int(what_is_calibrated)
            correct = True
            if what_is_calibrated == 1:
                while True:
                    what_variable_is_calibrated = input(
                        "\nType the number of what variable you want to calibrated:\nTemperature (1), Oxygen concentration (2), Chlorophylle (3) or all (4): ")
                    if what_variable_is_calibrated not in ("1", "2", "3", "4"):
                        print("%s is not a option, please select '1', '2', '3' or '4'.\n" % what_variable_is_calibrated)
                        continue
                    else:
                        # user enter correct lake name
                        what_variable_is_calibrated = int(what_variable_is_calibrated)
                        if what_variable_is_calibrated in [1,2,3,4]:
                            break
                        else:
                            print(
                                "For now, the calibration of other variable than Temperature is not supported, please select option '1'")
                            correct_variable = False
                        if correct_variable:
                            correct = True
                            break
                        else:
                            continue
            else:
                what_variable_is_calibrated = 2
            #     print("For now, the calibration of the sediment is not supported, please select option '1'\n")
            #     correct = False
            if correct:
                break
            else:
                continue

    return what_is_calibrated, what_variable_is_calibrated


def ask_what_to_save():
    """
    function to get in input the lake that will be calibrated. For now, options are "Nedre" or "Ovre".
    function is not case sensible.
    :return: the lake name ("Nedre" or "Ovre")
    """
    default = [True, True, False, False]
    while True:
        print(
            "Saving option are:\nSaving calibration report? %s\nSaving calibration figure? %s\nSaving all output data? %s\nSaving comparison data? %s\n" % (
            default[0], default[1], default[2], default[3]))
        saving_option = input("Do you want to continue with the default saving option (y or n)? ")
        if saving_option.upper() not in ("N", "Y"):
            print("answer giving is not an option, choose between 'y' and 'n'.\n")
            continue
        else:
            if saving_option.upper() == "N":
                while True:
                    report = input("Do you want to save a report of the calibration (y or n)? ")
                    if report.upper() not in ("N", "Y"):
                        print("answer giving is not an option, choose between 'y' and 'n'.\n")
                        continue
                    else:
                        if report.upper() == "N":
                            default[0] = False
                        else:
                            default[0] = True

                        break
                while True:
                    figure = input("Do you want to save the figures of the calibration (y or n)? ")
                    if figure.upper() not in ("N", "Y"):
                        print("answer giving is not an option, choose between 'y' and 'n'.\n")
                        continue
                    else:
                        if figure.upper() == "N":
                            default[1] = False
                        else:
                            default[1] = True

                        break
                while True:
                    output = input("Do you want to save the output the calibration (y or n)? ")
                    if output.upper() not in ("N", "Y"):
                        print("answer giving is not an option, choose between 'y' and 'n'.\n")
                        continue
                    else:
                        if output.upper() == "N":
                            default[2] = False
                        else:
                            default[2] = True

                        break
                while True:
                    comparison = input("Do you want to save the comparison files of the calibration (y or n)? ")
                    if comparison.upper() not in ("N", "Y"):
                        print("answer giving is not an option, choose between 'y' and 'n'.\n")
                        continue
                    else:
                        if comparison.upper() == "N":
                            default[3] = False
                        else:
                            default[3] = True

                        break

            break

    return default


### General Function needed by method of Lake CLass ###
def variables_by_depth(observation_folder, lakeName, output_folder, variable='T'):
    """
    Creates a new csv file with the observed temperatures separated in columns by depths.
    :param observation_folder: String
    :param lakeName: String
    :param output_folder: String
    :return: None
    """
    test = "{}/{}_bathymetry.csv".format(observation_folder, lakeName)

    with open("{}/{}_bathymetry.csv".format(observation_folder, lakeName)) as bathymetric_file:
        maxDepth = int(float(list(csv.reader(bathymetric_file))[-1][1]))
        depthLevels = list(range(maxDepth + 1))

    with open("{}/{}_observation_data.csv".format(observation_folder, lakeName)) as obs_file:
        reader = list(csv.reader(obs_file))[1:]
        depthlist = []
        loop = 0
        for depth in reader:
            if float(depth[3]) not in depthlist:
                depthlist.append(float(depth[3]))
            elif float(depth[3]) == depthlist[0]:
                loop += 1
            if loop == 1000:
                break
        outputdir2 = list(output_folder.split("/"))
        outputdir3 = 'Postproc_code'
        for i in outputdir2:
            if i == 'output':
                outputdir3 = i
            else:
                if not os.path.exists(outputdir3):
                    os.mkdir(outputdir3)
                outputdir3 = os.path.join(outputdir3, i)
        if not os.path.exists(outputdir3):
            os.mkdir(outputdir3)

        depthlist.sort()

        with open("{}/Observed_{}.csv".format(observation_folder, variable), "w", newline='') as csvfile:
            print("{}/Observed_{}.csv".format(observation_folder, variable))
            header = "{}, {}\n".format(np.nan, depthlist)
            csvfile.write(header.translate({ord(i): None for i in '[]'}))

            out = csv.writer(csvfile)
            rows = {}
            dates = []
            for i in depthlist:
                rows[i] = []

            for observation in reader:
                obs_date = observation[2]  # "%s%s%s"%(observation[2][0:3],observation[2][4:5],observation[2][5:])
                if obs_date not in dates:
                    dates.append(obs_date)
                    # print(obs_date)

            temp_list = []
            number = 0
            for date in dates:
                for observation in reader:
                    # print(int(observation[2]), date)
                    if str(observation[2]) == str(date):
                        temp_list.append(observation)

                for depth in depthlist:

                    missing_temp = True

                    for t in temp_list:

                        if float(t[3]) == depth:

                            if len(rows[depth]) <= number:
                                if t[variable_pos[variable]] == "":
                                    rows[depth].append("None")
                                else:
                                    rows[depth].append(float(t[variable_pos[variable]]))

                                missing_temp = False

                    if missing_temp:
                        rows[depth].append("None")
                number += 1
                temp_list.clear()
            temp_list.clear()
            for date in dates:
                temp_list.append(date)
                for x in depthlist:
                    temp_list.append(rows[x][dates.index(date)])
                out.writerow(temp_list)
                temp_list.clear()

    print("observation done ... ... ... ... ")


def get_dates_of_simulation(start_year, stop_year):
    """
    Finds the dates for each day of simulation. The dates are found by beginning at the first of January of the given
    start year and adding a day until the 31st of December of the stop year, accounting for leap years.
    :param start_year: An integer, the chosen year for the beginning of the simulation.
    :param stop_year: An integer, the chosen year for the end of the simulation.
    :return: A list of all the dates of simulation, in order from first to last. Dates are integers in the form YYYYMMDD.
    """
    date_list = []
    year = start_year
    nb_year = stop_year - start_year
    if nb_year == 0:
        nb_year = 1
    for i in range(0, nb_year):
        if i % 400 == 0 or (i % 4 == 0 and i % 100 != 0):
            for x in range(1, 367):
                date = 0
                str_x = str(x)

                if x <= 31:
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "01" + str_x)
                elif 31 < x <= 60:
                    str_x = str(int(str_x) - 31)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "02" + str_x)
                elif 60 < x <= 91:
                    str_x = str(int(year) - 60)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "03" + str_x)
                elif 91 < x <= 121:
                    str_x = str(int(year) - 91)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "04" + str_x)
                elif 121 < x <= 152:
                    str_x = str(int(str_x) - 121)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "05" + str_x)
                elif 152 < x <= 182:
                    str_x = str(int(str_x) - 152)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "06" + str_x)
                elif 182 < x <= 213:
                    str_x = str(int(str_x) - 182)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "07" + str_x)
                elif 213 < x <= 243:
                    str_x = str(int(str_x) - 213)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "08" + str_x)
                elif 243 < x <= 274:
                    str_x = str(int(str_x) - 243)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "09" + str_x)
                elif 274 < x <= 305:
                    str_x = str(int(str_x) - 274)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "10" + str_x)
                elif 305 < x <= 335:
                    str_x = str(int(str_x) - 305)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "11" + str_x)
                elif 335 < x <= 366:
                    str_x = str(int(str_x) - 335)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "12" + str_x)

                date_list.append(date)
        else:
            for x in range(1, 366):
                date = 0
                str_x = str(x)

                if x <= 31:
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "01" + str_x)
                elif 31 < x <= 59:
                    str_x = str(int(str_x) - 31)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "02" + str_x)
                elif 59 < x <= 90:
                    str_x = str(int(str_x) - 59)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "03" + str_x)
                elif 90 < x <= 120:
                    str_x = str(int(str_x) - 90)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "04" + str_x)
                elif 120 < x <= 151:
                    str_x = str(int(str_x) - 120)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "05" + str_x)
                elif 151 < x <= 181:
                    str_x = str(int(str_x) - 151)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "06" + str_x)
                elif 181 < x <= 212:
                    str_x = str(int(str_x) - 181)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "07" + str_x)
                elif 212 < x <= 242:
                    str_x = str(int(str_x) - 212)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "08" + str_x)
                elif 242 < x <= 273:
                    str_x = str(int(str_x) - 242)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "09" + str_x)
                elif 273 < x <= 304:
                    str_x = str(int(str_x) - 273)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "10" + str_x)
                elif 304 < x <= 334:
                    str_x = str(int(str_x) - 304)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "11" + str_x)
                elif 334 < x <= 365:
                    str_x = str(int(str_x) - 334)
                    if len(str_x) == 1:
                        str_x = "0" + str_x
                    date = int(str(year) + "12" + str_x)

                date_list.append(date)
        year += 1
    return date_list


def findYPoint(xa, xb, ya, yb, xc):
    m = (float(ya) - float(yb)) / (float(xa) - float(xb))
    yc = (float(xc) - (float(xb))) * m + float(ya)
    return yc


def date_range(start, end):
    delta = end - start  # as timedelta
    days = [start + timedelta(days=i) for i in range(delta.days + 1)]
    days_str = [int(day.strftime("%Y%m%d")) for day in days]
    return days_str


### General Function related to Figures (Used by Class Graphic) ###
def timeline_plot(modeleddata: list, modeleddates: list, observeddata: list = None, observeddates: list = None, ax=None,
                  ylimit=[-0.5, 30],
                  line_kwargs: dict = {},
                  sct_kwargs: dict = {}, linestyle="-"):
    """
    plot timeline with modeled data (line) and observed data (dot)

    :param modeleddata: modeled data for each day included in years presenting observed data (at the same depth)
    :param modeleddates: all days included in the period modeled
    :param observeddata: Observed data to be plotted
    :param observeddates: Dates Where observed data have been measured
    :param ax: Axe where you want to plot. If None, look for the last ax used in the current figure
    :param line_kwargs: Other arguments given to the line plot (measures' plot)
    :param sct_kwargs: Other arguments given to the scatterplot (Observations' plot)
    :return: None
    """
    if ax is None:
        ax = plt.gca()
    if modeleddates is None or modeleddata is None:
        raise TypeError

    sns.lineplot(x=modeleddates, y=modeleddata, ax=ax, **line_kwargs)
    if observeddata is not None:
        sns.scatterplot(x=observeddates, y=observeddata, ax=ax, **sct_kwargs)

    ax.set_xlabel("Dates")
    ax.set_xlim(min(modeleddates), max(modeleddates))
    ax.set_ylim(ylimit[0], ylimit[1])


def line_plot(lineStart=None, lineEnd=None, ax=None, linearg={}):
    if ax is None:
        ax = plt.gca()
    if lineStart is None:
        lineStart = 0
        lineEnd = 20

    ax.plot([lineStart, lineEnd], [lineStart, lineEnd], **linearg)
    ax.set_xlim(lineStart, lineEnd)
    ax.set_ylim(lineStart, lineEnd)


def linear_regression_plot(x2, y2, ax=None, linearregressionarg={}, confidentintervalarg={}, predictionintervalarg={}):
    if ax is None:
        ax = plt.gca()

    x = smodels.add_constant(x2)  # constant intercept term
    # Model: y ~ x + c
    model = smodels.OLS(y2, x)
    fitted = model.fit()
    x_pred = np.linspace(x.min(), x.max(), 50)
    x_pred2 = smodels.add_constant(x_pred)
    y_pred = fitted.predict(x_pred2)

    ax.plot(x_pred, y_pred, **linearregressionarg)

    # print(fitted.params)  # the estimated parameters for the regression line
    # print(fitted.summary())  # summary statistics for the regression

    y_hat = fitted.predict(x)  # x is an array from line 12 above
    y_err = y2 - y_hat
    mean_x = x.T[1].mean()
    n = len(x)
    dof = n - fitted.df_model - 1

    t = stats.t.ppf(1 - 0.025, df=dof)
    s_err = np.sum(np.power(y_err, 2))
    conf = t * np.sqrt((s_err / (n - 2)) * (1.0 / n + (
                np.power((x_pred - mean_x), 2) / ((np.sum(np.power(x_pred, 2))) - n * (np.power(mean_x, 2))))))
    upper = y_pred + abs(conf)
    lower = y_pred - abs(conf)

    # Prediction Interval
    sdev, lower, upper = wls_prediction_std(fitted, exog=x_pred2, alpha=0.025)
    ax.fill_between(x_pred, lower, upper, **predictionintervalarg)
    # ax.plot(x_pred, lower, **predictionintervalarg, alpha=1, linestyle='-.')
    # ax.plot(x_pred, upper, **predictionintervalarg, alpha=1, linestyle='-.')
    # ax.legend(bbox_to_anchor=(1,0), loc="lower right")

    # Confidence Interval
    # ax.fill_between(x_pred, lower, upper, **confidentintervalarg)
    # ax.plot(x_pred, lower, **confidentintervalarg, alpha=1, linestyle='-')
    # ax.plot(x_pred, upper, **confidentintervalarg, alpha=1, linestyle='-')
    sns.regplot(x2, y2, ci=95, scatter_kws={"color": "white", "alpha": 1, "zorder": -10}, line_kws=confidentintervalarg,
                ax=ax)


def error_bar_plot(x, y, xerr, yerr, ax=None, errorbararg={}, markerwidth=4):
    if ax is None:
        ax = plt.gca()

    (_, caps, _) = plt.errorbar(x, y, xerr=xerr, yerr=yerr, **errorbararg)
    for cap in caps:
        cap.set_markeredgewidth(markerwidth)


def base_plot_comparison(x, y, lineStart=None, lineEnd=None, ax=None, linecolor="k", lake=""):
    if ax is None:
        ax = plt.gca()
    linearg = {'color': linecolor, 'label': "y= x", 'linewidth': 4, 'linestyle': '--'}
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    rmse, rsr, nrmse = root_mean_square(x, y)
    linearregressionarg = {'color': 'k', 'linewidth': 1.5,
                           "label": "linear regression (y = %0.3f x + %0.3f) \n R\u00b2 : %0.3f RMSE: %0.3f" % (
                           slope, intercept, r_value, rmse)}
    confidenceintervalarg = {'color': 'k', 'linewidth': 1.5, 'label': "Confidence interval",
                             "zorder": 5}  # 'alpha': 0.4,
    predictionintervalarg = {'color': '#888888', 'alpha': 0.1, 'label': "Prediction interval",
                             "zorder": 0}  # 'alpha': 0.1,

    line_plot(lineStart=lineStart, lineEnd=lineEnd, ax=ax, linearg=linearg)

    linear_regression_plot(x, y, ax, linearregressionarg, confidenceintervalarg, predictionintervalarg)
    legendNames = ["Confidence interval", "Prediction interval"]
    # ax.legend(legendNames,loc='bottom right')

    ax.set_ylim()
    ax.set_xlim()
    ax.text(lineStart + 0.5, lineEnd - 1, "R\u00b2 : %0.3f RMSE: %0.3f" % (r_value, rmse),
            horizontalalignment='left',
            verticalalignment='top')
    return r_value, rmse


class Lake:
    """
    A class used to gather all information related to the selected lake

    ...

    Attributes
    ----------

    lake_name : str
        The name of the lake (For now, should be "Nedre" or "Ovre", error otherwise)


    Methods
    -------
    manual_calibration_loop(report = True, save_figures = False, save_output_data = False, save_comparison_data = False)
        Prints the animal's name and what sound it makes
    """

    def __init__(self, lake_name):
        """
        Parameters
        ----------
        lake_name : str
            The name of the lake (For now, should be "Nedre" or "Ovre", error otherwise)

        observation_folder : str
            folder directory for where the formatted observations are.

        input_folder : str
            folder directory for where the formatted inputs are.

        output_folder : str
            folder directory for where the data generated are.

        observation_folder : str
            folder directory for where the figures generated by the functions are.

        """

        self.name = lake_name
        self.observation_folder = r"obs/%s" % lake_name
        self.save_date = datetime.now().strftime('%Y%m%d')

        # Default value for parameters related to Temperature
        self.kz_N0 = 0.00007
        self.c_shelter = "NaN"
        self.i_scv = 1
        self.i_sct = 0
        self.swa_b0 = 2.5
        self.swa_b1 = 1

        # Default value for parameters related to Oxygen
        self.I_scDOC = 1
        self.I_scO = 1

        # Default value for parameter related to Chl_a
        self.I_scChl = 1
        self.k_Chl = 0.4

        # Default value for parameter related to sediment
        self.k_POP =  0.04
        self.k_POC =  0.02
        self.k_DOP = 0.04
        self.k_DOC = 0.02
        self.k_pdesorb_a = 100
        self.k_pdesorb_b = 100

        self.input_folder = r"IO/%s" % lake_name
        self.output_folder = r"Postproc_code/%s" % lake_name
        self.figures_folder = r"Postproc_code/%s/figures" % lake_name

        # generate observed data
        for variable in ["T", "O2", "Chl"]:
            variables_by_depth(self.observation_folder, self.name, self.output_folder, variable)

    def manual_calibration_loop(self, report=True, save_figures=False, save_output_data=False,
                                save_comparison_data=False, matlab=matlab_folder):
        """ Main function, loop through iteration of simulation of the temperature, oxygen, and chl_a for the selected lake.

        This function calls the function asking for parameters' values, launches the simulation with those values and
        generate a figure with those values.
        The function will be by default

        Parameters
        ----------
        sound : str, optional
            The sound the animal makes (default is None)

        Raises
        ------
        NotImplementedError
            If no sound is set for the animal or passed in as a
            parameter.
        """
        start_time = datetime.now()
        self.save_date = start_time.strftime('%Y%m%d_%H%M')

        dict_variable = [{1: "T", 2: "O2", 3: "Chl", 4: "O2"}, {1: "T", 2: "O2", 3: "Chl", 4: "O2"},
                         {1: "T", 2: "O2", 3: "Chl", 4: "O2"}]
        outputdir = os.path.join(self.output_folder, "save_output_all")
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)

        what_is_calibrated, what_variable_is_calibrated = aks_for_what_is_modeled()
        if what_is_calibrated in [2, 3]:
            enable_sediment = 1
        else:
            enable_sediment = 0
        dict_variable = dict_variable[what_is_calibrated - 1]

        if report:
            report_core_title = "manual_calibration_%s_%s_%s.csv" % (
            self.name, dict_variable[what_variable_is_calibrated], start_time.strftime('%Y%m%d_%H%M'))
            table_report = [
                ["Iteration", "Lake", "Variable_calibrated", "kz_N0", "c_shelter", "i_scv", "i_sct", "swa_b0",
                 "swa_b1", "I_scDOC", "RMSE_all", "R2_all", "RMSE", "NSE", "RSR", "Pbias", "R2",
                 "SOS", "nrmse", "M_score", "Note", "time"]]

        iteration_number = 0
        iteration_continue = True
        while True:
            if iteration_number != 0:
                while True:
                    stop_iteration = input("Continue to test value for manual calibration (y or n)? ")
                    if stop_iteration.upper() not in ('Y', 'N'):
                        print("answer giving is not an option, choose between 'y' and 'n'.\n")
                        continue
                    else:
                        if stop_iteration.upper() == 'N':
                            iteration_continue = False
                        break

            if iteration_continue:
                iteration_number += 1
                start_iteration = datetime.now()
                print("\n**** Start Iteration %s %s ****" % (iteration_number, start_iteration))

                if self.ask_parameters_value(what_is_calibrated, what_variable_is_calibrated) == 0:
                    return 0

                print("\nStart MyLake run with kz_N0 = %s, c_shelter = %s, i_scv = %s, i_sct = %s, swa_b0 = %s, "
                      "swa_b1 = %s, I_scDOC = %s, I_scO = %s, I_scChl  = %s, k_Chl = %s, k_POP = %s, k_POC = %s, "
                      "k_DOP = %s, k_DOC = %s, k_pdesorb_a = %s, k_pdesorb_b = %s" % (self.kz_N0, self.c_shelter,
                                                                                      self.i_scv, self.i_sct,
                                                                                      self.swa_b0, self.swa_b1,
                                                                                      self.I_scDOC,self.I_scO, self.I_scChl ,self.k_Chl,
                                                                                      self.k_POP, self.k_POC,
                                                                                      self.k_DOP, self.k_DOC,
                                                                                      self.k_pdesorb_a,
                                                                                      self.k_pdesorb_b))
                while True:
                    stop_iteration = input("Is that right (y or n)? ")
                    if stop_iteration.upper() not in ('Y', 'N'):
                        print("answer giving is not an option, choose between 'y' and 'n'.\n")
                        continue
                    else:
                        if stop_iteration.upper() == 'N':
                            if self.ask_parameters_value(what_is_calibrated, what_variable_is_calibrated) == 0:
                                return 0
                        break

                # Run MyLake
                if self.c_shelter == str(self.c_shelter):
                    cmd = r'%s -wait -r -nosplash -nodesktop MyLake_Milsbo_run(%d,%d,%s,%f,%s,%f,%f,%f,%f,%f,%f,%f,%f,' \
                          r'%f,%f,%f,%f,%f,%f,%d);quit' % ( '"%s"' % matlab, 2018, 2020, "'%s'" % self.name,
                                                            self.kz_N0, "'%s'" % self.c_shelter, self.i_scv,
                                                            self.i_sct, self.swa_b0, self.swa_b1, self.I_scDOC,self.I_scO,
                                                            self.I_scChl, self.k_Chl, self.k_POP, self.k_POC, self.k_DOP,
                                                            self.k_DOC, self.k_pdesorb_a, self.k_pdesorb_b,
                                                            enable_sediment)
                else:
                    cmd = r'%s -wait -r -nosplash -nodesktop MyLake_Milsbo_run(%d,%d,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,' \
                          r'%f,%f,%f,%f,%f,%f,%d);quit' % ('"%s"' % matlab, 2018, 2020, "'%s'" % self.name,
                                                           self.kz_N0, self.c_shelter, self.i_scv, self.i_sct,
                                                           self.swa_b0, self.swa_b1, self.I_scDOC,self.I_scO,self.I_scChl, self.k_Chl,
                                                           self.k_POP, self.k_POC, self.k_DOP, self.k_DOC,
                                                           self.k_pdesorb_a, self.k_pdesorb_b, enable_sediment)
                print("Run MyLake model with parameter\n" + cmd)
                try:
                    os.system(cmd)
                    print("run MyLake sucess")
                except:
                    print('error with matlab')
                    return 1

                if save_output_data:
                    outputdir2 = os.path.join(outputdir, "raw_output")
                    if not os.path.exists(outputdir2):
                        os.mkdir(outputdir2)
                    shutil.copy2(os.path.join(self.output_folder, "%s_result_run.mat" % self.name),
                                 os.path.join(outputdir2, "%s_%s_ite%s_result_run.mat" % (
                                 self.name, self.save_date, iteration_number)))

                print("Performance calcul")
                mat_data = load_data('%s/%s_result_run.mat' % (self.output_folder, lake_name), enable_sediment)
                for comp_variable in list(dict_variable.keys())[0:3]:
                    self.make_comparison_file_allobs(dict_variable[comp_variable], enable_sediment, mat_data=mat_data)

                if save_comparison_data:
                    outputdir2 = os.path.join(outputdir, "comparison_output")
                    if not os.path.exists(outputdir2):
                        os.mkdir(outputdir2)
                    for comp_variable in list(dict_variable.keys())[0:3]:
                        shutil.copy2(
                            os.path.join(self.output_folder, "%s_comparisonall.csv" % dict_variable[comp_variable]),
                            os.path.join(outputdir2, "%s_%s_ite%s_%s_comparisonall.csv" % (
                            self.name, self.save_date, iteration_number, dict_variable[comp_variable])))

                print("Figure generation")
                performances, final_performance, M_score, Note = Graphics(
                    Lake=self).figures_comparison_timeseries_and_profiles(save_figures=save_figures,
                                                                          iteration=iteration_number,
                                                                          variable_calibrated=dict_variable[
                                                                              what_variable_is_calibrated],
                                                                          mat_data=mat_data)

                if report:
                    RMSE_all = [round(float(x), 3) for x in final_performance['rmse']]
                    R2_all = [round(float(x), 3) for x in final_performance['r2']]
                    RMSE = [[performances[0][0][0], performances[0][1][0]],
                            [performances[1][0][0], performances[1][1][0]],
                            [performances[2][0][0], performances[2][1][0]]]
                    NSE = [[performances[0][0][1], performances[0][1][1]],
                           [performances[1][0][1], performances[1][1][1]],
                           [performances[2][0][1], performances[2][1][1]]]
                    RSR = [[performances[0][0][2], performances[0][1][2]],
                           [performances[1][0][2], performances[1][1][2]],
                           [performances[2][0][2], performances[2][1][2]]]
                    Pbias = [[performances[0][0][3], performances[0][1][3]],
                             [performances[1][0][3], performances[1][1][3]],
                             [performances[2][0][3], performances[2][1][3]]]
                    R2 = [[performances[0][0][4], performances[0][1][4]],
                          [performances[1][0][4], performances[1][1][4]],
                          [performances[2][0][4], performances[2][1][4]]]
                    SOS = [[performances[0][0][5], performances[0][1][5]],
                           [performances[1][0][5], performances[1][1][5]],
                           [performances[2][0][5], performances[2][1][5]]]
                    nrmse = [[performances[0][0][6], performances[0][1][6]],
                             [performances[1][0][6], performances[1][1][6]],
                             [performances[2][0][6], performances[2][1][6]]]

                    report_iteration_line = [iteration_number, self.name, what_variable_is_calibrated, self.kz_N0, self.c_shelter, self.i_scv, self.i_sct,
                                                           self.swa_b0, self.swa_b1, self.I_scDOC,self.I_scO,self.I_scChl, self.k_Chl,
                                                           self.k_POP, self.k_POC, self.k_DOP, self.k_DOC,
                                                           self.k_pdesorb_a, self.k_pdesorb_b, RMSE_all, R2_all,
                                             RMSE, NSE, RSR, Pbias, R2, SOS, nrmse, M_score, Note,
                                             str(datetime.now() - start_iteration)]
                    table_report.append(report_iteration_line)

                continue
            else:
                break

        if report:
            df = pd.DataFrame(table_report[1:], columns=table_report[0])
            df.to_csv(os.path.join(self.output_folder, report_core_title))
        return 1

    def ask_parameters_value(self, what_is_calibrated, what_variable_is_calibrated):
        """

        :return:
        """
        if what_is_calibrated == 1:
            if what_variable_is_calibrated == 1:
                print("The parameters value are set to:"
                      "\nTemperature parameters"
                      "\n kz_N0 = %s"
                      "\n c_shelter = %s"
                      "\n i_scv = %s"
                      "\n i_sct = %s"
                      "\n swa_b0 = %s"
                      "\n swa_b1 = %s" % (self.kz_N0, self.c_shelter, self.i_scv, self.i_sct, self.swa_b0, self.swa_b1))
            elif what_variable_is_calibrated == 2:
                print("The parameters value are set to:"
                      "\nOxygen parameters"
                      
                      "\n I_scDOC = %s"
                      "\n I_scO = %s"% ( self.I_scDOC,self.I_scO))

            elif what_is_calibrated == 3:
                print("The parameters value are set to:"
                      "\nChlorophyl parameters"
                      "\n I_scChl = %s"
                      "\n k_Chl = %s" % (self.k_Chl, self.I_scChl))

            else:
                print("The parameters value are set to:"
                      "\nTemperature, Oxygen and Chlorophyl parameters"
                      "\n kz_N0 = %s"
                      "\n c_shelter = %s"
                      "\n i_scv = %s"
                      "\n i_sct = %s"
                      "\n swa_b0 = %s"
                      "\n swa_b1 = %s"
                      "\nOxygen parameters"
                      
                      "\n I_scDOC = %s"
                      "\n I_scO = %s"
                      "\nChlorophyl parameters"
                      "\n I_scChl = %s"
                      "\n k_Chl = %s" % (
                      self.kz_N0, self.c_shelter, self.i_scv, self.i_sct, self.swa_b0, self.swa_b1, self.I_scDOC,self.I_scO,self.I_scChl, self.k_Chl))
                #return 0

        elif what_is_calibrated == 2:
            print("The parameters value are set to:"
                  "\n\nSediment parameters"
                  "\n k_POP = %s"
                  "\n k_POC = %s"
                  "\n k_DOP = %s"
                  "\n k_DOC = %s"
                  "\n k_pdesorb_a = %s"
                  "\n k_pdesorb_b = %s" % (
                  self.k_POP, self.k_POC, self.k_DOP, self.k_DOC, self.k_pdesorb_a, self.k_pdesorb_b))

        else:
            print("The parameters value are set to:"
                  "\n Water Column parameters\n"
                  "\nTemperature, Oxygen and Chlorophyl parameters"
                  "\n kz_N0 = %s"
                  "\n c_shelter = %s"
                  "\n i_scv = %s"
                  "\n i_sct = %s"
                  "\n swa_b0 = %s"
                  "\n swa_b1 = %s"
                  "\nOxygen parameters"
                 
                  "\n I_scDOC = %s"
                  "\n I_scO = %s"
                  "\nChlorophyl parameters"
                  "\n I_scChl = %s"
                  "\n k_Chl = %s"
                  "\nSediment parameters\n"
                  "\n k_POP = %s"
                  "\n k_POC = %s"
                  "\n k_DOP = %s"
                  "\n k_DOC = %s"
                  "\n k_pdesorb_a = %s"
                  "\n k_pdesorb_b = %s" % (
                  self.kz_N0, self.c_shelter, self.i_scv, self.i_sct, self.swa_b0, self.swa_b1,
                  self.I_scDOC,self.I_scO,self.I_scChl, self.k_Chl, self.k_POP, self.k_POC, self.k_DOP, self.k_DOC, self.k_pdesorb_a,
                  self.k_pdesorb_b))
            # return 0

        while True:
            change_parameters = input("\nKeep parameters value set (y or n)? ")
            if change_parameters.upper() not in ("N", "Y"):
                print("Answer giving is not an option, choose between 'y' and 'n'.\n")
                continue
            else:
                if change_parameters.upper() == "N":
                    if what_is_calibrated == 1:
                        if what_variable_is_calibrated == 1:
                            while True:
                                change_parameters = input("\nChange kz_N0 (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.kz_N0 = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new kz_N0: ")
                                            try:
                                                self.kz_N0 = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange c_shelter (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.c_shelter = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new c_shelter: ")
                                            try:
                                                self.c_shelter = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange i_scv (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.i_scv = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new i_scv: ")
                                            try:
                                                self.i_scv = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange i_sct (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.i_sct = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new i_sct: ")
                                            try:
                                                self.i_sct = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange swa_b0 (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.swa_b0 = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new swa_b0: ")
                                            try:
                                                self.swa_b0 = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange swa_b1 (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.swa_b1 = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new swa_b1: ")
                                            try:
                                                self.swa_b1 = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                        elif what_variable_is_calibrated == 2:

                            while True:
                                change_parameters = input("\nChange I_scDOC (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.I_scDOC = float(change_parameters)
                                        break
                                    except:
                                        print(
                                            "Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new I_scDOC: ")
                                            try:
                                                self.I_scDOC = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange I_scO (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.I_scO = float(change_parameters)
                                        break
                                    except:
                                        print(
                                            "Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new I_scDOC: ")
                                            try:
                                                self.I_scO = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                        elif what_variable_is_calibrated == 3:
                            while True:
                                change_parameters = input("\nChange I_scChl (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.I_scChl = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new I_scChl: ")
                                            try:
                                                self.I_scChl = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange k_Chl (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.k_Chl = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new k_Chl: ")
                                            try:
                                                self.k_Chl = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break

                        else:
                            while True:
                                change_parameters = input("\nChange kz_N0 (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.kz_N0 = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new kz_N0: ")
                                            try:
                                                self.kz_N0 = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange c_shelter (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.c_shelter = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new c_shelter: ")
                                            try:
                                                self.c_shelter = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange i_scv (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.i_scv = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new i_scv: ")
                                            try:
                                                self.i_scv = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange i_sct (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.i_sct = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new i_sct: ")
                                            try:
                                                self.i_sct = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange swa_b0 (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.swa_b0 = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new swa_b0: ")
                                            try:
                                                self.swa_b0 = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange swa_b1 (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.swa_b1 = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new swa_b1: ")
                                            try:
                                                self.swa_b1 = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange I_scDOC (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.I_scDOC = float(change_parameters)
                                        break
                                    except:
                                        print(
                                            "Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new I_scDOC: ")
                                            try:
                                                self.I_scDOC = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange I_scO (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.I_scO = float(change_parameters)
                                        break
                                    except:
                                        print(
                                            "Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new I_scDOC: ")
                                            try:
                                                self.I_scO = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange I_scChl (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.I_scChl = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new I_scChl: ")
                                            try:
                                                self.I_scChl = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break
                            while True:
                                change_parameters = input("\nChange k_Chl (y or n or value)? ")
                                if change_parameters.upper() not in ("N", "Y"):
                                    try:
                                        self.k_Chl = float(change_parameters)
                                        break
                                    except:
                                        print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                        continue
                                else:
                                    if change_parameters.upper() == "Y":
                                        while True:
                                            parameter_value = input("Enter new k_Chl: ")
                                            try:
                                                self.k_Chl = float(parameter_value)
                                                break
                                            except:
                                                continue
                                    break


                    elif what_is_calibrated == 2:
                        while True:
                            change_parameters = input("\nChange  k_POP (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.k_POP = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new  k_POP: ")
                                        try:
                                            self.k_POP = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange  k_POC (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.k_POC = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new  k_POC: ")
                                        try:
                                            self.k_POC = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange  k_DOP (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.k_DOP = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new  k_DOP: ")
                                        try:
                                            self.k_DOP = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange  k_DOC (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.k_DOC = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new  k_DOC: ")
                                        try:
                                            self.k_DOC = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange k_pdesorb_a (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.k_pdesorb_a = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new k_pdesorb_a: ")
                                        try:
                                            self.k_pdesorb_a = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange k_pdesorb_b (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.k_pdesorb_b = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new k_pdesorb_b: ")
                                        try:
                                            self.k_pdesorb_b = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                    else:
                        while True:
                            change_parameters = input("\nChange kz_N0 (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.kz_N0 = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new kz_N0: ")
                                        try:
                                            self.kz_N0 = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange c_shelter (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.c_shelter = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new c_shelter: ")
                                        try:
                                            self.c_shelter = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange i_scv (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.i_scv = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new i_scv: ")
                                        try:
                                            self.i_scv = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange i_sct (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.i_sct = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new i_sct: ")
                                        try:
                                            self.i_sct = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange swa_b0 (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.swa_b0 = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new swa_b0: ")
                                        try:
                                            self.swa_b0 = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange swa_b1 (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.swa_b1 = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new swa_b1: ")
                                        try:
                                            self.swa_b1 = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange I_scDOC (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.I_scDOC = float(change_parameters)
                                    break
                                except:
                                    print(
                                        "Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new I_scDOC: ")
                                        try:
                                            self.I_scDOC = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange I_scO (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.I_scO = float(change_parameters)
                                    break
                                except:
                                    print(
                                        "Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new I_scDOC: ")
                                        try:
                                            self.I_scO = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange I_scO (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.I_scO = float(change_parameters)
                                    break
                                except:
                                    print(
                                        "Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new I_scDOC: ")
                                        try:
                                            self.I_scO = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange I_scChl (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.I_scChl = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new I_scChl: ")
                                        try:
                                            self.I_scChl = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange k_Chl (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.k_Chl = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new k_Chl: ")
                                        try:
                                            self.k_Chl = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break

                        while True:
                            change_parameters = input("\nChange  k_POP (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.k_POP = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new  k_POP: ")
                                        try:
                                            self.k_POP = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange  k_POC (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.k_POC = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new  k_POC: ")
                                        try:
                                            self.k_POC = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange  k_DOP (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.k_DOP = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new  k_DOP: ")
                                        try:
                                            self.k_DOP = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange  k_DOC (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.k_DOC = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new  k_DOC: ")
                                        try:
                                            self.k_DOC = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange k_pdesorb_a (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.k_pdesorb_a = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new k_pdesorb_a: ")
                                        try:
                                            self.k_pdesorb_a = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break
                        while True:
                            change_parameters = input("\nChange k_pdesorb_b (y or n or value)? ")
                            if change_parameters.upper() not in ("N", "Y"):
                                try:
                                    self.k_pdesorb_b = float(change_parameters)
                                    break
                                except:
                                    print("Answer giving is not an option, choose between 'y', 'n' or a value.\n")
                                    continue
                            else:
                                if change_parameters.upper() == "Y":
                                    while True:
                                        parameter_value = input("Enter new k_pdesorb_b: ")
                                        try:
                                            self.k_pdesorb_b = float(parameter_value)
                                            break
                                        except:
                                            continue
                                break

                break

        return 1

    def make_comparison_file_allobs(self, variable='T', enable_sediment=0, start_year_comparison=2019, mat_data=[]):
        """
        Search a given output folder for an observation file, containing measured temperatures for a lake on a finite period,
        and a simulated temperatures file. Then writes corresponding observed and simulated temperatures to a CSV file,
        where each column is a list of temperatures for a given depth.
        :param output_folder: A string containing the folder to search and write to.
        :param mat_data     tuple or list of tuple
        :return: None
        """
        dt = 0.5
        if len(mat_data) == 2:
            water, sediment = mat_data[0], mat_data[1]
        else:
            water = mat_data
        # #exemple to create dataframe as it was csv
        # temp = water['T'][0, 0]
        # df = pd.DataFrame(temp).transpose()

        if variable in ['T', "O2", "Chl"]:  # Mean water column comparison
            with open("{}/{}_comparisonall.csv".format(self.output_folder, variable), "w", newline="\n") as file:
                observation_dict = {}
                simulation_dict = {}
                depth_levels = []
                col1, col2 = False, False
                with open("{}/Observed_{}.csv".format(self.observation_folder, variable), "r") as observation_file:
                    reader = list(csv.reader(observation_file))
                    depth_levels = reader[0][1:]
                    for obs in reader[1:]:
                        observation_dict[int(obs[0].replace('-', ''))] = list(obs[1:])

                sims_dates = date_range(datetime(2018, 1, 1), datetime(2021, 1, 1))

                date_format = "%m/%d/%Y"

                a = datetime.strptime('01/01/2018', date_format)
                b = datetime.strptime('01/01/%s' % start_year_comparison, date_format)
                number_skip_data = (b - a).days

                if variable == "T":
                    data_from_mat = pd.DataFrame(water['T'][0, 0]).transpose()
                    reader = data_from_mat.values.tolist()
                else:
                    data = water['concentrations'][0, 0][variable][0, 0]
                    data_from_mat = pd.DataFrame(data).transpose()
                    reader = data_from_mat.values.tolist()

                for y in range(0, len(reader)):
                    temp_at_date = []

                    sim = reader[y]

                    if 1 == 1:  # try:
                        for i in range(0, len(depth_levels)):

                            if float(depth_levels[i]) < 0 or float(depth_levels[i]) > (len(sim) * dt - 0.5):
                                temp_at_date.append(np.nan)

                            elif (float(depth_levels[i])) != int(float(depth_levels[i])) and (
                                    float(depth_levels[i])) != (int(float(depth_levels[i])) + 0.5):
                                if (float(depth_levels[i]) - int(float(depth_levels[i]))) < 0.5:
                                    xa = int(float(depth_levels[i]))
                                    xb = xa + dt
                                else:
                                    xa = int(float(depth_levels[i])) + 0.5
                                    xb = xa + dt

                                if sim[int(xa * (1 / dt))] != "None" and sim[int(xb * (1 / dt))] != "None":
                                    valueyc = findYPoint(xa, xb, sim[int(xa * (1 / dt))], sim[int(xb * (1 / dt))],
                                                         float(depth_levels[i]))
                                    # if valueyc < 0:
                                    #     print("here")
                                    if variable == 'O2':
                                        temp_at_date.append(valueyc*0.001)
                                    else:
                                        temp_at_date.append(valueyc)

                                else:
                                    temp_at_date.append(np.nan)

                            else:
                                depth1 = int(float(depth_levels[i]) * (1 / dt))
                                if sim[depth1] == "None":
                                    temp_at_date.append("None")
                                else:
                                    if variable == "O2":
                                        temp_at_date.append(float(sim[depth1])*0.001)
                                    else:
                                        temp_at_date.append(float(sim[depth1]))

                    simulation_dict[sims_dates[y]] = temp_at_date

                # with open("{}/{}zt.csv".format(self.output_folder, variable), "r") as simulation_file:
                #     reader = list(csv.reader(simulation_file))
                #

                y = 0

                csvFile = csv.writer(file)
                csvFile.writerow(
                    ["Date", "Depth", "Observations", "Simulations"])

                for date in sims_dates[number_skip_data:]:

                    if self.output_folder == "D:\output\PT\Alqueva\EWEMBI\historical":
                        date1 = date + 30000
                    else:
                        date1 = date

                    if date1 in observation_dict.keys():

                        if 1 == 1:  # try:
                            for i in range(0, len(observation_dict[date1])):
                                if str(observation_dict[date1][i])[0].isnumeric() and str(simulation_dict[date][i])[
                                    0].isnumeric():
                                    csvFile.writerow(
                                        [date1, depth_levels[i], observation_dict[date1][i], simulation_dict[date][i]])

                        # except KeyError:
                        #     continue
                return None

        if enable_sediment == 1:  # mean comparison sediment
            with open("{}/{}_comparisonall.csv".format(self.output_folder, variable), "w", newline="\n") as file:
                observation_dict = {}
                simulation_dict = {}
                depth_levels = []
                col1, col2 = False, False
                with open("{}/Observed_{}.csv".format(self.observation_folder, variable), "r") as observation_file:
                    reader = list(csv.reader(observation_file))
                    depth_levels = reader[0][1:]
                    for obs in reader[1:]:
                        observation_dict[int(obs[0].replace('-', ''))] = list(obs[1:])

                start_year = 2019
                end_year = 2021

                start_date = datetime(2019, 1, 1)
                end_date = datetime(2021, 1, 1)

                sims_dates = date_range(start_date, end_date)

                date_format = "%m/%d/%Y"

                a = datetime.strptime('01/01/2019', date_format)
                b = datetime.strptime('01/01/%s' % start_year, date_format)
                number_skip_data = (b - a).days
                with open("{}/{}zt.csv".format(self.output_folder, variable), "r") as simulation_file:
                    reader = list(csv.reader(simulation_file))

                    for y in range(0, len(reader)):
                        temp_at_date = []

                        sim = reader[y]

                        if 1 == 1:  # try:
                            for i in range(0, len(depth_levels)):

                                if float(depth_levels[i]) < 0 or float(depth_levels[i]) > (len(sim) * dt - 0.5):
                                    temp_at_date.append(np.nan)

                                elif (float(depth_levels[i])) != int(float(depth_levels[i])) and (
                                        float(depth_levels[i])) != (int(float(depth_levels[i])) + 0.5):
                                    if (float(depth_levels[i]) - int(float(depth_levels[i]))) < 0.5:
                                        xa = int(float(depth_levels[i]))
                                        xb = xa + dt
                                    else:
                                        xa = int(float(depth_levels[i])) + 0.5
                                        xb = xa + dt

                                    if sim[int(xa * (1 / dt))] != "None" and sim[int(xb * (1 / dt))] != "None":
                                        valueyc = findYPoint(xa, xb, sim[int(xa * (1 / dt))], sim[int(xb * (1 / dt))],
                                                             float(depth_levels[i]))
                                        # if valueyc < 0:
                                        #     print("here")
                                        if variable == 'O2':
                                            temp_at_date.append(valueyc*0.001)
                                        else:
                                            temp_at_date.append(valueyc)

                                    else:
                                        temp_at_date.append(np.nan)

                                else:
                                    depth1 = int(float(depth_levels[i]) * (1 / dt))
                                    if variable == 'O2':
                                        temp_at_date.append(float(sim[depth1])*0.001)
                                    else:
                                        temp_at_date.append(float(sim[depth1]))

                        simulation_dict[sims_dates[y]] = temp_at_date

                y = 0

                csvFile = csv.writer(file)
                csvFile.writerow(
                    ["Date", "Depth", "Observations", "Simulations"])

                for date in sims_dates[number_skip_data:]:

                    if self.output_folder == "D:\output\PT\Alqueva\EWEMBI\historical":
                        date1 = date + 30000
                    else:
                        date1 = date

                    if date1 in observation_dict.keys():

                        if 1 == 1:  # try:
                            for i in range(0, len(observation_dict[date1])):

                                if str(observation_dict[date1][i])[0].isnumeric() and str(simulation_dict[date][i])[
                                    0].isnumeric():
                                    csvFile.writerow(
                                        [date1, depth_levels[i], observation_dict[date1][i], simulation_dict[date][i]])

                        # except KeyError:
                        #     continue
                return None

    def stats_lake(self, Observation, Modelisation):
        final_stat = ["RMSE", "NSE", "RSR", "Pbias", "R2", "SOS", "nrmse"]
        final_stat[3] = round(percent_bias_pbias(Observation, Modelisation), 3)
        correlation_matrix = np.corrcoef(Observation, Modelisation)
        correlation_xy = correlation_matrix[0, 1]
        final_stat[4] = round(correlation_xy ** 2, 3)
        rmse = root_mean_square(Observation, Modelisation)
        final_stat[0] = round(rmse[0], 3)
        final_stat[2] = round(rmse[1], 3)
        final_stat[6] = round(rmse[2], 3)
        final_stat[1] = round(nash_sutcliffe_efficiency_nse(Observation, Modelisation), 3)
        final_stat[5] = round(sums_of_squares(Observation, Modelisation), 3)

        return final_stat


class Graphics:
    """
    The main class that regroups all functions to plot modeled data
    """

    def __init__(self, Lake, width=6.5, height=6.5, font_family="Arial", size=11):

        self.lake = Lake
        self.name = Lake.name
        self.observation_folder = Lake.observation_folder
        self.save_date = Lake.save_date

        # Default value for parameters related to Temperature
        self.kz_N0 = Lake.kz_N0
        self.c_shelter = Lake.c_shelter
        self.i_scv = Lake.i_scv
        self.i_sct = Lake.i_sct
        self.swa_b0 = Lake.swa_b0
        self.swa_b1 = Lake.swa_b1

        # Default value for parameters related to Oxygen

        self.I_scDOC = Lake.I_scDOC
        self.I_scO = Lake.I_scO

        # Default value for parameter related to Chl_a
        # ???

        self.input_folder = Lake.input_folder
        self.output_folder = Lake.output_folder
        self.figures_folder = Lake.figures_folder

        self.time = time.strftime("%Y%m%d-%H%M%S")
        # self.time = 4
        self.SMALL_SIZE = size - 2
        self.MEDIUM_SIZE = size - 1
        self.BIGGER_SIZE = size
        self.font_family = font_family
        self.scatterDotSize = 80
        self.lineswidth = {"surface": 1.8, "deepwater": 1.8}

        plt.style.context('seaborn-paper')
        plt.rc('font', size=self.MEDIUM_SIZE)  # controls default text sizes
        plt.rc('axes', titlesize=self.MEDIUM_SIZE)  # fontsize of the axes title
        plt.rc('axes', labelsize=self.BIGGER_SIZE)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=self.MEDIUM_SIZE)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=self.MEDIUM_SIZE)  # fontsize of the tick labels
        plt.rc('legend', fontsize=self.SMALL_SIZE)  # legend fontsize
        plt.rc('figure', titlesize=self.MEDIUM_SIZE)  # fontsize of the figure title
        plt.rcParams["font.family"] = self.font_family
        self.width, self.height = width, height

    def figures_comparison_timeseries_and_profiles(self, save_figures=False, iteration=0, overall=True, byvariable=True,
                                                   variable_calibrated='T', mat_data=[]):
        """

        :param variable:
        :param save_figures:
        :return:
        """
        all_performances = [['s', 'd'], ['s', 'd'], ['s', 'd']]
        variables_list = ['T', 'O2', 'Chl']
        variables_limit = {'T': [-0.2, 24], 'O2': [-0.2, 16], 'Chl': [-0.2, 120]}
        depth_levels = [{'surface': 0.5, 'deepwater': 5}, {'surface': 0.5}]
        if len(mat_data) == 2:
            water, sediment = mat_data[0], mat_data[1]
        else:
            water = mat_data

        # figure by variable
        if byvariable:
            fig1 = plt.figure(frameon=False)
            ax1 = plt.subplot(2, 3, 1)
            ax2 = plt.subplot(2, 3, 2)
            ax3 = plt.subplot(2, 3, 3)
            ax4 = plt.subplot(2, 1, 2)
            axsV = [ax1, ax2, ax3, ax4]
            manager = plt.get_current_fig_manager()
            manager.window.showMaximized()

            # plt.tight_layout(h_pad=0.5, rect=[0.05, 0.05, 0.05, 0.05])
            # print("\n%s" % variable_calibrated)

            # define data
            comparison_file = os.path.join(self.output_folder, "%s_comparisonall.csv" % (variable_calibrated))

            comparison = pd.read_csv(comparison_file, names=['Date', 'Depth', 'Observations', 'Modelisation'],
                                     header=None, skiprows=1)

            comparison['variable'] = variables_dict[variable_calibrated]
            comparison['Dates'] = pd.to_datetime(comparison["Date"].astype(str), format='%Y%m%d')
            comparison = comparison.set_index(comparison['Dates'])

            dates = list(comparison['Date'].unique())

            if variable_calibrated == "T":
                sim_file = pd.DataFrame(water['T'][0, 0]).transpose()
            else:
                data = water['concentrations'][0, 0][variable_calibrated][0, 0]
                sim_file = pd.DataFrame(data).transpose()


            # sim_file = pd.read_csv(os.path.join(self.output_folder, "%szt.csv" % variable_calibrated), header=None)
            modeldata = pd.DataFrame(columns=["Dates", "Model_surface", "Model_deepwater"])
            modeldata['Dates'] = pd.date_range(start='1/1/2019', end='31/12/2020')

            if variable_calibrated == "O2":
                modeldata['Model_surface'] = sim_file[0]*0.001
                modeldata['Model_deepwater'] = sim_file[len(sim_file.columns) - 1]*0.001
            else:
                modeldata['Model_surface'] = sim_file[0]
                modeldata['Model_deepwater'] = sim_file[len(sim_file.columns) - 1]

            results, performances = self.comparison_obs_sims_plot_V2(
                variable_analized=variables_dict[variable_calibrated], modeldata=modeldata, obsdata=comparison,
                depthlayers=depth_levels[0], ax=axsV[3], individual=True)

            # plt.tight_layout(h_pad=0.5, rect=[0.05, 0.05, 0.05, 0.05])
            # plt.show()

            datesprofile = [dates[0], dates[floor(len(dates) / 4)], dates[floor(len(dates) / 2) - 2],
                            dates[floor(len(dates) / 2)], dates[floor(len(dates) / 4) * 3 + 2], dates[-1]]

            final_performance = {'rmse': [], 'r2': []}
            for variable in variables_list:
                variables_index = {'T': 0, 'O2': 1, 'Chl': 2}
                # define data
                comparison_file = os.path.join(self.output_folder, "%s_comparisonall.csv" % (variable))

                comparison = pd.read_csv(comparison_file, names=['Date', 'Depth', 'Observations', 'Modelisation'],
                                         header=None, skiprows=1)
                comparison['Dates'] = pd.to_datetime(comparison["Date"].astype(str), format='%Y%m%d')
                comparison = comparison.set_index(comparison['Dates'])

                r_value, rmse = self.graphiqueTOC(comparison['Observations'], comparison['Modelisation'],
                                                  comparison['Depth'], variable=variable,
                                                  ax=axsV[variables_index[variable]])
                final_performance['rmse'].append(round(rmse, 3))
                final_performance['r2'].append(round(r_value, 3))
            plt.show(fig1)
            if save_figures:
                if not os.path.exists(self.figures_folder):
                    os.mkdir(self.figures_folder)
                fig1.savefig(os.path.join(self.figures_folder, "%s_%s_ite_%s_TOC_comparison.png" % (
                variable_calibrated, self.save_date, iteration)))
            fig1.canvas.mpl_connect('close_event', lambda _: fig1.canvas.manager.window.destroy())
            plt.close('all')

            # overall figure
        if overall:
            fig2, axs = plt.subplots(3, 3, gridspec_kw={'width_ratios': [3, 1, 1], 'height_ratios': [1, 1, 1]},
                                     frameon=False)

            manager = plt.get_current_fig_manager()
            manager.window.showMaximized()

            for variable in variables_list:
                print("\n%s" % variable)
                axis_index = variables_list.index(variable)

                # [0,0] [0,1] [0,2]
                # [1,0] [1,1] [1,2]
                # [2,0] [2,1] [2,2]

                # define data
                comparison_file = os.path.join(self.output_folder, "%s_comparisonall.csv" % (variable))

                comparison = pd.read_csv(comparison_file, names=['Date', 'Depth', 'Observations', 'Modelisation'],
                                         header=None, skiprows=1)

                comparison['variable'] = variables_dict[variable]
                comparison['Dates'] = pd.to_datetime(comparison["Date"].astype(str), format='%Y%m%d')
                comparison = comparison.set_index(comparison['Dates'])

                dates = list(comparison['Date'].unique())

                if variable == "T":
                    sim_file = pd.DataFrame(water['T'][0, 0]).transpose()
                else:
                    data = water['concentrations'][0, 0][variable][0, 0]
                    sim_file = pd.DataFrame(data).transpose()

                # sim_file = pd.read_csv(os.path.join(self.output_folder, "%szt.csv" % variable), header=None)
                modeldata = pd.DataFrame(columns=["Dates", "Model_surface", "Model_deepwater"])
                modeldata['Dates'] = pd.date_range(start='1/1/2019', end='31/12/2020')
                if variable == "O2":
                    modeldata['Model_surface'] = sim_file[0]*0.001
                    modeldata['Model_deepwater'] = sim_file[len(sim_file.columns) - 1]*0.001
                else:
                    modeldata['Model_surface'] = sim_file[0]
                    modeldata['Model_deepwater'] = sim_file[len(sim_file.columns) - 1]

                if variable == "Chl":
                    results, performances = self.comparison_obs_sims_plot_V2(variable_analized=variables_dict[variable],
                                                                             modeldata=modeldata,
                                                                             obsdata=comparison,
                                                                             depthlayers=depth_levels[1],
                                                                             ax=axs[axis_index][0])

                else:
                    results, performances = self.comparison_obs_sims_plot_V2(variable_analized=variables_dict[variable],
                                                                             modeldata=modeldata,
                                                                             obsdata=comparison,
                                                                             depthlayers=depth_levels[0],
                                                                             ax=axs[axis_index][0])
                all_performances[axis_index] = performances

                datesprofile = [dates[0], dates[floor(len(dates) / 4)], dates[floor(len(dates) / 2) - 2],
                                dates[floor(len(dates) / 2)], dates[floor(len(dates) / 4) * 3 + 2], dates[-1]]

                if variable == "T":
                    profile_pos = [[1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
                else:
                    profile_pos = [[1, 2], [1, 1], [2, 2], [2, 1], [1, 0], [2, 0]]
                for date in datesprofile:
                    data = comparison[comparison['Date'] == date]
                    Depth = data['Depth'].tolist()
                    variable_obs = data['Observations'].tolist()
                    variable_sim = data['Modelisation'].tolist()
                    self.profiles(variable_sim, variable_obs, Depth, variable, variables_limit, date,
                                  axs[profile_pos[datesprofile.index(date)][1]][
                                      profile_pos[datesprofile.index(date)][0]],
                                  profile_pos[datesprofile.index(date)][1], profile_pos[datesprofile.index(date)][0])

                    # plt.tight_layout(h_pad=0.5, rect=[0.05, 0.05, 0.05, 0.05])
                if variable == "T":
                    axs[2, 1].xaxis.set_visible(True)
                    axs[2, 2].xaxis.set_visible(True)
            plt.subplots_adjust(wspace=0.1, hspace=0.1)
            plt.show(fig2)
        if save_figures:
            if not os.path.exists(self.figures_folder):
                os.mkdir(self.figures_folder)
            fig2.savefig(os.path.join(self.figures_folder, "%s_%s_ite_%s_timeseries_and_profiles.png" % (
            variable_calibrated, self.save_date, iteration)))
        fig2.canvas.mpl_connect('close_event', lambda _: fig2.canvas.manager.window.destroy())
        plt.close('all')
        while True:
            M_score = input("Calibration Score (between 0 to 10): ")
            if M_score.upper() not in ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'):
                print("answer giving is not an option, choose between 0 and 10.\n")
                continue
            else:
                M_score = int(M_score)
                break
        Note = input("Calibration Note (press Enter is nothing to note): ")

        return all_performances, final_performance, M_score, Note

    def graphiqueTOC(self, x, y, z, variable, ax=None):

        if variable == "T":
            variable = "Temperature (°C)"
        # Esthetic Parameters
        colorpalette = sns.color_palette("dark", 10)
        sns.set_style("ticks", {"xtick.major.size": 100, "ytick.major.size": 100})
        plt.xticks(rotation=15)
        plt.rcParams.update({"font.family": self.font_family})

        if ax == None:
            ax = plt.gcf()
        # Figure
        if variable == "DO":
            ax.set_xticks(np.arange(-0, 20, 2))
            ax.set_yticks(np.arange(-0, 20, 2))
            ax.set_ylim(-0, 20)
            ax.set_xlim(-0, 20)
            lineStart = -0
            lineEnd = 20
        elif variable == "Chl":
            ax.set_xticks(np.arange(-0, 112, 10))
            ax.set_yticks(np.arange(-0, 112, 10))
            ax.set_ylim(-0, 112)
            ax.set_xlim(-0, 112)
            lineStart = -0
            lineEnd = 112
        else:
            ax.set_xticks(np.arange(-0, 27, 2))
            ax.set_yticks(np.arange(-0, 27, 2))
            ax.set_ylim(-0, 27)
            ax.set_xlim(-0, 27)
            lineStart = -0
            lineEnd = 27

        ccmap = 'seismic_r'
        xall, yall = [float(i) for i in x], [float(i) for i in y]

        if variable == "Temperature (°C)":
            linearg = {'color': colorpalette[0], 'label': "y= x", 'linewidth': 1, 'linestyle': '--'}
            r_value, rmse = base_plot_comparison(xall, yall, lineStart=lineStart, lineEnd=lineEnd, ax=ax, linecolor="k")
            line_plot(lineStart=lineStart, lineEnd=lineEnd, ax=ax, linearg=linearg)
            # ccmap = 'Blues'
            markers = ['s'] * 12
        elif variable == "O2":
            linearg = {'color': colorpalette[0], 'label': "y= x", 'linewidth': 1, 'linestyle': '--'}
            r_value, rmse = base_plot_comparison(xall, yall, lineStart=lineStart, lineEnd=lineEnd, ax=ax, linecolor="k")
            line_plot(lineStart=lineStart, lineEnd=lineEnd, ax=ax, linearg=linearg)
            # ccmap = 'Blues'
            markers = ['o'] * 12
        else:
            linearg = {'color': colorpalette[0], 'label': "y= x", 'linewidth': 1, 'linestyle': '--'}
            r_value, rmse = base_plot_comparison(xall, yall, lineStart=lineStart, lineEnd=lineEnd, ax=ax, linecolor="k")
            line_plot(lineStart=lineStart, lineEnd=lineEnd, ax=ax, linearg=linearg)
            # ccmap = 'Blues'
            markers = ['^'] * 12

        # markers = ["o", "v", "^", "s", "P", "*", ">", "X", "D", "<", "p", "d"]

        # print(markers[0])
        cs = ax.scatter(xall, yall, c=z, marker=markers[0], s=self.scatterDotSize, cmap=ccmap, linewidths=1,
                        edgecolors='k',
                        alpha=0.8)

        ax.set_xlabel("Observed %s" % variable, fontsize=self.BIGGER_SIZE)
        ax.set_ylabel("Modeled %s" % variable, fontsize=self.BIGGER_SIZE)
        if variable == "DO":
            ax.set_xticks(np.arange(-0, 20, 2))
            ax.set_yticks(np.arange(-0, 20, 2))
            ax.set_ylim(-0, 20)
            ax.set_xlim(-0, 20)
        elif variable == "Chl":
            ax.set_xticks(np.arange(-0, 112, 10))
            ax.set_yticks(np.arange(-0, 112, 10))
            ax.set_ylim(-0, 112)
            ax.set_xlim(-0, 112)
        else:
            ax.set_xticks(np.arange(-0, 27, 2))
            ax.set_yticks(np.arange(-0, 27, 2))
            ax.set_ylim(-0, 27)
            ax.set_xlim(-0, 27)

        return r_value, rmse

    def profiles(self, variable_sim, variable_obs, Depth, variable, variables_limit, date, ax=None, horizontal=0,
                 vertical=0):
        markerStyleByVariable = {"T": "s", "O2": "o", "Chl": "^"}
        markerColorByVariable = {"T": "b", "O2": "r", "Chl": "g"}
        if ax is None:
            ax = plt.gca()

        if variable != "T":
            axisV = ax.twiny()
        else:
            axisV = ax

        offset = 60
        if variable == "Chl":
            # Move twinned axis ticks and label from top to bottom
            axisV.xaxis.set_ticks_position("top")
            axisV.xaxis.set_label_position("top")

            axisV.set_xlabel("Chl")

        elif variable == "T":
            axisV.set_xlabel("T")
            axisV.set_ylabel("Depth")
            axisV.annotate("%s-%s-%s" % (str(date)[0:4], num2month[int(str(date)[4:6])], str(date)[6:]), xy=(1, 0),
                           xycoords='axes fraction', horizontalalignment='right',
                           verticalalignment='bottom')

        else:
            axisV.set_xlabel("DO")

        axisV.plot(variable_sim, Depth, 'k%s--' % (markerStyleByVariable[variable]), fillstyle='none',
                   markersize=int(self.scatterDotSize / 10))
        axisV.plot(variable_obs, Depth, '%s%s-' % (markerColorByVariable[variable], markerStyleByVariable[variable]),
                   markersize=int(self.scatterDotSize / 10))
        axisV.set_ylim(6, 0)
        axisV.set_xlim(variables_limit[variable][0], variables_limit[variable][1])

        if horizontal == 0 and variable != 'T':
            axisV.tick_params(axis='x', colors=markerColorByVariable[variable])
            if variable == "O2":
                axisV.set_xlabel(variable, horizontalalignment='right', x=1.0)

            else:
                # Offset the twin axis below the host
                axisV.spines["top"].set_position(("axes", 1.23))
                axisV.set_xlabel(variable, horizontalalignment='left', x=0.0)
                axisV.get_xaxis().set_label_coords(-0.03, 1.12)

            axisV.xaxis.set_visible(True)
            axisV.xaxis.label.set_color(markerColorByVariable[variable])

        else:
            axisV.xaxis.set_visible(False)

        if horizontal == 2 and variable == 'T':
            axisV.xaxis.set_visible(True)
            axisV.tick_params(axis='x', colors=markerColorByVariable[variable])
            axisV.set_xlabel(variable, horizontalalignment='right', x=0.0)
            axisV.xaxis.label.set_color(markerColorByVariable[variable])

    def comparison_obs_sims_plot_V2(self, variable_analized, modeldata, obsdata, depthlayers, ax=None,
                                    individual=False):
        """

        :type obsdata: pd.DataFrame
        :param variable_analized:
        :param calibration_methods:
        :param modeldata:
        :param obsdata:
        :param depthlayers:
        :param ice_cover:
        :param icecovertarea:
        :return:
        """

        markerStyleByVariable = {"Temperature": "s", "DO concentration": "o", "Chlorophyll a": "^"}
        markerColorByVariable = {"Temperature": "blue", 'DO concentration': "red", "Chlorophyll a": "green"}
        lineColorByWaterLevel = {"surface": 'black', "deepwater": markerColorByVariable[variable_analized]}
        lineStyleByWaterLevel = {"surface": "-.", "deepwater": "-"}
        markerfillingByWaterLevel = {"surface": 'none', "deepwater": markerColorByVariable[variable_analized]}
        # sns.set_style("ticks", {"xtick.major.size": 100, "ytick.major.size": 100})

        linewidthByWaterLevel = self.lineswidth
        orderPlotPresentedByWaterLevel = {"surface": 100, "deepwater": 10}
        scatterPlotDotSize = self.scatterDotSize
        transparenceLineplot = 0.8
        # iceCovertAreaColor = "grey"
        # iceCovertAreaAlpha = 0.2

        legendNames = ["simulated_surface", "simulated_deepwater", "Observations_surface", "Observations_deepwater"]

        # fig, axs = plt.subplots(2, 1, figsize=(15, 8))#, gridspec_kw={'height_ratios': [1, 1]})
        if ax is None:
            ax = plt.gca()

        plt.rcParams["font.family"] = self.font_family
        if variable_analized == "DO Concentration":
            limit = [-0.5, 20]
        elif variable_analized == "Temperature":
            limit = [-0.5, 25]
        else:
            limit = [-0.5, 25]
        results = pd.DataFrame(columns=[])
        stop = False

        all_performances = []
        for depth_level in ["surface", "deepwater"]:
            if depth_level == "deepwater1" and variable_analized == "Chlorophyll a":
                print("no deep data")
            else:
                try:
                    depthlayer = depthlayers[depth_level]

                except:
                    stop = True

                if not stop:
                    if not obsdata.loc[obsdata['Depth'] == depthlayer]["Observations"].empty:
                        result = pd.DataFrame()
                        result["Dates"] = obsdata.loc[obsdata['Depth'] == depthlayer]["Dates"]
                        result['Observations'] = obsdata.loc[obsdata['Depth'] == depthlayer]["Observations"]
                        result.set_index('Dates', inplace=True)

                        # Plot each method and observations
                        lineplotstyle = {"linewidth": linewidthByWaterLevel[depth_level],
                                         "color": lineColorByWaterLevel[depth_level],
                                         "zorder": orderPlotPresentedByWaterLevel[depth_level],
                                         "linestyle": lineStyleByWaterLevel[depth_level],
                                         "alpha": transparenceLineplot}

                        model_data = pd.DataFrame()
                        model_data["Model_%s" % (depth_level)] = modeldata["Model_%s" % (depth_level)]
                        model_data["Model_%s" % (depth_level)] = modeldata["Model_%s" % (depth_level)]
                        model_data['Dates'] = modeldata['Dates']
                        model_data.set_index('Dates')

                        concact = pd.concat([model_data, result], axis=1)
                        result = concact

                        scatterplotstyle = {'marker': markerStyleByVariable[variable_analized],
                                            's': scatterPlotDotSize, "edgecolor": 'k',
                                            "facecolors": markerfillingByWaterLevel[depth_level],

                                            "linewidth": linewidthByWaterLevel[depth_level],
                                            "zorder": orderPlotPresentedByWaterLevel[depth_level]}

                        if not individual:
                            performances = self.lake.stats_lake(
                                obsdata.loc[obsdata['Depth'] == depthlayer]["Observations"],
                                obsdata.loc[obsdata['Depth'] == depthlayer]["Modelisation"])
                            print(
                                "Performances Result for %s :\nRMSE : %s NSE : %s RSR : %s Pbias : %s \nR2 : %s SOS : %s nrmse : %s" % (
                                depth_level, performances[0], performances[1], performances[2], performances[3],
                                performances[4], performances[5], performances[6]))
                        timeline_plot(modeleddata=model_data["Model_%s" % (depth_level)],
                                      modeleddates=model_data['Dates'],
                                      observeddata=obsdata.loc[obsdata['Depth'] == depthlayer]["Observations"],
                                      observeddates=obsdata.loc[obsdata['Depth'] == depthlayer]["Dates"],
                                      ax=ax, ylimit=limit, sct_kwargs=scatterplotstyle, line_kwargs=lineplotstyle)

                        # linepos = ["surface", "deepwater"].index(depth_level)
                        # ax.lines[linepos].set_linestyle(lineStyleByWaterLevel[depth_level])

                        results = pd.concat([results, result], axis=1, sort=True)
                        results = results.loc[:, ~results.columns.duplicated()]

                        if variable_analized == "Temperature":
                            unity = "°C"
                            ax.set_ylabel("%s (°C)" % (variable_analized))
                        elif variable_analized == "DO concentration":
                            unity = "mg*$L^-1$"
                            variable = "DO concentration"
                            # axs.set( ylabel="%s \n at %s m (mg*$\mathregular{L^{\-1}}$)" % (variable_analized, depthlayer))
                            ax.set_ylabel("%s (mg*$\mathregular{L^{-1}}$)" % (variable),
                                          linespacing=0.8)
                        else:
                            unity = "ug*$L^-1$"
                            variable = "Chloropyll a"
                            # axs.set( ylabel="%s \n at %s m (mg*$\mathregular{L^{\-1}}$)" % (variable_analized, depthlayer))
                            ax.set_ylabel("%s (ug*$\mathregular{L^{-1}}$)" % (variable), linespacing=0.8)

                        # # legend
                        # ax.get_legend().remove()
                        # legends = ax.legend(legendNames,ncol= len(legendNames))
                        # legends.get_frame().set_facecolor('#FFFFFF')
                    else:

                        result = pd.DataFrame()
                        result["Dates"] = obsdata.loc[obsdata['Depth'] == depthlayer]["Dates"]
                        result['Observations'] = obsdata.loc[obsdata['Depth'] == depthlayer]["Observations"]
                        result.set_index('Dates', inplace=True)

                        # Plot each method and observations
                        lineplotstyle = {"linewidth": linewidthByWaterLevel[depth_level],
                                         "color": lineColorByWaterLevel[depth_level],
                                         "zorder": orderPlotPresentedByWaterLevel[depth_level],
                                         "alpha": transparenceLineplot}

                        model_data = pd.DataFrame()
                        if variable_analized == "O2":
                            model_data["Model_%s" % (depth_level)] = modeldata["Model_%s" % (depth_level)]

                        else:
                            model_data["Model_%s" % (depth_level)] = modeldata["Model_%s" % (depth_level)]
                        model_data['Dates'] = modeldata['Dates']
                        model_data.set_index('Dates')

                        concact = pd.concat([model_data, result], axis=1)
                        result = concact

                        scatterplotstyle = {'marker': markerStyleByVariable[variable_analized],
                                            's': scatterPlotDotSize, "edgecolor": 'k',
                                            "facecolors": markerfillingByWaterLevel[depth_level],
                                            "linewidth": linewidthByWaterLevel[depth_level],
                                            "linestyle": lineStyleByWaterLevel[depth_level],
                                            "zorder": orderPlotPresentedByWaterLevel[depth_level]}

                        if not individual:
                            performances = self.lake.stats_lake(
                                obsdata.loc[obsdata['Depth'] == depthlayer]["Observations"],
                                obsdata.loc[obsdata['Depth'] == depthlayer]["Modelisation"])
                            print(
                                "Performances Result for %s :\nRMSE : %s NSE : %s RSR : %s Pbias : %s \nR2 : %s SOS : %s nrmse : %s" % (
                                    depth_level, performances[0], performances[1], performances[2], performances[3],
                                    performances[4],
                                    performances[5], performances[6]))

                        timeline_plot(modeleddata=model_data["Model_%s" % (depth_level)],
                                      modeleddates=model_data['Dates'],
                                      observeddata=obsdata.loc[obsdata['Depth'] == depthlayer]["Observations"],
                                      observeddates=obsdata.loc[obsdata['Depth'] == depthlayer]["Dates"],
                                      ax=ax, ylimit=limit, sct_kwargs=scatterplotstyle, line_kwargs=lineplotstyle)
                        linepos = ["surface", "deepwater"].index(depth_level)
                        ax.lines[linepos].set_linestyle(lineStyleByWaterLevel[depth_level])

                        results = pd.concat([results, result], axis=1)
                        results = results.loc[:, ~results.columns.duplicated()]

                        if variable_analized == "Temperature":
                            unity = "°C"
                            ax.set_ylabel("%s \n at %s m (°C)" % (variable_analized, depthlayer))
                        elif variable_analized == "DO Concentration":
                            unity = "mg*$L^-1$"
                            variable = "DO concentration"
                            # axs.set( ylabel="%s \n at %s m (mg*$\mathregular{L^{\-1}}$)" % (variable_analized, depthlayer))
                            ax.set_ylabel("%s \n at %s m (mg*$\mathregular{L^{-1}}$)" % (variable, depthlayer),
                                          linespacing=0.8)
                        else:
                            unity = "mg*$L^-1$"
                            variable = "Chloropyll a"
                            # axs.set( ylabel="%s \n at %s m (mg*$\mathregular{L^{\-1}}$)" % (variable_analized, depthlayer))
                            ax.set_ylabel("%s \n at %s m (mg*$\mathregular{L^{-1}}$)" % (variable, depthlayer),
                                          linespacing=0.8)

                        # legend
                        ax.get_legend().remove()
                        # legends = ax.legend(legendNames,ncol= len(legendNames))
                        # legends.get_frame().set_facecolor('#FFFFFF')
                # ax.set_xticklabels([2019,2020,2021], rotation=(0), va='bottom', ha='center')
                # ax.set_yticklabels([" 0"," 0", "10", "20", "30"],va= 'center', ha='left')

                if not individual:
                    all_performances.append(performances)
                    if variable_analized == "Chlorophyll a":
                        ax.xaxis.set_visible(True)
                    else:
                        ax.xaxis.set_visible(False)

        return (results, all_performances)


if __name__ == "__main__":
    # ask for lake
    # print("Figure generation")
    # Graphics(Lake('Nedre')).figures_comparison_timeseries_and_profiles()
    # matlab_directory = matlab_folder
    matlab_directory = ask_for_matlab_directory()
    lake_name = aks_for_which_lake_wanted_to_calibrated()
    save_option = ask_what_to_save()

    print("\n------------------- Manual Calibration will Start for lake %s -------------------\n" % lake_name)
    lake = Lake(lake_name).manual_calibration_loop(save_option[0], save_option[1], save_option[2], save_option[3],
                                                   matlab=matlab_directory)
