#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : Marianne Cote
# Created Date: 2022-06-10
# version ='1.0'
# ---------------------------------------------------------------------------
""" Script formatting the Raw data into the correct format for MyLake"""
# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------

import os
import pandas as pd
import csv

# ---------------------------------------------------------------------------
# Global Variables
# ---------------------------------------------------------------------------

# Directory Path
where_I_am = "OI/"
raw_weather_directory = "../Raw_data/meteo"
raw_obs_directory = "../Raw_data/Donnees_mouillage_Bromont"
raw_geo_data_directory = "../Raw_data/geo_data"
observation_directory = r"obs"

# Information Dict
lakes_dict = {"Bromont": "Bromont"}
variables_dict = {"Datum": "Date",
                  "Nederbördsmängd": "Precipitation",
                  "Global Irradians (svenska stationer)": "Global radiation",
                  "Relativ Luftfuktighet": "Relative humidity",
                  "Vindhastighet": "Wind speed",
                  "Lufttemperatur": "Air temperature",
                  "Lufttryck reducerat havsytans nivå": "Air pressure"}

all_variable_needed = ["Date", "Precipitation", "GlobalRadiation", "Relative humidity", "Wind speed", "Air temperature",
                       "Air pressure", "Cloud cover"]


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

def get_all_files_in_directory(directory: str):
    """
    a function that walks in the directory and gets all file names. Return a list of all the files in the directory.
    :param directory:   directory path to the folder we want to get all the file names from.
    :type directory:    str
    :return:            list of all files (with the relative path to the directory) in the folder
    """
    final_list_of_all_files = []
    if os.path.exists(directory):
        print("Start walk into %s folder" % directory)
        for (root, dirs, files) in os.walk(directory):
            if files:
                final_list_of_all_files.extend([os.path.join(root, file) for file in files])
        print("End walk into %s folder" % directory)
    else:
        print("Folder given does not exist.")
    return final_list_of_all_files


def get_information_from_file_as_dataframe(file_path: str, skipuntil=None):
    """
    A function that opens the file and extracts the information
    :param skipuntil:   if part of the script is not needed, the script will remove lines of the data until this string
                        is found. the function will start to read data after that variable is found. If the variable is
                        None, the function considers that the data start at the first line.
    :type skipuntil:    str

    :param file_path:   file path + complete file name.
    :type file_path:    str

    :return:            if raw data are in matrix format (for txt file), return a DataFrame of the columns (column name
                        given if first data is a string, else used default number) and rows from the file.
                        Return None if the file does not contain data in matrix format (more than one column),
                        does not exist, or is empty.
    """

    if os.path.exists(file_path):
        # if file is a csv or xlsx
        if file_path.lower().endswith(('.xlsx', '.csv')):
            if os.path.splitext(file_path)[1] == ".csv":
                with open(file_path) as csvfile:
                    sniffer = csv.Sniffer()
                    dialect = sniffer.sniff(csvfile.read())

                dataframe_of_the_file = pd.read_csv(file_path, sep=dialect.delimiter)
                if skipuntil is not None:
                    dataframe_of_the_file = \
                        dataframe_of_the_file.iloc[
                            (dataframe_of_the_file.loc[dataframe_of_the_file[0] ==
                                                       skipuntil].index[0] + 1):, :].reset_index(drop=True)

            else:
                dataframe_of_the_file = pd.read_excel(file_path)
            if dataframe_of_the_file.empty:
                print("file given is empty")
                return None
            else:
                return dataframe_of_the_file
        elif file_path.lower().endswith('.txt'):
            print()
        else:
            print("File format %s is not supported by this script" % os.path.splitext(file_path)[1])
            return None

    else:
        print("File %s does not exist or the relative path is not in function of this script" % file_path)
        return None


def extract_bathymetry(lakes: dict, rawdata_path: str = raw_geo_data_directory,
                       rawnamefile: str = "area_export_corrected.csv", observation_path: str = observation_directory):
    """
    Generate bathymetry files for each lake. This function takes the information from the raw data to create formatted
    files in the observation folder.

    :param lakes:               Lake names used in the bathymetry files. The keys are the shorter lake name that
                                will name the folders and the values are the full names of the lakes used by
                                the raw files.
    :type lakes:                dict

    :param rawdata_path:        Relative path where the raw data are.
    :type rawdata_path          str

    :param rawnamefile:         Name of the file with the area at each depth.
    :type rawnamefile:          str

    :param observation_path:    Relative path where the formatted files of the observations are.
    :type observation_path:     str

    :return:                    1 if the function was successful, 0 if an error happens during the extraction.

    """
    problem_during_loop = False
    for lake in lakes.keys():
        try:
            # Extract data from Raw file
            raw_data = pd.read_csv(os.path.join(rawdata_path, rawnamefile))
            data_for_lake = raw_data.loc[raw_data["lake_name"] == lakes[lake]]

            # Generate bathymetry file
            if not os.path.exists(os.path.join(observation_path, lake)):
                os.makedirs(os.path.join(observation_path, lake))
            bathymetry_file = os.path.join(observation_path, lake, "%s_bathymetry.csv" % lake)
            bathymetry_data = pd.DataFrame(columns=["Lake", "Depth", "Area"])
            bathymetry_data['Depth'] = list(data_for_lake["depth m"])
            bathymetry_data['Area'] = list(data_for_lake["area m2"])
            bathymetry_data['Lake'] = lake
            bathymetry_data = bathymetry_data.sort_values("Depth")
            bathymetry_data.to_csv(bathymetry_file, index=False)

        except:
            print("Error with lake %s" % lake)
            problem_during_loop = True

    if problem_during_loop:
        print("The function could not extract bathymetry for one or more lakes")
        return 0
    else:
        print("bathymetry files generated")
        return 1


if __name__ == "__main__":
    print("Hello")
    print(extract_bathymetry(lakes=lakes_dict, rawnamefile="Bromont_Bathymetry.csv"))
