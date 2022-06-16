#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : Marianne Cote
# Created Date: 2022-06-10
# version ='1.0'
# ---------------------------------------------------------------------------
""" Script formatting the Raw data into correct format for MyLake"""
# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------

import os
import pandas as pd

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
