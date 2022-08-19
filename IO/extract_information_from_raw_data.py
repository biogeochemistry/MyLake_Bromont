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
import numpy as np
from itertools import *
import xlrd
import numbers

# ---------------------------------------------------------------------------
# Global Variables
# ---------------------------------------------------------------------------

# Directory Path
raw_weather_directory = "../Raw_data/meteo"
raw_obs_directory = "../Raw_data/Donnees_mouillage_Bromont"
raw_geo_data_directory = "../Raw_data/geo_data"
observation_directory = r"../obs"
input_directory = r"../IO"

# Information Dict
lakes_dict = {"Bromont": "Bromont"}
variables_dict = {"Date": "Date",
                  "Date de mesure": "Date",
                  "Total precip": "Precipitation",
                  "Préc, (mm)": "Precipitation",
                  " Pluie (mm)/jour": "Precipitation",
                  "Global Irradians (svenska stationer)": "Global radiation",
                  "HR, % ": "Relative humidity",
                  "Vitesse du vent, km/hr ": "Wind speed",
                  "Vitesse du vent, m/s": "Wind speed",
                  "Température de l'air ": "Air temperature",
                  "Temp, °C": "Air temperature",
                  "Temp., °C ":"Air temperature",
                  "Pression, Kpa": "Air pressure"}

all_variable_needed = ["Date","Precipitation","Global radiation","Relative humidity","Wind speed","Air temperature",
                       "Air pressure","Cloud cover"]



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

    :return:            if raw data are in matrix format (for txt file), return a dict of DataFrame of the columns
                        (column name given if first data is a string, else used default number) and rows from the file.
                        if excel file has multiple sheets, return a dictionnary of Dataframe with the sheet name as keys.
                        Return None if the file does not contain data in matrix format (more than one column),
                        does not exist, or is empty.
    """

    def isBadLine(line, skipuntil=skipuntil):
        return (skipuntil not in line)

    if os.path.exists(file_path):
        # if file is a csv or xlsx
        if file_path.lower().endswith(('.xlsx', '.csv')):
            if os.path.splitext(file_path)[1] == ".csv":
                with open(file_path) as csvfile:
                    sniffer = csv.Sniffer()
                    dialect = sniffer.sniff(csvfile.read())

                dataframe_of_the_file = pd.read_csv(file_path, sep=dialect.delimiter)
                if skipuntil is not None:
                    start_line = dataframe_of_the_file.loc[dataframe_of_the_file[0] == skipuntil].index[0] + 1
                    dataframe_of_the_file = dataframe_of_the_file.iloc[start_line:, :].reset_index(drop=True)

            else:
                xls = xlrd.open_workbook(file_path, on_demand=True)
                keys = xls.sheet_names()
                dataframe_of_the_file = {}
                for key in keys:
                    data_of_the_file = pd.read_excel(file_path, sheet_name=key)
                    dataframe_of_the_file[key] = data_of_the_file
                # data_of_the_file = pd.read_excel(file_path, sheet_name=0)
                # dataframe_of_the_file["depth_%s"%( float(''.join(c for c in file_path if (c.isdigit() or c == '.'))))] = data_of_the_file


            if type(dataframe_of_the_file) == dict:
                return dataframe_of_the_file
            else:
                if dataframe_of_the_file.empty:
                    print("file given is empty")
                    return None
                else:
                    return {"data": dataframe_of_the_file}



        elif file_path.lower().endswith('.txt'):

            with open(file_path) as f:
                number_of_line_skip = 0
                alllines = []
                for line in dropwhile(isBadLine, f):
                    line = line[:-1]
                    line = line.replace("            ", "")
                    line = line.replace("   ", "")
                    alllines.append(line.split(","))
                dataframe_of_the_file = pd.DataFrame(alllines[2:], columns=alllines[0])

                if "top" in file_path:
                    dataframe_of_the_file["depth"] = -6
                elif "bottom" in file_path:
                    dataframe_of_the_file["depth"] = -1
                return {"data": dataframe_of_the_file}
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
            data_for_lake = raw_data[raw_data["lake_name"] == lakes[lake]]

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

def extract_observations(lakes: dict, rawdata_path: str = raw_obs_directory,
                         observation_path: str = observation_directory):
    """
    Extract and merge observations from the lakes. This function takes the information from the raw data to create
    formatted files in the observation folder.

    :param lakes:               Lake names used in the bathymetry files. The keys are the shorter lake name that
                                will name the folders and the values are the full names of the lakes used by
                                the raw files.
    :type lakes:                dict

    :param rawdata_path:        Relative path where the raw data are.
    :type rawdata_path:         str

    :param observation_path:    Relative path where the formatted files of the observations are.
    :type observation_path:     str

    :return:                    1 if the function was successful, 0 if an error happens during the extraction.

    """

    files_with_information = get_all_files_in_directory(raw_obs_directory)

    for lake in lakes.keys():
        # Create file with all years
        observation_data = pd.DataFrame(
            columns=["Lake", "Full_Name", "Date", "Depth", 'Date_Depth', "Temp", "DO"])
        for file in files_with_information:
            print(file)
            # Create file with all years
            table_data = pd.DataFrame(
                columns=["Lake", 'Full_Name', "Date","Time", "Depth", 'Date_Depth', "Temp", "DO"])
            # Extract data from both Raw file
            if "Cat" in file:
                base_data = get_information_from_file_as_dataframe(file, "            Unix Timestamp,")['data']
                base_data["Unix Timestamp"] = base_data["Unix Timestamp"].replace(" ", "")
                raw_data = pd.DataFrame(columns=['Short_Name', 'Full_Name', 'Date', 'Depth', "Temp", "DO", "Secchi"])
                raw_data["Date"] = pd.to_datetime(base_data['Unix Timestamp'], unit='s')
                raw_data["Short_Name"] = lake
                raw_data['Depth'] = base_data["depth"].astype(float)
                raw_data['Full_Name'] = lakes[lake]
                raw_data["Temp"] = (base_data["Temperature"].replace(" ", "")).astype(float)
                raw_data["DO"] = (base_data["  Dissolved Oxygen"].replace(" ", "")).astype(float)

                # raw_data = pd.read_csv(os.path.join(rawdata_path, files_with_information[file]), encoding='iso-8859-1')
                data_for_lake = raw_data[raw_data["Short_Name"] == lake]
                table_data[['Date','Time']] = data_for_lake['Date'].astype(str).str.split(' ', 1, expand=True)
                table_data['Date'] = pd.to_datetime(table_data['Date'])
                table_data['Date'] = pd.to_datetime(table_data['Date'].astype(str), format='%Y-%m-%d')
                table_data['Depth'] = data_for_lake['Depth'].astype(float)
                table_data['Date_Depth'] = table_data['Date'].astype(str) + "_" + table_data['Depth'].astype(str)
                # table_data['Date'] = table_data['Date'].values.astype(np.int64) // 10 ** 9
                # table_data['Date'] = table_data['Date'].map(pd.Timestamp.timestamp)
                # table_data['Date'] = (table_data['Date'] - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s')

                table_data['Lake'] = lake
                table_data['Full_Name'] = lakes[lake]

                table_data["Temp"] = data_for_lake["Temp"]
                if "DO" in data_for_lake:
                    table_data["DO"] = data_for_lake["DO"]
                else:
                    table_data["DO"] = np.nan

                if file == files_with_information[0]:
                    observation_data = table_data
                else:
                    table_data = table_data.reset_index(drop=True)
                    observation_data = observation_data.append(table_data)
                    observation_data = observation_data.reset_index(drop=True)

            elif "_Détails" not in file:
                base_data_dict = get_information_from_file_as_dataframe(file, "1")
                for key in base_data_dict.keys():
                    base_data = base_data_dict[key]
                    base_data["depth"] = float(''.join(c for c in key if (c.isdigit() or c == '.')))*-1
                    raw_data = pd.DataFrame(columns=['Short_Name', 'Full_Name', 'Date', 'Depth', "Temp", "DO", "Secchi"])
                    raw_data["Date"] = base_data['Date Heure, GMT-04:00']
                    raw_data["Date"] = pd.to_datetime(raw_data["Date"], format='%y/%m/%d %H:%M:%S')
                    raw_data["Short_Name"] = lake
                    raw_data['Depth'] = base_data["depth"].astype(float)
                    raw_data['Full_Name'] = lakes[lake]
                    raw_data["Temp"] = base_data.filter(regex='Temp').astype(float)


                    # raw_data = pd.read_csv(os.path.join(rawdata_path, files_with_information[file]), encoding='iso-8859-1')
                    data_for_lake = raw_data[raw_data["Short_Name"] == lake]
                    table_data[['Date','Time']] = data_for_lake['Date'].astype(str).str.split(' ', 1, expand=True)
                    table_data['Date'] = pd.to_datetime(table_data['Date'])
                    table_data['Depth'] = data_for_lake['Depth'].astype(float)
                    table_data['Date_Depth'] = table_data['Date'].astype(str) + "_" + table_data['Depth'].astype(str)
                    # table_data['Date'] = table_data['Date'].values.astype(np.int64) // 10 ** 9
                    # table_data['Date'] = table_data['Date'].map(pd.Timestamp.timestamp)
                    # table_data['Date'] = (table_data['Date'] - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s')
                    table_data['Date'] = pd.to_datetime(table_data['Date'].astype(str), format='%Y-%m-%d')


                    table_data['Lake'] = lake
                    table_data['Full_Name'] = lakes[lake]

                    table_data["Temp"] = data_for_lake["Temp"]
                    if "DO" in data_for_lake:
                        table_data["DO"] = data_for_lake["DO"]
                    else:
                        table_data["DO"] = np.nan

                    if file == files_with_information[0]:
                        observation_data = table_data
                    else:
                        table_data = table_data.reset_index(drop=True)
                        observation_data = observation_data.append(table_data)
                        observation_data = observation_data.reset_index(drop=True)

        # Add data to DataFrame

        observation_data = observation_data.sort_values("Date_Depth")

        # Generate bathymetry file
        observation_file = os.path.join(observation_path, lake, "%s_observation_data_all.csv" % lake)
        observation_data.to_csv(observation_file, index=False)

        # Add data to DataFrame

        observation_data_names = list(observation_data.columns).remove("Time")
        observation_data['DO'] = observation_data['DO'].astype(float)
        observation_data1 = observation_data[["Full_Name","Lake","Date","Depth","Date_Depth",'DO']].groupby(["Full_Name","Lake","Date","Depth","Date_Depth"], dropna=False).mean()
        observation_data = observation_data.groupby(["Full_Name","Lake","Date","Depth","Date_Depth"],dropna=False).mean()
        observation_data.reset_index(inplace=True)


        # Generate bathymetry file
        observation_file = os.path.join(observation_path, lake, "%s_observation_data.csv" % lake)
        observation_data.to_csv(observation_file, index=False)


    return "Files generated"

def extract_climate(climat_directory: str = raw_weather_directory, observation_path: str = observation_directory):
    all_files = get_all_files_in_directory(climat_directory)
    all_data = pd.DataFrame()
    n = 1
    for file in all_files:
        if "original_with_issue" not in file:
            print("##### File %s of %s : %s #####"%(n, len(all_files),file))
            n+=1
            data_dict = get_information_from_file_as_dataframe(file)
            wind_ajustment = True
            if type(data_dict) is dict:
                for key in data_dict:
                    data = data_dict[key]
                    if 'Date' in data.columns:
                        data = data[data['Date'].notnull()]
                        data.loc[:,'Date'] = data.loc[:,'Date'].dt.strftime('%Y-%m-%d')
                        data.loc[:,'Date'] = pd.to_datetime(data.loc[:,'Date']).dt.floor('d')
                        data = data[data['Date'].notnull()]

                    elif 'Dates' in data.columns:
                        data = data[data['Dates'].notnull()]
                        # data['Dates'] = data['Dates'].dt.strftime('%Y-%m-%d')
                        data.loc[:,'Date'] = pd.to_datetime(data.loc[:,'Dates']).dt.floor('d')
                        data = data[data['Date'].notnull()]

                    elif 'Date de mesure'in data.columns:
                        data = data[data['Date de mesure'].notnull()]
                        # data['Date de mesure'] = data['Date de mesure'].dt.strftime('%Y-%m-%d')
                        data.loc[:,'Date'] = pd.to_datetime(data.loc[:,'Date de mesure']).dt.floor('d')
                        data = data[data['Date'].notnull()]

                    if all_data.empty:
                        for col in data.columns:
                            if col in variables_dict:
                                all_data[variables_dict[col]] = data[col]
                                if col == "Vitesse du vent, m/s":
                                    wind_ajustment = False
                        if len(all_data.columns) > 1:
                            all_data["Date"] = pd.to_datetime(all_data['Date'].astype(str), format='%Y-%m-%d')
                            organized_data = pd.DataFrame()
                            for column in all_data.columns:
                                if column == "Precipitation":
                                    subselect_data = all_data.groupby(["Date"]).sum()
                                    subselect_data.reset_index(inplace=True)
                                    subselect_data = subselect_data[["Date",column]]
                                    organized_data = pd.concat([organized_data, subselect_data], sort=False, ignore_index=True)
                                elif column == "Global radiation":
                                    subselect_data = all_data.groupby(["Date"]).mean() * (86000 / 1000000)
                                    subselect_data.reset_index(inplace=True)
                                    subselect_data = subselect_data[["Date", column]]
                                    organized_data = pd.concat([organized_data, subselect_data], sort=False, ignore_index=True)
                                elif column == "Wind speed" and wind_ajustment:
                                    subselect_data = all_data.groupby(["Date"]).mean() * (1000 / 3600)
                                    subselect_data.reset_index(inplace=True)
                                    subselect_data = subselect_data[["Date", column]]
                                    organized_data = pd.concat([organized_data, subselect_data], sort=False, ignore_index=True)
                                elif column != "Date":
                                    subselect_data = all_data.groupby(["Date"]).mean()
                                    subselect_data.reset_index(inplace=True)
                                    subselect_data = subselect_data[["Date", column]]
                                    organized_data = pd.concat([organized_data, subselect_data], sort=False, ignore_index=True)


                            all_data = organized_data

                            # if "Precipitation" in all_data.columns:
                            #     all_data = all_data.groupby(["Date"]).sum()
                            #     all_data.reset_index(inplace=True)
                            # elif "Global radiation" in all_data.columns:
                            #     all_data = all_data.groupby(["Date"]).mean() * (86000 / 1000000)
                            # elif "Wind speed" in all_data.columns:
                            #     all_data = all_data.groupby(["Date"]).mean() * (1000 / 3600)
                            # else:
                            #     all_data = all_data.groupby(["Date"]).mean()
                        else:
                            all_data = pd.DataFrame()
                    else:
                        select_data = pd.DataFrame()
                        for col in data.columns:
                            if col in variables_dict:
                                select_data[variables_dict[col]] = data[col]
                                if col == "Vitesse du vent, m/s":
                                    wind_ajustment = False
                        if len(select_data.columns) > 1:
                            select_data["Date"] = pd.to_datetime(select_data['Date'].astype(str), format='%Y-%m-%d')
                            for column in select_data.columns:
                                if column == "Precipitation":
                                    subselect_data = select_data.groupby(["Date"]).sum()
                                    subselect_data.reset_index(inplace=True)
                                    subselect_data = subselect_data[["Date",column]]
                                    all_data = pd.concat([all_data, subselect_data], sort=False, ignore_index=True)
                                elif column == "Global radiation":
                                    subselect_data = select_data.groupby(["Date"]).mean() * (86000 / 1000000)
                                    subselect_data.reset_index(inplace=True)
                                    subselect_data = subselect_data[["Date", column]]
                                    all_data = pd.concat([all_data, subselect_data], sort=False, ignore_index=True)
                                elif column == "Wind speed" and wind_ajustment:
                                    subselect_data = select_data.groupby(["Date"]).mean() * (1000 / 3600)
                                    subselect_data.reset_index(inplace=True)
                                    subselect_data = subselect_data[["Date", column]]
                                    all_data = pd.concat([all_data, subselect_data], sort=False, ignore_index=True)
                                elif column != "Date":
                                    subselect_data = select_data[["Date", column]]
                                    subselect_data.loc[:,column] = select_data.loc[:,column].astype(float)
                                    subselect_data = subselect_data.groupby(["Date"]).mean()
                                    subselect_data.reset_index(inplace=True)
                                    # subselect_data = subselect_data[["Date", column]]
                                    all_data = pd.concat([all_data, subselect_data], sort=False, ignore_index=True)

                            # elif "Global radiation" in select_data.columns:
                            #     select_data = select_data.groupby(["Date"]).mean() * (86000 / 1000000)
                            # elif "Wind speed" in select_data.columns:
                            #     select_data = select_data.groupby(["Date"]).mean() * (1000 / 3600)
                            # else:
                            #     select_data = select_data.groupby(["Date"]).mean()

            else:
                data = data_dict
                if all_data.empty:
                    for col in data.columns:
                        if col in variables_dict:
                            all_data[variables_dict[col]] = data[col]
                            if col == "Vitesse du vent, m/s":
                                wind_ajustment = False
                    if len(all_data.columns)> 1:
                        all_data["Date"] = pd.to_datetime(all_data['Date'].astype(str), format='%Y-%m-%d')
                        organized_data = pd.DataFrame()
                        for column in all_data.columns:
                            if column == "Precipitation":
                                subselect_data = all_data.groupby(["Date"]).sum()
                                subselect_data.reset_index(inplace=True)
                                subselect_data = subselect_data[["Date", column]]
                                organized_data = pd.concat([organized_data, subselect_data], sort=False, ignore_index=True)
                            elif column == "Global radiation":
                                subselect_data = all_data.groupby(["Date"]).mean() * (86000 / 1000000)
                                subselect_data.reset_index(inplace=True)
                                subselect_data = subselect_data[["Date", column]]
                                organized_data = pd.concat([organized_data, subselect_data], sort=False, ignore_index=True)
                            elif column == "Wind speed" and wind_ajustment:
                                subselect_data = all_data.groupby(["Date"]).mean() * (1000 / 3600)
                                subselect_data.reset_index(inplace=True)
                                subselect_data = subselect_data[["Date", column]]
                                organized_data = pd.concat([organized_data, subselect_data], sort=False, ignore_index=True)
                            elif column != "Date":
                                subselect_data = all_data.groupby(["Date"]).mean()
                                subselect_data.reset_index(inplace=True)
                                subselect_data = subselect_data[["Date", column]]
                                organized_data = pd.concat([organized_data, subselect_data], sort=False, ignore_index=True)
                        all_data = organized_data
                        # if "Precipitation" in all_data.columns:
                        #     all_data = all_data.groupby(["Date"]).sum()
                        #     all_data.reset_index(inplace=True)
                        # elif "Global radiation" in all_data.columns:
                        #     all_data = all_data.groupby(["Date"]).mean() * (86000/1000000)
                        # elif "Wind speed" in all_data.columns:
                        #     all_data = all_data.groupby(["Date"]).mean() * (1000/3600)
                        # else:
                        #     all_data = all_data.groupby(["Date"]).mean()
                    else:
                        all_data = pd.DataFrame()
                else:
                    select_data = pd.DataFrame()
                    for col in data.columns:
                        if col in variables_dict:
                            select_data[variables_dict[col]] = data[col]
                            if col == "Vitesse du vent, m/s":
                                wind_ajustment = False
                    if len(select_data.columns) > 1:
                        select_data["Date"] = pd.to_datetime(select_data['Date'].astype(str), format='%Y-%m-%d')

                        for column in select_data.columns:
                            if column == "Precipitation":
                                subselect_data = select_data.groupby(["Date"]).sum()
                                subselect_data.reset_index(inplace=True)
                                subselect_data = subselect_data[["Date", column]]
                                all_data = pd.concat([all_data, subselect_data], sort=False, ignore_index=True)
                            elif column == "Global radiation":
                                subselect_data = select_data.groupby(["Date"]).mean() * (86000 / 1000000)
                                subselect_data.reset_index(inplace=True)
                                subselect_data = subselect_data[["Date", column]]
                                all_data = pd.concat([all_data, subselect_data], sort=False, ignore_index=True)
                            elif column == "Wind speed" and wind_ajustment:
                                subselect_data = select_data.groupby(["Date"]).mean() * (1000 / 3600)
                                subselect_data.reset_index(inplace=True)
                                subselect_data = subselect_data[["Date", column]]
                                all_data = pd.concat([all_data, subselect_data], sort=False, ignore_index=True)
                            elif column != "Date":
                                subselect_data = select_data.groupby(["Date"]).mean()
                                subselect_data.reset_index(inplace=True)
                                subselect_data = subselect_data[["Date", column]]
                                all_data = pd.concat([all_data, subselect_data], sort=False, ignore_index=True)

                        # if "Precipitation" in select_data.columns:
                        #     select_data = select_data.groupby(["Date"]).sum()
                        #     select_data.reset_index(inplace=True)
                        # elif "Global radiation" in select_data.columns:
                        #     select_data = select_data.groupby(["Date"]).mean() * (86000 / 1000000)
                        # elif "Wind speed" in select_data.columns:
                        #     select_data = select_data.groupby(["Date"]).mean() * (1000 / 3600)
                        # else:
                        #     select_data = select_data.groupby(["Date"]).mean()
                        # all_data = pd.concat([all_data, select_data])
        else:
            print("##### File %s of %s Skipped: %s #####" % (n, len(all_files), file))
            n += 1

    print("Data extraction finished, start creating climate file")
    # all_data = all_data.dropna()
    all_data.reset_index(inplace=True, drop=True)
    df2 = all_data

    df2 = df2.groupby(by=df2.index, axis=0).apply(lambda g: g.mean() if isinstance(g.iloc[0, 0], numbers.Number) else g.iloc[0])

    all_data = df2

    # all_data = all_data.groupby(level=0, axis=1).apply(lambda x: x.apply(sjoin, axis=1))



    for variable in all_variable_needed:
        if variable not in all_data.columns:
            if variable =="Cloud cover":
                all_data[variable] = 0.65
            elif variable == "Relative humidity":
                all_data[variable] = 50
            elif variable == "Global radiation":
                all_data[variable] = np.nan


    climate_file = os.path.join(observation_path,  "climate_data.csv")
    all_data['Date'] = all_data['Date'].astype(str)

    first_column = all_data.pop('Date')
    all_data.insert(0, 'Date', first_column)

    all_data.to_csv(climate_file, index=False)

    # return all_data
    return "Climate Files generated"

if __name__ == "__main__":
    print("Hello")
    # print(extract_bathymetry(lakes=lakes_dict, rawnamefile="Bromont_Bathymetry.csv"))
    print(extract_climate())
    #extract_observations(lakes=lakes_dict)
