
Lake Bromont Project
=====================

This project adapted the scripts from the [MyLake Milsbo lakes Project](https://github.com/biogeochemistry/MyLake_Milsbo_lakes.git) 
to simulate the distribution of phosphorus (P) in sediments. This lake is treated with lanthanum (Phoslock) to improve 
the water quality by reducing the P concentration in the water column and the algae bloom issues.
This project aims to simulate the present state (2019-2022) of the lake (water column and sediment) using the model 
MyLake-Sediment ([Markelov et al.,2019](https://doi.org/10.1029/2019JG005254)). For this project, we may add the impact 
of the lanthanum on the P transportation from the water column to the sediment to MyLake-Sediment. 

As done for the [MyLake Milsbo lakes Project](https://github.com/biogeochemistry/MyLake_Milsbo_lakes.git), 
I first simulate the water column by calibrating Temperature, Oxygen Concentration and Chlorophyll a. 
Second, I will calibrate MyLake-Sediment to simulate the P distribution.

I divided this project into three main parts:

- **Formatting (Part 1):**
    
    Production of the scripts formatting the raw data for the lake in files that the model can use.

- **Manual Calibration (Part 2):** 

    Production of scripts manually calibrating the water column and the sediment module for the lake. The goal is to 
    find the general values ​​of the model parameters to use when selecting the range and get the "guessed" values ​​
    that will be used when calibrating with an algorithm.

- **Calibration with algorithm (Part 3):**
    
    Production of the scripts leading to the calibration using an algorithm and generating the final results. 
These scripts will run the simulation, run the post-processing and generate the numbers for each lake analyzed.


## Repository layout ##

    .
    ├── IO                              # Formatted files for MyLake model and script to format information from Raw Data
    ├── MyLake-v2.0                     # Submodule of MyLake-v2.0, used to run MyLake model
    ├── Postproc_code                   # Output files from the model, Script to analyse/visualized data and final dataset. 
    ├── Raw Data                        # All information given for the simulations
    │   ├── Donnees_mouillage_Bromont   # observations files for the water column
    │   ├── meteo                       # Weather data
    │   ├── geo_data                    # GIS information for Bromont lake
    │── Sediment-v2.0                   # Submodule of Sediment-v2.0, used to run the sediment simulation for MyLake
    │── img                             # Exemple of figures and other pictures used in the readme files 
    │── obs                             # Files merging information from raw data to be ready to be formatted. 
    │                                     Included the formatted files of observations, ready to be used to evaluate the performance of the model.
    └── README.md


## Description of scripts ##

The following section describes each script used in this project. This includes the purposes of these scripts, 
the information needed, the functions included in the script, and how these scripts have been used.

The scripts included in this project are:

##### **Part 1**

- [x] extract\_information\_from\_raw\_data.py
- [x] create\_inputs.py

##### **Part 2**

- [x] manual\_calibration.py

##### **Part 3**

- [ ] final\_calibration.py (TODO)
- [ ] post\_processing.py (TODO)

### Part 1 ###

#### [**extract\_information\_from\_raw\_data.py**](OI/extract_information_from_raw_data.py) ####

Script extracting raw data for easier manipulation by create_input.py and verifying the information we have.
It extracts information from hypsometry and observations from the water column of the lake. 
It took the data from the folder Raw_Data and placed the output files in the obs folder.
This script includes three main functions: extract_hypsometry(), extract_observations(), and extract_climate() 
and the auxiliary functions:  get_all_files_in_directory(directory() and get_information_from_file_as_dataframe())

the main function call for the specific file (in the hypsometry's case) using the filename and all files in specific folders 
(for observations and climate) using the auxiliary functions. This script only compiles the information in one file for each main function, 
changes the unit if needed (see Raw data conversion into wanted units by MyLake), and saves those files in the obs folder. 

For more information on those functions, use *function.\_\_doc\_\_* :
``` {.}
# example for extract_hypsometry()
$ cd IO
$ python -c "from extract_information_from_raw_data import *; print(extract_hypsometry.__doc__)"
```

##### **Raw data conversion into wanted units by MyLake**

Mylake model uses daily input with specific units.  The raw data given for those simulations present multiple measures 
per day and different units than the ones used by the model. To resolve this issue, the below solutions have been used 
to convert the raw data into usable inputs.

**Table. Weather variables and the equations used to convert the raw data into the formatted input** 

| Variable          | Raw data Unit (multiple times a day)   | MyLake Unit (daily)   | Solution    |
| ----------------- | -------------                         | -------------------   | --------------------------------------------------------------------- |
| Global radiation  | no data given                         | MJ/(m^2 day)          | let the model estimate GR                                             |
| Cloud cover       | no data given                         | -                     | default value: 0.65                                                   |
| Air temperature   | deg. C                                |deg. C                 | $\overline{T}_{day}$                                                  |
| Relative humidity | %                                     | %                     | $\overline{H}_{day}$                                                  |
| Air pressure      | kpa                                   | mbar                  | $\overline{P}_{day}$ * 10                                             |
| Wind speed        | km/hr                                 | m/s                   | $\overline{W}_{day}$ * $\frac{1000 m * km^{-1}}{3600 s * hr^{-1}}$    |
| Precipitation     | mm                                    | mm/day                | $\sum_{time = 0}^{24} Prep(time)$                                     |



#### [**IO/create\_inputs\_model\_from\_obs.py**](IO/create_inputs_model_from_obs.py) ####

Script formatting the raw data in the correct format for MyLake. If called directly, it will run both lakes, giving 
initialization and input files from hypsometry, initial concentration, inflow, and meteorological data.

This script includes two Classes: Lake and Simulation_informations.

**Class Lake**: This class contains all information specific to the lake analysed(if multiple lakes are simulated). 
The attributes are the lake name (ask when the class is called) and the specific path to the input, observation, and output folders 
(generated from the attribute from the class Simulation_informations also ask when the class Lake is called). 
It also contains the class Simulation_informations as an attribute that can be shared with multiple lakes. 
This class has two methods: create_file.py and mylakeinit.py. 
The first one calls the functions to create two of the three files needed for the simulation: input_file and init_file 
(par_file is generated by the model when specific values are given during the manual calibration).
The second is the function call to create the file with the initial concentrations (init_file). 

**Simulation_informations**: This class contains all information specific to the simulation. 
The attributes are the specific path to the input, observation, and output folders (given by default for the present structure)
and the information for the simulations(start and end years, the depth resolution for the water column simulation and 
the switch to enable river inflow (turn off for this project since not given for the lake Bromont, if this data is added, 
we may turn the switch on)). It also contains the class simulation information as an attribute that can be shared with multiple lakes. 
This class has one method: mylakeinput.py. 
This function is called to create the file with the input information(climate variables and inflows if the switch for 
the inflow is turned on (generate inflow columns of 0 if not)). This function formats the data and estimates if needed, 
the missing data to generate complete time series of the desired period.  

Again, uses *function.\_\_doc\_\_* for more information on those functions.

##### **Estimation of missing data needed for simulation**

Meteorological data used for the simulation were from multiple stations and the dataset is missing daily values for 
some periods and variables. To estimate those values, two approaches were used: interpolation and the use of average daily value. 
I did the first interpolation considering the dates of the missing data and days with data before and after the days skipped.
It was first done with the restriction of interpolating a maximum of 3 consecutive days, to remove the possibility of a 
long period with constant values. The script after identified the days still without data and replace those values, 
when available, with an average daily value calculated from the raw data. Finally, for the days that are still missing 
(with days in a long period without data and never measure over the 5 years of data given), larger interpolation is 
done to fill the gaps and ensure a completed dataset.


### Part 2 ###

#### [**script\_manual\_calibratin.py**](script_manual_calibration.py)

This script allows testing different values ​​for the selected parameters for a specific variable fpr the water column (Temperature or Oxygen)
and for the sediments (). 

When the main script is run, the user will be prompted for input to determine what is calibrated, which lake is
calibrated, and what the backup options will be. Once these questions are answered, the main loop starts, and,
for each iteration, the script asks the user to select the value of each parameter. 
These values ​​are used to launch the MyLake model. Model output is saved if \"Save all output data\" is True and used
to generate a file comparing modeled and observed data. This file, saved if the \"Save comparison data\" option 
is True, is used to produce the figures (time series, profiles, and observed vs simulated comparison). 
It also calculates the calibration performance and prints the result to the console. 
It saved the numbers if the save numbers option was selected, and the iteration information 
(value of each parameter, performance analysis, score, and grades) is added to the report.

Just call the script in the console to start the manual calibration.

``` {.}
$ python script_manual_calibration.py
```

Or [**Use the Jupiter notebook script**](run_manual_calibration_from_jupyter_notebook.ipynb) to test run. 
Note that because of incompatibility, this version cannot launch the Matlab session to run MyLake but does not break,
 therefore the script continues using the last run data to calculate performances and generate the figures. 
 It will be resolved, but for the moment it can verify if everything has been installed correctly or showcase the script.


Example of figures generated (for lac Bromont, using the initial values
for parameters related to Temperature):

![The first digit generated. Give a comparison of modeled data to observations.](img/example_TOC_comparison_for_document.png)

Figure 1. Observed vs simulated temperature (squares), oxygen (circles), and chlorophyll a (triangles) for all depths 
(from red (0; surface) to blue (7; max depth)) of Lake Bromont. the fourth figure is a comparison of observations (squares)
 and simulated (lines) of the calibrated variable (here, Temperature) for the surface (black) and deep water (blue).

![Second digit generated. Give a time series of the three variables and profiles of these variables for different dates.](img/example_timeseries_and_profils.png)

Figure 2. Comparison of observed (Square: Temperature, Circle: Oxygen, Triangle: Chlorophyll a) and simulated (lines) 
values ​​for surface (black) and deep water (blue: Temperature, red: Oxygen, green: Chlorophyll a).
The first three figures present this comparison over time, while the other six figures present the water profiles 
(colored profiles: Observations, black profiles: Simulations) for spring, summer, and autumn withdrawals in 2020 and 2021.
