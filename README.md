
Lake Bromont Project
=====================

This project adapted the scripts from the [MyLake Milsbo lakes Project](https://github.com/biogeochemistry/MyLake_Milsbo_lakes.git) 
to simulate the distribution of phosphorus (P) in sediments. This lake is treated with lanthanum (Phoslock) to improve 
the water quality by reducing the P concentration in the water column and the algae blooms issues.
This project aims to simulate the present state (2019-2022) of the lake (water column and sediment) using the model 
MyLake-Sediment ([Markelov et al.,2019](https://doi.org/10.1029/2019JG005254)). For this project, we may add the impact 
of the lanthanum on the P transportation from the water column to the sediment to MyLake-Sediment. 

As done for the [MyLake Milsbo lakes Project](https://github.com/biogeochemistry/MyLake_Milsbo_lakes.git), 
I first simulate the water column by calibrating Temperature, Oxygen Concentration and Chlorophyll a. 
Second, I will calibrate MyLake-Sediment to simulate the P distribution.



**Raw data conversion into wanted units by MyLake**

| Variable          | Raw data Unit (multiple time a day)   | MyLake Unit (daily)   | Solution    |
| ----------------- | -------------                         | -------------------   | --------------------------------------------------------------------- |
| Global radiation  |  no data given                                     | MJ/(m^2 day)          |    let the model estimate GR    |
| Cloud cover       | no data given                                     | -                     | default value: 0.65      |
|  Air temperature  | deg. C                                |deg. C                 | $\overline{T}_{day}$                                                  |
| Relative humidity | %                                     | %                     | $\overline{H}_{day}$                                                  |
| Air pressure      | kpa                                   | mbar                  | $\overline{P}_{day}$ * 10                                             |
| Wind speed        | km/hr                                 | m/s                   | $\overline{W}_{day}$ * $\frac{1000 m * km^{-1}}{3600 s * hr^{-1}}$    |
| Precipitation     | mm                                    | mm/day                | $\sum_{time = 0}^{24} Prep(time)$                                     |



<!--- 
I divided this project into three main parts:

- **Formatting (Part 1):**
    
    Production of the scripts formatting the raw data for the two lakes in files that the model can use.

- **Manual Calibration (Part 2):** 

    Production of scripts manually calibrating the water column and the sediment module for each lake. The goal is to 
    find the general values ​​of the model parameters to use when selecting the range and to get the "guessed" values ​​
    that will be used when calibrating with an algorithm.

- **Calibration with algorithm (Part 3):**
    
    Production of the scripts leading to the calibration using an algorithm and generating the final results. 
These scripts will run the simulation, run the post-processing and generate the numbers for each lake analyzed.


## Repository layout ##

    .
    ├── IO                      # Formatted files for MyLake model and script to format information from Raw Data
    ├── MyLake-v2.0             # Submodule of MyLake-v2.0, used to run MyLake model
    ├── Postproc_code           # Ooutput files from the model, Script to analyse/visualized data and final dataset. 
    ├── Raw Data                # All information given for the simulations
    │   ├── Data                # observations files for the water column and the inflows.
    │   ├── Data descrition     # Description of the information given in the Raw data folder
    │   ├── GIS_milsbo_sweden   # GIS information for the two lakes analysed
    │── Sediment-v2.0           # Submodule of Sediment-v2.0, used to run the sediment simulation for MyLake
    │── img                     # Exemple of figures and other pictures used in the readme files 
    │── obs                     # Files merging information from raw data to be ready to be formatted. Included the formatted files of observations, ready to be used to evaluation the performance of the model.
    └── README.md




## Description of scripts ##

The following section describes each script used in this project. This includes the purposes of these scripts, 
the information needed, the functions included in the script and how these scripts have been used.

The script included in this project are:

##### **Part 1**

- [x] structure\_observations_data.py
- [x] create\_inputs.py

##### **Part 2**

- [x] manual\_calibration.py

##### **Part 3**

- [ ] final\_calibration.py (TODO)
- [ ] post\_processing.py (TODO)

### Part 1 ###

#### [**structure\_observations_data.py**](IO/structure_observations_data.py) ####

Script extracting raw data for easier manipulation by create_input.py and verifying the information we have.
It extracts information from hypsometry and observations from the water column and the inflows for the two lakes. 
It took the data from the folder Raw_Data and placed the output files in the obs folder.
This script includes three functions: extract_hypsometry(), extract_observations(), and extract_inflows(). 

For more informations on those functions, use *function.\_\_doc\_\_* :
``` {.}
# example for extract_hypsometry()
$ cd IO
$ python -c "from structure_observations_data import *; print(extract_hypsometry.__doc__)"
```

**Solution for the multiple inlets and missing inflow data**

MyLake model uses one value by variable for each day analyzed, therefore the values from multiple inlets need to be 
summarised to be employed by the model.

To calculate the overall velocity of the inflow given to the model, I did the sum of the velocity of all available 
inlets for each sampling day. 
For the overall temperature, pH, and for the other variables' concentration, I calculated a weighted average of all 
available inlets, considering each inlet velocity for each sampling day.
Also, with those two lakes, the outlet of the lake Övre is one inlet of the lake Nedre. 
To consider this connection, I used the simulation of the surface layer of Övre to estimate the inflow concentration of 
this Nedre's inlet, and I weighted those concentrations in the inflows' average using the velocity information of this 
inlet.

The figure below presents the equations used for each variable estimated for both lakes.   

![eexemple for the equations used to stimated inflows with map of lakes](img/Lake_inflow_solutions.jpg)


As mentioned in the figure, missing data between sampling campaigns were interpolated using the information from the 
sample dates. 
I estimated dates outside the range of sampling campaigns using the information from the closest sample date. 
Finally, to consider the smaller flow during winter (season not sampled), I divided the inflow information estimated 
for those periods by two.
I used [structure\_observations_data.py](IO/structure_observations_data.py) to calculate the inflow average, 
and [create\_inputs.py](IO/create_inputs.py) to estimate of the missing data.

#### [**create\_inputs.py**](IO/create_inputs.py)

Script formatting the raw data in the correct format for MyLake. If called directly, it will run both lakes, giving 
initialization and input files from hypsometry, initial concentration, inflow, and meteorological data.

This script include two Class: Lake and simulation_information.

**Class Lake**: This class 


### Part 2 ###

#### [**script\_manual\_calibratin.py**](script_manual_calibration.py)

This script allows to test different values ​​for the selected parameters for a specific variable (Temperature, Oxygen 
or Chlorophyll a). At the moment I only wrote this script to calibrate the temperature in the water column of Lake Nedre.
 I will add the other options in the future.

When the main script is run, the user will be prompted for input to determine what is calibrated, which lake is
 calibrated, and what the backup options will be. Once these questions are answered, the main loop starts and,
  for each iteration, the script asks the user to select the value of each parameter. 
  These values ​​are used to launch the MyLake model. Model output is saved if \"Save all output data\" is True and used
   to generate a file comparing modeled and observed data. This file, saved if the \"Save comparison data\" option 
   is True, is used to produce the figures (time series, profiles and observed vs simulated comparison). 
   It also calculates the calibration performance and prints the result to the console. 
   It saved the numbers if the save numbers option was selected, and the iteration information 
   (value of each parameter, performance analysis, score, and grades) is added to the report.

Just call the script in the console to start the manual calibration.

``` {.}
$ python script_manual_calibration.py
```

Or [**Use the Jupiter notebook script**](run_manual_calibration_from_jupyter_notebook.ipynb) to test run. 
Note that because of incompatibility, this version cannot launch the Matlab session to run MyLake but does not break,
 therefore the script continues using the last run data to calculate performances and to generate the figures. 
 It will be resolved, but for the moment it can verify if everything has been installed correctly or showcase the script.


Example of figures generated (for Nedre, using the initial values
for parameters related to Temperature):

![First digit generated. Give a comparison of modeled data to observations.](img/example_TOC_comparison_for_document.png)

Figure 1. Observed vs simulated temperature (squares), oxygen (circles) and chlorophyll a (triangles) for all depths 
(from red (0; surface) to blue (6;max depth)) of Lake Nedre. the fourth figure is a comparison of observations (squares)
 and simulated (lines) of the calibrated variable (here, Temperature) for the surface (black) and deep water (blue).

![Second digit generated. Give a time series of the three variables and profiles of these variables for different dates.](img/example_timeseries_and_profils.png)

Figure 2. Comparison of observed (Square: Temperature, Circle: Oxygen, Triangle: Chlorophyll a) and simulated (lines) 
values ​​for surface (black) and deep water (blue: Temperature, red: Oxygen, green: Chlorophyll a).
The first three figures present this comparison over time, while the other six figures present the water profiles 
(colored profiles: Observations, black profiles: Simulations) for spring, summer and autumn withdrawals in 2019 and 2020.

--->