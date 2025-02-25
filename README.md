# HydroBoost
HydroBoost is an optimization solver for hydro power and battery storage systems. Its primary objective is to optimize revenue based on hourly dispatch of the hydro and battery assets with forecast energy market prices. It can be applied to both run-of-river and reservoir storage sites.


## Initial Setup to run HydroBoost
Currently, HydroBoost is ran using both Python and Julia. You must setup both environments to run the tool. If you do not have either installed you can find downloads for  Python here https://www.anaconda.com/download/success, and Julia can be found here https://julialang.org/downloads/.

This will walk you through setting up the environment and running the tool. There are four steps:

* Create an input spreadsheet for you project (examples provided)
* Generating forecast data needed for optimization (Python)
* Running the optimization (Julia)
* Generating interactive figures (Python)


## Step 1: Input Spreadsheet
You need to fill out an input spreadsheet with the characteristics of your system. The examples provided will guide you in filling out the input spreadsheet.


## Step 2: Generating Forecast Price Data
First, you need to open the Generate_forecast_data.py file and change the file_name to the name of the input spreadsheet for your system.

It is recommended that you use a virtual environment but you can run using the global environment. To create a virtual environment and install the needed packages:
```
python3 -m venv <environment_name>
```
Then you need to activate the virtual environment you just created.<br>
On Linux/Mac
```
source <environment_name>/bin/activate
```
On Windows
```
<environment_name>\Scripts\activate
```
Once the virtual environment is active you can install all the needed Python packages with
```
pip install -r requirements.txt
```

Now you can generate the forecasting price data needed for the optimization. This is done with the command line using 
```
python Generate_forecast_price_data.py
```
The forecast data will be generated and the results will be located in the Core_Models/HydroBoost/generated_data/ directory. The HydroBoost Julia optimization solver will use these files.


## Step 3: Optimization Solver

Now you will run the solver using the data you generated above. First, you need to install the needed packages and add the HydroBoost_Sim package. We will use the Julia REPL for this. In the terminal command line type the command julia and hit enter to start the REPL.

Now add the needed packages with the following lines:
```julia
using Pkg
Pkg.add(path=pwd(), subdir="Core_Models/HydroBoost")
```

### Execute HydroBoost

To run HydroBoost, you simply call the "execute_HydroBoost_model" function, which will execute the simulation process:

```julia
using Pkg 
using HydroBoost_Sim
HydroBoost_Sim.execute_HydroBoost_model(project_id)

```
where project_id is the name of your input file.

The results will be save in the Simulation_Results directory.

## Step 4: Generating Figures
This is not necessary but serves as a convenience to generate many interactive plots to view the results data. 

You will need to activate the virtual environment again as in step 2 and run the Generate_results_figures.py file. Again, ensure that the project_name is set to the correct name. The figures will be added to the Simulation_Results directory.