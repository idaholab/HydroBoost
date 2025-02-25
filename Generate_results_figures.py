# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import pandas as pd
import plotly.express as px

file_name = "Project_1"

df_hydro = pd.read_csv(f"Simulation_Results/{file_name}/ALEAF_HydroBoost_Project_1__hydro_dispatch.csv")

print(df_hydro)