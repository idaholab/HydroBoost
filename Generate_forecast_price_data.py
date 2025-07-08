# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

# Name of the input file located in the Input_Spreadsheets directory. Do not include the .xlsx extension here
file_name = "Project_2"

# Functions used to generate the forecast
from Core_Models.Price_Forecasting.import_data import read_price_and_forecasting, read_water_flow, read_daily_water_flow, copy_Excel
from Core_Models.Price_Forecasting.perfect_foresight_model import create_perfect_foresight_forecast
from Core_Models.Price_Forecasting.mean_persistence_model import create_mean_persistence_forecast
from Core_Models.Price_Forecasting.additive_models import additive_model_no_regressors, additive_model_with_regressors
from Core_Models.Price_Forecasting.generate_figures import plot_forecast_model_errors

# Price data and forecast features (if provided)
price_df, forecast_features = read_price_and_forecasting(file_name)

# Perfect foresight and mean persistence forecast
# Dicts with keys: DA-LMP, Regulation Up, Regulation Down, Spinning Reserve
perfect_foresight_forecast_dict = create_perfect_foresight_forecast(price_df, file_name)
mean_persistence_forecast_dict = create_mean_persistence_forecast(perfect_foresight_forecast_dict, file_name, num_days=7)

# # Add the additive forecast models. Note: you need regressors to run the second model
additive_no_regressor_forecast = additive_model_no_regressors(price_df, perfect_foresight_forecast_dict, file_name)
additive_with_regressors_forecast = additive_model_with_regressors(price_df, perfect_foresight_forecast_dict, forecast_features, file_name)

# Currently adding Autoregressive Neural Network models

# Add water flow data to the generated_data directory
hourly_water = read_water_flow(file_name)
daily_water = read_daily_water_flow(file_name)

# The Julia optimization need additional information from the input Excel file
copy_Excel(file_name)

# Generate forecast errors figure
fig_errors = plot_forecast_model_errors(file_name, perfect_foresight_forecast_dict, mean_persistence_forecast_dict, additive_no_regressor_forecast, additive_with_regressors_forecast)