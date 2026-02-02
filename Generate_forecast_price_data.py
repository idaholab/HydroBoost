# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

# ============================================================
# SELECT MODELS TO RUN
# 
# NOTE: Perfect Foresight and Mean Persistence models are 
#       always run by default as they are required baselines.
# ============================================================
RUN_ADDITIVE_NO_REGRESSORS = True
RUN_ADDITIVE_WITH_REGRESSORS = False
RUN_AUTOREGRESSIVE_NN_NO_REGRESSORS = True
RUN_AUTOREGRESSIVE_NN_WITH_REGRESSORS = False
RUN_RANDOM_FOREST_NO_REGRESSORS = False
RUN_RANDOM_FOREST_WITH_REGRESSORS = False

# Functions used to generate the forecast
from Core_Models.Price_Forecasting.helper_functions import parse_CLI_arguments
from Core_Models.Price_Forecasting.import_data import read_price_and_forecasting, read_water_flow, read_daily_water_flow, read_RES_profiles, copy_Excel
from Core_Models.Price_Forecasting.perfect_foresight_model import create_perfect_foresight_forecast
from Core_Models.Price_Forecasting.mean_persistence_model import create_mean_persistence_forecast
from Core_Models.Price_Forecasting.additive_models import additive_model_no_regressors, additive_model_with_regressors
from Core_Models.Price_Forecasting.autoregressive_NN_model import autoregressive_NN_no_regressors, autoregressive_NN_with_regressors
from Core_Models.Price_Forecasting.random_forest_model import random_forest_no_regressors, random_forest_with_regressors
from Core_Models.Price_Forecasting.generate_figures import plot_forecast_model_errors


# Parse command line arguments to get the input filename
file_name = parse_CLI_arguments()

# Price data and forecast features (if provided)
price_df, forecast_features = read_price_and_forecasting(file_name)

# Perfect foresight and mean persistence forecast
# Dicts with keys: DA-LMP, Regulation Up, Regulation Down, Spinning Reserve
perfect_foresight_forecast_dict = create_perfect_foresight_forecast(price_df, file_name)
mean_persistence_forecast_dict = create_mean_persistence_forecast(perfect_foresight_forecast_dict, file_name, num_mean_persistence_days=7)

# Add the additive forecast models. Note: you need regressors in your inputs spreadsheet to run the second model
additive_no_regressors_forecast = (
    additive_model_no_regressors(price_df, perfect_foresight_forecast_dict, file_name) 
    if RUN_ADDITIVE_NO_REGRESSORS else None
)
additive_with_regressors_forecast = (
    additive_model_with_regressors(price_df, perfect_foresight_forecast_dict, forecast_features, file_name) 
    if RUN_ADDITIVE_WITH_REGRESSORS else None
)

# Add the autoregressive neural network models. Note: you need regressors in your inputs spreadsheet to run the second model
autoregressive_NN_no_regressors_forecast = (
    autoregressive_NN_no_regressors(price_df, perfect_foresight_forecast_dict, file_name) 
    if RUN_AUTOREGRESSIVE_NN_NO_REGRESSORS else None
)
autoregressive_NN_with_regressors_forecast = (
    autoregressive_NN_with_regressors(price_df, perfect_foresight_forecast_dict, forecast_features, file_name) 
    if RUN_AUTOREGRESSIVE_NN_WITH_REGRESSORS else None
)

# Add the random forest models. Note: you need regressors in your inputs spreadsheet to run the second model
random_forest_no_regressors_forecast = (
    random_forest_no_regressors(price_df, perfect_foresight_forecast_dict, file_name) 
    if RUN_RANDOM_FOREST_NO_REGRESSORS else None
)
random_forest_with_regressors_forecast = (
    random_forest_with_regressors(price_df, perfect_foresight_forecast_dict, forecast_features, file_name) 
    if RUN_RANDOM_FOREST_WITH_REGRESSORS else None
)

# Add water flow data to the generated_data directory
hourly_water = read_water_flow(file_name)
daily_water = read_daily_water_flow(file_name)

# Add RES profiles to the generated_data directory
res_profiles = read_RES_profiles(file_name)

# The Julia optimization needs additional information from the input Excel file so we copy it into the dirctory
copy_Excel(file_name)

# Generate forecast errors figure
fig_errors = plot_forecast_model_errors(
    file_name, 
    perfect_foresight_forecast_dict, 
    mean_persistence_forecast_dict, 
    additive_no_regressors_forecast, 
    additive_with_regressors_forecast,
    autoregressive_NN_no_regressors_forecast,
    autoregressive_NN_with_regressors_forecast,
    random_forest_no_regressors_forecast,
    random_forest_with_regressors_forecast
)
