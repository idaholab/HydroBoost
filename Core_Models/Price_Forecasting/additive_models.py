# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import os
import pandas as pd
import numpy as np

# This can fix numpy issue based on version being used
np.float_ = np.float64

from datetime import timedelta
from prophet import Prophet

import logging
logging.getLogger('cmdstanpy').disabled = True

def additive_model_no_regressors(price_df, perfect_foresight_forecast_dict, file_name):
    print("Started additive model no regressors.")

    df_prophet = price_df[["DA-LMP"]].copy()

    df_prophet = df_prophet.reset_index()
    df_prophet = df_prophet.rename(columns={"Date-Time": "ds", "DA-LMP": "y"})
    m = Prophet(changepoint_prior_scale=1, 
                changepoint_range=1)
    m.fit(df_prophet)

    forecast = perfect_foresight_forecast_dict["DA-LMP"].T.copy()

    results = m.predict(m.make_future_dataframe(periods=24*8))["yhat"].to_numpy()[25:]

    for i, date in enumerate(forecast.index):
        temp = results[(i)*24:(i+6)*24]
        forecast.iloc[i, 24:] = temp

    # Transpose to get in format for HydroBoost
    forecast = forecast.T

    # We know the directory exists
    forecast.to_csv(f"Core_Models/HydroBoost/generated_data/{file_name}/Market/Additive_model_no_regressors.csv")

    print("Finished additive model no regressors.\n")

    return forecast



def additive_model_with_regressors(price_df, perfect_foresight_forecast_dict, forecast_features, file_name):

    print(f"Started additive model with regressors: {', '.join(forecast_features)}.")

    # Make sure there are forecast regressors to add
    if len(forecast_features) == 0:
        print("You did not provide any forecast features so there are no regressors to add to the model.")
        return None

    # Get the dataframe needed
    df_prophet = price_df.copy()
    df_prophet = df_prophet.drop(columns=["Regulation Up", "Regulation Down", "Spinning Reserve"])
    df_prophet = df_prophet.reset_index()
    df_prophet = df_prophet.rename(columns={"Date-Time": "ds", "DA-LMP": "y"})

    # Make the model
    m = Prophet(changepoint_prior_scale=1, changepoint_range=1)

    # Find and add the regressors to the model
    regressors = df_prophet.columns.to_list()[2:]

    for regressor in regressors:
        m.add_regressor(regressor)

    # Train the model
    m.fit(df_prophet)

    # Need to update our dataframe, need to extend dataframe at the end of the year 6 days 
    last_day = df_prophet.iloc[-24:, :]

    for i in range(7):
        temp = last_day.copy()
        temp['ds'] = temp["ds"] + timedelta(days=(i+1))
        df_prophet = pd.concat([df_prophet, temp])

    predictions = m.predict(df_prophet)["yhat"].to_numpy()[25:]

    forecast = perfect_foresight_forecast_dict["DA-LMP"].T.copy()

    for i, date in enumerate(forecast.index):
        temp = predictions[i*24:(i+6)*24]
        forecast.iloc[i, 24:] = temp

    # Transpose for format used in HydroBoost
    forecast = forecast.T

    # We know the directory exists
    forecast.to_csv(f"Core_Models/HydroBoost/generated_data/{file_name}/Market/Additive_model_with_regressors.csv")

    print(f"Finished additive model with regressors.\n")
    return forecast