# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import os
import pandas as pd
import numpy as np
from datetime import timedelta

from Core_Models.Price_Forecasting.helper_functions import generate_path

def create_mean_persistence_forecast(perfect_foresight_forecast_dict, file_name, num_days=7):
    """
    Creates an Excel sheet with 7-day forecast using mean persistence
    Input: resulting dataframe form perfect foresight model, number of days to use for mean persistance 
    Note: Writing to Excel file is slow.
    """

    print("Started mean persistence forecast")

    def get_forecast(forecast_df):
        """
        Helper function to run each price: LMP, Reg up/down, Spinning reserve
        """
        
        for idx, date in enumerate(forecast_df.index):
            
            # Use the previous days mean to forecast days 2-7, with day-ahead being correct values
            if idx == 0:
                hour_means = forecast_df.iloc[0, 0:24].to_numpy()
            elif idx < num_days:
                hour_means = forecast_df.iloc[0:idx, 0:24].mean().to_numpy()
            else:
                hour_means = forecast_df.iloc[idx-num_days:idx, 0:24].mean().to_numpy()

            hour_means = np.tile(hour_means, 6)

            forecast_df.iloc[idx, 24:] = hour_means

        return forecast_df.T
    
    # Get mean persistance for LMP, reg up/down and spinning reserve
    LMP_forecast = get_forecast(perfect_foresight_forecast_dict["DA-LMP"].T.copy())
    Reg_up_forecast = get_forecast(perfect_foresight_forecast_dict["Regulation Up"].T.copy())
    Reg_down_forecast = get_forecast(perfect_foresight_forecast_dict["Regulation Down"].T.copy())
    Spin_forecast = get_forecast(perfect_foresight_forecast_dict["Spinning Reserve"].T.copy())
    
    # Now we will save the mean persistence DA-LMP, Reg_up, Reg_down, Spin
    # Make sure the directory structure is there
    generate_path(["generated_data", file_name, "Market", "Mean_persistence"])
    LMP_forecast.to_csv(f"Core_Models/HydroBoost/generated_data/{file_name}/Market/Mean_persistence/DA_LMP.csv")
    Reg_up_forecast.to_csv(f"Core_Models/HydroBoost/generated_data/{file_name}/Market/Mean_persistence/Regulation_up.csv")
    Reg_down_forecast.to_csv(f"Core_Models/HydroBoost/generated_data/{file_name}/Market/Mean_persistence/Regulation_down.csv")
    Spin_forecast.to_csv(f"Core_Models/HydroBoost/generated_data/{file_name}/Market/Mean_persistence/Spin.csv")

    print("Finished mean persistence forecast.\n")

    # Return the results as a dict
    return {"DA-LMP": LMP_forecast, 
            "Regulation Up": Reg_up_forecast, 
            "Regulation Down": Reg_down_forecast,
            "Spinning Reserve": Spin_forecast,}