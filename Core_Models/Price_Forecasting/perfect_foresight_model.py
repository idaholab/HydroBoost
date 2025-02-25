# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import os
import pandas as pd
import numpy as np
from datetime import timedelta

from Core_Models.Price_Forecasting.helper_functions import generate_path


def initialize_empty_dataframe_for_forecast(df):
    """
    Initialize an empty dataframe for the 7-day forecast. Index of dates and columns are hours.
    """
    # Need values for the dataframe index(forecast hours) and columns (dates)
    dates = list(set(df.index.date))
    dates = pd.to_datetime(dates)
    dates = dates.sort_values()

    hours = [i for i in range(0, 7*24)]

    x = pd.DataFrame(index=dates, columns=hours)

    return x


def create_perfect_foresight_forecast(price_df, file_name):
    """
    Creates a dict with the following: ['DA-LMP', 'Regulation Up', 'Regulation Down', 'Spinning Reserve']
    Input: dataframe constructed using read_ALEAF_input_spreadsheet function.
    Note: Writing to Excel file is slow.
    """

    print("Started perfect foresight model.")

    # Get the dataframe of Data from the input_data dict
    df = price_df.copy()

    def get_forecast(forecast, price_column):
        """
        Helper function to run each price: LMP, Reg up/down, Spinning reserve
        """

        for date in forecast.index:

            # Start with an empty array
            x = np.empty(7*24)

            # Find what the actual values are and add them to empty array
            temp = df.loc[date:date + timedelta(days=7) - timedelta(hours=1), price_column].to_numpy()  

            # The last week doesn't have a week of future data, so we will duplicate the last day
            if len(temp) < 168:
                # Get the last day of values and duplicate it
                vals = temp[-24:]
                vals = np.tile(vals, 5)
    
                # Append them and take the number of values needed
                temp = np.append(temp, vals)
                temp = temp[0: 168]

            x[:len(temp)] = temp      
            
            forecast.loc[date] = x

        return forecast
    
    # Initialize an empty dataframe for the forecast
    empty_forecast = initialize_empty_dataframe_for_forecast(df)

    # Get all the forecast using helper function 
    LMP_forecast = get_forecast(empty_forecast.copy(), "DA-LMP")
    Reg_up_forecast = get_forecast(empty_forecast.copy(), "Regulation Up")
    Reg_down_forecast = get_forecast(empty_forecast.copy(), "Regulation Down")
    Spin_forecast = get_forecast(empty_forecast.copy(), "Spinning Reserve")

    # Lets drop the timestamp from the column names M-d-Y H:mm:ss --> M-d-Y 
    dates_to_string = LMP_forecast.index.astype(str)

    # Want the index as strings not datetimes
    LMP_forecast.index = dates_to_string
    Reg_up_forecast.index = dates_to_string
    Reg_down_forecast.index = dates_to_string
    Spin_forecast.index = dates_to_string

    # HydroBoost optimization solver is set with date columns and hour index, so we will transpose
    LMP_forecast = LMP_forecast.T
    Reg_up_forecast = Reg_up_forecast.T
    Reg_down_forecast = Reg_down_forecast.T
    Spin_forecast = Spin_forecast.T

    # Now we will save the perfect foresight DA-LMP, Reg_up, Reg_down, Spin
    # Make sure the directory structure is there
    generate_path(["generated_data", file_name, "Market", "Perfect_foresight"])
    LMP_forecast.to_csv(f"Core_Models/HydroBoost/generated_data/{file_name}/Market/Perfect_foresight/DA_LMP.csv")
    Reg_up_forecast.to_csv(f"Core_Models/HydroBoost/generated_data/{file_name}/Market/Perfect_foresight/Regulation_up.csv")
    Reg_down_forecast.to_csv(f"Core_Models/HydroBoost/generated_data/{file_name}/Market/Perfect_foresight/Regulation_down.csv")
    Spin_forecast.to_csv(f"Core_Models/HydroBoost/generated_data/{file_name}/Market/Perfect_foresight/Spin.csv")

    print("Finished perfect foresight model.\n")

    # Return the results as a dict
    return {"DA-LMP": LMP_forecast, 
            "Regulation Up": Reg_up_forecast, 
            "Regulation Down": Reg_down_forecast,
            "Spinning Reserve": Spin_forecast}