# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

from io import BytesIO
import shutil
import pandas as pd
from Core_Models.Price_Forecasting.helper_functions import generate_path

def read_price_and_forecasting(file_name):

    print("Started data import.")
    # Needs to handle calls from python as well as the application
    if isinstance(file_name, str):
        df = pd.read_excel(f"Input_Spreadsheets/{file_name}.xlsx", 
                        sheet_name="Hourly Price and Forecasting", 
                        skiprows=[0], 
                        index_col="Date-Time", 
                        )
    elif hasattr(file_name, 'read'):
        # If the input has a 'read' method, treat it as a file-like object
        df = pd.read_excel(BytesIO(file_name.read()),
                           sheet_name="Hourly Price and Forecasting", 
                           skiprows=[0], 
                           index_col="Date-Time")
    else:
        raise ValueError("Invalid input type. Expected string or file-like object.")


    df = df.bfill()

    # Add check to make sure required columns are present:
    # DA-LMP, Regulation Up, Regulation Down, Spinning Reserve

    # Drop any column headers with no data
    df = df.dropna(axis=1, how='all')

    # Find if forecasting features were provided
    if len(df.columns) > 4:
        forecast_features = df.columns.to_list()[4:]
    else:
        forecast_features = None

    print("Finished importing data.\n")

    return df, forecast_features


def read_water_flow(file_name):
    df = pd.read_excel(f"Input_Spreadsheets/{file_name}.xlsx", 
                    sheet_name="Water Flow and Availability", 
                    skiprows=[0], 
                    index_col="Date-Time",
                    )

    # Drop generation columns that are empty
    df = df.dropna(axis=1, how='all')

    # Make sure the directory structure is there
    generate_path(["generated_data", file_name, "Hydro"])
    df.to_csv(f"Core_Models/HydroBoost/generated_data/{file_name}/Hydro/Hourly_flow.csv")

    return df

def read_daily_water_flow(file_name):
    df = pd.read_excel(f"Input_Spreadsheets/{file_name}.xlsx", 
                    sheet_name="Daily Flow Constraints", 
                    skiprows=[0], 
                    )

    df["Date-time"] = pd.to_datetime(df["Date-time"])
    df = df.set_index("Date-time")

    # Make sure the directory structure is there
    generate_path(["generated_data", file_name, "Hydro"])
    df.to_csv(f"Core_Models/HydroBoost/generated_data/{file_name}/Hydro/Daily_flow.csv")

    return df


def copy_Excel(file_name):
    shutil.copy(f"Input_Spreadsheets/{file_name}.xlsx", f"Core_Models/HydroBoost/generated_data/{file_name}")
    return None
