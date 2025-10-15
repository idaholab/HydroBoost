# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED


import os
import math
import pandas as pd
import numpy as np
from datetime import timedelta
from typing import Dict
from .helper_functions import generate_path

HOURS_PER_DAY = 24


def initialize_empty_dataframe_for_forecast(df: pd.DataFrame, num_forecast_days: int = 7) -> pd.DataFrame:
    """
    Initialize an empty DataFrame for a forecast of specified length.
    """

    dates = sorted(pd.to_datetime(list({d.date() for d in df.index})))
    forecast_hours = num_forecast_days * HOURS_PER_DAY
    hours = list(range(forecast_hours))
    return pd.DataFrame(index=dates, columns=hours)


def create_perfect_foresight_forecast(
    price_df: pd.DataFrame, 
    file_name: str, 
    num_forecast_days: int = 7
) -> Dict[str, pd.DataFrame]:
    """
    Generates perfect-foresight forecasts for electricity market prices.
    
    This function creates forecasts by extracting actual historical values for each
    forecast window. If insufficient data exists at the end of the dataset, the last
    available day is repeated to fill the forecast window.
    
    Parameters
    ----------
    price_df : pd.DataFrame
        Historical price data with datetime index and columns:
        - 'DA-LMP': Day-ahead locational marginal pricing
        - 'Regulation Up': Regulation up reserve prices
        - 'Regulation Down': Regulation down reserve prices
        - 'Spinning Reserve': Spinning reserve prices
    file_name : str
        Project identifier (e.g., 'Model_A', 'Model_B') used for output directory
        structure under .../generated_data/{file_name}/Market/Perfect_foresight/
    num_forecast_days : int, optional
        Number of days in forecast window (default: 7, resulting in 168 hours)
    
    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary mapping price type names to forecast DataFrames. Each forecast
        DataFrame has:
        - Index: Hour numbers (0 to num_forecast_days*24-1) as strings
        - Columns: Forecast start dates as strings
        - Values: Forecasted prices for each hour
    
    Raises
    ------
    ValueError
        If required price columns are missing from price_df
    IOError
        If there are issues writing forecast files to disk
    
    Notes
    -----
    - Output CSV files are saved to:
      .../generated_data/{file_name}/Market/Perfect_foresight/
        * DA_LMP.csv
        * Regulation_up.csv
        * Regulation_down.csv
        * Spin.csv
    - Default forecast horizon is 7 days (168 hours)
    - Custom forecast lengths can be specified via num_forecast_days parameter
    
    Examples
    --------
    >>> # Generate default 7-day forecast
    >>> forecasts = create_perfect_foresight_forecast(price_df, 'Model_A')
    >>> # Generate custom 5-day forecast
    >>> forecasts = create_perfect_foresight_forecast(price_df, 'Model_A', num_forecast_days=5)
    """
    # Calculate forecast hours based on num_forecast_days
    forecast_hours = num_forecast_days * HOURS_PER_DAY

    # Validate input DataFrame has required columns
    required_columns = ["DA-LMP", "Regulation Up", "Regulation Down", "Spinning Reserve"]
    missing_cols = set(required_columns) - set(price_df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    df = price_df.copy()

    def get_forecast(template: pd.DataFrame, col: str) -> pd.DataFrame:
        """
        Extract actual values for each forecast window.
        """

        for date in template.index:
            # Extract actual values for forecast_hours from df
            end_date = date + timedelta(days=num_forecast_days, hours=-1)
            arr = df.loc[date:end_date, col].to_numpy()
            
            # Handle case where we don't have enough data (end of dataset)
            if arr.size < forecast_hours:
                # Determine how many hours to pad
                hours_needed = forecast_hours - arr.size
                
                # Get the last available day for padding
                if arr.size >= HOURS_PER_DAY:
                    last_day = arr[-HOURS_PER_DAY:]
                else:
                    # Edge case: less than 24 hours available
                    # Tile the available data to create a full day
                    repeats_for_day = math.ceil(HOURS_PER_DAY / arr.size)
                    last_day = np.tile(arr, repeats_for_day)[:HOURS_PER_DAY]
                
                # Calculate how many times to repeat the last day
                repeats = math.ceil(hours_needed / HOURS_PER_DAY)
                arr = np.concatenate([arr, np.tile(last_day, repeats)])[:forecast_hours]
            
            template.loc[date] = arr
        
        return template.T

    # Build forecasts for each price type
    empty = initialize_empty_dataframe_for_forecast(df, num_forecast_days)
    
    forecasts = {
        "DA-LMP": get_forecast(empty.copy(), "DA-LMP"),
        "Regulation Up": get_forecast(empty.copy(), "Regulation Up"),
        "Regulation Down": get_forecast(empty.copy(), "Regulation Down"),
        "Spinning Reserve": get_forecast(empty.copy(), "Spinning Reserve"),
    }

    # Convert index to string dates for CSV output
    for fc in forecasts.values():
        fc.index = fc.index.astype(str)

    # Ensure output folder exists
    out_dir = generate_path([file_name, "Market", "Perfect_foresight"])

    # Save CSVs with error handling
    output_files = {
        "DA-LMP": "DA_LMP.csv",
        "Regulation Up": "Regulation_up.csv",
        "Regulation Down": "Regulation_down.csv",
        "Spinning Reserve": "Spin.csv",
    }
    
    try:
        for price_type, filename in output_files.items():
            forecasts[price_type].to_csv(os.path.join(out_dir, filename))
        print(f"→ Perfect-foresight files saved to: {out_dir}")
    except IOError as e:
        print(f"Error saving forecast files to {out_dir}: {e}")
        raise

    return forecasts


if __name__ == "__main__":
    """
    Command-line interface for testing perfect foresight models.
    
    Usage:
        python perfect_foresight_models.py Model_A
        python perfect_foresight_models.py Model_B
        python perfect_foresight_models.py Model_C
    
    This will:
    1. Load price data from Input_Spreadsheets/{project}.xlsx
    2. Generate perfect foresight forecast
    """
    
    from .import_data import read_price_and_forecasting
    from .helper_functions import parse_CLI_arguments
    
    # Use CLI to get project name
    project = parse_CLI_arguments()
    
    # Load data and generate forecasts
    print(f"\n=== Generating Perfect Foresight Forecasts for {project} ===")
    price_df, _ = read_price_and_forecasting(project)
    
    # Generate with default 7-day forecast
    forecasts = create_perfect_foresight_forecast(price_df, project, num_forecast_days=7)
    
    # Print results to console for testing
    print("\n" + "="*70)
    print("PERFECT FORESIGHT FORECAST SUMMARY")
    print("="*70)
    
    for name, df in forecasts.items():
        print(f"\n{name}:")
        print(f"  Shape: {df.shape[0]} hours × {df.shape[1]} forecasts")
        print(f"  Forecast dates: {df.columns[0]} to {df.columns[-1]}")
        print(f"  First forecast statistics:")
        print(f"    Mean: ${df.iloc[:, 0].mean():>8.2f}")
        print(f"    Min:  ${df.iloc[:, 0].min():>8.2f}")
        print(f"    Max:  ${df.iloc[:, 0].max():>8.2f}")
    
    # Verify last forecast completeness
    last_forecast_col = forecasts["DA-LMP"].columns[-1]
    last_forecast = forecasts["DA-LMP"][last_forecast_col]
    expected_hours = 7 * HOURS_PER_DAY
    print(f"\n" + "-"*70)
    print(f"Last Forecast Validation ({last_forecast_col}):")
    print(f"  Total hours: {len(last_forecast)}/{expected_hours}")
    print(f"  Non-null values: {last_forecast.notna().sum()}")
    print(f"  Data completeness: {'✓ Complete' if len(last_forecast) == expected_hours else '✗ Incomplete'}")
    print("="*70)
    print("Forecast generation complete\n")
