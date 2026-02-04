# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED


import os
import pandas as pd
import numpy as np
from typing import Dict

# Make it runnable with the IDE "Run" button (works both as package and as script)
try:
    from .helper_functions import generate_path
except ImportError:
    from helper_functions import generate_path


HOURS_PER_DAY = 24


def create_mean_persistence_forecast(
    perfect_foresight_forecast_dict: Dict[str, pd.DataFrame], 
    file_name: str,
    num_forecast_days: int = 7,
    num_mean_persistence_days: int = 7,
    use_perfect_foresight_first_day: bool = True
) -> Dict[str, pd.DataFrame]:
    """
    Creates forecasts using mean persistence with configurable forecast horizon (default: 7 days).
    
    This function generates forecasts by computing the mean of historical values
    from the previous num_mean_persistence_days and repeating that daily pattern for the 
    forecast horizon. Optionally, the first 24 hours can use perfect foresight values.
    
    Parameters
    ----------
    perfect_foresight_forecast_dict : Dict[str, pd.DataFrame]
        Dictionary of DataFrames containing perfect foresight forecasts with keys:
        - 'DA-LMP': Day-ahead locational marginal pricing
        - 'RT-LMP': Real-time locational marginal pricing
        - 'Regulation Up': Regulation up reserve prices
        - 'Regulation Down': Regulation down reserve prices
        - 'Spinning Reserve': Spinning reserve prices
        Each DataFrame has:
        - Index: Hour numbers (0 to num_forecast_days*24-1) as strings
        - Columns: Forecast start dates
        - Values: Perfect foresight prices
    file_name : str
        Project identifier (e.g., 'Model_A', 'Model_B') used for output directory
        structure under .../generated_data/{file_name}/Market/Mean_persistence/
    num_forecast_days : int, optional
        Number of days in forecast window (default: 7, resulting in 168 hours).
        This should match the forecast window used in perfect_foresight_forecast_dict.
    num_mean_persistence_days : int, optional
        Number of days to look back for computing mean persistence (default: 7).
        For the first forecast, uses available history up to num_mean_persistence_days.
    use_perfect_foresight_first_day : bool, optional
        If True (default), the first 24 hours (day 1) of each forecast use perfect 
        foresight values, and remaining hours use mean persistence.
        If False, all forecast hours use mean persistence.
        Default is True for more realistic day-ahead forecasting.
    
    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary mapping price type names to forecast DataFrames. Each forecast
        DataFrame has:
        - Index: Hour numbers (0 to num_forecast_days*24-1) as strings
        - Columns: Forecast start dates as strings
        - Values: Forecasted prices using mean persistence
    
    Raises
    ------
    ValueError
        If required keys are missing from perfect_foresight_forecast_dict or
        if num_days is less than 1
    IOError
        If there are issues writing forecast files to disk
    
    Notes
    -----
    Output CSV files are saved to:
    .../generated_data/{file_name}/Market/Mean_persistence/
        - DA_LMP.csv
        - RT_LMP.csv
        - Regulation_up.csv
        - Regulation_down.csv
        - Spin.csv
        - Delta_RU.csv
        - Delta_RD.csv
    """
    # Validate input parameters
    required_keys = ["DA-LMP", "Regulation Up", "Regulation Down", "Spinning Reserve"]
    missing_keys = set(required_keys) - set(perfect_foresight_forecast_dict.keys())
    if missing_keys:
        raise ValueError(f"Missing required keys in perfect_foresight_forecast_dict: {missing_keys}")
    
    # Check if RT-LMP is available (needed for Model C)
    has_rt_lmp = "RT-LMP" in perfect_foresight_forecast_dict
    if not has_rt_lmp:
        print("WARNING: RT-LMP not found in perfect foresight dict. Using DA-LMP as proxy.")
    
    if num_forecast_days < 1:
        raise ValueError(f"num_forecast_days must be at least 1, got {num_forecast_days}")
    
    if num_mean_persistence_days < 1:
        raise ValueError(f"num_mean_persistence_days must be at least 1, got {num_mean_persistence_days}")

    pf_dict = dict(perfect_foresight_forecast_dict)
    if "Delta RU" not in pf_dict:
        print("WARNING: Delta RU not found in perfect foresight dict. Filling with zeros.")
        pf_dict["Delta RU"] = pf_dict["DA-LMP"].copy()
        pf_dict["Delta RU"].loc[:, :] = 0.0
    if "Delta RD" not in pf_dict:
        print("WARNING: Delta RD not found in perfect foresight dict. Filling with zeros.")
        pf_dict["Delta RD"] = pf_dict["DA-LMP"].copy()
        pf_dict["Delta RD"].loc[:, :] = 0.0

    def get_forecast(df_tpl: pd.DataFrame) -> pd.DataFrame:
        """
        Apply mean persistence to create forecast.
        """
        df_result = df_tpl.copy()
        forecast_days = df_tpl.shape[1] // HOURS_PER_DAY  # Calculate from template shape
        
        for idx in range(len(df_result.index)):
            # Determine which previous forecasts to average
            if idx == 0:
                # First forecast: use only the first day's values
                means = df_result.iloc[0, :HOURS_PER_DAY].to_numpy()
            elif idx < num_mean_persistence_days:
                # Not enough history: use all available previous forecasts
                means = df_result.iloc[0:idx, :HOURS_PER_DAY].mean().to_numpy()
            else:
                # Full history available: use last num_mean_persistence_days forecasts
                means = df_result.iloc[idx-num_mean_persistence_days:idx, :HOURS_PER_DAY].mean().to_numpy()
            
            # Apply mean persistence based on mode
            if use_perfect_foresight_first_day:
                # Keep perfect foresight for first 24 hours, use persistence for rest
                df_result.iloc[idx, HOURS_PER_DAY:] = np.tile(means, forecast_days - 1)
            else:
                # Use mean persistence for all hours
                df_result.iloc[idx, :] = np.tile(means, forecast_days)
        
        return df_result.T

    # Build forecasts for each price type
    forecasts = {
        "DA-LMP": get_forecast(pf_dict["DA-LMP"].T.copy()),
        "Regulation Up": get_forecast(pf_dict["Regulation Up"].T.copy()),
        "Regulation Down": get_forecast(pf_dict["Regulation Down"].T.copy()),
        "Spinning Reserve": get_forecast(pf_dict["Spinning Reserve"].T.copy()),
        "Delta RU": get_forecast(pf_dict["Delta RU"].T.copy()),
        "Delta RD": get_forecast(pf_dict["Delta RD"].T.copy()),
    }
    
    # Add RT-LMP (use actual if available, otherwise use DA-LMP as proxy)
    if has_rt_lmp:
        forecasts["RT-LMP"] = get_forecast(pf_dict["RT-LMP"].T.copy())
    else:
        # Fallback: Use DA-LMP as proxy for RT-LMP (Model C requires this)
        forecasts["RT-LMP"] = forecasts["DA-LMP"].copy()

    # Prepare output directory
    out_dir = generate_path([file_name, "Market", "Mean_persistence"])

    # Save CSVs with error handling
    output_files = {
        "DA-LMP": "DA_LMP.csv",
        "RT-LMP": "RT_LMP.csv",  # Now always generated (Model C requires this)
        "Regulation Up": "Regulation_up.csv",
        "Regulation Down": "Regulation_down.csv",
        "Spinning Reserve": "Spin.csv",
        "Delta RU": "Delta_RU.csv",
        "Delta RD": "Delta_RD.csv",
    }
    
    try:
        for price_type, filename in output_files.items():
            forecasts[price_type].to_csv(os.path.join(out_dir, filename))
        print(f"→ Mean persistence files saved to: {out_dir}")
    except IOError as e:
        print(f"Error saving mean persistence files to {out_dir}: {e}")
        raise

    return forecasts


if __name__ == "__main__":
    """
    Click Run to execute (no CLI required).

    This will:
    1. Load price data from Input_Spreadsheets/{project}.xlsx
    2. Generate perfect foresight forecast
    3. Generate mean persistence forecast
    """

    # Imports that work both as package and as script
    try:
        from .import_data import read_price_and_forecasting
    except ImportError:
        from import_data import read_price_and_forecasting

    try:
        from .perfect_foresight_model import create_perfect_foresight_forecast
    except ImportError:
        from perfect_foresight_model import create_perfect_foresight_forecast

    # Choose the project here (uncomment one)
    # project = "Model_A"
    # project = "Model_B"
    project = "Model_C"

    # Load data and generate perfect foresight forecasts first
    print(f"\n=== Generating Mean Persistence Forecasts for {project} ===")
    price_df, _ = read_price_and_forecasting(project)
    perfect_dict = create_perfect_foresight_forecast(price_df, project)

    # Generate mean persistence forecasts
    mean_persist_dict = create_mean_persistence_forecast(
        perfect_dict,
        project,
        num_forecast_days=7,
        num_mean_persistence_days=7,
        use_perfect_foresight_first_day=True
    )

    # Print results to console for testing
    print("\n" + "="*70)
    print("MEAN PERSISTENCE FORECAST SUMMARY")
    print("="*70)
    print(f"Configuration:")
    print(f"  Forecast horizon: {7} days")
    print(f"  Look-back window: {7} days")
    print(f"  First day mode: Perfect foresight (hours 0-23)")
    print(f"  Remaining days: Mean persistence")

    for name, df in mean_persist_dict.items():
        print(f"\n{name}:")
        print(f"  Shape: {df.shape[0]} hours × {df.shape[1]} forecasts")
        print(f"  Forecast dates: {df.columns[0]} to {df.columns[-1]}")
        print(f"  First forecast statistics:")
        print(f"    Mean: ${df.iloc[:, 0].mean():>8.2f}")
        print(f"    Min:  ${df.iloc[:, 0].min():>8.2f}")
        print(f"    Max:  ${df.iloc[:, 0].max():>8.2f}")

    # Show comparison between perfect foresight and mean persistence for first forecast
    print(f"\n" + "-"*70)
    print("First Forecast Comparison (DA-LMP):")
    pf_first = perfect_dict["DA-LMP"].iloc[:, 0]
    mp_first = mean_persist_dict["DA-LMP"].iloc[:, 0]

    diff = np.abs(pf_first.iloc[HOURS_PER_DAY:] - mp_first.iloc[HOURS_PER_DAY:])
    print(f"  Day 1 (hours 0-23): Perfect foresight values used")
    print(f"  Remaining days: Mean persistence applied")
    print(f"    Mean absolute error: ${diff.mean():.2f}")
    print(f"    Max absolute error: ${diff.max():.2f}")
    print("="*70)
    print("Forecast generation complete\n")
