# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

"""
Additive forecasting models using Facebook Prophet for DA-LMP price forecasting.

This module provides two Prophet-based forecasting approaches:
1. No regressors: Time series decomposition using only historical DA-LMP values
2. With regressors: Time series decomposition including external features (e.g., temperature, demand)

Both models generate forecasts for a configurable number of days (default: 7 days, 168 hours).
When use_perfect_foresight_first_day=True (default):
- First day (hours 0-23): Perfect foresight (actual values)
- Remaining days: Prophet-based predictions
"""

import math
import os
from datetime import timedelta
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

# NumPy compatibility for older versions
np.float_ = np.float64

from prophet import Prophet

from .helper_functions import generate_path, parse_CLI_arguments
from .import_data import read_price_and_forecasting
from .perfect_foresight_model import create_perfect_foresight_forecast

# Module-level constants
HOURS_PER_DAY = 24


def additive_model_no_regressors(
    price_df: pd.DataFrame,
    perfect_dict: Dict[str, pd.DataFrame],
    project: str,
    num_forecast_days: int = 7,
    use_perfect_foresight_first_day: bool = True
) -> pd.DataFrame:
    """
    Create DA-LMP forecast using Prophet time series decomposition without external regressors.
    
    Uses Facebook Prophet to model the additive components of the time series (trend, 
    seasonality, holidays) based solely on historical DA-LMP values. Forecasts can optionally
    use perfect foresight for the first day and Prophet predictions for remaining days.
    
    Parameters
    ----------
    price_df : pd.DataFrame
        Historical price data with DatetimeIndex and columns including 'DA-LMP'.
        Index must be named 'Date-Time' and contain hourly datetime values.
    perfect_dict : Dict[str, pd.DataFrame]
        Dictionary of perfect foresight forecasts with keys including 'DA-LMP'.
        Used as a template for output structure (num_forecast_days*24 hours × N dates).
    project : str
        Project name (e.g., 'Model_A') for determining output file path.
    num_forecast_days : int, optional
        Number of days in forecast window (default: 7, resulting in 168 hours).
    use_perfect_foresight_first_day : bool, optional
        If True (default), hours 0-23 use perfect foresight values, hours 24+ use Prophet.
        If False, all hours use Prophet predictions.
    
    Returns
    -------
    pd.DataFrame
        Forecast DataFrame with:
        - Index: Hour numbers as strings ('0' to 'num_forecast_days*24-1')
        - Columns: Forecast start dates
        - Values: Forecasted DA-LMP prices ($/MWh)
        
        If use_perfect_foresight_first_day is True, hours 0-23 contain perfect 
        foresight values and remaining hours contain Prophet predictions.
        Otherwise, all hours contain Prophet predictions.
    
    Raises
    ------
    KeyError
        If 'DA-LMP' column is missing from price_df or perfect_dict.
    ValueError
        If price_df has insufficient data for Prophet fitting.
    RuntimeError
        If Prophet model fails to fit or predict.
    OSError
        If output file cannot be written.
    
    Notes
    -----
    - Prophet predicts with an extra buffer day beyond the forecast window for safety
    - If predictions run short for any date, they are padded with the last
      predicted value to maintain the forecast window
    - Output is saved to: generated_data/<project>/Market/Additive_Models/
      Additive_model_no_regressors.csv
    - Model uses default Prophet hyperparameters:
      * changepoint_prior_scale=1 (moderate trend flexibility)
      * changepoint_range=1 (changepoints throughout entire history)
    
    Examples
    --------
    >>> price_df, _ = read_price_and_forecasting('Model_A')
    >>> perfect_dict = create_perfect_foresight_forecast(price_df, 'Model_A')
    >>> # Default: 7-day forecast with perfect foresight for first day
    >>> forecast_df = additive_model_no_regressors(price_df, perfect_dict, 'Model_A')
    >>> # Custom: 5-day forecast with all Prophet predictions
    >>> forecast_df = additive_model_no_regressors(
    ...     price_df, perfect_dict, 'Model_A', num_forecast_days=5, 
    ...     use_perfect_foresight_first_day=False
    ... )
    """
    # Calculate forecast parameters
    forecast_hours = num_forecast_days * HOURS_PER_DAY
    
    # Input validation
    if 'DA-LMP' not in price_df.columns:
        raise KeyError("price_df must contain 'DA-LMP' column")
    
    if 'DA-LMP' not in perfect_dict:
        raise KeyError("perfect_dict must contain 'DA-LMP' key")
    
    if len(price_df) < HOURS_PER_DAY:
        raise ValueError(
            f"price_df must contain at least {HOURS_PER_DAY} hours of data for Prophet fitting"
        )
    
    try:
        # Prepare data for Prophet (requires 'ds' for datetime and 'y' for target)
        df_prophet = price_df[['DA-LMP']].reset_index().rename(
            columns={'Date-Time': 'ds', 'DA-LMP': 'y'}
        )
        
        # Initialize and fit Prophet model
        model = Prophet(
            changepoint_prior_scale=1,  # Controls trend flexibility
            changepoint_range=1  # Use full history for changepoints
        )
        model.fit(df_prophet)
        
        # Create future dataframe for predictions with buffer
        # Predict extra days beyond forecast window for safety
        prediction_periods = HOURS_PER_DAY * (num_forecast_days + 1)
        future = model.make_future_dataframe(periods=prediction_periods, freq='h')
        
        # Generate predictions and extract relevant portion
        predictions = model.predict(future)
        # Skip only first hour to get all forecast predictions
        predicted_values = predictions['yhat'].iloc[1:].to_numpy()
        
    except Exception as e:
        raise RuntimeError(f"Prophet model fitting or prediction failed: {str(e)}") from e
    
    # Use perfect foresight template structure (will be transposed for filling)
    template = perfect_dict['DA-LMP'].T.copy()
    num_forecasts = template.shape[0]
    
    # Fill each forecast date with corresponding predictions
    for i in range(num_forecasts):
        start_idx = i * HOURS_PER_DAY
        forecast_block = predicted_values[start_idx:start_idx + forecast_hours]
        
        # Handle edge case: pad if predictions are insufficient
        if forecast_block.shape[0] < forecast_hours:
            pad_length = forecast_hours - forecast_block.shape[0]
            last_value = forecast_block[-1] if forecast_block.size > 0 else np.nan
            forecast_block = np.concatenate([
                forecast_block,
                np.full(pad_length, last_value)
            ])
        
        # Apply predictions based on mode
        if use_perfect_foresight_first_day:
            # Keep perfect foresight for hours 0-23, use Prophet for hours 24+
            template.iloc[i, HOURS_PER_DAY:] = forecast_block[HOURS_PER_DAY:]
        else:
            # Use Prophet predictions for all hours
            template.iloc[i, :] = forecast_block
    
    # Transpose back to standard format (hours as rows, dates as columns)
    result = template.T
    
    # Save output
    try:
        out_dir = generate_path([project, 'Market', 'Additive_Models'])
        output_path = os.path.join(out_dir, 'Additive_model_no_regressors.csv')
        result.to_csv(output_path, index=True)
        print(f"✓ Saved additive model (no regressors) to: {output_path}")
    except OSError as e:
        raise OSError(f"Failed to write output file: {str(e)}") from e
    
    return result


def additive_model_with_regressors(
    price_df: pd.DataFrame,
    perfect_dict: Dict[str, pd.DataFrame],
    features: Optional[List[str]],
    project: str,
    num_forecast_days: int = 7,
    use_perfect_foresight_first_day: bool = True
) -> Optional[pd.DataFrame]:
    """
    Create DA-LMP forecast using Prophet with external regressors (e.g., temperature, demand).
    
    Uses Facebook Prophet to model the time series including external features as regressors.
    This captures relationships between DA-LMP and other variables. Future regressor values
    are generated by repeating the last 24 hours cyclically for each forecast day.
    
    Parameters
    ----------
    price_df : pd.DataFrame
        Historical price data with DatetimeIndex and columns including 'DA-LMP' and
        all features specified in the features list. Index must be named 'Date-Time'.
    perfect_dict : Dict[str, pd.DataFrame]
        Dictionary of perfect foresight forecasts with keys including 'DA-LMP'.
        Used as a template for output structure (num_forecast_days*24 hours * N dates).
    features : Optional[List[str]]
        List of column names in price_df to use as external regressors.
        If None or empty, the function returns None without generating forecasts.
    project : str
        Project name (e.g., 'Model_A') for determining output file path.
    num_forecast_days : int, optional
        Number of days in forecast window (default: 7, resulting in 168 hours).
    use_perfect_foresight_first_day : bool, optional
        If True (default), hours 0-23 use perfect foresight values, hours 24+ use Prophet.
        If False, all hours use Prophet predictions.
    
    Returns
    -------
    Optional[pd.DataFrame]
        Forecast DataFrame with same structure as additive_model_no_regressors,
        or None if no features are provided.
        
        If returned:
        - Index: Hour numbers as strings ('0' to 'num_forecast_days*24-1')
        - Columns: Forecast start dates
        - Values: Forecasted DA-LMP prices ($/MWh)
    
    Raises
    ------
    KeyError
        If 'DA-LMP' column or any feature column is missing from price_df,
        or if 'DA-LMP' key is missing from perfect_dict.
    ValueError
        If price_df has fewer than 24 hours of data (needed for regressor extension).
    RuntimeError
        If Prophet model fails to fit or predict.
    OSError
        If output file cannot be written.
    
    Examples
    --------
    >>> price_df, features = read_price_and_forecasting('Model_A')
    >>> perfect_dict = create_perfect_foresight_forecast(price_df, 'Model_A')
    >>> if features:
    ...     # Default: 7-day forecast with perfect foresight for first day
    ...     forecast_df = additive_model_with_regressors(
    ...         price_df, perfect_dict, features, 'Model_A'
    ...     )
    ...     print(f"Used regressors: {features}")
    """
    # Calculate forecast parameters
    forecast_hours = num_forecast_days * HOURS_PER_DAY
    
    # Handle case with no regressors
    if not features:
        print("ℹ No regressors provided; skipping additive model with regressors.")
        return None
    
    # Input validation
    if 'DA-LMP' not in price_df.columns:
        raise KeyError("price_df must contain 'DA-LMP' column")
    
    missing_features = [f for f in features if f not in price_df.columns]
    if missing_features:
        raise KeyError(
            f"price_df is missing required feature columns: {missing_features}"
        )
    
    if 'DA-LMP' not in perfect_dict:
        raise KeyError("perfect_dict must contain 'DA-LMP' key")
    
    if len(price_df) < HOURS_PER_DAY:
        raise ValueError(
            f"price_df must contain at least {HOURS_PER_DAY} hours of data "
            f"for regressor extension"
        )
    
    try:
        # Prepare data including all regressors
        df_prophet = price_df.reset_index().rename(
            columns={'Date-Time': 'ds', 'DA-LMP': 'y'}
        )
        columns_to_keep = ['ds', 'y'] + features
        df_prophet = df_prophet[columns_to_keep]
        
        # Initialize Prophet and add each regressor
        model = Prophet(
            changepoint_prior_scale=1,
            changepoint_range=1
        )
        for feature in features:
            model.add_regressor(feature)
        
        model.fit(df_prophet)
        
        # Extend dataframe with future dates and cyclically repeated regressor values
        # Strategy: Repeat the last 24 hours of regressor values for each future day
        last_24_hours = df_prophet.iloc[-HOURS_PER_DAY:].copy()
        extended_df = df_prophet.copy()
        
        # Extend for all forecast days
        for day_offset in range(1, num_forecast_days):
            future_block = last_24_hours.copy()
            future_block['ds'] = future_block['ds'] + timedelta(days=day_offset)
            extended_df = pd.concat([extended_df, future_block], ignore_index=True)
        
        # Generate predictions using extended dataframe
        predictions = model.predict(extended_df)
        # Skip only first hour to get all forecast predictions
        predicted_values = predictions['yhat'].iloc[1:].to_numpy()
        
    except Exception as e:
        raise RuntimeError(
            f"Prophet model with regressors fitting or prediction failed: {str(e)}"
        ) from e
    
    # Use perfect foresight template structure
    template = perfect_dict['DA-LMP'].T.copy()
    num_forecasts = template.shape[0]
    
    # Fill each forecast date with corresponding predictions
    for i in range(num_forecasts):
        start_idx = i * HOURS_PER_DAY
        forecast_block = predicted_values[start_idx:start_idx + forecast_hours]
        
        # Handle edge case: pad if predictions are insufficient
        if forecast_block.shape[0] < forecast_hours:
            pad_length = forecast_hours - forecast_block.shape[0]
            last_value = forecast_block[-1] if forecast_block.size > 0 else np.nan
            forecast_block = np.concatenate([
                forecast_block,
                np.full(pad_length, last_value)
            ])
        
        # Apply predictions based on mode
        if use_perfect_foresight_first_day:
            # Keep perfect foresight for hours 0-23, use Prophet for hours 24+
            template.iloc[i, HOURS_PER_DAY:] = forecast_block[HOURS_PER_DAY:]
        else:
            # Use Prophet predictions for all hours
            template.iloc[i, :] = forecast_block
    
    # Transpose back to standard format
    result = template.T
    
    # Save output
    try:
        out_dir = generate_path([project, 'Market', 'Additive_Models'])
        output_path = os.path.join(out_dir, 'Additive_model_with_regressors.csv')
        result.to_csv(output_path, index=True)
        print(f"✓ Saved additive model (with regressors) to: {output_path}")
        print(f"  Regressors used: {', '.join(features)}")
    except OSError as e:
        raise OSError(f"Failed to write output file: {str(e)}") from e
    
    return result


if __name__ == '__main__':
    """
    Command-line interface for testing additive models.
    
    Usage:
        python additive_models.py Model_A
        python additive_models.py Model_B
        python additive_models.py Model_C
    
    This will:
    1. Load price data and features from Input_Spreadsheets/{project}.xlsx
    2. Generate perfect foresight forecast
    3. Generate additive model without regressors
    4. Generate additive model with regressors (if features are available)
    """
    # Parse command line arguments for project name
    project = parse_CLI_arguments()
    
    print("=" * 70)
    print(f"ADDITIVE MODELS FORECAST GENERATION - {project}")
    print("=" * 70)
    
    try:
        # Load data
        print(f"\nLoading data for {project}...")
        price_df, features = read_price_and_forecasting(project)
        print(f"✓ Loaded {len(price_df)} hours of price data")
        if features:
            print(f"✓ Found {len(features)} regressor features: {', '.join(features)}")
        else:
            print("ℹ No regressor features found")
        
        # Generate perfect foresight forecast (needed as template)
        print(f"\nGenerating perfect foresight forecast...")
        perfect_dict = create_perfect_foresight_forecast(price_df, project)
        num_forecasts = perfect_dict['DA-LMP'].shape[1]
        print(f"✓ Generated {num_forecasts} forecast dates")
        
        # Generate additive model without regressors
        print(f"\nGenerating additive model (no regressors)...")
        result_no_reg = additive_model_no_regressors(
            price_df, perfect_dict, project, 
            num_forecast_days=7, 
            use_perfect_foresight_first_day=True
        )
        print(f"✓ Generated forecast shape: {result_no_reg.shape}")
        
        # Generate additive model with regressors (if available)
        print(f"\nGenerating additive model (with regressors)...")
        result_with_reg = additive_model_with_regressors(
            price_df, perfect_dict, features, project,
            num_forecast_days=7,
            use_perfect_foresight_first_day=True
        )
        if result_with_reg is not None:
            print(f"✓ Generated forecast shape: {result_with_reg.shape}")
            
            # Show comparison statistics
            print(f"\n" + "=" * 70)
            print("FORECAST COMPARISON (Hours 24-167)")
            print("=" * 70)
            
            # Compare predictions for the first forecast date
            first_forecast = result_no_reg.columns[0]
            no_reg_pred = result_no_reg.loc['24':'167', first_forecast].astype(float)
            with_reg_pred = result_with_reg.loc['24':'167', first_forecast].astype(float)
            
            print(f"Sample forecast date: {first_forecast}")
            print(f"  No regressors  - Mean: ${no_reg_pred.mean():.2f}, "
                  f"Std: ${no_reg_pred.std():.2f}")
            print(f"  With regressors - Mean: ${with_reg_pred.mean():.2f}, "
                  f"Std: ${with_reg_pred.std():.2f}")
            print(f"  Mean absolute difference: ${abs(no_reg_pred - with_reg_pred).mean():.2f}")
        
        print(f"\n" + "=" * 70)
        print("FORECAST GENERATION COMPLETE")
        print("=" * 70)
        
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        raise
