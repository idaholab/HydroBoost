# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

"""
Autoregressive neural network forecasting models using NeuralProphet for DA-LMP price forecasting.

This module provides two NeuralProphet-based forecasting approaches:
1. No regressors: Autoregressive neural network using only historical DA-LMP values
2. With regressors: Autoregressive neural network with external features (e.g., temperature, demand)

Both models generate forecasts for a configurable number of days (default: 7 days, 168 hours).
When use_perfect_foresight_first_day=True (default):
- First day (hours 0-23): Perfect foresight (actual values)
- Remaining days: NeuralProphet-based predictions
"""

import os
import warnings
from datetime import timedelta
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from neuralprophet import NeuralProphet

from helper_functions import generate_path, parse_CLI_arguments
from import_data import read_price_and_forecasting
from perfect_foresight_model import create_perfect_foresight_forecast

# Module-level constants
HOURS_PER_DAY = 24

# Suppress NeuralProphet warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=FutureWarning)


def autoregressive_NN_no_regressors(
    price_df: pd.DataFrame,
    perfect_dict: Dict[str, pd.DataFrame],
    project: str,
    num_forecast_days: int = 7,
    use_perfect_foresight_first_day: bool = True,
    n_lags: int = 24,
    epochs: int = 50,
    learning_rate: float = 0.01,
    batch_size: int = 32
) -> pd.DataFrame:
    """
    Create DA-LMP forecast using NeuralProphet autoregressive neural network without external regressors.
    
    Uses NeuralProphet to model the time series with an autoregressive neural network (AR-Net)
    that learns complex patterns from historical DA-LMP values. The model uses past values 
    as features to predict future values.
    
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
        If True (default), hours 0-23 use perfect foresight values, hours 24+ use NeuralProphet.
        If False, all hours use NeuralProphet predictions.
    n_lags : int, optional
        Number of lagged observations to use as autoregressive features (default: 24).
        E.g., n_lags=24 means the model uses the previous 24 hours to predict.
    epochs : int, optional
        Number of training epochs (default: 50).
    learning_rate : float, optional
        Learning rate for neural network training (default: 0.01).
    batch_size : int, optional
        Batch size for training (default: 32).
    
    Returns
    -------
    pd.DataFrame
        Forecast DataFrame with:
        - Index: Hour numbers as strings ('0' to 'num_forecast_days*24-1')
        - Columns: Forecast start dates
        - Values: Forecasted DA-LMP prices ($/MWh)
        
        If use_perfect_foresight_first_day is True, hours 0-23 contain perfect 
        foresight values and remaining hours contain NeuralProphet predictions.
        Otherwise, all hours contain NeuralProphet predictions.
    
    Raises
    ------
    KeyError
        If 'DA-LMP' column is missing from price_df or perfect_dict.
    ValueError
        If price_df has insufficient data for NeuralProphet fitting.
    RuntimeError
        If NeuralProphet model fails to fit or predict.
    OSError
        If output file cannot be written.
    
    Notes
    -----
    - NeuralProphet uses autoregressive features (n_lags previous values) for prediction
    - Training is performed silently (progress output suppressed)
    - Model uses AR-Net architecture with configurable hyperparameters
    - Output is saved to: generated_data/<project>/Market/Autoregressive_NN/
      Autoregressive_NN_no_regressors.csv
    
    Examples
    --------
    >>> price_df, _ = read_price_and_forecasting('Model_A')
    >>> perfect_dict = create_perfect_foresight_forecast(price_df, 'Model_A')
    >>> # Default: 7-day forecast with 24-hour lags
    >>> forecast_df = autoregressive_NN_no_regressors(price_df, perfect_dict, 'Model_A')
    >>> # Custom: 5-day forecast with 48-hour lags, 100 epochs
    >>> forecast_df = autoregressive_NN_no_regressors(
    ...     price_df, perfect_dict, 'Model_A', 
    ...     num_forecast_days=5,
    ...     n_lags=48,
    ...     epochs=100
    ... )
    """
    # Calculate forecast parameters
    forecast_hours = num_forecast_days * HOURS_PER_DAY
    
    # Input validation
    if 'DA-LMP' not in price_df.columns:
        raise KeyError("price_df must contain 'DA-LMP' column")
    
    if 'DA-LMP' not in perfect_dict:
        raise KeyError("perfect_dict must contain 'DA-LMP' key")
    
    if len(price_df) < n_lags + HOURS_PER_DAY:
        raise ValueError(
            f"price_df must contain at least {n_lags + HOURS_PER_DAY} hours of data "
            f"for NeuralProphet fitting with n_lags={n_lags}"
        )
    
    try:
        # Prepare data for NeuralProphet (requires 'ds' for datetime and 'y' for target)
        df_neural = price_df[['DA-LMP']].reset_index().rename(
            columns={'Date-Time': 'ds', 'DA-LMP': 'y'}
        )
        
        # Initialize NeuralProphet model with autoregressive components
        model = NeuralProphet(
            n_lags=n_lags,
            n_forecasts=forecast_hours,  # Forecast all hours at once
            epochs=epochs,
            learning_rate=learning_rate,
            batch_size=batch_size,
            yearly_seasonality=True,
            weekly_seasonality=True,
            daily_seasonality=True
        )
        
        # Fit the model
        metrics = model.fit(df_neural, freq='h')
        
        # Generate predictions
        # Create future dataframe with buffer for safety (similar to Prophet approach)
        prediction_periods = HOURS_PER_DAY * (num_forecast_days + 1)
        future = model.make_future_dataframe(
            df_neural,
            periods=prediction_periods,
            n_historic_predictions=len(df_neural)
        )
        
        # Generate forecast
        forecast = model.predict(future)
        
        # Extract predictions - NeuralProphet with n_forecasts creates columns yhat1, yhat2, ..., yhat{n_forecasts}
        # For rolling forecasts: each row contains a full forecast horizon
        # We want to extract 365 forecasts of 168 hours each (one per day)
        predicted_values = []
        
        # Get future predictions only (skip historical predictions)
        future_forecast = forecast.iloc[len(df_neural):].copy()
        
        # Extract forecasts - for each forecast origin, get the 168-hour forecast
        for i in range(len(future_forecast)):
            forecast_row = []
            # Extract yhat1 through yhat{forecast_hours} for this forecast origin
            for h in range(1, forecast_hours + 1):
                col_name = f'yhat{h}'
                if col_name in future_forecast.columns:
                    forecast_row.append(future_forecast[col_name].iloc[i])
                else:
                    # If column doesn't exist, pad with NaN
                    forecast_row.append(np.nan)
            
            predicted_values.extend(forecast_row)
            
            # Only take first 365 forecasts (one per day of year)
            if i >= 364:
                break
        
        predicted_values = np.array(predicted_values)
        
    except Exception as e:
        raise RuntimeError(f"NeuralProphet model fitting or prediction failed: {str(e)}") from e
    
    # Use perfect foresight template structure (will be transposed for filling)
    template = perfect_dict['DA-LMP'].T.copy()
    num_forecasts = template.shape[0]
    
    # Fill each forecast date with corresponding predictions
    for i in range(min(num_forecasts, len(predicted_values) // forecast_hours)):
        start_idx = i * forecast_hours
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
            # Keep perfect foresight for hours 0-23, use NeuralProphet for hours 24+
            template.iloc[i, HOURS_PER_DAY:] = forecast_block[HOURS_PER_DAY:]
        else:
            # Use NeuralProphet predictions for all hours
            template.iloc[i, :] = forecast_block
    
    # Transpose back to standard format (hours as rows, dates as columns)
    result = template.T
    
    # Save output
    try:
        out_dir = generate_path([project, 'Market', 'Autoregressive_NN'])
        output_path = os.path.join(out_dir, 'Autoregressive_NN_no_regressors.csv')
        result.to_csv(output_path, index=True)
        print(f"✓ Saved autoregressive NN model (no regressors) to: {output_path}")
    except OSError as e:
        raise OSError(f"Failed to write output file: {str(e)}") from e
    
    return result


def autoregressive_NN_with_regressors(
    price_df: pd.DataFrame,
    perfect_dict: Dict[str, pd.DataFrame],
    features: Optional[List[str]],
    project: str,
    num_forecast_days: int = 7,
    use_perfect_foresight_first_day: bool = True,
    n_lags: int = 24,
    epochs: int = 50,
    learning_rate: float = 0.01,
    batch_size: int = 32
) -> Optional[pd.DataFrame]:
    """
    Create DA-LMP forecast using NeuralProphet AR-Net with external regressors (e.g., temperature, demand).
    
    Uses NeuralProphet to model the time series with an autoregressive neural network that
    includes external features as additional inputs. Future regressor values are generated
    by repeating the last 24 hours cyclically for each forecast day.
    
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
        If True (default), hours 0-23 use perfect foresight values, hours 24+ use NeuralProphet.
        If False, all hours use NeuralProphet predictions.
    n_lags : int, optional
        Number of lagged observations for autoregressive features (default: 24).
    epochs : int, optional
        Number of training epochs (default: 50).
    learning_rate : float, optional
        Learning rate for neural network training (default: 0.01).
    batch_size : int, optional
        Batch size for training (default: 32).
    
    Returns
    -------
    Optional[pd.DataFrame]
        Forecast DataFrame with same structure as autoregressive_NN_no_regressors,
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
        If price_df has insufficient data for NeuralProphet fitting.
    RuntimeError
        If NeuralProphet model fails to fit or predict.
    OSError
        If output file cannot be written.
    
    Examples
    --------
    >>> price_df, features = read_price_and_forecasting('Model_A')
    >>> perfect_dict = create_perfect_foresight_forecast(price_df, 'Model_A')
    >>> if features:
    ...     # Default: 7-day forecast with features
    ...     forecast_df = autoregressive_NN_with_regressors(
    ...         price_df, perfect_dict, features, 'Model_A'
    ...     )
    ...     print(f"Used regressors: {features}")
    """
    # Calculate forecast parameters
    forecast_hours = num_forecast_days * HOURS_PER_DAY
    
    # Handle case with no regressors
    if not features:
        print("ℹ No regressors provided; skipping autoregressive NN with regressors.")
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
    
    if len(price_df) < n_lags + HOURS_PER_DAY:
        raise ValueError(
            f"price_df must contain at least {n_lags + HOURS_PER_DAY} hours of data "
            f"for NeuralProphet fitting"
        )
    
    try:
        # Prepare data including all regressors
        df_neural = price_df.reset_index().rename(
            columns={'Date-Time': 'ds', 'DA-LMP': 'y'}
        )
        columns_to_keep = ['ds', 'y'] + features
        df_neural = df_neural[columns_to_keep]
        
        # Initialize NeuralProphet model
        model = NeuralProphet(
            n_lags=n_lags,
            n_forecasts=1,
            epochs=epochs,
            learning_rate=learning_rate,
            batch_size=batch_size,
            yearly_seasonality=True,
            weekly_seasonality=True,
            daily_seasonality=True
        )
        
        # Add each regressor
        for feature in features:
            model.add_future_regressor(feature)
        
        # Fit the model
        metrics = model.fit(df_neural, freq='h')
        
        # Extend dataframe with future dates and cyclically repeated regressor values
        # Strategy: Repeat the last 24 hours of regressor values for each future day
        last_24_hours = df_neural.iloc[-HOURS_PER_DAY:].copy()
        extended_df = df_neural.copy()
        
        # Extend for all forecast days
        for day_offset in range(1, num_forecast_days + 1):
            future_block = last_24_hours.copy()
            future_block['ds'] = future_block['ds'] + timedelta(days=day_offset)
            extended_df = pd.concat([extended_df, future_block], ignore_index=True)
        
        # Generate predictions using extended dataframe
        predictions = model.predict(extended_df)
        # Skip only first hour to get all forecast predictions
        predicted_values = predictions['yhat1'].iloc[1:].to_numpy()
        
    except Exception as e:
        raise RuntimeError(
            f"NeuralProphet model with regressors fitting or prediction failed: {str(e)}"
        ) from e
    
    # Use perfect foresight template structure
    template = perfect_dict['DA-LMP'].T.copy()
    num_forecasts = template.shape[0]
    
    # Fill each forecast date with corresponding predictions
    for i in range(min(num_forecasts, len(predicted_values) // forecast_hours)):
        start_idx = i * forecast_hours
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
            # Keep perfect foresight for hours 0-23, use NeuralProphet for hours 24+
            template.iloc[i, HOURS_PER_DAY:] = forecast_block[HOURS_PER_DAY:]
        else:
            # Use NeuralProphet predictions for all hours
            template.iloc[i, :] = forecast_block
    
    # Transpose back to standard format
    result = template.T
    
    # Save output
    try:
        out_dir = generate_path([project, 'Market', 'Autoregressive_NN'])
        output_path = os.path.join(out_dir, 'Autoregressive_NN_with_regressors.csv')
        result.to_csv(output_path, index=True)
        print(f"✓ Saved autoregressive NN model (with regressors) to: {output_path}")
        print(f"  Regressors used: {', '.join(features)}")
    except OSError as e:
        raise OSError(f"Failed to write output file: {str(e)}") from e
    
    return result


if __name__ == '__main__':
    """
    Command-line interface for testing autoregressive NN models.
    
    Usage:
        python autoregressive_NN_model.py Model_A
        python autoregressive_NN_model.py Model_B
        python autoregressive_NN_model.py Model_C
    
    This will:
    1. Load price data and features from Input_Spreadsheets/{project}.xlsx
    2. Generate perfect foresight forecast
    3. Generate autoregressive NN model without regressors
    4. Generate autoregressive NN model with regressors (if features are available)
    """
    # Parse command line arguments for project name
    project = parse_CLI_arguments()
    
    print("=" * 70)
    print(f"AUTOREGRESSIVE NEURAL NETWORK FORECAST GENERATION - {project}")
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
        
        # Generate autoregressive NN model without regressors
        print(f"\nGenerating autoregressive NN model (no regressors)...")
        print(f"  Training neural network (n_lags=24, epochs=50)...")
        result_no_reg = autoregressive_NN_no_regressors(
            price_df, perfect_dict, project, 
            num_forecast_days=7, 
            use_perfect_foresight_first_day=True,
            n_lags=24,
            epochs=50
        )
        print(f"✓ Generated forecast shape: {result_no_reg.shape}")
        
        # Generate autoregressive NN model with regressors (if available)
        print(f"\nGenerating autoregressive NN model (with regressors)...")
        if features:
            print(f"  Training neural network with regressors...")
        result_with_reg = autoregressive_NN_with_regressors(
            price_df, perfect_dict, features, project,
            num_forecast_days=7,
            use_perfect_foresight_first_day=True,
            n_lags=24,
            epochs=50
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
        print(f"\nError: {str(e)}")
        import traceback
        traceback.print_exc()
        raise
