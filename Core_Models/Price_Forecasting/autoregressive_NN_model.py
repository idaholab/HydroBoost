# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

"""
NeuralProphet-based autoregressive forecasts for DA-LMP.
Includes a no-regressor model and a regressor model.
"""

import os
import warnings
from datetime import timedelta
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from neuralprophet import NeuralProphet

from .helper_functions import generate_path, parse_CLI_arguments
from .import_data import read_price_and_forecasting
from .perfect_foresight_model import create_perfect_foresight_forecast

# NumPy 2.0 compatibility: np.NaN was removed, use np.nan instead
if not hasattr(np, 'NaN'):
    np.NaN = np.nan

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
    NeuralProphet AR forecast without regressors.

    Returns a DataFrame with hours as rows and forecast start dates as columns.
    If use_perfect_foresight_first_day is True, hours 0-23 are perfect foresight.
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
        
        # Initialize NeuralProphet model (same as with_regressors, just without regressors)
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
        
        # Fit the model
        metrics = model.fit(df_neural, freq='h')
        
        # Extend dataframe with future dates
        # Repeat the last 24 hours for each forecast day
        last_24_hours = df_neural.iloc[-HOURS_PER_DAY:][['ds', 'y']].copy()
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
        raise RuntimeError(f"NeuralProphet model fitting or prediction failed: {str(e)}") from e
    
    # Use perfect foresight template structure (will be transposed for filling)
    template = perfect_dict['DA-LMP'].T.copy()
    num_forecasts = template.shape[0]
    
    # Fill each forecast date with corresponding predictions
    for i in range(num_forecasts):
        start_idx = i * HOURS_PER_DAY
        end_idx = start_idx + forecast_hours
        forecast_block = predicted_values[start_idx:end_idx]
        
        # Handle edge case: pad if predictions are insufficient
        if len(forecast_block) < forecast_hours:
            pad_length = forecast_hours - len(forecast_block)
            last_value = forecast_block[-1] if forecast_block.size > 0 else np.nan
            forecast_block = np.concatenate([
                forecast_block,
                np.full(pad_length, last_value)
            ])
        
        # Apply predictions based on mode
        if use_perfect_foresight_first_day:
            # Keep perfect foresight for hours 0-23, use NeuralProphet for hours 24+
            template.iloc[i, HOURS_PER_DAY:] = forecast_block[HOURS_PER_DAY:].astype(float)
        else:
            # Use NeuralProphet predictions for all hours
            template.iloc[i, :] = forecast_block.astype(float)
    
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
    NeuralProphet AR forecast with regressors.

    Returns a DataFrame with hours as rows and forecast start dates as columns,
    or None if no features are provided.
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
        
        # Extend dataframe with future dates
        # Repeat the last 24 hours of regressor values for each forecast day
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
    CLI entry point for quick testing.
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
