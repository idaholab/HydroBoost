# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

"""
Random Forest forecasting models for DA-LMP price forecasting.

This module provides two Random Forest-based forecasting approaches:
1. No regressors: Random Forest using only engineered time series features
2. With regressors: Random Forest with additional external features (e.g., temperature, demand)

Both models use direct multi-output prediction (entire horizon in one call) for maximum speed.

Both models generate forecasts for a configurable number of days (default: 7 days, 168 hours).
When use_perfect_foresight_first_day=True (default):
- First day (hours 0-23): Perfect foresight (actual values)
- Remaining days: Random Forest predictions
"""

import os
import warnings
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

from helper_functions import generate_path, parse_CLI_arguments
from import_data import read_price_and_forecasting
from perfect_foresight_model import create_perfect_foresight_forecast

# Module-level constants
HOURS_PER_DAY = 24

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=FutureWarning)


# ========================================
# HELPER FUNCTIONS FOR VECTORIZED FEATURES
# ========================================

def _time_features(index: pd.DatetimeIndex) -> pd.DataFrame:
    """Create time-based features from datetime index."""
    return pd.DataFrame({
        "hour": index.hour.astype(np.int16),
        "dow": index.dayofweek.astype(np.int16),
        "month": index.month.astype(np.int16),
        "doy": index.dayofyear.astype(np.int16),
    }, index=index)


def _lag_matrix(series: pd.Series, n_lags: int) -> pd.DataFrame:
    """Create lagged features (lag_1 to lag_n)."""
    df = pd.concat(
        {f"lag_{i}": series.shift(i) for i in range(1, n_lags + 1)},
        axis=1
    )
    return df


def _rolling_stats(series: pd.Series, window: int = 24) -> pd.DataFrame:
    """Create rolling statistics features."""
    return pd.DataFrame({
        "roll_mean": series.rolling(window, min_periods=window).mean(),
        "roll_std": series.rolling(window, min_periods=window).std()
    }, index=series.index)


def _build_direct_targets(series: pd.Series, start: int, horizon: int) -> pd.DataFrame:
    """
    Build target matrix for direct multi-output forecasting.
    
    Parameters
    ----------
    series : pd.Series
        Target time series
    start : int
        Starting hour offset (e.g., 24 for perfect foresight first day)
    horizon : int
        Number of hours to forecast
        
    Returns
    -------
    pd.DataFrame
        Target matrix with columns for each forecast hour
    """
    cols = {}
    for h in range(start, start + horizon):
        cols[f"y+h{h}"] = series.shift(-h)
    Y = pd.concat(cols, axis=1)
    return Y


def _origin_rows_mask(df: pd.DataFrame) -> pd.Series:
    """Return mask for rows without NaN values."""
    return ~df.isna().any(axis=1)


def _astype_float32(df: pd.DataFrame) -> pd.DataFrame:
    """Convert DataFrame to float32 for memory efficiency."""
    return df.astype(np.float32)


# ========================================
# MAIN FORECASTING FUNCTIONS
# ========================================

def random_forest_no_regressors(
    price_df: pd.DataFrame,
    perfect_dict: Dict[str, pd.DataFrame],
    project: str,
    num_forecast_days: int = 7,
    use_perfect_foresight_first_day: bool = True,
    n_lags: int = 24,
    n_estimators: int = 200,
    max_depth: Optional[int] = 24
) -> pd.DataFrame:
    """
    Create DA-LMP forecast using direct multi-output Random Forest.
    
    This optimized version trains a single model that predicts the entire forecast
    horizon in one call, eliminating the need for recursive per-hour prediction.
    
    Features used:
    - Lagged values (previous n_lags hours)
    - Time features (hour, day of week, month, day of year)
    - Rolling statistics (24-hour mean and std)
    
    Parameters
    ----------
    price_df : pd.DataFrame
        Historical price data with DatetimeIndex and columns including 'DA-LMP'.
    perfect_dict : Dict[str, pd.DataFrame]
        Dictionary of perfect foresight forecasts with keys including 'DA-LMP'.
        Used as a template for output structure.
    project : str
        Project name (e.g., 'Model_A') for determining output file path.
    num_forecast_days : int, optional
        Number of days in forecast window (default: 7, resulting in 168 hours).
    use_perfect_foresight_first_day : bool, optional
        If True (default), hours 0-23 use perfect foresight values.
    n_lags : int, optional
        Number of lagged observations to use as features (default: 24).
    n_estimators : int, optional
        Number of trees in the Random Forest (default: 200).
    max_depth : Optional[int], optional
        Maximum depth of trees (default: 24 for speed/accuracy balance).
    
    Returns
    -------
    pd.DataFrame
        Forecast DataFrame with hours as rows and dates as columns.
    
    Raises
    ------
    KeyError
        If 'DA-LMP' column is missing from price_df or perfect_dict.
    ValueError
        If price_df has insufficient data.
    RuntimeError
        If Random Forest model fails to fit or predict.
    OSError
        If output file cannot be written.
    """
    # Input validation
    if 'DA-LMP' not in price_df.columns:
        raise KeyError("price_df must contain 'DA-LMP' column")
    if 'DA-LMP' not in perfect_dict:
        raise KeyError("perfect_dict must contain 'DA-LMP' key")
    
    target = price_df['DA-LMP']
    start = HOURS_PER_DAY if use_perfect_foresight_first_day else 0
    horizon = num_forecast_days * HOURS_PER_DAY - start
    
    if len(price_df) < n_lags + HOURS_PER_DAY + horizon:
        raise ValueError(
            f"Insufficient data: need at least {n_lags + HOURS_PER_DAY + horizon} hours, "
            f"but only have {len(price_df)}"
        )
    
    try:
        # ---- BUILD TRAINING MATRIX (VECTORIZED) ----
        print(f"  Creating features (n_lags={n_lags}, horizon={horizon})...")
        
        # Time features
        tf = _time_features(price_df.index)
        
        # Lag features
        lags = _lag_matrix(target, n_lags=n_lags)
        
        # Rolling statistics
        rolls = _rolling_stats(target, window=HOURS_PER_DAY)
        
        # Combine all features
        X_all = pd.concat([tf, lags, rolls], axis=1)
        
        # Build target matrix for direct multi-output
        Y_all = _build_direct_targets(target, start=start, horizon=horizon)
        
        # Combine and filter valid rows
        XY = pd.concat([X_all, Y_all], axis=1)
        mask = _origin_rows_mask(XY)
        XY = XY.loc[mask]
        
        X = _astype_float32(XY[X_all.columns])
        Y = _astype_float32(XY[Y_all.columns])
        
        print(f"  Training data shape: X={X.shape}, Y={Y.shape}")
        
        # ---- TRAIN RANDOM FOREST (MULTI-OUTPUT) ----
        print(f"  Training Random Forest (n_estimators={n_estimators}, max_depth={max_depth})...")
        rf = RandomForestRegressor(
            n_estimators=n_estimators,
            max_depth=max_depth,
            n_jobs=-1,
            random_state=42,
            min_samples_split=2,
            min_samples_leaf=1,
            verbose=0
        )
        rf.fit(X.values, Y.values)
        print(f"  ✓ Training complete")
        
        # ---- INFERENCE (ONE PREDICT PER FORECAST DATE) ----
        print(f"  Generating forecasts for each date...")
        template = perfect_dict['DA-LMP'].T.copy()  # rows=dates, cols=hours
        result = template.copy()
        
        num_forecasts = len(template)
        for i, forecast_date in enumerate(pd.to_datetime(template.index)):
            if (i + 1) % 50 == 0:
                print(f"    Progress: {i+1}/{num_forecasts} forecasts...")
            
            if forecast_date not in X_all.index:
                # Forecast date not in historical data - skip or use fallback
                continue
            
            # Get feature row at origin time
            x_row = X_all.loc[[forecast_date]].copy()
            
            if x_row.isna().any(axis=1).iloc[0]:
                # Insufficient history - use simple persistence fallback
                if forecast_date in price_df.index:
                    last_val = price_df.loc[:forecast_date, 'DA-LMP'].iloc[-1]
                    if use_perfect_foresight_first_day:
                        result.iloc[i, HOURS_PER_DAY:] = last_val
                    else:
                        result.iloc[i, :] = last_val
                continue
            
            # Predict entire horizon in one call
            yhat_vec = rf.predict(_astype_float32(x_row).values)[0]
            
            # Write predictions to template
            if use_perfect_foresight_first_day:
                result.iloc[i, HOURS_PER_DAY:] = yhat_vec
            else:
                result.iloc[i, :] = yhat_vec
        
        print(f"  ✓ Generated {num_forecasts} forecasts")
        
        # Transpose back to standard format (hours as rows, dates as columns)
        result = result.T
        
    except Exception as e:
        raise RuntimeError(f"Random Forest model failed: {str(e)}") from e
    
    # Save output
    try:
        out_dir = generate_path([project, 'Market', 'Random_Forest'])
        output_path = os.path.join(out_dir, 'Random_Forest_no_regressors.csv')
        result.to_csv(output_path, index=True)
        print(f"✓ Saved Random Forest model (no regressors) to: {output_path}")
    except OSError as e:
        raise OSError(f"Failed to write output file: {str(e)}") from e
    
    return result


def random_forest_with_regressors(
    price_df: pd.DataFrame,
    perfect_dict: Dict[str, pd.DataFrame],
    features: Optional[List[str]],
    project: str,
    num_forecast_days: int = 7,
    use_perfect_foresight_first_day: bool = True,
    n_lags: int = 24,
    n_estimators: int = 200,
    max_depth: Optional[int] = 24
) -> Optional[pd.DataFrame]:
    """
    Create DA-LMP forecast using direct multi-output Random Forest with external regressors.
    
    This optimized version builds horizon-specific regressor features by cycling the
    last 24 hours of regressor values, allowing the model to learn hour-specific effects
    while maintaining fast single-call prediction per forecast date.
    
    Features used:
    - Lagged price values (previous n_lags hours)
    - Time features (hour, day of week, month, day of year)
    - Rolling statistics (24-hour mean and std)
    - Horizon-specific regressor values (cycled from last 24 hours)
    
    Parameters
    ----------
    price_df : pd.DataFrame
        Historical price data with DatetimeIndex and columns including 'DA-LMP'
        and all features specified in the features list.
    perfect_dict : Dict[str, pd.DataFrame]
        Dictionary of perfect foresight forecasts.
    features : Optional[List[str]]
        List of column names in price_df to use as external regressors.
        If None or empty, the function returns None.
    project : str
        Project name (e.g., 'Model_A') for determining output file path.
    num_forecast_days : int, optional
        Number of days in forecast window (default: 7).
    use_perfect_foresight_first_day : bool, optional
        If True (default), hours 0-23 use perfect foresight values.
    n_lags : int, optional
        Number of lagged observations for features (default: 24).
    n_estimators : int, optional
        Number of trees in the Random Forest (default: 200).
    max_depth : Optional[int], optional
        Maximum depth of trees (default: 24).
    
    Returns
    -------
    Optional[pd.DataFrame]
        Forecast DataFrame, or None if no features provided.
    
    Raises
    ------
    KeyError
        If required columns are missing.
    ValueError
        If price_df has insufficient data.
    RuntimeError
        If Random Forest model fails.
    OSError
        If output file cannot be written.
    """
    # Handle case with no regressors
    if not features:
        print("ℹ No regressors provided; skipping Random Forest with regressors.")
        return None
    
    # Input validation
    if 'DA-LMP' not in price_df.columns:
        raise KeyError("price_df must contain 'DA-LMP' column")
    
    missing_features = [f for f in features if f not in price_df.columns]
    if missing_features:
        raise KeyError(f"price_df is missing required feature columns: {missing_features}")
    
    if 'DA-LMP' not in perfect_dict:
        raise KeyError("perfect_dict must contain 'DA-LMP' key")
    
    target = price_df['DA-LMP']
    start = HOURS_PER_DAY if use_perfect_foresight_first_day else 0
    horizon = num_forecast_days * HOURS_PER_DAY - start
    
    if len(price_df) < n_lags + HOURS_PER_DAY + horizon:
        raise ValueError(
            f"Insufficient data: need at least {n_lags + HOURS_PER_DAY + horizon} hours"
        )
    
    try:
        # ---- BUILD TRAINING MATRIX ----
        print(f"  Creating features with regressors (n_lags={n_lags}, horizon={horizon})...")
        
        # Base features
        tf = _time_features(price_df.index)
        lags = _lag_matrix(target, n_lags=n_lags)
        rolls = _rolling_stats(target, window=HOURS_PER_DAY)
        X_base = pd.concat([tf, lags, rolls], axis=1)
        
        # Build horizon regressor features
        # For each regressor and each horizon step, create feature using
        # last 24h pattern indexed by hour-of-day
        def make_horizon_reg_block(df: pd.DataFrame, reg_name: str) -> pd.DataFrame:
            """
            Create horizon-specific regressor features using 24h cyclic pattern.
            
            For each origin time t and horizon step h, the regressor value is taken
            from the last 24 hours based on the hour-of-day of (t+h).
            """
            reg = df[reg_name]
            idx = df.index
            origin_hour = idx.hour.values  # 0..23
            
            # Create lag columns for last 24 hours (most recent first)
            last24_cols = []
            for j in range(1, HOURS_PER_DAY + 1):
                last24_cols.append(reg.shift(j))
            last24_df = pd.concat(last24_cols, axis=1)
            last24_vals = last24_df.to_numpy(dtype=np.float32)  # shape: [t, 24]
            
            # For each horizon step, compute hour-of-day and map to last24 index
            # The last24 array is ordered [t-1, t-2, ..., t-24]
            # We want to index by hour-of-day, so we need to map correctly
            
            # Build a 24-hour template from last 24 values
            # Index 0 corresponds to hour (current_hour - 1) % 24
            # We'll index based on target hour-of-day
            
            horizon_features = []
            for h in range(start, start + horizon):
                # Hour of day for time (t + h)
                target_hour = (origin_hour + h) % HOURS_PER_DAY
                
                # Map to last24 index: we want value from last occurrence of this hour
                # If current hour is 10 and we want hour 11, that's 1 hour ago (index 0 in last24)
                # If current hour is 10 and we want hour 9, that's 23 hours ago (index 22 in last24)
                hours_back = (origin_hour - target_hour) % HOURS_PER_DAY
                hours_back[hours_back == 0] = HOURS_PER_DAY  # If same hour, use 24h ago
                last24_idx = hours_back - 1  # Convert to 0-based index
                
                # Gather values
                t_idx = np.arange(len(idx))
                feature_vals = last24_vals[t_idx, last24_idx]
                horizon_features.append(feature_vals)
            
            # Create DataFrame
            horizon_arr = np.column_stack(horizon_features)
            cols = [f'{reg_name}_h{h}' for h in range(start, start + horizon)]
            return pd.DataFrame(horizon_arr, index=idx, columns=cols)
        
        # Build regressor blocks for all features
        X_reg_blocks = [make_horizon_reg_block(price_df, r) for r in features]
        X_all = pd.concat([X_base] + X_reg_blocks, axis=1)
        
        # Build targets
        Y_all = _build_direct_targets(target, start=start, horizon=horizon)
        
        # Combine and filter valid rows
        XY = pd.concat([X_all, Y_all], axis=1)
        mask = _origin_rows_mask(XY)
        XY = XY.loc[mask]
        
        X = _astype_float32(XY[X_all.columns])
        Y = _astype_float32(XY[Y_all.columns])
        
        print(f"  Training data shape: X={X.shape}, Y={Y.shape}")
        print(f"  Regressors used: {', '.join(features)}")
        
        # ---- TRAIN RANDOM FOREST ----
        print(f"  Training Random Forest (n_estimators={n_estimators}, max_depth={max_depth})...")
        rf = RandomForestRegressor(
            n_estimators=n_estimators,
            max_depth=max_depth,
            n_jobs=-1,
            random_state=42,
            min_samples_split=2,
            min_samples_leaf=1,
            verbose=0
        )
        rf.fit(X.values, Y.values)
        print(f"  ✓ Training complete")
        
        # ---- INFERENCE ----
        print(f"  Generating forecasts for each date...")
        template = perfect_dict['DA-LMP'].T.copy()
        result = template.copy()
        
        num_forecasts = len(template)
        for i, forecast_date in enumerate(pd.to_datetime(template.index)):
            if (i + 1) % 50 == 0:
                print(f"    Progress: {i+1}/{num_forecasts} forecasts...")
            
            if forecast_date not in X_all.index:
                continue
            
            x_row = X_all.loc[[forecast_date]].copy()
            
            if x_row.isna().any(axis=1).iloc[0]:
                # Fallback to persistence
                if forecast_date in price_df.index:
                    last_val = price_df.loc[:forecast_date, 'DA-LMP'].iloc[-1]
                    if use_perfect_foresight_first_day:
                        result.iloc[i, HOURS_PER_DAY:] = last_val
                    else:
                        result.iloc[i, :] = last_val
                continue
            
            # Predict entire horizon
            yhat_vec = rf.predict(_astype_float32(x_row).values)[0]
            
            if use_perfect_foresight_first_day:
                result.iloc[i, HOURS_PER_DAY:] = yhat_vec
            else:
                result.iloc[i, :] = yhat_vec
        
        print(f"  ✓ Generated {num_forecasts} forecasts")
        
        # Transpose back
        result = result.T
        
    except Exception as e:
        raise RuntimeError(
            f"Random Forest model with regressors failed: {str(e)}"
        ) from e
    
    # Save output
    try:
        out_dir = generate_path([project, 'Market', 'Random_Forest'])
        output_path = os.path.join(out_dir, 'Random_Forest_with_regressors.csv')
        result.to_csv(output_path, index=True)
        print(f"✓ Saved Random Forest model (with regressors) to: {output_path}")
        print(f"  Regressors used: {', '.join(features)}")
    except OSError as e:
        raise OSError(f"Failed to write output file: {str(e)}") from e
    
    return result


if __name__ == '__main__':
    """
    Command-line interface for testing Random Forest models.
    
    Usage:
        python random_forest_model.py Model_A
        python random_forest_model.py Model_B
        python random_forest_model.py Model_C
    
    This will:
    1. Load price data and features from Input_Spreadsheets/{project}.xlsx
    2. Generate perfect foresight forecast
    3. Generate Random Forest model without regressors
    4. Generate Random Forest model with regressors (if features are available)
    """
    # Parse command line arguments for project name
    project = parse_CLI_arguments()
    
    print("=" * 70)
    print(f"RANDOM FOREST FORECAST GENERATION - {project}")
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
        
        # Generate Random Forest model without regressors
        print(f"\nGenerating Random Forest model (no regressors)...")
        result_no_reg = random_forest_no_regressors(
            price_df, perfect_dict, project, 
            num_forecast_days=7, 
            use_perfect_foresight_first_day=True,
            n_lags=24,
            n_estimators=200,
            max_depth=24
        )
        print(f"✓ Generated forecast shape: {result_no_reg.shape}")
        
        # Generate Random Forest model with regressors (if available)
        print(f"\nGenerating Random Forest model (with regressors)...")
        result_with_reg = random_forest_with_regressors(
            price_df, perfect_dict, features, project,
            num_forecast_days=7,
            use_perfect_foresight_first_day=True,
            n_lags=24,
            n_estimators=200,
            max_depth=24
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
