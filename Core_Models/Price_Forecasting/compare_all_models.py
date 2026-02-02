# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

"""
Comprehensive comparison script for all DA-LMP forecasting models.

This script loads forecasts from all models and compares them against perfect foresight
to evaluate forecast accuracy using Mean Absolute Error (MAE) and other statistics.

Models compared:
1. Mean Persistence
2. Additive Model (no regressors)
3. Additive Model (with regressors)
4. Autoregressive NN (no regressors)
5. Autoregressive NN (with regressors)

Usage:
    python compare_all_models.py Model_A
    python compare_all_models.py Model_B
    python compare_all_models.py Model_C
"""

import os
import pandas as pd
import numpy as np
from typing import Dict, Tuple
from .helper_functions import parse_CLI_arguments

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir, os.pardir))


def repo_path(*parts: str) -> str:
    return os.path.join(REPO_ROOT, *parts)

# Module-level constants
HOURS_PER_DAY = 24


def load_forecast(project: str, model_path: str, model_name: str) -> pd.DataFrame:
    """
    Load a forecast CSV file.
    
    Parameters
    ----------
    project : str
        Project name (e.g., 'Model_A')
    model_path : str
        Relative path to the forecast file within generated_data/{project}/Market/
    model_name : str
        Display name for the model (used in error messages)
    
    Returns
    -------
    pd.DataFrame
        Forecast dataframe with hours as index and dates as columns
        
    Raises
    ------
    FileNotFoundError
        If the forecast file doesn't exist
    """
    base_dir = repo_path("Core_Models", "HydroBoost", "generated_data", project, "Market")
    full_path = os.path.join(base_dir, model_path)
    
    if not os.path.exists(full_path):
        raise FileNotFoundError(
            f"{model_name} forecast not found at: {full_path}\n"
            f"Please run the model to generate forecasts first."
        )
    
    df = pd.read_csv(full_path, index_col=0)
    return df


def compute_statistics(
    forecast: pd.DataFrame,
    perfect_foresight: pd.DataFrame,
    start_hour: int = 24,
    end_hour: int = 168
) -> Dict[str, float]:
    """
    Compute forecast statistics vs perfect foresight.
    
    Parameters
    ----------
    forecast : pd.DataFrame
        Forecast values
    perfect_foresight : pd.DataFrame
        Perfect foresight (ground truth) values
    start_hour : int, optional
        Starting hour for comparison (default: 24, after first day)
    end_hour : int, optional
        Ending hour for comparison (default: 168, exclusive end for hours 24-167)
    
    Returns
    -------
    Dict[str, float]
        Dictionary containing:
        - 'mean': Mean of forecast values
        - 'std': Standard deviation of forecast values
        - 'min': Minimum forecast value
        - 'max': Maximum forecast value
        - 'mae': Mean Absolute Error vs perfect foresight
        - 'rmse': Root Mean Squared Error vs perfect foresight
        - 'mape': Mean Absolute Percentage Error vs perfect foresight
    """
    # Extract comparison window (hours 24-167 across all dates)
    # Use iloc for integer-based indexing (row 24 to row 167 inclusive)
    forecast_values = forecast.iloc[start_hour:end_hour].values.flatten()
    perfect_values = perfect_foresight.iloc[start_hour:end_hour].values.flatten()
    
    # Remove any NaN values
    valid_mask = ~(np.isnan(forecast_values) | np.isnan(perfect_values))
    forecast_values = forecast_values[valid_mask]
    perfect_values = perfect_values[valid_mask]
    
    # Compute errors
    errors = forecast_values - perfect_values
    abs_errors = np.abs(errors)
    squared_errors = errors ** 2
    pct_errors = abs_errors / np.abs(perfect_values) * 100
    
    return {
        'mean': float(np.mean(forecast_values)),
        'std': float(np.std(forecast_values)),
        'min': float(np.min(forecast_values)),
        'max': float(np.max(forecast_values)),
        'mae': float(np.mean(abs_errors)),
        'rmse': float(np.sqrt(np.mean(squared_errors))),
        'mape': float(np.mean(pct_errors))
    }


def print_comparison_table(results: Dict[str, Dict[str, float]]) -> None:
    """
    Print a formatted comparison table of all model statistics.
    
    Parameters
    ----------
    results : Dict[str, Dict[str, float]]
        Dictionary mapping model names to their statistics
    """
    print("\n" + "=" * 100)
    print("FORECAST MODEL COMPARISON - Hours 24-167 (Days 2-7)")
    print("=" * 100)
    print(f"\n{'Model':<35} {'Mean':>10} {'Std':>10} {'MAE':>10} {'RMSE':>10} {'MAPE':>10}")
    print("-" * 100)
    
    # Print perfect foresight first
    if 'Perfect Foresight' in results:
        stats = results['Perfect Foresight']
        print(f"{'Perfect Foresight (Ground Truth)':<35} "
              f"${stats['mean']:>9.2f} ${stats['std']:>9.2f} "
              f"{'N/A':>10} {'N/A':>10} {'N/A':>10}")
        print("-" * 100)
    
    # Print other models sorted by MAE
    other_models = {k: v for k, v in results.items() if k != 'Perfect Foresight'}
    sorted_models = sorted(other_models.items(), key=lambda x: x[1]['mae'])
    
    for model_name, stats in sorted_models:
        print(f"{model_name:<35} "
              f"${stats['mean']:>9.2f} ${stats['std']:>9.2f} "
              f"${stats['mae']:>9.2f} ${stats['rmse']:>9.2f} "
              f"{stats['mape']:>9.2f}%")
    
    print("=" * 100)
    
    # Print summary
    best_model = sorted_models[0][0]
    best_mae = sorted_models[0][1]['mae']
    worst_model = sorted_models[-1][0]
    worst_mae = sorted_models[-1][1]['mae']
    
    print(f"\n📊 SUMMARY:")
    print(f"  Best Model (Lowest MAE):  {best_model} (${best_mae:.2f})")
    print(f"  Worst Model (Highest MAE): {worst_model} (${worst_mae:.2f})")
    print(f"  MAE Range: ${best_mae:.2f} - ${worst_mae:.2f}")
    print(f"  Improvement: {((worst_mae - best_mae) / worst_mae * 100):.1f}% better")
    print()


def print_detailed_comparison(
    model1_name: str,
    model1: pd.DataFrame,
    model2_name: str,
    model2: pd.DataFrame,
    perfect_foresight: pd.DataFrame
) -> None:
    """
    Print detailed comparison between two specific models.
    
    Parameters
    ----------
    model1_name : str
        Name of first model
    model1 : pd.DataFrame
        First model forecast
    model2_name : str
        Name of second model
    model2 : pd.DataFrame
        Second model forecast
    perfect_foresight : pd.DataFrame
        Perfect foresight (ground truth)
    """
    print(f"\n" + "=" * 100)
    print(f"DETAILED COMPARISON: {model1_name} vs {model2_name}")
    print("=" * 100)
    
    # First forecast date comparison
    first_date = model1.columns[0]
    
    # Extract hours 24-167 for first forecast (use iloc for integer indexing)
    m1_values = model1.iloc[24:168][first_date].astype(float)
    m2_values = model2.iloc[24:168][first_date].astype(float)
    pf_values = perfect_foresight.iloc[24:168][first_date].astype(float)
    
    # Compute errors
    m1_errors = np.abs(m1_values - pf_values)
    m2_errors = np.abs(m2_values - pf_values)
    
    print(f"\nFirst Forecast Date: {first_date}")
    print(f"  {model1_name}:")
    print(f"    Mean: ${m1_values.mean():.2f}, Std: ${m1_values.std():.2f}, MAE: ${m1_errors.mean():.2f}")
    print(f"  {model2_name}:")
    print(f"    Mean: ${m2_values.mean():.2f}, Std: ${m2_values.std():.2f}, MAE: ${m2_errors.mean():.2f}")
    print(f"  Perfect Foresight:")
    print(f"    Mean: ${pf_values.mean():.2f}, Std: ${pf_values.std():.2f}")
    
    # Direct comparison
    diff = np.abs(m1_values - m2_values)
    print(f"\n  Direct Comparison:")
    print(f"    Mean absolute difference: ${diff.mean():.2f}")
    print(f"    Max absolute difference: ${diff.max():.2f}")
    
    if m1_errors.mean() < m2_errors.mean():
        winner = model1_name
        improvement = ((m2_errors.mean() - m1_errors.mean()) / m2_errors.mean() * 100)
    else:
        winner = model2_name
        improvement = ((m1_errors.mean() - m2_errors.mean()) / m1_errors.mean() * 100)
    
    print(f"    Better model: {winner} ({improvement:.1f}% lower MAE)")
    print("=" * 100)


if __name__ == '__main__':
    """
    Compare all forecasting models for a given project.
    """
    # Parse command line arguments
    project = parse_CLI_arguments()
    
    print("=" * 100)
    print(f"COMPREHENSIVE FORECAST MODEL COMPARISON - {project}")
    print("=" * 100)
    print("\nComparing all models against perfect foresight (ground truth)")
    print("Evaluation window: Hours 24-167 (Days 2-7 of forecast)")
    
    try:
        # Load all forecasts
        print(f"\nLoading forecasts for {project}...")
        
        forecasts = {}
        
        # Perfect foresight (ground truth)
        forecasts['Perfect Foresight'] = load_forecast(
            project, 'Perfect_foresight/DA_LMP.csv', 'Perfect Foresight'
        )
        print(f"✓ Loaded Perfect Foresight")
        
        # Mean Persistence
        forecasts['Mean Persistence'] = load_forecast(
            project, 'Mean_persistence/DA_LMP.csv', 'Mean Persistence'
        )
        print(f"✓ Loaded Mean Persistence")
        
        # Additive Models
        forecasts['Additive (no regressors)'] = load_forecast(
            project, 'Additive_Models/Additive_model_no_regressors.csv', 
            'Additive (no regressors)'
        )
        print(f"✓ Loaded Additive Model (no regressors)")
        
        try:
            forecasts['Additive (with regressors)'] = load_forecast(
                project, 'Additive_Models/Additive_model_with_regressors.csv',
                'Additive (with regressors)'
            )
            print(f"✓ Loaded Additive Model (with regressors)")
        except FileNotFoundError:
            print(f"⚠ Additive Model (with regressors) not found - skipping")
        
        # Autoregressive NN Models
        try:
            forecasts['Autoregressive NN (no regressors)'] = load_forecast(
                project, 'Autoregressive_NN/Autoregressive_NN_no_regressors.csv',
                'Autoregressive NN (no regressors)'
            )
            print(f"✓ Loaded Autoregressive NN (no regressors)")
        except FileNotFoundError:
            print(f"⚠ Autoregressive NN (no regressors) not found - skipping")
        
        try:
            forecasts['Autoregressive NN (with regressors)'] = load_forecast(
                project, 'Autoregressive_NN/Autoregressive_NN_with_regressors.csv',
                'Autoregressive NN (with regressors)'
            )
            print(f"✓ Loaded Autoregressive NN (with regressors)")
        except FileNotFoundError:
            print(f"⚠ Autoregressive NN (with regressors) not found - skipping")
        
        # Random Forest Models
        try:
            forecasts['Random Forest (no regressors)'] = load_forecast(
                project, 'Random_Forest/Random_Forest_no_regressors.csv',
                'Random Forest (no regressors)'
            )
            print(f"✓ Loaded Random Forest (no regressors)")
        except FileNotFoundError:
            print(f"⚠ Random Forest (no regressors) not found - skipping")
        
        try:
            forecasts['Random Forest (with regressors)'] = load_forecast(
                project, 'Random_Forest/Random_Forest_with_regressors.csv',
                'Random Forest (with regressors)'
            )
            print(f"✓ Loaded Random Forest (with regressors)")
        except FileNotFoundError:
            print(f"⚠ Random Forest (with regressors) not found - skipping")
        
        # Compute statistics for all models
        print(f"\nComputing statistics...")
        perfect_foresight = forecasts['Perfect Foresight']
        results = {}
        
        for model_name, forecast_df in forecasts.items():
            if model_name == 'Perfect Foresight':
                # For perfect foresight, just compute descriptive stats
                # Use iloc for integer-based indexing
                try:
                    values = forecast_df.iloc[24:168].values.flatten()
                    # Remove NaN values
                    values = values[~np.isnan(values)]
                    if len(values) == 0:
                        print(f"⚠ Warning: No valid data for {model_name}, skipping statistics")
                        continue
                    results[model_name] = {
                        'mean': float(np.mean(values)),
                        'std': float(np.std(values)),
                        'min': float(np.min(values)),
                        'max': float(np.max(values)),
                        'mae': 0.0,
                        'rmse': 0.0,
                        'mape': 0.0
                    }
                except Exception as e:
                    print(f"⚠ Error computing stats for {model_name}: {e}")
                    continue
            else:
                try:
                    results[model_name] = compute_statistics(forecast_df, perfect_foresight)
                except Exception as e:
                    print(f"⚠ Error computing stats for {model_name}: {e}")
                    continue
        
        # Print comparison table
        print_comparison_table(results)
        
        # Print detailed comparisons if multiple models available
        model_names = [k for k in forecasts.keys() if k != 'Perfect Foresight']
        
        if len(model_names) >= 2:
            # Compare no-regressor models
            if 'Additive (no regressors)' in forecasts and 'Autoregressive NN (no regressors)' in forecasts:
                print_detailed_comparison(
                    'Additive (no regressors)',
                    forecasts['Additive (no regressors)'],
                    'Autoregressive NN (no regressors)',
                    forecasts['Autoregressive NN (no regressors)'],
                    perfect_foresight
                )
            
            # Compare with-regressor models
            if 'Additive (with regressors)' in forecasts and 'Autoregressive NN (with regressors)' in forecasts:
                print_detailed_comparison(
                    'Additive (with regressors)',
                    forecasts['Additive (with regressors)'],
                    'Autoregressive NN (with regressors)',
                    forecasts['Autoregressive NN (with regressors)'],
                    perfect_foresight
                )
        
        print("\n" + "=" * 100)
        print("COMPARISON COMPLETE")
        print("=" * 100)
        print()
        
    except FileNotFoundError as e:
        print(f"\nError: {e}")
        print("\nPlease ensure all models have been run to generate forecasts.")
        print("Run the following scripts in order:")
        print("  1. python perfect_foresight_model.py {project}")
        print("  2. python mean_persistence_model.py {project}")
        print("  3. python additive_models.py {project}")
        print("  4. python autoregressive_NN_model.py {project}")
    except Exception as e:
        print(f"\nUnexpected error: {e}")
        import traceback
        traceback.print_exc()
        raise
