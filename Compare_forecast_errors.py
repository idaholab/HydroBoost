# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

"""
Compare forecast errors across different price forecasting models.

This script compares the performance of multiple forecasting models:
- Mean Persistence
- Additive Model (NeuralProphet-based)
- Autoregressive Neural Network
- Random Forest

The comparison uses Perfect Foresight as ground truth and calculates:
- Mean Absolute Error (MAE)
- Root Mean Squared Error (RMSE)

Results are visualized as line plots showing error vs forecast hour (0 to 168 hours).
"""

import os
import sys
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def load_forecast_data(project: str, forecast_type: str, market_product: str = 'DA_LMP') -> Optional[pd.DataFrame]:
    """
    Load forecast data for a specific model type.
    
    Parameters
    ----------
    project : str
        Project name (e.g., 'Model_A', 'Model_B', 'Model_C')
    forecast_type : str
        Type of forecast ('Perfect_foresight', 'Mean_persistence', 'Additive_Models', 
        'Autoregressive_NN', 'Random_Forest')
    market_product : str, optional
        Market product to analyze (default: 'DA_LMP')
    
    Returns
    -------
    Optional[pd.DataFrame]
        Forecast DataFrame with hours as index and dates as columns, or None if file not found
    """
    # Build file path based on forecast type
    base_path = f'Core_Models/HydroBoost/generated_data/{project}/Market/{forecast_type}'
    
    if forecast_type == 'Perfect_foresight' or forecast_type == 'Mean_persistence':
        file_path = os.path.join(base_path, f'{market_product}.csv')
    elif forecast_type == 'Additive_Models':
        file_path = os.path.join(base_path, 'Additive_model_no_regressors.csv')
    elif forecast_type == 'Autoregressive_NN':
        file_path = os.path.join(base_path, 'Autoregressive_NN_no_regressors.csv')
    elif forecast_type == 'Random_Forest':
        file_path = os.path.join(base_path, 'Random_Forest_no_regressors.csv')
    else:
        print(f"Warning: Unknown forecast type '{forecast_type}'")
        return None
    
    if not os.path.exists(file_path):
        print(f"Warning: File not found: {file_path}")
        return None
    
    try:
        df = pd.read_csv(file_path, index_col=0)
        return df
    except Exception as e:
        print(f"Error loading {file_path}: {str(e)}")
        return None


def calculate_error_by_hour(
    forecast_df: pd.DataFrame,
    perfect_foresight_df: pd.DataFrame,
    metric: str = 'MAE'
) -> np.ndarray:
    """
    Calculate forecast error for each hour across all forecast dates.
    
    Parameters
    ----------
    forecast_df : pd.DataFrame
        Forecast data with hours as index (0-167) and dates as columns
    perfect_foresight_df : pd.DataFrame
        Perfect foresight (ground truth) with same structure
    metric : str, optional
        Error metric to calculate: 'MAE' (Mean Absolute Error) or 
        'RMSE' (Root Mean Squared Error). Default: 'MAE'
    
    Returns
    -------
    np.ndarray
        Array of errors for each forecast hour (length = number of hours)
    """
    # Ensure both dataframes have the same columns (dates)
    common_dates = forecast_df.columns.intersection(perfect_foresight_df.columns)
    
    if len(common_dates) == 0:
        raise ValueError("No common dates found between forecast and perfect foresight")
    
    forecast_aligned = forecast_df[common_dates].astype(float)
    perfect_aligned = perfect_foresight_df[common_dates].astype(float)
    
    # Calculate error for each hour (row)
    errors_by_hour = []
    
    for hour_idx in forecast_aligned.index:
        forecast_values = forecast_aligned.loc[hour_idx].values
        perfect_values = perfect_aligned.loc[hour_idx].values
        
        # Remove NaN values
        mask = ~(np.isnan(forecast_values) | np.isnan(perfect_values))
        forecast_clean = forecast_values[mask]
        perfect_clean = perfect_values[mask]
        
        if len(forecast_clean) == 0:
            errors_by_hour.append(np.nan)
            continue
        
        # Calculate error based on metric
        if metric == 'MAE':
            error = np.mean(np.abs(forecast_clean - perfect_clean))
        elif metric == 'RMSE':
            error = np.sqrt(np.mean((forecast_clean - perfect_clean) ** 2))
        else:
            raise ValueError(f"Unknown metric: {metric}. Use 'MAE' or 'RMSE'")
        
        errors_by_hour.append(error)
    
    return np.array(errors_by_hour)


def plot_forecast_errors(
    errors_dict: Dict[str, np.ndarray],
    metric: str,
    project: str,
    market_product: str = 'DA_LMP',
    save_path: Optional[str] = None
) -> None:
    """
    Create line plot comparing forecast errors across models using Plotly.
    
    Parameters
    ----------
    errors_dict : Dict[str, np.ndarray]
        Dictionary mapping model names to error arrays
    metric : str
        Error metric name ('MAE' or 'RMSE')
    project : str
        Project name for plot title
    market_product : str, optional
        Market product name for plot title
    save_path : Optional[str], optional
        Path to save the plot. If None, plot is displayed but not saved
    """
    # Create figure
    fig = go.Figure()
    
    # Define colors and line styles for each model (matching paper style)
    colors = {
        'Mean Persistence': '#636EFA',  # Blue
        'Additive Model': '#EF553B',    # Red
        'Autoregressive NN': '#00CC96', # Green
        'Random Forest': '#AB63FA'      # Purple
    }
    
    dash_styles = {
        'Mean Persistence': 'solid',
        'Additive Model': 'dash',
        'Autoregressive NN': 'dot',
        'Random Forest': 'dashdot'
    }
    
    # Plot each model
    for model_name, errors in errors_dict.items():
        hours = np.arange(len(errors))
        
        fig.add_trace(go.Scatter(
            x=hours,
            y=errors,
            mode='lines',
            name=model_name,
            line=dict(
                color=colors.get(model_name, '#000000'),
                width=2,
                dash=dash_styles.get(model_name, 'solid')
            ),
            hovertemplate='<b>%{fullData.name}</b><br>' +
                         'Hour: %{x}<br>' +
                         f'{metric}: $%{{y:.2f}}/MWh<br>' +
                         '<extra></extra>'
        ))
    
    # Add vertical lines to mark days
    for day in range(1, 8):
        fig.add_vline(
            x=day*24,
            line_dash="dash",
            line_color="lightgray",
            line_width=1,
            opacity=0.5
        )
    
    # Update layout to match paper style
    fig.update_layout(
        title={
            'text': f'Forecast Error Comparison - {project}<br><sub>{metric} by Forecast Hour</sub>',
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 16}
        },
        xaxis=dict(
            title='Forecast Hour',
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(128, 128, 128, 0.2)',
            zeroline=False,
            dtick=24,  # Tick every day
            range=[0, 168]
        ),
        yaxis=dict(
            title=f'{metric} ($/MWh)',
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(128, 128, 128, 0.2)',
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor='rgba(128, 128, 128, 0.3)'
        ),
        hovermode='x unified',
        legend=dict(
            orientation='v',
            yanchor='top',
            y=0.99,
            xanchor='right',
            x=0.99,
            bgcolor='rgba(255, 255, 255, 0.8)',
            bordercolor='rgba(128, 128, 128, 0.5)',
            borderwidth=1
        ),
        plot_bgcolor='white',
        width=1200,
        height=600,
        margin=dict(l=80, r=80, t=100, b=80)
    )
    
    # Save or show the plot
    if save_path:
        # Save as HTML for interactive viewing
        html_path = save_path.replace('.png', '.html')
        fig.write_html(html_path)
        print(f"✓ Interactive plot saved to: {html_path}")
        
        # Also save as static PNG
        try:
            fig.write_image(save_path, width=1200, height=600, scale=2)
            print(f"✓ Static plot saved to: {save_path}")
        except Exception as e:
            print(f"Note: Could not save PNG (kaleido may not be installed). HTML version available.")
            print(f"  To enable PNG export: pip install kaleido")
    else:
        fig.show()


def print_error_summary(errors_dict: Dict[str, np.ndarray], metric: str) -> None:
    """
    Print summary statistics for forecast errors.
    
    Parameters
    ----------
    errors_dict : Dict[str, np.ndarray]
        Dictionary mapping model names to error arrays
    metric : str
        Error metric name
    """
    print(f"\n{'='*80}")
    print(f"FORECAST ERROR SUMMARY - {metric}")
    print(f"{'='*80}")
    print(f"{'Model':<25} {'Overall':>12} {'Day 1':>12} {'Day 2-7':>12} {'Last Day':>12}")
    print(f"{'-'*80}")
    
    for model_name, errors in errors_dict.items():
        # Calculate statistics for different time periods
        overall_error = np.nanmean(errors)
        day1_error = np.nanmean(errors[0:24])  # Hours 0-23
        day2_7_error = np.nanmean(errors[24:])  # Hours 24+
        last_day_error = np.nanmean(errors[-24:])  # Last 24 hours
        
        print(f"{model_name:<25} {overall_error:>11.2f}$ {day1_error:>11.2f}$ "
              f"{day2_7_error:>11.2f}$ {last_day_error:>11.2f}$")
    
    print(f"{'='*80}\n")


def compare_forecast_models(
    project: str,
    market_product: str = 'DA_LMP',
    metric: str = 'MAE',
    save_plots: bool = True
) -> Dict[str, np.ndarray]:
    """
    Main function to compare all forecast models.
    
    Parameters
    ----------
    project : str
        Project name (e.g., 'Model_A', 'Model_B', 'Model_C')
    market_product : str, optional
        Market product to analyze (default: 'DA_LMP')
    metric : str, optional
        Error metric to use: 'MAE' or 'RMSE' (default: 'MAE')
    save_plots : bool, optional
        Whether to save plots to file (default: True)
    
    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary of errors for each model
    """
    print(f"\n{'='*80}")
    print(f"FORECAST MODEL COMPARISON - {project} - {market_product}")
    print(f"{'='*80}\n")
    
    # Load perfect foresight (ground truth)
    print("Loading forecast data...")
    perfect_foresight = load_forecast_data(project, 'Perfect_foresight', market_product)
    
    if perfect_foresight is None:
        raise FileNotFoundError(f"Perfect foresight data not found for {project}")
    
    print(f"✓ Perfect foresight loaded: {perfect_foresight.shape}")
    
    # Load all forecast models
    models = {
        'Mean Persistence': ('Mean_persistence', market_product),
        'Additive Model': ('Additive_Models', None),
        'Autoregressive NN': ('Autoregressive_NN', None),
        'Random Forest': ('Random_Forest', None)
    }
    
    forecast_data = {}
    for model_name, (forecast_type, product) in models.items():
        df = load_forecast_data(project, forecast_type, product or market_product)
        if df is not None:
            forecast_data[model_name] = df
            print(f"✓ {model_name} loaded: {df.shape}")
        else:
            print(f"✗ {model_name} not found - skipping")
    
    if len(forecast_data) == 0:
        raise ValueError("No forecast data found to compare")
    
    # Calculate errors for each model
    print(f"\nCalculating {metric} for each model...")
    errors_dict = {}
    
    for model_name, forecast_df in forecast_data.items():
        try:
            errors = calculate_error_by_hour(forecast_df, perfect_foresight, metric)
            errors_dict[model_name] = errors
            print(f"✓ {model_name}: Average {metric} = ${np.nanmean(errors):.2f}/MWh")
        except Exception as e:
            print(f"✗ Error calculating {metric} for {model_name}: {str(e)}")
    
    # Print summary statistics
    print_error_summary(errors_dict, metric)
    
    # Create plot
    print("Generating comparison plot...")
    if save_plots:
        output_dir = f'Core_Models/HydroBoost/generated_data/{project}/Figures'
        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f'Forecast_error_comparison_{metric}.png')
    else:
        save_path = None
    
    plot_forecast_errors(errors_dict, metric, project, market_product, save_path)
    
    return errors_dict


def main():
    """
    Command-line interface for forecast model comparison.
    
    Usage:
        python Compare_forecast_errors.py Model_A
        python Compare_forecast_errors.py Model_B MAE
        python Compare_forecast_errors.py Model_C RMSE
    """
    # Parse command line arguments
    if len(sys.argv) < 2:
        print("Usage: python Compare_forecast_errors.py <project> [metric]")
        print("  project: Model_A, Model_B, or Model_C")
        print("  metric: MAE (default) or RMSE")
        sys.exit(1)
    
    project = sys.argv[1]
    metric = sys.argv[2].upper() if len(sys.argv) > 2 else 'MAE'
    
    if metric not in ['MAE', 'RMSE']:
        print(f"Error: Invalid metric '{metric}'. Use 'MAE' or 'RMSE'")
        sys.exit(1)
    
    try:
        # Run comparison for both MAE and RMSE
        for error_metric in ['MAE', 'RMSE']:
            errors_dict = compare_forecast_models(
                project=project,
                market_product='DA_LMP',
                metric=error_metric,
                save_plots=True
            )
        
        print(f"\n{'='*80}")
        print("COMPARISON COMPLETE")
        print(f"{'='*80}")
        print(f"Plots saved to: Core_Models/HydroBoost/generated_data/{project}/Figures/")
        
    except Exception as e:
        print(f"\nError: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
