# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED


import os
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def plot_forecast_model_errors(file_name,
                               perfect_foresight_forecast_dict,
                               mean_persistence_forecast_dict,
                               additive_no_regressor_forecast=None,
                               additive_with_regressors_forecast=None,
                               autoregressive_NN_no_regressors_forecast=None,
                               autoregressive_NN_with_regressors_forecast=None,
                               random_forest_no_regressors_forecast=None,
                               random_forest_with_regressors_forecast=None):
    """Generate heatmap comparison of forecast errors for all models."""
    error = dict()
    max_error = 0

    
    e = mean_persistence_forecast_dict["DA-LMP"] - perfect_foresight_forecast_dict["DA-LMP"]
    error["Mean Persistence"] = e
    if e.abs().max().max():
        max_error = e.abs().max().max()

    # Check to see what forecast we have generated
    if additive_no_regressor_forecast is not None:
        e = additive_no_regressor_forecast - perfect_foresight_forecast_dict["DA-LMP"]
        error["Additive No Regressors"] = e
        if e.abs().max().max():
            max_error = e.abs().max().max()

    if additive_with_regressors_forecast is not None:
        e = additive_with_regressors_forecast - perfect_foresight_forecast_dict["DA-LMP"]
        error["Additive With Regressors"] = e
        if e.abs().max().max():
            max_error = e.abs().max().max()

    if autoregressive_NN_no_regressors_forecast is not None:
        e = autoregressive_NN_no_regressors_forecast - perfect_foresight_forecast_dict["DA-LMP"]
        error["Autoregressive NN No Regressors"] = e
        if e.abs().max().max():
            max_error = e.abs().max().max()

    if autoregressive_NN_with_regressors_forecast is not None:
        e = autoregressive_NN_with_regressors_forecast - perfect_foresight_forecast_dict["DA-LMP"]
        error["Autoregressive NN With Regressors"] = e
        if e.abs().max().max():
            max_error = e.abs().max().max()

    if random_forest_no_regressors_forecast is not None:
        e = random_forest_no_regressors_forecast - perfect_foresight_forecast_dict["DA-LMP"]
        error["Random Forest No Regressors"] = e
        if e.abs().max().max():
            max_error = e.abs().max().max()

    if random_forest_with_regressors_forecast is not None:
        e = random_forest_with_regressors_forecast - perfect_foresight_forecast_dict["DA-LMP"]
        error["Random Forest With Regressors"] = e
        if e.abs().max().max():
            max_error = e.abs().max().max()


    # Set up the subplot figure
    fig = make_subplots(rows=int((len(error)+1)/2), 
                        cols=2, 
                        subplot_titles=list(error.keys()))

    # Create and add each heatmap to the figure
    for i, (method, df) in enumerate(error.items()):

        heatmap = px.imshow(df)
        
        # Correct row, col for heatmap
        row = i // 2 + 1
        col = i % 2 + 1
        
        for trace in heatmap.data:
            fig.add_trace(trace, row=row, col=col)

        fig.update_xaxes(title_text="Date", row=row, col=col)
        fig.update_yaxes(title_text="Forecast Hours", row=row, col=col)

        fig.update_layout(coloraxis=dict(colorscale='RdBu_r',cmid=0))
        
    # Update layout for better appearance
    fig.update_layout(
        title=dict(text="Comparison of Forecasting (model-actual)", x=.5, xanchor='center', font=dict(size=25)),
        height=int((len(error)+1)/2)*400,
        width=1200
    )

    # Save and show the results. First make sure there is a directory to save in
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, os.pardir, os.pardir))
    figure_directory = os.path.join(repo_root, "Core_Models", "HydroBoost", "generated_data", file_name, "Figures")
    
    if not os.path.exists(figure_directory):
        os.makedirs(figure_directory, exist_ok=True)
        print(figure_directory)
    
    # Save as HTML for interactive viewing
    html_path = f"{figure_directory}/Forecast_error.html"
    fig.write_html(html_path)
    print(f"✓ Forecast error heatmap saved to: {html_path}")
    
    # Also save as static PNG
    try:
        png_path = f"{figure_directory}/Forecast_error.png"
        fig.write_image(png_path, width=1200, height=int((len(error)+1)/2)*400, scale=2)
        print(f"✓ Forecast error heatmap saved to: {png_path}")
    except Exception as e:
        print(f"Note: Could not save PNG (kaleido may not be installed). HTML version available.")

    return fig


def calculate_error_by_hour(forecast_df, perfect_foresight_df, metric='MAE'):
    """Compute error by forecast hour across all dates."""
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


def plot_forecast_model_errors_MAE(file_name,
                                   perfect_foresight_forecast_dict,
                                   mean_persistence_forecast_dict,
                                   additive_no_regressor_forecast=None,
                                   additive_with_regressors_forecast=None,
                                   autoregressive_NN_no_regressors_forecast=None,
                                   autoregressive_NN_with_regressors_forecast=None,
                                   random_forest_no_regressors_forecast=None,
                                   random_forest_with_regressors_forecast=None,
                                   metric='MAE'):
    """Generate line plot of MAE/RMSE errors by forecast hour for all models."""
    
    # Define colors and line styles for each model (matching publication style)
    colors = {
        'Mean Persistence': 'rgb(64, 164, 223)',    # Blue
        'Additive No Regressors': 'rgb(255, 127, 14)',      # Orange
        'Additive With Regressors': 'rgb(255, 187, 120)',   # Light Orange
        'Autoregressive NN No Regressors': 'rgb(44, 160, 44)',    # Green
        'Autoregressive NN With Regressors': 'rgb(152, 223, 138)', # Light Green
        'Random Forest No Regressors': 'rgb(214, 39, 40)',        # Red
        'Random Forest With Regressors': 'rgb(255, 152, 150)'     # Light Red
    }
    
    dash_styles = {
        'Mean Persistence': 'solid',
        'Additive No Regressors': 'dash',
        'Additive With Regressors': 'dot',
        'Autoregressive NN No Regressors': 'dashdot',
        'Autoregressive NN With Regressors': 'dash',
        'Random Forest No Regressors': 'dot',
        'Random Forest With Regressors': 'dashdot'
    }
    
    # Calculate errors for each model
    errors_dict = {}
    perfect_foresight = perfect_foresight_forecast_dict["DA-LMP"]
    
    # Mean Persistence (always present)
    mean_persistence = mean_persistence_forecast_dict["DA-LMP"]
    errors_dict['Mean Persistence'] = calculate_error_by_hour(mean_persistence, perfect_foresight, metric)
    
    # Check to see what forecasts we have generated
    if additive_no_regressor_forecast is not None:
        errors_dict['Additive No Regressors'] = calculate_error_by_hour(
            additive_no_regressor_forecast, perfect_foresight, metric)
    
    if additive_with_regressors_forecast is not None:
        errors_dict['Additive With Regressors'] = calculate_error_by_hour(
            additive_with_regressors_forecast, perfect_foresight, metric)
    
    if autoregressive_NN_no_regressors_forecast is not None:
        errors_dict['Autoregressive NN No Regressors'] = calculate_error_by_hour(
            autoregressive_NN_no_regressors_forecast, perfect_foresight, metric)
    
    if autoregressive_NN_with_regressors_forecast is not None:
        errors_dict['Autoregressive NN With Regressors'] = calculate_error_by_hour(
            autoregressive_NN_with_regressors_forecast, perfect_foresight, metric)
    
    if random_forest_no_regressors_forecast is not None:
        errors_dict['Random Forest No Regressors'] = calculate_error_by_hour(
            random_forest_no_regressors_forecast, perfect_foresight, metric)
    
    if random_forest_with_regressors_forecast is not None:
        errors_dict['Random Forest With Regressors'] = calculate_error_by_hour(
            random_forest_with_regressors_forecast, perfect_foresight, metric)
    
    # Create figure
    fig = go.Figure()
    
    # Plot each model with publication styling
    for model_name, errors in errors_dict.items():
        hours = np.arange(len(errors))
        
        fig.add_trace(go.Scatter(
            x=hours,
            y=errors,
            mode='lines',
            name=model_name,
            line=dict(
                color=colors.get(model_name, '#000000'),
                width=3,  # Thicker lines for publication
                dash=dash_styles.get(model_name, 'solid')
            ),
            hovertemplate='<b>%{fullData.name}</b><br>' +
                         'Hour: %{x}<br>' +
                         f'{metric}: $%{{y:.2f}}/MWh<br>' +
                         '<extra></extra>'
        ))
    
    # Update layout to match publication style
    fig.update_layout(
        xaxis_title="Forecast Hour",
        yaxis_title=f"{metric} ($/MWh)",
        xaxis=dict(
            title_font=dict(size=28, family='Arial, sans-serif'),  # Larger for publication
            tickfont=dict(size=24, family='Arial, sans-serif'),    # Larger for publication
            showgrid=True,
            zeroline=False,
            linecolor='black',
            linewidth=2,
            gridcolor='lightgrey',
            dtick=24,  # Tick every day
            range=[0, 168]
        ),
        yaxis=dict(
            title_font=dict(size=28, family='Arial, sans-serif'),  # Larger for publication
            tickfont=dict(size=24, family='Arial, sans-serif'),    # Larger for publication
            showgrid=True,
            zeroline=False,
            rangemode='tozero',  # Set y-axis minimum to 0
            linecolor='black',
            linewidth=2,
            gridcolor='lightgrey'
        ),
        title_font=dict(size=22, family='Arial, sans-serif'),
        legend=dict(
            title="",
            orientation='h',
            yanchor='top',
            y=1.12,
            xanchor='center',
            x=0.5,
            font=dict(size=20, family='Arial, sans-serif'),       # Larger legend font
            title_font=dict(size=22, family='Arial, sans-serif')
        ),
        plot_bgcolor='white',
        margin=dict(l=100, r=40, t=60, b=80),  # More margin for larger fonts
        showlegend=True,
        hovermode='x unified',
        height=450,  # Smaller height makes text appear larger relatively
        width=800    # Smaller width makes text appear larger relatively
    )
    
    # Save and show the results. First make sure there is a directory to save in
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, os.pardir, os.pardir))
    figure_directory = os.path.join(repo_root, "Core_Models", "HydroBoost", "generated_data", file_name, "Figures")
    
    if not os.path.exists(figure_directory):
        os.makedirs(figure_directory, exist_ok=True)
        print(figure_directory)
    
    # Save as HTML for interactive viewing
    html_path = os.path.join(figure_directory, f'Forecast_error_{metric}.html')
    fig.write_html(html_path)
    print(f"✓ Forecast {metric} error plot saved to: {html_path}")
    
    # Also save as static PNG
    try:
        png_path = os.path.join(figure_directory, f'Forecast_error_{metric}.png')
        fig.write_image(png_path, width=800, height=450, scale=2)
        print(f"✓ Forecast {metric} error plot saved to: {png_path}")
    except Exception as e:
        print(f"Note: Could not save PNG (kaleido may not be installed). HTML version available.")
    
    return fig
