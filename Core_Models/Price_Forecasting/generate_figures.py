# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED


import os
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
    fig.write_html(f"{figure_directory}/Forecast_error.html")

    return fig
