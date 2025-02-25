# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import os
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def input_data_plot(df):

    fig = px.line(df)

    fig.update_layout(title=dict(text="Input Data from Spreadsheet", xanchor="center", x=.5, font=dict(size=20)),
                      margin=dict(t=80),
                      legend=dict(title="Data",
                                  font=dict(size=12),
                                  xanchor='right',
                                  x=.99,
                                  yanchor='top',
                                  y=.99,
                                  bgcolor='rgba(255, 255, 255, 0.4)',
                                  borderwidth=1,
                                  ),
                        modebar=dict(bgcolor='rgba(255, 255, 255, 0.3)'),
                        font=dict(size=14),
                        )
    
    # If no regressors are give all values are price
    if len(df.columns) == 4:
        fig.update_yaxes(title="Price ($)")

    return fig


def plot_forecast_model_errors(file_name,
                               perfect_foresight_forecast_dict,
                               mean_persistence_forecast_dict,
                               additive_no_regressor_forecast=None,
                               additive_with_regressors_forecast=None,):
    
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
        title=dict(text="Comparison of Forecasting Errors (model-actual)", x=.5, xanchor='center', font=dict(size=25)),
        # height=int((len(error)+1)/2)*400,
        # width=1200
    )

    # Save and show the results. First make sure there is a directory to save in
    figure_directory = os.path.join(os.getcwd(), f"Core_Models/HydroBoost/generated_data/{file_name}/Figures")
    if os.path.exists(figure_directory) == False:
        os.mkdir(figure_directory)
    fig.write_html(f"{figure_directory}/Forecast_error.html")

    print("Generated error figure.")

    return fig


