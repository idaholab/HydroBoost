# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import os
import argparse

# Select the project folder in teh Simulation_Results directory
#project = "Project_1_parallel"


# Create color scales for cell backgrounds
def get_cell_color(value, is_total_row=False, column=None):
    
    # Don't color the total row if specified
    if is_total_row:
        return '#F0F0F0'  # Light gray background for Total row
    
    # Remove any non-numeric characters and convert to float
    try:
        if 'M' in value:
            num_value = float(value.replace('M', '')) * 1000000
        elif 'K' in value:
            num_value = float(value.replace('K', '')) * 1000
        else:
            num_value = float(value)
        
        # Get max value for this column (for scaling)
        max_val = max_values.get(column, 1000000)
        
        # Define color based on value
        if num_value > 0:
            # Green scale for positive: lighter green -> darker green
            intensity = min(abs(num_value) / max_val, 1) 
            # Make the green more intense overall
            r = int(255 * (1 - intensity * 0.9))
            g = 255
            b = int(255 * (1 - intensity * 0.9))
            return f'rgb({r},{g},{b})'
        elif num_value < 0:
            # Red scale for negative: lighter red -> darker red
            intensity = min(abs(num_value) / max_val, 1)
            # Make the red more intense overall
            r = 255
            g = int(255 * (1 - intensity * 0.9))
            b = int(255 * (1 - intensity * 0.9))
            return f'rgb({r},{g},{b})'
        else:
            # White for zero
            return 'rgb(255,255,255)'
    except:
        return 'rgb(255,255,255)'  # Default value, should not happen

# Format the data based on value
    def format_value(x):
        if abs(x) >= 1000000:
            return f"{x/1000000:.1f}M"
        elif abs(x) >= 1000:
            return f"{x/1000:.1f}K"
        else:
            return f"{x:.1f}"

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Please specify the name of the Project")

    parser.add_argument("--project", required=True, type=str)
    args = parser.parse_args()
    project = args.project

    fig_path = os.path.join("Simulation_Results",project,"Figures")

    if not os.path.exists(fig_path):
        os.makedirs(fig_path)
    

    ##############################################
    # Import and format the data
    ##############################################

    df_hydro = pd.read_csv(os.path.join("Simulation_Results", project, f"ALEAF_HydroBoost_{project}__hydro_dispatch.csv"))
    df_plant = pd.read_csv(os.path.join("Simulation_Results", project, f"ALEAF_HydroBoost_{project}__plant_dispatch.csv"))
    df_storage = pd.read_csv(os.path.join("Simulation_Results", project, f"ALEAF_HydroBoost_{project}__storage_dispatch.csv"))

    dates = pd.date_range(start='2025-01-01', periods=8760, freq='h')

    df_hydro = df_hydro.drop(columns=['day', 'hour', 'time',])
    df_plant = df_plant.drop(columns=['day', 'hour', 'time',])
    df_storage = df_storage.drop(columns=['day', 'hour', 'time',])

    df_hydro.index =  dates
    df_plant.index =  dates
    df_storage.index =  dates

    df_hydro = df_hydro.rename(columns={"u_H_jht": "Hydro Plant Commitment",
                                        "a_H_jht": "Hydro Plant Startup",
                                        "z_H_jht": "Hydro Plant Shutdown",
                                        "p_H_jht": "Hydro Power Output (MW)",
                                        "u_jht": "Hydro Water Discharge (ft^3/s)",
                                        "u_0_jht": "Minimum Water Discharge (ft^3/s)",
                                        "u_1_jht": "Water Discharge from Piecewise Block 1 (ft^3/s)",
                                        "u_2_jht": "Water Discharge from Piecewise Block 2 (ft^3/s)",
                                        "u_3_jht": "Water Discharge from Piecewise Block 3 (ft^3/s)",
                                        "u_4_jht": "Water Discharge from Piecewise Block 4 (ft^3/s)",
                                        })

    df_plant = df_plant.rename(columns={'p_B_DT_ht': 'Battery Discharged (MW)', 
                                        'p_B_CT_ht': 'Battery Charging (MW)', 
                                        'p_GB_ht': 'Grid to Battery (MW)', 
                                        'p_HT_ht': 'Total Hydropower (MW)', 
                                        'p_HB_ht': 'Hydropower to Battery (MW)', 
                                        'p_HG_ht': 'Hydropower to Grid',
                                        's_ht': 'Spillage (ft^3/s)', 
                                        'u_ht': 'Total Water Discharge (ft^3/s)', 
                                        'e_H_ht': 'Volume of Reservoir (Acre*feet)', 
                                        'I_ht': 'Water Inflow to Reservoir (ft^3/s)',
                                        })

    df_storage = df_storage.rename(columns={'u_B_iht': 'Binary Battery Charging Mode',
                                            'p_B_D_iht': 'Battery Discharge (MW)', 
                                            'p_B_C_iht': 'Battery Charge (MW)', 
                                            'e_B_iht': 'State of Charge (MWh)', 
                                            'r_RU_D_iht': 'Reserve for Regulation Up Discharge (MW)', 
                                            'r_RD_D_iht': 'Reserve for Regulation Down Discharge (MW)',
                                            'r_RU_C_iht': 'Reserve for Regulation Up Charge (MW)', 
                                            'r_RD_C_iht': 'Reserve for Regulation Down Charge (MW)', 
                                            'r_RU_iht': 'Reserve Regulation Up Market Sell (MW)', 
                                            'r_RD_iht': 'Reserve Regulation Down Market Sell (MW)', 
                                            'r_SR_D_iht': 'Reserve Spinning Discharge (MW)',
                                            'r_SR_C_iht': 'Reserve Spinning Charge (MW)', 
                                            'r_SR_iht': 'Reserve Spinning Market Sell (MW)'
                                            })

    df = df_hydro.join(df_plant.join(df_storage.drop(columns=['unit_id', 'UnitGroup', 'Unit_Category', 'Unit_Type'])))
    df.index.name = 'Datetime'

    # Create revenue columns based on results and price
    df["Hydro Revenue"] = df["Hydropower to Grid"] * df["LMP"]
    df["Battery Revenue"] = df["Battery Discharged (MW)"] * df["LMP"]
    df["Battery Cost"] = df["Grid to Battery (MW)"] * - df["LMP"]
    df["Total Battery Revenue"] = df["Battery Revenue"] + df["Battery Cost"]
    df["Hydro and Battery Revenue"] = df["Hydro Revenue"] + df["Total Battery Revenue"]

    df["Regulation Up Revenue"] = df["Reserve Regulation Up Market Sell (MW)"] * df["Reg Up Price"]
    df["Regulation Down Revenue"] = df["Reserve Regulation Down Market Sell (MW)"] * df["Reg Dn Price"]
    df["Spinning Reserve Revenue"] = df["Reserve Spinning Market Sell (MW)"] * df["Spin Price"]



    # Define the revenue columns we want to plot
    revenue_columns = ["Hydro Revenue", "Battery Revenue", "Battery Cost", 
                      "Total Battery Revenue", "Hydro and Battery Revenue"]

    # Define the colors for the plots
    colors = {
        "Hydro Revenue": "#1E88E5",             # Blue
        "Battery Revenue": "#43A047",           # Green
        "Battery Cost": "#E53935",              # Red
        "Total Battery Revenue": "#8E24AA",     # Purple
        "Hydro and Battery Revenue": "#FB8C00"  # Orange
    }


    ###########################################
    # Line plot
    ###########################################

    print(df.columns)

    fig = px.line(df.drop(columns=["unit_id", "UnitGroup", "Unit_Category", "Unit_Type"]))
    fig.write_html(f"Simulation_Results/{project}/Figures/All_data.html")


    ###########################################
    # Overview of all the monthly revenue
    ###########################################

    # Create a figure
    fig = go.Figure()

    monthly_data = df[revenue_columns].resample('M').sum()

    # Create month labels
    month_names = ['January', 'February', 'March', 'April', 'May', 'June', 
                   'July', 'August', 'September', 'October', 'November', 'December']

    x_labels = [f"{m[:3]}" for m in month_names]

    # Add yearly data traces (initially visible)
    for col in revenue_columns:
        fig.add_trace(
            go.Bar(
                x=x_labels,
                y=monthly_data[col],
                name=col,
                marker_color=colors[col],
                visible=True
            )
        )

    ###########################################
    # Each month plot for hour of the day
    ###########################################

    for month_num in range(1, 13):
        # Filter data for the current month
        month_data = df[df.index.month == month_num]
        
        # Calculate hourly averages
        monthly_hourly_data = []
        for hour in range(24):
            hour_data = month_data[month_data.index.hour == hour]
            if not hour_data.empty:
                hour_avg = hour_data[revenue_columns].mean()
                monthly_hourly_data.append(hour_avg)
            else:
                hour_avg = pd.Series({col: 0 for col in revenue_columns})
                monthly_hourly_data.append(hour_avg)
        
        hourly_avg = pd.DataFrame(monthly_hourly_data)
        hour_labels = [f"{h:02d}:00" for h in range(24)]
        
        # Add traces for this month (initially hidden)
        for col in revenue_columns:
            fig.add_trace(
                go.Bar(
                    x=hour_labels,
                    y=hourly_avg[col],
                    name=col,
                    marker_color=colors[col],
                    visible=False
                )
            )


    #############################################
    # Dropdown menu
    #############################################

    buttons = []

    # Button for yearly overview
    buttons.append(
        dict(
            label="Yearly Overview",
            method="update",
            args=[
                {"visible": [True, True, True, True, True] + [False]*60}, 
                {
                    "title": "Monthly Revenue and Cost Breakdown",
                    "xaxis": {"title": "Month", "visible": True, "domain": [0.1, 0.9]},
                    "yaxis": {"title": "Revenue ($)", "visible": True, "domain": [0.1, 0.9]},
                    "showlegend": True
                }
            ]
        )
    )

    # Add buttons for each month
    for i, month_name in enumerate(month_names):
        month_visible = [False]*5  # Hide yearly traces
        
        # Set visibility for monthly traces
        for j in range(12):
            if j == i:  # This is the month we want to show
                month_visible.extend([True, True, True, True, True])
            else:
                month_visible.extend([False, False, False, False, False])
        
        # Add button for this month
        buttons.append(
            dict(
                label=month_name,
                method="update",
                args=[
                    {"visible": month_visible},
                    {
                        "title": f"Average Hourly Revenue and Cost for {month_name}",
                        "xaxis": {"title": "Hour of Day", "visible": True, "domain": [0.1, 0.9]},
                        "yaxis": {"title": "Revenue ($)", "visible": True, "domain": [0.1, 0.9]},
                        "showlegend": True
                    }
                ]
            )
        )


    updatemenus = [
        dict(
            active=0,
            buttons=buttons,
            direction="down",
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0.99,
            xanchor="right",
            y=1.17,
            yanchor="top",
            bgcolor="#F9F9F9",
            bordercolor="#CCCCCC",
            font=dict(color="#333333", size=14)
        )
    ]


    fig.update_layout(
        updatemenus=updatemenus,
        title=dict(
            text="Monthly Revenue and Cost Breakdown",
            x=0.5,
            font=dict(size=22) 
        ),
        template="plotly_white",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
            font=dict(size=14)
        ),
        barmode="group",
        height=700,
        width=1000,
        margin=dict(l=50, r=50, t=100, b=50),
        font=dict(family="Arial, sans-serif", size=14, color="#333333")
    )


    fig.update_xaxes(
        title_text="Month",
        title_font=dict(size=16),
        tickangle=-45,
        showgrid=True,
        gridcolor="#EEEEEE",
        tickfont=dict(size=14),
        range=[-0.5, len(x_labels)-0.5]
    )

    fig.update_yaxes(
        title_text="Revenue ($)",
        title_font=dict(size=16),
        showgrid=True,
        gridcolor="#EEEEEE",
        tickfont=dict(size=14)
    )


    ###########################################
    # Summary table
    ###########################################

    yearly_total = df[revenue_columns].sum()
    table_data = pd.DataFrame(index=revenue_columns)
    table_data['Total'] = yearly_total

    # Add monthly data
    for i, month in enumerate(month_names):
        # Get the data for this month
        month_data = df[df.index.month == i+1][revenue_columns].sum()
        table_data[month[:3]] = month_data

    #Format the data based on value
    def format_value(x):
        if abs(x) >= 1000000:
            return f"{x/1000000:.1f}M"
        elif abs(x) >= 1000:
            return f"{x/1000:.1f}K"
        else:
            return f"{x:.1f}"    
    #print(table_data)
    table_data_formatted = table_data.map(format_value)
    table_data_formatted = table_data_formatted.T

    # Move the 'Total' row from top to bottom
    total_row = table_data_formatted.iloc[0:1]
    table_data_formatted = pd.concat([table_data_formatted.iloc[1:], total_row])

    # Calculate maximum values for color scaling
    max_values = {}
    for col in revenue_columns:
        # Get all values except the Total row
        values = []
        for idx, val_str in enumerate(table_data_formatted[col]):
            if idx == len(table_data_formatted) - 1:  # Skip the Total row
                continue
            try:
                if 'M' in val_str:
                    val = float(val_str.replace('M', '')) * 1000000
                elif 'K' in val_str:
                    val = float(val_str.replace('K', '')) * 1000
                else:
                    val = float(val_str)
                values.append(abs(val))
            except:
                pass
        
        # Store max value for this column
        max_values[col] = max(values) if values else 1000000  # Default if no values (shouldn't be an issue)



    table_fig = go.Figure()

    # Create cell color arrays
    cell_colors = []
    for col in revenue_columns:
        col_colors = []
        for idx, value in enumerate(table_data_formatted[col]):
            is_total_row = (idx == len(table_data_formatted) - 1)
            col_colors.append(get_cell_color(value, is_total_row, col))
        cell_colors.append(col_colors)

    # Convert the Total row to bold format
    bold_values = []
    bold_values.append([f"<b>{val}</b>" if val == "Total" else val for val in table_data_formatted.index])
    for col in revenue_columns:
        bold_values.append([f"<b>{val}</b>" if idx == (len(table_data_formatted) - 1) else val 
                          for idx, val in enumerate(table_data_formatted[col])])


    table_fig.add_trace(go.Table(
        header=dict(
            values=['<b>Period</b>'] + [f'<b>{col}</b>' for col in revenue_columns],
            line_color='darkslategray',
            fill_color='#F9F9F9',
            align='left',
            font=dict(color='black', size=14)
        ),
        cells=dict(
            values=bold_values,
            line_color='darkslategray',
            fill_color=['white'] + cell_colors,
            align='left',
            font=dict(
                color='black', 
                size=13,
                family="Arial, sans-serif"
            ),
            height=30
        )
    ))

    table_fig.update_layout(
        title=dict(
            text="Revenue and Cost Summary Table",
            x=0.5,  # Center the title
            font=dict(size=22)  # Larger title font
        ),
        height=600,
        width=1000,
        margin=dict(l=50, r=50, t=100, b=50),
        font=dict(family="Arial, sans-serif", size=14, color="#333333")
    )





    directory = f"Simulation_Results/{project}/Figures"

    os.makedirs(directory, exist_ok=True)

    fig.write_html(f"Simulation_Results/{project}/Figures/Monthly_Revenue.html")
    table_fig.write_html(f"Simulation_Results/{project}/Figures/Monthly_Summary_Table.html")
