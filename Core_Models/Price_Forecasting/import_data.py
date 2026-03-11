# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED


import os
import shutil
import pandas as pd
import openpyxl

# Make it runnable with the IDE "Run" button (works both as package and as script)
try:
    from .helper_functions import generate_path
except ImportError:
    from helper_functions import generate_path

# 1) Directory of this script:
#    C:\Users\agrimaldi\a-leaf-dev_HydroBoost_AG_version\Core_Models\Price_Forecasting
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# 2) Ascend two levels to reach project root:
#    C:\Users\agrimaldi\a-leaf-dev_HydroBoost_AG_version
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir, os.pardir))

# 3) Directory where input Excel files live:
INPUT_DIR = os.path.join(PROJECT_ROOT, "Input_Spreadsheets")


def read_price_and_forecasting(file_name):
    """
    Read the 'Hourly Price and Forecasting' sheet from:
      <PROJECT_ROOT>/Input_Spreadsheets/{file_name}.xlsx

    Returns:
        df               Pandas DataFrame of prices and forecasts, with Date-Time as index
        forecast_features  List of column names used as forecasting features (if any)
    """
    xlsx_path = os.path.join(INPUT_DIR, f"{file_name}.xlsx")
    if not os.path.isfile(xlsx_path):
        raise FileNotFoundError(f"Excel file not found: {xlsx_path}")

    df = pd.read_excel(
        xlsx_path,
        sheet_name="Hourly Price and Forecasting",
        skiprows=[0],
        index_col="Date-Time"
    ).bfill()  # backfill any missing values

    # If there are more than 4 columns, the extra columns are forecasting features
    if len(df.columns) > 4:
        forecast_features = df.columns.to_list()[4:]
    else:
        forecast_features = None

    return df, forecast_features


def read_water_flow(file_name):
    """
    Read the 'Water Flow and Availability' sheet from:
      <PROJECT_ROOT>/Input_Spreadsheets/{file_name}.xlsx

    Saves the hourly hydro flow data to:
      generated_data/{file_name}/Hydro/Hourly_flow.csv

    Returns:
        df  Pandas DataFrame of hourly water flows
    """
    xlsx_path = os.path.join(INPUT_DIR, f"{file_name}.xlsx")
    if not os.path.isfile(xlsx_path):
        raise FileNotFoundError(f"Excel file not found: {xlsx_path}")

    df = pd.read_excel(
        xlsx_path,
        sheet_name="Water Flow and Availability",
        skiprows=[0],
        index_col="Date-Time"
    ).dropna(axis=1, how="all")

    out_dir = generate_path([file_name, "Hydro"])
    df.to_csv(os.path.join(out_dir, "Hourly_flow.csv"))
    print(f"→ Hydro hourly flow saved to: {out_dir}")

    return df


def read_daily_water_flow(file_name):
    """
    Read the 'Daily Flow Constraints' sheet from:
      <PROJECT_ROOT>/Input_Spreadsheets/{file_name}.xlsx

    Saves the daily hydro flow constraints to:
      generated_data/{file_name}/Hydro/Daily_flow.csv

    Returns:
        df  Pandas DataFrame of daily water flow constraints
    """
    xlsx_path = os.path.join(INPUT_DIR, f"{file_name}.xlsx")
    if not os.path.isfile(xlsx_path):
        raise FileNotFoundError(f"Excel file not found: {xlsx_path}")

    df = pd.read_excel(
        xlsx_path,
        sheet_name="Daily Flow Constraints",
        skiprows=[0]
    )
    # Convert date column to datetime and set as index
    df["Date-time"] = pd.to_datetime(df["Date-time"])
    df = df.set_index("Date-time")

    out_dir = generate_path([file_name, "Hydro"])
    df.to_csv(os.path.join(out_dir, "Daily_flow.csv"))
    print(f"→ Hydro daily flow saved to: {out_dir}")

    return df


def read_RES_profiles(file_name):
    """
    Read the 'RES' sheet from:
      <PROJECT_ROOT>/Input_Spreadsheets/{file_name}.xlsx

    Saves the renewable profiles to:
      generated_data/{file_name}/RES/RES_profiles.csv

    Returns:
        df  Pandas DataFrame of RES generation profiles
    """
    xlsx_path = os.path.join(INPUT_DIR, f"{file_name}.xlsx")
    if not os.path.isfile(xlsx_path):
        raise FileNotFoundError(f"Excel file not found: {xlsx_path}")

    df = pd.read_excel(
        xlsx_path,
        sheet_name="RES",
        skiprows=[0],
        index_col="Date-Time"
    ).dropna(axis=1, how="all")

    out_dir = generate_path([file_name, "RES"])
    df.to_csv(os.path.join(out_dir, "RES_profiles.csv"))
    print(f"→ RES profiles saved to: {out_dir}")

    return df


def copy_Excel(file_name):
    """
    Copy the input Excel file and modify it without corrupting macros.
    """

    import os
    import shutil
    import xlwings as xw

    xlsx_path = os.path.join(INPUT_DIR, f"{file_name}.xlsx")  # Ensure correct extension
    if not os.path.isfile(xlsx_path):
        raise FileNotFoundError(f"Excel file not found: {xlsx_path}")

    out_dir = generate_path([file_name])
    dest_path = os.path.join(out_dir, f"{file_name}.xlsx")  # Save as .xlsx always
    
    # Copy the file first
    shutil.copy2(xlsx_path, dest_path)
    print(f"→ Excel file copied to: {dest_path}")

    # Open the copied file using xlwings to modify it
    try:
        app = xw.App(visible=False)  # Start Excel in the background
        wb = app.books.open(dest_path)
        
        if 'Simulation Configuration' in [sheet.name for sheet in wb.sheets]:
            ws = wb.sheets['Simulation Configuration']
            # Check rows 3-6 (assuming header is in row 2)
            for row in range(3, 7):
                cell_value = ws.range(f"B{row}").value  # Column B is the Run_Flag column
                if cell_value == 'true':
                    ws.range(f"B{row}").value = True
                elif cell_value == 'false':
                    ws.range(f"B{row}").value = False
        
        wb.save()  # Save as .xlsx
        wb.close()
        app.quit()
        print(f"→ Run_Flag values fixed in {dest_path}")
    except Exception as e:
        print(f"Warning: Could not modify {dest_path}. Error: {e}")
        if app:
            app.quit()


if __name__ == "__main__":
    # project_name = "Model_A"
    # project_name = "Model_B"
    project_name = "Model_C"
    # project_name = "Model_D"

    # 1) Load price and forecasting data
    price_df, features = read_price_and_forecasting(project_name)

    # 2) Load and save hydro flow data
    _ = read_water_flow(project_name)
    _ = read_daily_water_flow(project_name)

    # 3) Load and save renewable profiles
    _ = read_RES_profiles(project_name)
