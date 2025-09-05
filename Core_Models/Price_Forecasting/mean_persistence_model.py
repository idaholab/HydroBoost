#

# HydroBoost Model

# Main Authors:
# Jonghwan Kwon; Argonne National Laboratory; kwonj@anl.gov
# Carlos Josue Lopez; Argonne National Laboratory; clopezsalgado@anl.gov
# Alberto Grimaldi; Argonne National Laboratory; agrimaldi@anl.gov

# Current version: 2.0
# Last update: 07.31.2025

#

import os
import sys
import pandas as pd
import numpy as np
from datetime import timedelta

# Ensure we can import helper_functions from this directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from helper_functions import generate_path


def create_mean_persistence_forecast(perfect_foresight_forecast_dict, file_name, num_days=7):
    """
    Creates a 7-day forecast using mean persistence and saves CSVs.

    perfect_foresight_forecast_dict: dict of DataFrames for keys:
      'DA-LMP', 'Regulation Up', 'Regulation Down', 'Spinning Reserve'
    file_name: project name under generated_data (e.g. 'Model_A')
    num_days: look-back window for mean persistence
    Returns a dict of the forecast DataFrames.
    """
    def get_forecast(df_tpl):
        # df_tpl: rows = forecast dates, cols = 0..167 hours
        for idx in range(len(df_tpl.index)):
            if idx == 0:
                means = df_tpl.iloc[0, :24].to_numpy()
            elif idx < num_days:
                means = df_tpl.iloc[0:idx, :24].mean().to_numpy()
            else:
                means = df_tpl.iloc[idx-num_days:idx, :24].mean().to_numpy()
            df_tpl.iloc[idx, 24:] = np.tile(means, 6)
        return df_tpl.T

    # Build forecasts for each category
    LMP_fc = get_forecast(perfect_foresight_forecast_dict["DA-LMP"].T.copy())
    RU_fc  = get_forecast(perfect_foresight_forecast_dict["Regulation Up"].T.copy())
    RD_fc  = get_forecast(perfect_foresight_forecast_dict["Regulation Down"].T.copy())
    SR_fc  = get_forecast(perfect_foresight_forecast_dict["Spinning Reserve"].T.copy())

    # Prepare output directory
    out_dir = generate_path([file_name, "Market", "Mean_persistence"])

    # Save CSVs under generated_data/<file_name>/Market/Mean_persistence
    LMP_fc.to_csv(os.path.join(out_dir, "DA_LMP.csv"))
    RU_fc.to_csv(os.path.join(out_dir, "Regulation_up.csv"))
    RD_fc.to_csv(os.path.join(out_dir, "Regulation_down.csv"))
    SR_fc.to_csv(os.path.join(out_dir, "Spin.csv"))

    print(f"→ Mean persistence CSVs saved in: {out_dir}")

    return {
        "DA-LMP": LMP_fc,
        "Regulation Up": RU_fc,
        "Regulation Down": RD_fc,
        "Spinning Reserve": SR_fc,
    }

if __name__ == "__main__":
    # Example invocation
    # project = "Model_A"
    # project = "Model_B"
    project = "Model_C"
    # project = "Model_D"
    from import_data import read_price_and_forecasting
    from perfect_foresight_model import create_perfect_foresight_forecast

    price_df, _ = read_price_and_forecasting(project)
    perfect_dict = create_perfect_foresight_forecast(price_df, project)
    create_mean_persistence_forecast(perfect_dict, project)
