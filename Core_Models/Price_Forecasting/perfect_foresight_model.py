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
import pandas as pd
import numpy as np
from datetime import timedelta
from helper_functions import generate_path


def initialize_empty_dataframe_for_forecast(df):
    """
    Initialize an empty DataFrame for a 7-day forecast.
    Rows = forecast start dates, Columns = hours 0..167.
    """
    dates = sorted(pd.to_datetime(list({d.date() for d in df.index})))
    hours = list(range(7 * 24))
    return pd.DataFrame(index=dates, columns=hours)

def create_perfect_foresight_forecast(price_df, file_name):
    """
    Create forecasts perfect-foresight (7 days) for:
      - 'DA-LMP'
      - 'RT-LMP'
      - 'Regulation Up'
      - 'Regulation Down'
      - 'Spinning Reserve'
      - 'Delta RU'
      - 'Delta RD'
    """
    df = price_df.copy()

    def get_forecast(template, col):
        for date in template.index:
            arr = df.loc[date: date + timedelta(days=7) - timedelta(hours=1), col].to_numpy()
            if arr.size < 168:
                last_day = arr[-24:]
                arr = np.concatenate([arr, np.tile(last_day, (168 - arr.size + 23) // 24)])[:168]
            template.loc[date] = arr
        return template.T

    empty = initialize_empty_dataframe_for_forecast(df)
    LMP_fc = get_forecast(empty.copy(), "DA-LMP")
    RT_LMP_fc = get_forecast(empty.copy(), "RT-LMP")
    RU_fc  = get_forecast(empty.copy(), "Regulation Up")
    RD_fc  = get_forecast(empty.copy(), "Regulation Down")
    SR_fc  = get_forecast(empty.copy(), "Spinning Reserve")
    RU_sig_fc = get_forecast(empty.copy(), "Delta RU")
    RD_sig_fc = get_forecast(empty.copy(), "Delta RD")

    for fc in (LMP_fc, RT_LMP_fc, RU_fc, RD_fc, SR_fc, RU_sig_fc, RD_sig_fc):
        fc.index = fc.index.astype(str)

    out_dir = generate_path([file_name, "Market", "Perfect_foresight"])

    LMP_fc.to_csv(os.path.join(out_dir, "DA_LMP.csv"))
    RT_LMP_fc.to_csv(os.path.join(out_dir, "RT_LMP.csv"))
    RU_fc.to_csv(os.path.join(out_dir, "Regulation_up.csv"))
    RD_fc.to_csv(os.path.join(out_dir, "Regulation_down.csv"))
    SR_fc.to_csv(os.path.join(out_dir, "Spin.csv"))
    RU_sig_fc.to_csv(os.path.join(out_dir, "Delta_RU.csv"))
    RD_sig_fc.to_csv(os.path.join(out_dir, "Delta_RD.csv"))

    print(f"→ Perfect-foresight files saved to: {out_dir}")

    return {
        "DA-LMP": LMP_fc,
        "RT-LMP": RT_LMP_fc,
        "Regulation Up": RU_fc,
        "Regulation Down": RD_fc,
        "Spinning Reserve": SR_fc,
        "Delta RU": RU_sig_fc,
        "Delta RD": RD_sig_fc,
    }

if __name__ == "__main__":
    # project = "Model_A"
    # project = "Model_B"
    project = "Model_C"
    # project = "Model_D"
    from import_data import read_price_and_forecasting
    price_df, _ = read_price_and_forecasting(project)
    create_perfect_foresight_forecast(price_df, project)
