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
from prophet import Prophet # run using hboost_py310 in Anaconda Prompt

# Ensure current directory on path for local imports
dir_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, dir_path)
from helper_functions import generate_path
from import_data import read_price_and_forecasting
from perfect_foresight_model import create_perfect_foresight_forecast

# numpy compatibility
np.float_ = np.float64


def additive_model_no_regressors(price_df, perfect_dict, project):
    """
    Forecast DA-LMP with Prophet (no regressors).
    Saves:
      generated_data/<project>/Market/Additive_Models/Additive_model_no_regressors.csv
    """
    # Prepare data for Prophet
    df_p = price_df[['DA-LMP']].reset_index().rename(columns={'Date-Time':'ds','DA-LMP':'y'})
    m = Prophet(changepoint_prior_scale=1, changepoint_range=1)
    m.fit(df_p)

    # Predict 8 days ahead, skip first 25 hours
    future = m.make_future_dataframe(periods=24*8)
    preds = m.predict(future)['yhat'].iloc[25:].to_numpy()

    # Template for 6-day, 24h blocks
    template = perfect_dict['DA-LMP'].T.copy()
    for i in range(template.shape[0]):
        start = i * 24
        block = preds[start:start+144]
        if block.shape[0] < 144:
            # pad with last value
            pad_len = 144 - block.shape[0]
            block = np.concatenate([block, np.full(pad_len, block[-1] if block.size>0 else np.nan)])
        template.iloc[i, 24:] = block

    result = template.T
    out_dir = generate_path([project, 'Market', 'Additive_Models'])
    result.to_csv(os.path.join(out_dir, 'Additive_model_no_regressors.csv'), index=True)
    print(f"→ Saved no_regressors to: {out_dir}")
    return result


def additive_model_with_regressors(price_df, perfect_dict, features, project):
    """
    Forecast DA-LMP with Prophet (with regressors).
    Saves:
      generated_data/<project>/Market/Additive_Models/Additive_model_with_regressors.csv
    """
    if not features:
        print("No regressors to add; skipping model.")
        return None

    # Prepare data including regressors
    df_p = price_df.reset_index().rename(columns={'Date-Time':'ds','DA-LMP':'y'})
    df_p = df_p[['ds','y'] + features]
    m = Prophet(changepoint_prior_scale=1, changepoint_range=1)
    for f in features:
        m.add_regressor(f)
    m.fit(df_p)

    # Extend last 24h for next 6 days
    last = df_p.iloc[-24:]
    ext = df_p.copy()
    for d in range(1,7):
        tmp = last.copy()
        tmp['ds'] = tmp['ds'] + timedelta(days=d)
        ext = pd.concat([ext, tmp], ignore_index=True)

    preds = m.predict(ext)['yhat'].iloc[25:].to_numpy()

    template = perfect_dict['DA-LMP'].T.copy()
    for i in range(template.shape[0]):
        start = i * 24
        block = preds[start:start+144]
        if block.shape[0] < 144:
            pad_len = 144 - block.shape[0]
            block = np.concatenate([block, np.full(pad_len, block[-1] if block.size>0 else np.nan)])
        template.iloc[i, 24:] = block

    result = template.T
    out_dir = generate_path([project, 'Market', 'Additive_Models'])
    result.to_csv(os.path.join(out_dir, 'Additive_model_with_regressors.csv'), index=True)
    print(f"→ Saved with_regressors to: {out_dir}")
    return result


if __name__ == '__main__':
    # project = "Model_A"
    # project = "Model_B"
    project = "Model_C"
    # project = "Model_D"
    price_df, features = read_price_and_forecasting(project)
    perfect = create_perfect_foresight_forecast(price_df, project)
    additive_model_no_regressors(price_df, perfect, project)
    additive_model_with_regressors(price_df, perfect, features, project)
