"""
Synchronize Observational and Model Data

This module provides functionality to synchronize observational data with corresponding model data. 
It is particularly useful in atmospheric and environmental science, where it's common to compare 
observational datasets that contain missing values (NaNs) with continuous model outputs.

The main function, `sync_obs_with_model`, takes two pandas DataFrames: one representing observational 
data that may contain missing values (NaNs), and another representing model data without missing values. 
It aligns these datasets by removing days from both where observations are missing, facilitating 
direct comparison and analysis.

Functions:
- sync_obs_with_model: Aligns observational and model data by removing days with missing observations.
"""

def sync_obs_with_model(obs_data, mod_data):
    """
    Synchronize observational data with model data by removing days with missing observations.

    Parameters
    ----------
    obs_data : pd.DataFrame
        Observation DataFrame containing NaNs.
    mod_data : pd.DataFrame
        Model DataFrame without NaNs.

    Returns
    -------
    tuple
        A tuple containing two DataFrames: 
        1. The observational data without NaNs.
        2. The model data corresponding to the same days as the observational data.
    """
    obs_data_no_nans = obs_data.dropna()
    mod_data_selected = mod_data.loc[obs_data_no_nans.index]

    return obs_data_no_nans, mod_data_selected
