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
import pandas as pd
'''
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
    mod_data = mod_data.dropna()
    mod_data_selected = mod_data.loc[obs_data_no_nans.index]

    return obs_data_no_nans, mod_data_selected
'''
def sync_obs_with_model(obs_data, model_data):
    """
    Sincroniza el conjunto de datos del modelo con el conjunto de datos de observación,
    seleccionando solo los días presentes en el modelo y que no son NaN.

    Parameters
    ----------
    model_data : pd.DataFrame
        DataFrame del modelo.
    observation_data : pd.DataFrame
        DataFrame de observación.

    Returns
    -------
    pd.DataFrame
        DataFrame de observación sincronizado con el modelo.
    """
    #eliminar valores menores o iguales a cero
    model_data = model_data[model_data['AbsBrC370'] > 0]
    model_data = model_data.dropna()

    # Eliminar duplicados en los índices, si los hay
    model_data = model_data[~model_data.index.duplicated(keep='first')]
    obs_data = obs_data[~obs_data.index.duplicated(keep='first')]

    obs = obs_data.loc[obs_data.index.isin(model_data.index)]
    model = model_data.loc[model_data.index.isin(obs.index)]

    return obs, model
