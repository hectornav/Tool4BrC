"""
Aerosol Absorption Calculator

This module is designed for calculating the absorption of light by various types
of aerosols, including Brown Carbon (BrC), using their optical properties. It provides
functions to calculate absorption based on the extinction coefficient, scattering coefficient, 
phase function, density, and effective radius of the aerosols. These calculations are crucial 
for atmospheric science research, particularly in the study of aerosol impacts on climate 
and air quality.

The module supports calculations from both pandas DataFrame and xarray Dataset, accommodating
different data structures commonly used in atmospheric science.

Functions:
- calculate_absorption: Computes absorption for aerosols from pandas DataFrame.
- calculate_absorption4brc: Similar to calculate_absorption but specifically tuned for BrC aerosols.
- calculate_absorption4brc_from_netcdf: Computes absorption for aerosols from xarray Dataset.
"""

import pandas as pd
import xarray as xr

def calculate_absorption(model, opt_params, delz=1):
    """
    Calculate the absorption for different aerosols using their optical properties.

    This function calculates the absorption of light for different aerosols 
    based on their optical properties (extinction coefficient, scattering coefficient, 
    phase function, density, and effective radius). These properties are provided 
    through the opt_params dictionary.

    Parameters
    ----------
    model : pd.DataFrame
        A DataFrame where each column represents an aerosol and each row represents 
        a different measurement or simulation. The values in the DataFrame are 
        aerosol concentrations in [g/m3].
    opt_params : dict
        A dictionary containing the optical properties for each aerosol. Each 
        key should be an aerosol name matching the columns in the model DataFrame. 
        Each value should be another dictionary with keys: 'qext' (extinction coefficient),
        'qsca' (scattering coefficient), 'vphi' (phase function), 'dens' (density, [g/cm³]),
        'reff' (effective radius).
    delz : float, optional
        The thickness of the atmospheric layer for which the absorption is calculated, 
        in the same units as the aerosol concentrations in the model DataFrame. 
        By default, this is set to 1, implying the concentrations are per unit thickness, [m].

    Returns
    -------
    df_abs : pd.DataFrame
        A DataFrame with the same structure as the model DataFrame, but with values 
        representing the calculated absorption for each aerosol.
    """
    aero = model.columns.tolist()
    #convert model from [ug/m3] to [g/m3]
    #model = model / 1e6 
    # get the optical properties
    abs_values = {}
    for i in aero:
        absorption = ((3.0 * (opt_params[i]['qext'] - opt_params[i]['qsca']) * opt_params[i]['vphi']) / 
                      (4.0 * opt_params[i]['dens'] * opt_params[i]['reff'])) * model[i] * delz
        abs_values[i] = absorption

    # create a dataframe with the absorption values
    df_abs = pd.DataFrame(abs_values, index=model.index)

    return df_abs


def calculate_absorption4brc(model, opt_params, delz=1):
    """
    Calculate the absorption for different aerosols using their optical properties.

    This function calculates the absorption of light for different aerosols 
    based on their optical properties (extinction coefficient, scattering coefficient, 
    phase function, density, and effective radius). These properties are provided 
    through the opt_params dictionary.

    Parameters
    ----------
    model : pd.DataFrame
        A DataFrame where each column represents an aerosol and each row represents 
        a different measurement or simulation. The values in the DataFrame are 
        aerosol concentrations in [g/m3].
    opt_params : dict
        A dictionary containing the optical properties for each aerosol. Each 
        key should be an aerosol name matching the columns in the model DataFrame. 
        Each value should be another dictionary with keys: 'qext' (extinction coefficient),
        'qsca' (scattering coefficient), 'vphi' (phase function), 'dens' (density, [g/cm³]),
        'reff' (effective radius).
    delz : float, optional
        The thickness of the atmospheric layer for which the absorption is calculated, 
        in the same units as the aerosol concentrations in the model DataFrame. 
        By default, this is set to 1, implying the concentrations are per unit thickness, [m].

    Returns
    -------
    df_abs : pd.DataFrame
        A DataFrame with the same structure as the model DataFrame, but with values 
        representing the calculated absorption for each aerosol in [Mm-1].
    """
    #aero = model.columns.tolist()
    aero = ['poabrc', 'soabrc']
    #print(model['poabrc'])
    #convert model from [ug/m3] to [g/m3]
    #model = model / 1e6
    # get the optical properties
    abs_values = {}
    for i in aero:
        absorption = ((3.0 * (opt_params[i]['qext'] - opt_params[i]['qsca']) * opt_params[i]['vphi']) / 
                      (4.0 * opt_params[i]['dens'] * opt_params[i]['reff'])) * model['pm2p5pbrc' if i=='poabrc' else 'pm2p5sbrc'] * delz
        abs_values[i] = absorption

    # create a dataframe with the absorption values
    df_abs = pd.DataFrame(abs_values, index=model.index)

    return df_abs

def calculate_absorption4brc_from_netcdf(model, opt_params, delz=1):
    """
    Calculate the absorption for different aerosols using their optical properties.

    This function calculates the absorption of light for different aerosols 
    based on their optical properties (extinction coefficient, scattering coefficient, 
    phase function, density, and effective radius). These properties are provided 
    through the opt_params dictionary.

    Parameters
    ----------
    model : xr.Dataset
        A Dataset where each DataArray represents an aerosol and each data point represents 
        a different measurement or simulation. The values in the Dataset are 
        aerosol concentrations in [g/m3].
    opt_params : dict
        A dictionary containing the optical properties for each aerosol. Each 
        key should be an aerosol name matching the DataArrays in the model Dataset. 
        Each value should be another dictionary with keys: 'qext' (extinction coefficient),
        'qsca' (scattering coefficient), 'vphi' (phase function), 'dens' (density, [g/cm³]),
        'reff' (effective radius).
    delz : float, optional
        The thickness of the atmospheric layer for which the absorption is calculated, 
        in the same units as the aerosol concentrations in the model Dataset. 
        By default, this is set to 1, implying the concentrations are per unit thickness, [m].

    Returns
    -------
    ds_abs : xr.Dataset
        A Dataset with the same structure as the model Dataset, but with values 
        representing the calculated absorption for each aerosol.
    """
    aero = list(model.data_vars.keys())

    # get the optical properties
    abs_values = {}
    for i in aero:
        absorption = ((3.0 * (opt_params[i]['qext'] - opt_params[i]['qsca']) * opt_params[i]['vphi']) / 
                      (4.0 * opt_params[i]['dens'] * opt_params[i]['reff'])) * model[i] * delz
        abs_values[i] = absorption

    # create a Dataset with the absorption values
    ds_abs = xr.Dataset(abs_values)

    return ds_abs

def calculateabs4oa(model, opt_params, delz=1):
    """
    Calculate the absorption for OA using their optical properties.

    This function calculates the absorption of light for OA 
    based on their optical properties (extinction coefficient, scattering coefficient, 
    phase function, density, and effective radius). These properties are provided 
    through the opt_params dictionary.

    Parameters
    ----------
    model : pd.DataFrame
        A DataFrame with OA and each row represents 
        a different measurement or simulation. The values in the DataFrame are 
        aerosol concentrations in [g/m3].
    opt_params : dict
        A dictionary containing the optical properties for OA. Each 
        key should be an aerosol name matching the columns in the model DataFrame. 
        Each value should be another dictionary with keys: 'qext' (extinction coefficient),
        'qsca' (scattering coefficient), 'vphi' (phase function), 'dens' (density, [g/cm³]),
        'reff' (effective radius).
    delz : float, optional
        The thickness of the atmospheric layer for which the absorption is calculated, 
        in the same units as the aerosol concentrations in the model DataFrame. 
        By default, this is set to 1, implying the concentrations are per unit thickness, [m].

    Returns
    -------
    df_abs : pd.DataFrame
        A DataFrame with the same structure as the model DataFrame, but with values 
        representing the calculated absorption for each aerosol.
    """
    aero = model.columns.tolist()
    #convert model from [ug/m3] to [g/m3]
    #model = model / 1e6 
    # get the optical properties
    abs_values = {}
    for i in aero:
        absorption = ((3.0 * (opt_params[i]['qext'] - opt_params[i]['qsca']) * opt_params[i]['vphi']) / 
                      (4.0 * opt_params[i]['dens'] * opt_params[i]['reff'])) * model[i] * delz
        abs_values[i] = absorption

    # create a dataframe with the absorption values
    df_abs = pd.DataFrame(abs_values, index=model.index)

    return df_abs