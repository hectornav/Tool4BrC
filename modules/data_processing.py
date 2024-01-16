"""
Data Processing Module

This module focuses on data manipulation and preparation for aerosol absorption models. 
It includes functions for specific data processing tasks, such as extracting data for 
particular stations and calculating absorption based on model concentrations.

Functions:
- get_data_for_station: Filters data for a specific station from a DataFrame.
- calculate_absorption: Computes absorption using model concentrations.
- calculate_abs_modeled: Calculates modeled absorption for a specific station.
"""
import pandas as pd

from . import aerosol_absorption_calculator as aac
from . import atmospheric_aerosol_optics as aao
from . import constants as const

def get_data_for_station(dataframe, station, columns):
    """
    Retrieve data for a specific station from a DataFrame.

    Parameters:
    dataframe (DataFrame): The DataFrame to filter.
    station (str): The name of the station.
    columns (list of str): The list of columns to include in the result.

    Returns:
    DataFrame: A DataFrame filtered for the specified station and columns.
    """
    # Check if 'station_name' or 'time' are not in the dataframe columns, return original dataframe.
    if 'station_name' not in dataframe.columns or 'time' not in dataframe.columns:
        return dataframe
    else:
        # Filter the data for the given station and select the specified columns.
        station_data = dataframe[dataframe['station_name'] == station][columns + ['time']]
        #print(station_data)
        # Set the 'time' column as the index of the DataFrame.
        return station_data.set_index('time')


def calculate_absorption(model_conc, optical_parameters):
    """
    Calculate absorption using model concentrations and optical parameters.

    Parameters:
    model_conc (DataFrame): DataFrame of model concentrations.
    optical_parameters (dict): Dictionary of optical parameters.

    Returns:
    DataFrame: DataFrame of calculated absorption.
    """
    model_abs = aac.calculate_absorption(model_conc, optical_parameters).sum(axis=1)
    return pd.DataFrame(model_abs, index=model_conc.index, columns=['AbsBrC370'])

def calculate_abs_modeled(station, model, ri_values):
    """
    Calculates the modeled absorption for a specific station based on provided refractive index values.

    Parameters:
    station (str): Station name used to filter concentration data.
    model (DataFrame): DataFrame containing model data.
    ri_values (list): List of refractive index values for different substances.

    Returns:
    Series: Initial modeled absorption values.
    """
    
    # Unpacking ri_values
    ri_gfs_poa, ri_gfs_soa, ri_res_poa, ri_res_soa, ri_shp_poa, ri_shp_soa,\
    ri_trf_poa, ri_trf_soa, ri_oth_poa, ri_oth_soa = ri_values
    
    # Convert SPECIES to uppercase
    upper_species = [i.upper() for i in const.SPECIES]
    
    # Calculate optical properties
    optical_parameters = aao.calculate_optical_properties(const.RELATIVE_HUMIDITY, upper_species, const.WAVELENGTH, 
                                                      ri_gfs_poa=ri_gfs_poa,
                                                      ri_gfs_soa=ri_gfs_soa,
                                                      ri_res_poa=ri_res_poa,
                                                      ri_res_soa=ri_res_soa,
                                                      ri_shp_poa=ri_shp_poa,
                                                      ri_shp_soa=ri_shp_soa,
                                                      ri_trf_poa=ri_trf_poa,
                                                      ri_trf_soa=ri_trf_soa,
                                                      ri_oth_poa=ri_oth_poa,
                                                      ri_oth_soa=ri_oth_soa)
    
    # Extract concentration data for the station from the model
    model_conc = get_data_for_station(model, station, const.SPECIES)
    
    #passing absorption in Mm-1
    calc_absorption = calculate_absorption(model_conc*1e-6, optical_parameters)
    
    return calc_absorption*1e6

def calculate_abs_modeled_oa(station, model, ri_values):
    """
    Calculates the modeled absorption for a specific station based on provided refractive index values.

    Parameters:
    station (str): Station name used to filter concentration data.
    model (DataFrame): DataFrame containing model data.
    ri_values (list): List of refractive index values for different substances.

    Returns:
    Series: Initial modeled absorption values.
    """
    
    # Unpacking ri_values
    ri_oa = ri_values
    
    # Convert SPECIES to uppercase
    upper_species = [i.upper() for i in const.SPECIES_OA]
    
    # Calculate optical properties
    optical_parameters = aao.calculate_optical_properties4oa(const.RELATIVE_HUMIDITY, upper_species, const.WAVELENGTH, 
                                                      ri_oa=ri_oa)
    
    # Extract concentration data for the station from the model
    model_conc = get_data_for_station(model, station, const.SPECIES_OA)
    
    #passing absorption in Mm-1
    calc_absorption = calculate_absorption(model_conc*1e-6, optical_parameters)
    
    return calc_absorption*1e6

