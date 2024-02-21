"""
Optimization Module

This module encapsulates the core logic for the optimization of refractive index values 
in aerosol absorption models. It includes functions for calculating error, defining cost 
functions, and running the optimization process using various methods and constraints.

Functions:
- calculate_error: Computes the error between measured and calculated absorption.
- cost_function: Defines the cost function for optimization.
- optimize_stations: Performs optimization for a list of stations.
- optimize_ri: Main function to initiate the optimization process for different cases.
"""
import pandas as pd
import numpy as np
import os
import time
from scipy.optimize import minimize

#for pso
#from pyswarm import pso

from . import sync_obs_model_data as sync
from . import data_processing as dp
from . import seasonal_data_grouper as sbys

#hide warning
import warnings
warnings.filterwarnings("ignore")

def calculate_error(measured_abs, calc_absorption):
    """
    Calcula el error entre la absorción medida y la calculada.

    Parámetros:
    measured_abs (DataFrame): DataFrame de absorción medida.
    calc_absorption (DataFrame): DataFrame de absorción calculada.

    Devoluciones:
    float: Error entre la absorción medida y la calculada.
    """
    # Convertir índice a datetime
    measured_abs.index = pd.to_datetime(measured_abs.index)
    calc_absorption.index = pd.to_datetime(calc_absorption.index)
    obs, mod = sync.sync_obs_with_model(measured_abs, calc_absorption)
    #print(len(obs), len(mod))
    # Calcular las diferencias absolutas
    differences = np.abs(obs - mod)

    # Calcular el error absoluto medio
    error = differences.mean()
    
    #print(f"Error: {error}")
    return error

def count_values(mod, obs):
    """
    Counts the number of values in the model and observation dataframes.

    Parameters:
    mod (DataFrame): The model data.
    obs (DataFrame): The observation data.

    Returns:
    int: The number of values in the model data.
    int: The number of values in the observation data.
    """
    # Count the number of values in the model and observation dataframes which are not NaN
    o,m = sync.sync_obs_with_model(obs, mod)
    mod_values = m.count().sum()
    obs_values = o.count().sum()
    
    if mod_values != obs_values:
        print(f"Warning: The number of values in the model and observation dataframes are not equal. {mod_values} != {obs_values}")
    else:
        #print(f"Number of values: {mod_values}")
        return mod_values


def cost_function(station, model, observation, ri_values):
    """
    Calculates the error between modelled and measured absorption for a specific station.

    Parameters:
    station (str): Station name used to filter concentration and absorption data.
    model (DataFrame): DataFrame containing model data.
    observation (DataFrame): DataFrame containing observation data.
    ri_values (list): List of refractive index values for different substances.

    Returns:
    float: Error between modelled and measured absorption.
    """
    calc_absorption = dp.calculate_abs_modeled(station, model, ri_values)
    
    measured_abs = dp.get_data_for_station(observation, station, ['AbsBrC370'])
    #conver index to datetime
    measured_abs.index = pd.to_datetime(measured_abs.index)

    error = calculate_error(measured_abs, calc_absorption)
    #print a df with number of values treated
    num_val = count_values(calc_absorption, measured_abs)
    df = pd.DataFrame({'station': station, 'num_val': num_val}, index=[0])
    #save df to csv in a folder called num_val with its station name
    if not os.path.exists('num_val'):
        os.makedirs('num_val')
    df.to_csv(f'num_val/{station}.csv', index=False)

    #rmse = calculate_rmse(measured_abs, calc_absorption*1e6)
    
    return error

def optimize_stations(stations, model_data, observed_data, method, bounds, constraints, initial_refractive_indices, **kwargs):
    """
    Perform optimization for stations.

    Parameters:
    stations (list): List of station names.
    model_data (DataFrame): The model data as a DataFrame.
    observed_data (DataFrame): The observed data as a DataFrame.
    method (str): The optimization method.
    bounds (Bounds): The optimization bounds.
    constraints (dict): The optimization constraints.
    initial_refractive_indices (list): The initial guess for refractive indices.

    Returns:
    Result: The result of the optimization process.
    """
    # consider only positive values from observed_data
    observed_data = observed_data[observed_data['AbsBrC370'] > 0]

    if kwargs.get('by_season') == 'yes':
        start_time = time.time()

        # Avoid stations with empty seasons in observed_data
        empty_stations = [
            station
            for station in stations
            if sbys.season(observed_data, station_name=station, kind='obs')[kwargs['season']].empty
        ]

        # Exclude stations that are in empty_stations
        stations = [station for station in stations if station not in empty_stations]

        # Define the objective function for optimization
        objective = lambda ri_values: sum(
            cost_function(
                station,
                sbys.season(model_data, station_name=station, kind='model')[kwargs['season']],
                sbys.season(observed_data, station_name=station, kind='obs')[kwargs['season']],
                ri_values
            ) for station in stations
        )

        result = minimize(objective, initial_refractive_indices, method=method, bounds=bounds, constraints=constraints)
        
        print(f"Optimization successful: {result.success}")
        print(f"Time: {time.time() - start_time}")
    
        return result
    else:
        max_attempts = 5  # Define el número máximo de intentos
        attempt = 0
        success = False

        while attempt < max_attempts and not success:
            start_time = time.time()
            # Esta es la función objetivo que se minimizará
            objective = lambda ri_values: sum(
                cost_function(
                    station,
                    model_data,
                    observed_data[observed_data['station_name'] == station],
                    ri_values
                ) for station in stations
            )

            result = minimize(objective, initial_refractive_indices, method=method, bounds=bounds, constraints=constraints)

            success = result.success
            #counting number of values in model and observation

            if success:
                print(f"Optimization successful: {result.success}")
            else:
                print(f"Optimization attempt {attempt + 1} failed")
                print(f"Optimization successful: {result.message}")


            print(f"Time: {time.time() - start_time}")
            attempt += 1

        return result     

def cost_function4oa(station, model, observation, ri_values):
    """
    Calculates the error between modelled and measured absorption for a specific station.

    Parameters:
    station (str): Station name used to filter concentration and absorption data.
    model (DataFrame): DataFrame containing model data.
    observation (DataFrame): DataFrame containing observation data.
    ri_values (list): List of refractive index values for different substances.

    Returns:
    float: Error between modelled and measured absorption.
    """
    calc_absorption = dp.calculate_abs_modeled_oa(station, model, ri_values)
    
    measured_abs = dp.get_data_for_station(observation, station, ['AbsBrC370'])
    #conver index to datetime
    measured_abs.index = pd.to_datetime(measured_abs.index)

    error = calculate_error(measured_abs, calc_absorption)
    #rmse = calculate_rmse(measured_abs, calc_absorption*1e6)
  
    return error

def opt4oa(stations, model_data, observed_data, method, bounds, constraints, initial_refractive_indices, **kwargs):
    """
    Perform optimization for stations for OA total.

    Parameters:
    stations (list): List of station names.
    model_data (DataFrame): The model data as a DataFrame.
    observed_data (DataFrame): The observed data as a DataFrame.
    method (str): The optimization method.
    bounds (Bounds): The optimization bounds.
    constraints (dict): The optimization constraints.
    initial_refractive_indices (list): The initial guess for refractive indices.

    Returns:
    Result: The result of the optimization process.
    """
    # Remove negative values from observed_data
    observed_data = observed_data[observed_data['AbsBrC370'] > 0]

    if kwargs.get('by_season') == 'yes':
        start_time = time.time()

        # Avoid stations with empty seasons in observed_data
        empty_stations = [
            station
            for station in stations
            if sbys.season(observed_data, station_name=station, kind='obs')[kwargs['season']].empty
        ]

        # Exclude stations that are in empty_stations
        stations = [station for station in stations if station not in empty_stations]

        # Define the objective function for optimization
        objective = lambda ri_values: sum(
            cost_function4oa(
                station,
                sbys.season(model_data, station_name=station, kind='model')[kwargs['season']],
                sbys.season(observed_data, station_name=station, kind='obs')[kwargs['season']],
                ri_values
            ) for station in stations
        )

        result = minimize(objective, initial_refractive_indices, method=method, bounds=bounds, constraints=constraints)
        
        print(f"Optimization successful: {result.success}")
        print(f"Time: {time.time() - start_time}")
    
        return result
    else:
        start_time = time.time()

        # This is the objective function that will be minimized
        objective = lambda ri_values: sum(
            cost_function4oa(
                station,
                model_data,
                observed_data[observed_data['station_name'] == station],
                ri_values
            ) for station in stations
        )

        result = minimize(objective, initial_refractive_indices, method=method, bounds=bounds, constraints=constraints)
        
        print(f"Optimization successful: {result.success}")
        print(f"Time: {time.time() - start_time}")

        return result