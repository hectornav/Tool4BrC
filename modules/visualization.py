"""
Visualization Module

This module provides functionality for visualizing the results of the aerosol absorption 
model optimization. It includes functions to plot and compare optimized refractive index 
values against initial values, enhancing the interpretability of the optimization results.

Functions:
- plot_optimized_vs_initial: Plots a comparison of optimized refractive indices vs. initial values.
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

from . import data_processing as dp


def plot_optimized_vs_initial(stations, observed_data, method, mode, initial_refractive_indices, optimization_result):
    """
    Plot the comparison of optimized refractive indices versus initial values for each station.
    
    Parameters:
    stations (list): List of station names.
    observed_data (DataFrame): The observed data as a DataFrame.
    method (str): The optimization method used.
    mode (str): The mode of optimization.
    initial_refractive_indices (list): The initial guess for refractive indices.
    optimization_result (Result): The result of the optimization process.
    """
    directory = "optimized_vs_initial"
    os.makedirs(directory, exist_ok=True)

    # Remove negative values from observed_data
    observed_data = observed_data[observed_data['AbsBrC370'] > 0]

    for station in stations:
        # Extract observed data with date index
        observed = observed_data[observed_data['station_name'] == station][['AbsBrC370', 'time']]
        observed.set_index('time', inplace=True)
        observed.index = pd.to_datetime(observed.index)

        # Extract initial modeled data (this might need adjustment depending on your model dataframe structure)
        initial_modeled = dp.calculate_abs_modeled(station, initial_refractive_indices)
        
        # Calculate optimized values (this will depend on how you get the optimized data using the optimized refractive indices)
        optimized_modeled = dp.calculate_abs_modeled(station, optimization_result.x)

        plt.figure(figsize=(10, 6))
        # Plot observed and modeled data
        plt.plot(observed.index, observed['AbsBrC370'], label='Observed')
        plt.plot(initial_modeled.index, initial_modeled['AbsBrC370'], label='Initial Modeled', marker='o')
        plt.plot(optimized_modeled.index, optimized_modeled['AbsBrC370'], label='Optimized')

        plt.title(f"Model vs Observed: {station}")
        plt.xlabel('Time')
        plt.ylabel('Absorption Coefficient (Mm$^{-1}$)')
        plt.legend()
        plt.savefig(os.path.join(directory, f"{station}-{method}-{mode}.png"))
        plt.close()
    