"""
Utilities Module

This module contains various utility functions that support the main functionalities of 
the aerosol absorption model optimization. These include general-purpose functions for 
printing results, handling data conversions, and other supportive tasks.

Functions:
- print_table: Prints a table of refractive index values post-optimization.
"""
import os
import pandas as pd
import plotly.figure_factory as ff

from . import constants as const

def print_table(optimization_result, mode, method, **kwargs):
    try:
        # Prepare the DataFrame
        ri_values = optimization_result.x
        df_result = pd.DataFrame({
            'RI_name': const.SPECIES,  # This assumes SPECIES is a predefined list of species names
            'ri_values': ri_values
        })
        df_result['ri_values'] = df_result['ri_values'].round(4)

        # Determine the directory path
        model_name = kwargs.get("model", "default_model")
        directory_path = os.path.join('ri_tables', model_name, f'RI_{mode}')
        os.makedirs(directory_path, exist_ok=True)

        # Construct the file path based on the mode
        file_parts = ['RI', method]
        if mode == 'all':
            file_parts.append(kwargs.get("mass_data"))
        elif mode in ['by_category', 'by_season', 'by_station', 'by_station_season']:
            file_parts.extend(kwargs.get(part) for part in mode.split('_')[1:] if kwargs.get(part))
            file_parts.append(kwargs.get("mass_data"))
        
        file_name = '_'.join(filter(None, file_parts))
        file_path = os.path.join(directory_path, file_name)
        
        print(f"Saving table to {file_path}.png and {file_path}.csv")

        # Save the table as an image and CSV
        fig = ff.create_table(df_result)
        fig.update_layout(autosize=False, width=500, height=300)
        fig.write_image(f'{file_path}.png', scale=2)
        df_result.to_csv(f'{file_path}.csv', index=False)

        print("Files saved successfully.")
        return df_result
    except Exception as e:
        print(f"An error occurred: {e}")
        return None