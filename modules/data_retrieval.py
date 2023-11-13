"""
Data Retrieval Module

This module is responsible for fetching and preprocessing atmospheric aerosol data.
It includes functions to load model data from different sources and formats, facilitating
the analysis and optimization of aerosol absorption models in atmospheric studies.

Functions:
- get_model_data_4monarch: Retrieves and preprocesses data from MONARCH model.
- get_model_data_4emac: Fetches and processes EMAC model data.
"""
import pandas as pd
import numpy as np

# Load the observed data from the CSV file.
#observed_data_clean = pd.read_csv('../../absorption/NInventory/obs/absorption/absorption_brc370.csv')

# Get a list of unique stations from the observed data.
#unique_stations = observed_data_clean[['station_name', 'category_station']].drop_duplicates()

# Get a list of station names from the unique stations dataframe.
#list_of_stations = unique_stations['station_name'].unique()

def get_model_data_4monarch(mass_data='all'):
    """
    Read the model data and return a DataFrame with the model data.

    Additionally, add the station name and category to the DataFrame.

    Parameters:
    mass_data (str): Specifies the type of model data to retrieve.

    Returns:
    DataFrame: A DataFrame with the model data, station name, and category.
    """
    # Define the path based on mass_data.
    path = '../absorption/NInventory/mod/4brc/2018/'
    if mass_data == 'best':
        path = '../absorption/NInventory/mod/4brc/2018/BestModObs/'

    # Get the list of stations from the observed data.
    stations_data = pd.read_csv(
        '../absorption/NInventory/obs/absorption/absorption_brc370.csv',
        usecols=['station_name', 'category_station']
    )
    unique_stations = stations_data['station_name'].unique()

    # Read data for each station and compile into a list of dataframes.
    model_dataframes = []
    for station in unique_stations:
        # Read station data.
        station_data = pd.read_csv(
            f'{path}{mass_data}_{station}.csv' if mass_data == 'best' else f'{path}4brc_{station}.csv',
            index_col=0,
            parse_dates=True
        )

        # Add the time column from the index.
        station_data['time'] = station_data.index

        # Add category for the station.
        station_category = stations_data[
            stations_data['station_name'] == station
        ]['category_station'].unique()[0]
        station_data['category_station'] = station_category

        # Add station name.
        station_data['station_name'] = station

        model_dataframes.append(station_data)

    # Concatenate all dataframes into one.
    combined_model_data = pd.concat(model_dataframes, ignore_index=True)

    return combined_model_data


def get_model_data_4emac(mass_data='all'):
    """
    Read the EMAC model data and return a DataFrame with the model data.

    Additionally, add the station name and category to the DataFrame.

    Parameters:
    mass_data (str): Specifies the type of model data to retrieve, either 'all' or 'best'.

    Returns:
    DataFrame: A DataFrame with the model data, station names, and categories.
    """
    # Define the path based on mass_data.
    path = '../absorption/scpy_apprch/EMAC_DATA/data/'
    if mass_data == 'best':
        path = '../absorption/scpy_apprch/EMAC_DATA/data/BestModObs/'

    # Get the list of stations from the observed data.
    observed_stations_data = pd.read_csv(
        '../absorption/NInventory/obs/absorption/absorption_brc370.csv',
        usecols=['station_name', 'category_station']
    )
    # Define a fixed list of stations.
    stations = [
        'Montseny', 'Rigi', 'Athens_Demokritos', 'Ispra', 'OPE',
        'Payerne', 'Paris_SIRTA', 'Barcelona_PalauReial', 'Hyytiala', 'Helsinki'
    ]

    # Read data for each station and compile into a list of dataframes.
    model_dataframes = []
    for station in stations:
        # Read station data.
        station_data = pd.read_csv(
            f'{path}best_{station}.csv',
            index_col=0,
            parse_dates=True
        )

        # Add the time column from the index.
        station_data['time'] = station_data.index

        # Add category for the station.
        station_category = observed_stations_data[
            observed_stations_data['station_name'] == station
        ]['category_station'].unique()[0]
        station_data['category_station'] = station_category

        # Add station name.
        station_data['station_name'] = station

        model_dataframes.append(station_data)

    # Concatenate all dataframes into one.
    combined_model_data = pd.concat(model_dataframes, ignore_index=True)

    return combined_model_data

def get_obsabs370(station, **kwargs):
    """
    Get absorption from observation for a specific station.
    Parameters:
    station (str): Station name.
    **kwargs: Keyword arguments passed to calculate_optical_properties.
    Returns:
    absorption (float): Absorption for the given mass concentrations.
    """

    df_abs_obs = pd.read_csv('../absorption/NInventory/obs/absorption/absorption_brc370.csv')
    #get observation for specific station
    stn_o = df_abs_obs[df_abs_obs['station_name'] == station][['time', 'AbsBrC370']]
    #set index as time and datetime format
    stn_o = stn_o.set_index('time')
    stn_o.index = pd.to_datetime(stn_o.index)

    if kwargs.get('remove_negatives'):
        #set negative values to nan
        stn_o[stn_o < 0] = np.nan
        return stn_o
    else:
        return stn_o

if __name__ == '__main__':
    print(get_model_data_4emac('best'))