"""
Seasonal Data Grouper

This module is designed to group atmospheric or environmental data into distinct 
seasonal categories: DJF (December, January, February), MAM (March, April, May), 
JJA (June, July, August), and SON (September, October, November). It is particularly 
useful in climate studies and atmospheric science, where analyzing data on a seasonal 
basis is crucial.

The module provides functions to ensure data completeness for each month, and to group 
data according to the specified seasons. It handles both observed data and model output, 
allowing for flexible application in various research contexts.

Functions:
- ensure_month_in_data: Ensures that a given month is present in the dataset, filling missing months with NaNs.
- get_season_df: Retrieves a DataFrame filtered by specific months.
- season: Groups data into seasonal categories, handling different types of data structures.
"""
import pandas as pd
import numpy as np



def ensure_month_in_data(data, month, year):
    """
    Asegura que el mes dado esté presente en los datos. Si no está, agrega un mes lleno de NaN.
    """
    month_data = data[(data.index.month == month) & (data.index.year == year)]
    
    if month_data.empty:
        # Crea un rango de fechas para ese mes específico
        date_range = pd.date_range(start=f"{year}-{month}-01", end=f"{year}-{month}-{pd.Timestamp(year, month, 1).days_in_month}")
        
        # Crea un DataFrame con NaN para ese rango de fechas
        month_df = pd.DataFrame(index=date_range, columns=data.columns).assign(**{col: np.nan for col in data.columns})
        
        # Concatena el DataFrame original con el mes NaN
        data = pd.concat([data, month_df], axis=0).sort_index()
    
    return data

def get_season_df(data, months, year):
    """ 
    Devuelve un DataFrame filtrado por los meses específicos.
    """
    for month in months:
        data = ensure_month_in_data(data, month, year)
    return data[data.index.month.isin(months)]

def season(data, **kwargs):
    '''
    This method groups passed data into DJF, MAM, JJA, SON months according to seasons.

    Args:
        data (pd.DataFrame): Data to be grouped by season.

    Returns:
        dict: A dictionary containing data grouped by seasons.
    '''
    #if kwargs is present, then we need to threat data different
    if kwargs:
        if kwargs['kind']=='model':
            #selct just the station
            data = data[data.station_name == kwargs['station_name']]
            #remove station_name, and category_station as columns
            data = data.drop(columns=['station_name', 'category_station'])
            #set time column as index
            data = data.set_index('time')
            #convert index to datetime
            if not isinstance(data.index, pd.DatetimeIndex):
                data.index = pd.to_datetime(data.index)
            #get the season
            dec= data[data.index.month.isin([12])]
            jan= data[data.index.month.isin([1])]
            feb= data[data.index.month.isin([2])]
            df_djf = pd.concat([dec, jan, feb], ignore_index=False, sort=True)
            trim_dict ={
                #iloc to have DEC as first month
                'DJF': df_djf, 
                'MAM': data[data.index.month.isin([3,4,5])],
                'JJA': data[data.index.month.isin([6,7,8])],
                'SON': data[data.index.month.isin([9,10,11])]
            }
            return trim_dict
        elif kwargs['kind']=='obs':
            #select just the station
            data = data[data.station_name == kwargs['station_name']]
            #in the new data select just columns time and AbsBrC370
            data = data[['time', 'AbsBrC370']]
            #set time column as index
            data = data.set_index('time')
            #convert index to datetime
            if not isinstance(data.index, pd.DatetimeIndex):
                data.index = pd.to_datetime(data.index)
            #get the season
            dec= data[data.index.month.isin([12])]
            jan= data[data.index.month.isin([1])]
            feb= data[data.index.month.isin([2])]
            df_djf = pd.concat([dec, jan, feb], ignore_index=False, sort=True)
            trim_dict ={
                #iloc to have DEC as first month
                'DJF': df_djf, 
                'MAM': data[data.index.month.isin([3,4,5])],
                'JJA': data[data.index.month.isin([6,7,8])],
                'SON': data[data.index.month.isin([9,10,11])]
            }
            return trim_dict
    else:
        if not isinstance(data.index, pd.DatetimeIndex):
            data.index = pd.to_datetime(data.index)

        #get the months in the data
        months = data.index.month.unique()

        #create a list with all months
        all_months = np.arange(1,13)
        if not np.array_equal(months, all_months):
            #get the missing months
            missing_months = np.setdiff1d(all_months, months)
            #create a dataframe with the missing months at dayly frequency
            missing_months_df = pd.DataFrame(index=pd.date_range(start=data.index.min(), end=data.index.max(), freq='D'))
            #concatenate the missing months to the original data
            data = pd.concat([data, missing_months_df], ignore_index=False, sort=True)
            #sort the index
            data = data.sort_index()
            #fill the missing months with NaN
            data = data.fillna(np.nan)
            
        #._._._._._._._._._._._._._._._._._._._
        dec= data[data.index.month.isin([12])]
        jan= data[data.index.month.isin([1])]
        feb= data[data.index.month.isin([2])]
        df_djf = pd.concat([dec, jan, feb], ignore_index=False, sort=True)
        trim_dict ={
            #iloc to have DEC as first month
            'DJF': df_djf, 
            'MAM': data[data.index.month.isin([3,4,5])],
            'JJA': data[data.index.month.isin([6,7,8])],
            'SON': data[data.index.month.isin([9,10,11])]
        }

        return trim_dict