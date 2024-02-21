
import pandas as pd
import inspect
from scipy.optimize import Bounds

from modules import data_retrieval, optimization, utils, constants

def optimize_ri(df_mod, df_obs, method, mode, bounds, constraints, initial_ri_values, mass_data, **kwargs):
    # Remove negative values from df_obs
    df_obs = df_obs[df_obs['AbsBrC370'] > 0]
    df_mod['oa'] = df_mod.filter(regex='^poa|^soa').sum(axis=1, skipna=False)
    #now remove poa and soa columns
    df_mod = df_mod.drop(columns=['poagfs','soagfs','poares','soares','poashp','soashp','poatrf','soatrf','poaoth','soaoth'])

    if mode == 'all':
        # Get a list of all station names
        stations = df_mod['station_name'].unique().tolist()
        #sum columns with start wit poa and soa and put in oa column mantain nan values as nan
        
        result = optimization.opt4oa(
            stations, df_mod, df_obs, method, bounds, constraints,
            initial_ri_values, by_season='no', model=kwargs['model']
        )
        utils.print_table4oa(
            result, mode, method, mass_data=mass_data, model=kwargs['model']
        )
        print(f'Optimization for mode: {mode} and method: {method} finished')
    
    elif mode == 'by_category':
        # Change categories
        df_mod['category_station'] = df_mod['category_station'].replace({
            'peri-urban background': 'suburban background',
            'rural background': 'regional background'
        })
        # Get a list of categories
        categories = df_mod['category_station'].unique().tolist()

        # Optimize for each category
        for category in categories:
            stations = df_mod[df_mod['category_station'] == category]['station_name'].unique().tolist()
            #sum columns with start wit poa and soa and put in oa column mantain nan values as nan
            df_mod['oa'] = df_mod.filter(regex='^poa|^soa').sum(axis=1, skipna=False)
            #now remove poa and soa columns
            df_mod = df_mod.drop(columns=['poagfs','soagfs','poares','soares','poashp','soashp','poatrf','soatrf','poaoth','soaoth'])

            result = optimization.opt4oa(
                stations, df_mod, df_obs, method, bounds, constraints,
                initial_ri_values, by_season='no', model=kwargs['model']
            )
            utils.print_table4oa(
                result, mode, method, category=category, mass_data=mass_data, model=kwargs['model']
            )
            print(f'Optimization for mode: {mode} and method: {method}: {category} finished')
    
    elif mode == 'by_station':
            # Get a list of all station names
            stations = df_mod['station_name'].unique().tolist()

            # Optimize for each station
            for station in stations:
                result = optimization.opt4oa(
                    [station], df_mod, df_obs, method, bounds, constraints,
                    initial_ri_values, by_season='no', model=kwargs['model']
                )
                utils.print_table4oa(
                    result, mode, method, station=station, mass_data=mass_data, model=kwargs['model']
                )
                print(f'Optimization for mode: {mode} and method: {method}, {station} finished')

    

if __name__ == '__main__':
    bounds = Bounds(0.001, 0.1)
    #set void constraints
    constraints = []
    # Get the list of stations
    mass_data = 'best' #or 'all'
    monarch_data = data_retrieval.get_model_data_4monarch(mass_data=mass_data)

    methods = ['SLSQP'] #SLSQP, COBYLA, trust-constr
    cases = ['all'] #all,  by_station
    observed_data_clean = pd.read_csv('../absorption/NInventory/obs/absorption/absorption_brc370.csv')

    # Get a list of unique stations from the observed data.
    unique_stations = observed_data_clean[['station_name', 'category_station']].drop_duplicates()

    # Get a list of station names from the unique stations dataframe.
    list_of_stations = unique_stations['station_name'].unique()
    scenario = 'oa'
    for method in methods:
        for case in cases:
            optimize_ri(
                monarch_data, observed_data_clean, method, case,
                bounds, constraints, constants.INITIAL_RI_VALUES_OA, mass_data=mass_data, model=f'monarch_{mass_data}_{scenario}'
            )

    print('Optimization finished')