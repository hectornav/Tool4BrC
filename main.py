import pandas as pd
import inspect
from scipy.optimize import Bounds
from modules import data_retrieval, data_processing, optimization, visualization, utils, constants, brown_carbon_ri_boundaries, ri_optimization_constraints


def optimize_ri(df_mod, df_obs, method, mode, bounds, constraints, initial_ri_values, mass_data, **kwargs):
    try:
        # Remove negative values from df_obs
        df_obs = df_obs[df_obs['AbsBrC370'] > 0]

        if mode == 'all':
            # Get a list of all station names
            stations = df_mod['station_name'].unique().tolist()

            result = optimization.optimize_stations(
                stations, df_mod, df_obs, method, bounds, constraints,
                initial_ri_values, by_season='no', model=kwargs['model']
            )
            utils.print_table(
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
                result = optimization.optimize_stations(
                    stations, df_mod, df_obs, method, bounds, constraints,
                    initial_ri_values, by_season='no', model=kwargs['model']
                )
                utils.print_table(
                    result, mode, method, category=category, mass_data=mass_data, model=kwargs['model']
                )
                print(f'Optimization for mode: {mode} and method: {method}: {category} finished')

        elif mode == 'by_season':
            # Get a list of all station names
            stations = df_mod['station_name'].unique().tolist()

            # Optimize for each season
            for season in ['DJF', 'MAM', 'JJA', 'SON']:
                result = optimization.optimize_stations(
                    stations, df_mod, df_obs, method, bounds, constraints,
                    initial_ri_values, by_season='yes', season=season, model=kwargs['model']
                )
                utils.print_table(
                    result, mode, method, season=season, mass_data=mass_data, model=kwargs['model']
                )
                print(f'Optimization for mode: {mode} and method: {method}, {season} finished')

        elif mode == 'by_station':
            # Get a list of all station names
            stations = df_mod['station_name'].unique().tolist()

            # Optimize for each station
            for station in stations:
                result = optimization.optimize_stations(
                    [station], df_mod, df_obs, method, bounds, constraints,
                    initial_ri_values, by_season='no', model=kwargs['model']
                )
                utils.print_table(
                    result, mode, method, station=station, mass_data=mass_data, model=kwargs['model']
                )
                print(f'Optimization for mode: {mode} and method: {method}, {station} finished')

        # Another by station and by season
        elif mode == 'by_station_season':
            stations = df_mod['station_name'].unique().tolist()
            for station in stations:
                for season in ['DJF', 'MAM', 'JJA', 'SON']:
                    result = optimization.optimize_stations(
                        [station], df_mod, df_obs, method, bounds, constraints,
                        initial_ri_values, by_season='yes', season=season, model=kwargs['model']
                    )
                    utils.print_table(
                        result, mode, method, station=station, season=season, mass_data=mass_data, model=kwargs['model']
                    )
                    print(f'Optimization for mode: {mode} and method: {method}, {station}, {season} finished')

    except Exception as e:
        print(f"Error: {e}")
        return None

def is_constraint_function(obj):
    return inspect.isfunction(obj) and obj.__name__.startswith("constraint")

if __name__ == '__main__':
    # Set bounds according to the RI bounds obtained for the wavelength (nm).
    brc_ri_bounds, tags = brown_carbon_ri_boundaries.get_ri_bounds(370)
    # save table of ri boundries
    visualization.plot_boundaries(brc_ri_bounds, tags)
    #para la ejecucion del codigo
  
  
    bounds = Bounds(
        [bound['start'] for bound in brc_ri_bounds.values()],
        [bound['end'] for bound in brc_ri_bounds.values()]
    )

    constraint_functions = inspect.getmembers(ri_optimization_constraints, is_constraint_function)
    #constraints = [
    #    {'type': 'ineq', 'fun': ri_optimization_constraints.constraint2},
    #]
    constraints = [{'type': 'ineq', 'fun': func} for _, func in constraint_functions]
    # Load the observed data from the CSV file.
    observed_data_clean = pd.read_csv('../absorption/NInventory/obs/absorption/absorption_brc370.csv')

    # Get a list of unique stations from the observed data.
    unique_stations = observed_data_clean[['station_name', 'category_station']].drop_duplicates()

    # Get a list of station names from the unique stations dataframe.
    list_of_stations = unique_stations['station_name'].unique()

    mass_data = 'best' #or 'all'
    monarch_data = data_retrieval.get_model_data_4monarch(mass_data=mass_data)
    methods = ['SLSQP'] #SLSQP, COBYLA, trust-constr
    cases = ['all'] #['all', 'by_category', 'by_season', 'by_station', 'by_station_season']
    for method in methods:
        for case in cases:
            optimize_ri(
                monarch_data, observed_data_clean, method, case,
                bounds, constraints, constants.INITIAL_RI_VALUES, mass_data=mass_data, model=f'monarch_{mass_data}'
            )