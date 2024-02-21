import pandas as pd
import inspect
from scipy.optimize import Bounds

from modules import data_retrieval, optimization, utils, constants, brown_carbon_ri_boundaries, ri_optimization_constraints, visualization_boundaries


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
            #stations = stations[6:]
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
    
def cases_bounds(scenario):
    #...................................................
    #give values to bounds to the ri for all species
    #...................................................
    # m : moderately
    # w : weakly
    # vw: very weakly
    # s: strongly
    #...................................................
    #case: strongly absorbing
    if scenario == 'strongly':
        poa_gfas = 's'
        soa_gfas = 'm'
        poa_res = 'm'
        soa_res = 'w'
        poa_shp = 'm'
        soa_shp = 'w'
        poa_trf = 'vw'
        soa_trf = 'vw'
        poa_oth = 'vw'
        soa_oth = 'vw'
    #case: moderately absorbing
    elif scenario == 'moderately':
        poa_gfas = 'm'
        soa_gfas = 'w'
        poa_res = 'm'
        soa_res = 'w'
        poa_shp = 'm'
        soa_shp = 'w'
        poa_trf = 'vw'
        soa_trf = 'vw'
        poa_oth = 'vw'
        soa_oth = 'vw'
    #case: weakly absorbing
    elif scenario == 'weakly':
        poa_gfas = 'w'
        soa_gfas = 'vw'
        poa_res = 'w'
        soa_res = 'vw'
        poa_shp = 'w'
        soa_shp = 'vw'
        poa_trf = 'vw'
        soa_trf = 'vw'
        poa_oth = 'vw'
        soa_oth = 'vw'
    #case: very weakly absorbing
    elif scenario == 'random':
        poa_gfas = 'vw'
        soa_gfas = 'vw'
        poa_res = 'vw'
        soa_res = 'vw'
        poa_shp = 'vw'
        soa_shp = 'vw'
        poa_trf = 'vw'
        soa_trf = 'vw'
        poa_oth = 'vw'
        soa_oth = 'vw'
    return poa_gfas, soa_gfas, poa_res, soa_res, poa_shp, soa_shp, poa_trf, soa_trf, poa_oth, soa_oth


def is_constraint_function(obj):
    return inspect.isfunction(obj) and obj.__name__.startswith("constraint")

if __name__ == '__main__':
    #...................................................
    #give values to bounds to the ri for all species
    #...................................................
    scenario = 'strongly' #or 'moderately' or 'strongly' or 'random
    poa_gfas, soa_gfas, poa_res, soa_res, poa_shp, soa_shp, poa_trf, soa_trf, poa_oth, soa_oth = cases_bounds(scenario)
    wavelength = 370
    brc_ri_bounds, tags = brown_carbon_ri_boundaries.get_ri_bounds(wavelength=wavelength, 
                                                                   poa_gfas=poa_gfas, 
                                                                   soa_gfas=soa_gfas, 
                                                                   poa_res=poa_res, 
                                                                   soa_res=soa_res, 
                                                                   poa_shp=poa_shp, 
                                                                   soa_shp=soa_shp, 
                                                                   poa_trf=poa_trf, 
                                                                   soa_trf=soa_trf, 
                                                                   poa_oth=poa_oth, 
                                                                   soa_oth=soa_oth)
    # save table of ri boundries
    visualization_boundaries.plot_boundaries(brc_ri_bounds, tags)
    print("Boundaries plotted")
    #................................................... 
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
    #monarch_data.to_csv('model_mass_data.csv')
    methods = ['SLSQP'] #SLSQP, COBYLA, trust-constr
    cases = ['by_station', 'all'] #['all', 'by_category', 'by_season', 'by_station', 'by_station_season']
    for method in methods:
        for case in cases:
            optimize_ri(
                monarch_data, observed_data_clean, method, case,
                bounds, constraints, constants.INITIAL_RI_VALUES, mass_data=mass_data, model=f'monarch_{mass_data}_{scenario}'
            )