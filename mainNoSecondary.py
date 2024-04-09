import pandas as pd
import inspect
from scipy.optimize import Bounds

from modules import data_retrieval, optimization, utils, constants, \
                    brown_carbon_ri_boundaries, ri_constraints_nosecondary, visualization_boundaries


def optimize_ri(df_mod, df_obs, method, mode, bounds, constraints, initial_ri_values, mass_data, **kwargs):
    try:
        # Remove negative values from df_obs
        df_obs = df_obs[df_obs['AbsBrC370'] > 0]

        if mode == 'all':
            # Get a list of all station names
            stations = df_mod['station_name'].unique().tolist()
            result = optimization.optimizeStationsNoSecondary(
                stations, df_mod, df_obs, method, bounds, constraints,
                initial_ri_values, by_season='no', model=kwargs['model']
            )
            utils.TableNoSecondary(
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
                result = optimization.optimizeStationsNoSecondary(
                    stations, df_mod, df_obs, method, bounds, constraints,
                    initial_ri_values, by_season='no', model=kwargs['model']
                )
                utils.TableNoSecondary(
                    result, mode, method, category=category, mass_data=mass_data, model=kwargs['model']
                )
                print(f'Optimization for mode: {mode} and method: {method}: {category} finished')

        elif mode == 'by_season':
            # Get a list of all station names
            stations = df_mod['station_name'].unique().tolist()

            # Optimize for each season
            for season in ['DJF', 'MAM', 'JJA', 'SON']:
                result = optimization.optimizeStationsNoSecondary(
                    stations, df_mod, df_obs, method, bounds, constraints,
                    initial_ri_values, by_season='yes', season=season, model=kwargs['model']
                )
                utils.TableNoSecondary(
                    result, mode, method, season=season, mass_data=mass_data, model=kwargs['model']
                )
                print(f'Optimization for mode: {mode} and method: {method}, {season} finished')

        elif mode == 'by_station':
            # Get a list of all station names
            stations = df_mod['station_name'].unique().tolist()
            #stations = stations[6:]
            # Optimize for each station
            for station in stations:
                result = optimization.optimizeStationsNoSecondary(
                    [station], df_mod, df_obs, method, bounds, constraints,
                    initial_ri_values, by_season='no', model=kwargs['model']
                )
                utils.TableNoSecondary(
                    result, mode, method, station=station, mass_data=mass_data, model=kwargs['model']
                )
                print(f'Optimization for mode: {mode} and method: {method}, {station} finished')

        # Another by station and by season
        elif mode == 'by_station_season':
            stations = df_mod['station_name'].unique().tolist()
            for station in stations:
                for season in ['DJF', 'MAM', 'JJA', 'SON']:
                    result = optimization.optimizeStationsNoSecondary(
                        [station], df_mod, df_obs, method, bounds, constraints,
                        initial_ri_values, by_season='yes', season=season, model=kwargs['model']
                    )
                    utils.TableNoSecondary(
                        result, mode, method, station=station, season=season, mass_data=mass_data, model=kwargs['model']
                    )
                    print(f'Optimization for mode: {mode} and method: {method}, {station}, {season} finished')
        elif mode == 'by_emissioninventory_nsoa':
            result = optimization.optimize_stations4bcnandmsy_ns(
                        [kwargs.get('station')], df_mod, df_obs, method, bounds, constraints,
                        initial_ri_values, by_season='no',  model=kwargs['model']
                    )
            utils.TableNoSecondary(
                    result, mode, method, station=kwargs.get('station'), mass_data=mass_data, model=kwargs['model']
                )
            print(f'Optimization for mode: {mode} and method: {method}, {kwargs.get("station")}, finished')


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
        gfas = 's'
        resi = 'm'
        ship = 'm'
        traf = 'vw'
        othr = 'vw'

    #case: moderately absorbing
    elif scenario == 'moderately':
        gfas = 'm'
        resi = 'm'
        ship = 'm'
        traf = 'vw'
        othr = 'vw'

    #case: weakly absorbing
    elif scenario == 'weakly':
        gfas = 'w'
        resi = 'w'
        ship = 'w'
        traf = 'vw'
        othr = 'vw'
    #case: very weakly absorbing
    elif scenario == 'random':
        gfas = 'vw'
        resi = 'vw'
        ship = 'vw'
        traf = 'vw'
        othr = 'vw'

    return gfas, resi, ship, traf, othr


def is_constraint_function(obj):
    return inspect.isfunction(obj) and obj.__name__.startswith("constraint")


def reduceColumnNames(df):
    df['gfas'] = df['poagfs'] + df['soagfs']
    df['resi'] = df['poares'] + df['soares']
    df['ship'] = df['poashp'] + df['soashp']
    df['traf'] = df['poatrf'] + df['soatrf']
    df['othr'] = df['poaoth'] + df['soaoth']
    df.drop(columns=['poagfs', 'soagfs', 'poares', 'soares', 'poashp', 'soashp', 'poatrf', 'soatrf', 'poaoth', 'soaoth'], inplace=True)

 
    return df

if __name__ == '__main__':
    #...................................................
    #give values to bounds to the ri for all species
    #...................................................
    scenario = 'weakly' #or 'moderately' or 'strongly' or 'random
    gfas, resi, ship, traf, othr = cases_bounds(scenario)
    wavelength = 370
    brc_ri_bounds, tags = brown_carbon_ri_boundaries.getRiBoundsNoSecondary(wavelength=wavelength, 
                                                                   gfas=gfas, 
                                                                   resi=resi, 
                                                                   ship=ship, 
                                                                   traf=traf, 
                                                                   othr=othr)
    # save table of ri boundries
    visualization_boundaries.plotBoundariesNoSecondary(brc_ri_bounds, tags)
    print("Boundaries plotted")
    #................................................... 
    bounds = Bounds(
        [bound['start'] for bound in brc_ri_bounds.values()],
        [bound['end'] for bound in brc_ri_bounds.values()]
    )


    constraint_functions = inspect.getmembers(ri_constraints_nosecondary, is_constraint_function)
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
    methods = ['SLSQP'] #SLSQP, COBYLA, trust-constr
    '''
    monarch_data = data_retrieval.get_model_data_4monarch(mass_data=mass_data)
    #monarch_data.to_csv('model_mass_data.csv')
    monarch_data = reduceColumnNames(monarch_data)
    
    cases = ['by_station', 'all'] #['all', 'by_category', 'by_season', 'by_station', 'by_station_season']
    for method in methods:
        for case in cases:
            optimize_ri(
                monarch_data, observed_data_clean, method, case,
                bounds, constraints, constants.INITIAL_RI_VALUES_NOSECONDARY, mass_data=mass_data, model=f'monarch_{mass_data}_{scenario}'
            )
    '''
    ### Calculating ri optimized for the Hermees emission inventory
    bcn_data = pd.read_csv('/home/hnavarro/Desktop/PHD_BSC/GIT/absorption/NInventory/mod/4brc/2018/BestModObs4bcnandmsy/best_Barcelona_PalauReial.csv')
    bcn_data = reduceColumnNames(bcn_data)
    msy_data = pd.read_csv('/home/hnavarro/Desktop/PHD_BSC/GIT/absorption/NInventory/mod/4brc/2018/BestModObs4bcnandmsy/best_Montseny.csv')
    msy_data = reduceColumnNames(msy_data)
    cases = ['by_emissioninventory_nsoa'] #['all', 'by_category', 'by_season', 'by_station', 'by_station_season']
    for method in methods:
        for case in cases:
            optimize_ri(
                msy_data, observed_data_clean, method, case,
                bounds, constraints, constants.INITIAL_RI_VALUES_NOSECONDARY, mass_data=mass_data, model=f'monarch_{mass_data}_{scenario}', 
                station='Montseny'
            )
    