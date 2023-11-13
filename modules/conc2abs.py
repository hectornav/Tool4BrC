import pandas as pd 
import numpy as np
import os

from . import constants
from . import atmospheric_aerosol_optics as aao
from . import aerosol_absorption_calculator as aac


def get_absorption(massconc, WAVELENGTH=constants.WAVELENGTH, REL_HUM=constants.RELATIVE_HUMIDITY, **kwargs):
    """
    Calculates the absorption for a given mass concentration.

    Parameters:
    mass_concentrations (dict): Dictionary containing mass concentrations for each species should be in ug/m3.
    WAVELENGTH (int): Wavelength used to calculate the absorption.
    REL_HUM (array): Relative humidity used to calculate the absorption.
    **kwargs: Keyword arguments passed to calculate_optical_properties.
    mass_mode (str): Mass mode for model (all or best)

    Returns:
    absorption (float): Absorption for the given mass concentrations.
    """
    if kwargs.get('mode') == 'all':
        ri_cat = pd.read_csv(f'ri_tables/{kwargs.get("model")}/RI_{kwargs.get("mode")}/RI_{kwargs.get("method")}_{kwargs.get("mass_mode")}.csv', index_col=0)
    elif kwargs.get('mode') == 'by_category':
        ri_cat = pd.read_csv(f'ri_tables/{kwargs.get("model")}/RI_{kwargs.get("mode")}/RI_{kwargs.get("method")}_{kwargs.get("mass_mode")}_{kwargs.get("category")}.csv', index_col=0)
    elif kwargs.get('mode') == 'by_season':
        ri_cat = pd.read_csv(f'ri_tables/{kwargs.get("model")}/RI_{kwargs.get("mode")}/RI_{kwargs.get("method")}_{kwargs.get("mass_mode")}_{kwargs.get("season")}.csv', index_col=0)
    elif kwargs.get('mode') == 'by_station':
        ri_cat = pd.read_csv(f'ri_tables/{kwargs.get("model")}/RI_{kwargs.get("mode")}/RI_{kwargs.get("method")}_{kwargs.get("station")}_{kwargs.get("mass_mode")}.csv', index_col=0)
    elif kwargs.get('mode') == 'apriori':
        ri_cat = pd.DataFrame(constants.INITIAL_RI_VALUES, index=constants.SPECIES, columns=['ri'])
    
    upper_species = [i.upper() for i in constants.SPECIES]


    #calculate optical parameters
    optical_parameters = aao.calculate_optical_properties(
                        REL_HUM,
                        upper_species,
                        WAVELENGTH,
                        ri_gfs_poa=ri_cat.loc['poagfs'].values[0],
                        ri_gfs_soa=ri_cat.loc['soagfs'].values[0],
                        ri_res_poa=ri_cat.loc['poares'].values[0],
                        ri_res_soa=ri_cat.loc['soares'].values[0],
                        ri_shp_poa=ri_cat.loc['poashp'].values[0],
                        ri_shp_soa=ri_cat.loc['soashp'].values[0],
                        ri_trf_poa=ri_cat.loc['poatrf'].values[0],
                        ri_trf_soa=ri_cat.loc['soatrf'].values[0],
                        ri_oth_poa=ri_cat.loc['poaoth'].values[0],
                        ri_oth_soa=ri_cat.loc['soaoth'].values[0],
                    )
    if kwargs.get('SA'):
        calculated_absorption = aac.calculate_absorption(
                                                        massconc*1e-6, #appling 1e-6 to convert from ug/m³ to g/m³
                                                        optical_parameters,
                                                        ) 
        columns = massconc.columns
        final_absorption = pd.DataFrame(calculated_absorption, columns=columns)
        #return absorption in Mm-1
        return final_absorption*1e6
        
    else:
        calculated_absorption = aac.calculate_absorption(
                                                        massconc*1e-6, #appling 1e-6 to convert from ug/m³ to g/m³
                                                        optical_parameters,
                                                        ).sum(axis=1)
        columns = ['AbsBrC370']
        final_absorption = pd.DataFrame(calculated_absorption, columns=columns)
        #return absorption in Mm-1
        return final_absorption*1e6

def get_absorption_by_stn(massconc, WAVELENGTH=constants.WAVELENGTH, REL_HUM=constants.RELATIVE_HUMIDITY, **kwargs):
    """
    Calculates the absorption for a given mass concentration.

    Parameters:
    mass_concentrations (dict): Dictionary containing mass concentrations for each species should be in ug/m3.
    WAVELENGTH (int): Wavelength used to calculate the absorption.
    REL_HUM (array): Relative humidity used to calculate the absorption.
    **kwargs: Keyword arguments passed to calculate_optical_properties.
    mass_mode (str): Mass mode for model (all or best)

    Returns:
    absorption (float): Absorption for the given mass concentrations.
    """
    if kwargs.get('mode') == 'all':
        ri_cat = pd.read_csv(f'ri_tables/{kwargs.get("model")}/RI_{kwargs.get("mode")}/RI_{kwargs.get("method")}_{kwargs.get("mass_mode")}.csv', index_col=0)
    elif kwargs.get('mode') == 'by_category':
        ri_cat = pd.read_csv(f'ri_tables/{kwargs.get("model")}/RI_{kwargs.get("mode")}/RI_{kwargs.get("method")}_{kwargs.get("mass_mode")}_{kwargs.get("category")}.csv', index_col=0)
    elif kwargs.get('mode') == 'by_season':
        ri_cat = pd.read_csv(f'ri_tables/{kwargs.get("model")}/RI_{kwargs.get("mode")}/RI_{kwargs.get("method")}_{kwargs.get("mass_mode")}_{kwargs.get("season")}.csv', index_col=0)
    elif kwargs.get('mode') == 'by_station':
        ri_cat = pd.read_csv(f'ri_tables/{kwargs.get("model")}/RI_{kwargs.get("mode")}/RI_{kwargs.get("method")}_{kwargs.get("mass_mode")}_{kwargs.get("station")}.csv', index_col=0)
    elif kwargs.get('mode') == 'apriori':
        ri_cat = pd.DataFrame(constants.INITIAL_RI_VALUES, index=constants.SPECIES, columns=['ri'])
    
    upper_species = [i.upper() for i in constants.SPECIES]

    #calculate optical parameters
    optical_parameters = aao.calculate_optical_properties(
                        REL_HUM,
                        upper_species,
                        WAVELENGTH,
                        ri_gfs_poa=ri_cat.loc['poagfs'].values[0],
                        ri_gfs_soa=ri_cat.loc['soagfs'].values[0],
                        ri_res_poa=ri_cat.loc['poares'].values[0],
                        ri_res_soa=ri_cat.loc['soares'].values[0],
                        ri_shp_poa=ri_cat.loc['poashp'].values[0],
                        ri_shp_soa=ri_cat.loc['soashp'].values[0],
                        ri_trf_poa=ri_cat.loc['poatrf'].values[0],
                        ri_trf_soa=ri_cat.loc['soatrf'].values[0],
                        ri_oth_poa=ri_cat.loc['poaoth'].values[0],
                        ri_oth_soa=ri_cat.loc['soaoth'].values[0],
                    )
    if kwargs.get('SA'):
        calculated_absorption = aac.calculate_absorption(
                                                        massconc*1e-6, #appling 1e-6 to convert from ug/m³ to g/m³
                                                        optical_parameters,
                                                        ) 
        columns = massconc.columns

    else:
        calculated_absorption = aac.calculate_absorption(
                                                        massconc*1e-6, #appling 1e-6 to convert from ug/m³ to g/m³
                                                        optical_parameters,
                                                        ) 
        print(calculated_absorption)
        columns = ['AbsBrC370'] 

    final_absorption = pd.DataFrame(calculated_absorption, columns=columns)
    #return absorption in Mm-1
    return final_absorption*1e6