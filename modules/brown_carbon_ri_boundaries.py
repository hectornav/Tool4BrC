"""
Brown Carbon Refractive Index Boundary Calculator

This script calculates the refractive index bounds for various types of aerosols 
using the Saleh methodology. It is designed to provide accurate refractive index 
estimations for different aerosol sources like GFAS (Global Fire Assimilation System), 
residential, shipping, traffic, and others. 

The script defines parameter ranges for absorption coefficients and wavelength 
exponents, and applies these to calculate the refractive index bounds at a given 
wavelength. The absorption properties are determined based on the aerosol type 
and source, with the order of absorption being GFAS > residential > shipping > 
traffic > other.

Functions:
    get_ri_bounds: Calculates and returns refractive index bounds for a specified wavelength.
"""
import numpy as np

params_saleh = {
    'vw_brc_k': {'start': 1e-4, 'end': 1e-3},
    'w_brc_k': {'start': 1e-3, 'end': 1e-2},
    'm_brc_k': {'start': 1e-2, 'end': 1e-1},
    's_brc_k': {'start': 1e-1, 'end': 0.38}, #upper limit should be revised
    'vw_brc_w': {'start': 6, 'end': 9},
    'w_brc_w': {'start': 4, 'end': 7},
    'm_brc_w': {'start': 1.5, 'end': 4},
    's_brc_w': {'start': 0.5, 'end': 1.5}
}

dict_of_params = {
    's': {'k': params_saleh['s_brc_k'], 'w': params_saleh['s_brc_w']},
    'm': {'k': params_saleh['m_brc_k'], 'w': params_saleh['m_brc_w']},
    'w':{'k': params_saleh['w_brc_k'], 'w': params_saleh['w_brc_w']},
    'vw': {'k': params_saleh['vw_brc_k'], 'w': params_saleh['vw_brc_w']},
}

names = {
    's': 'strongly',
    'm': 'moderately',
    'w': 'weakly',
    'vw': 'very weakly',
    }


def get_ri_bounds(wavelength, **kwargs):
    """
    Calculate refraction index bounds for different aerosol types and sources.

    Function uses Saleh methodology to create bounds for different aerosol types
    (Primary Organic Aerosol - POA, Secondary Organic Aerosol - SOA) and different
    sources (GFAS, residential, shipping, traffic, other). These bounds are based
    on the absorption properties of each aerosol type and source. The absorption
    order is: GFAS > residential > shipping > traffic > other.

    Parameters
    ----------
    wavelength : float
        The wavelength for which to calculate the refraction index bounds.
        This parameter should be in the nm units.

    Returns
    -------
    dict
        A dictionary containing the calculated refraction index bounds for each
        aerosol type and source.
    """
    
    def k_lambda(k, w, wavelength):
        """Calculates absorption of light at specific wavelength based on Saleh methodology."""
        return np.round(k * ((550/wavelength)**w), 4)

    def _create_ri_bounds_entry(k_start, k_end, w_start, w_end, wavelength):
        """
        Create a refraction index bounds entry for a given wavelength and parameters.

        Parameters
        ----------
        k_start : float
            The starting absorption coefficient.
        k_end : float
            The ending absorption coefficient.
        w_start : float
            The starting wavelength exponent.
        w_end : float
            The ending wavelength exponent.
        wavelength : float
            The wavelength for which to calculate the refraction index bounds.

        Returns
        -------
        dict
            A dictionary containing the calculated refraction index bounds.
        """
        return {'start': k_lambda(k_start, w_start, wavelength), 'end': k_lambda(k_end, w_end, wavelength)}

        """    ri_bounds = {
        'ri_gfs_poa': _create_ri_bounds_entry(params_saleh['m_brc_k']['start'], params_saleh['m_brc_k']['end'], params_saleh['m_brc_w']['start'], params_saleh['m_brc_w']['end'], wavelength),
        'ri_gfs_soa': _create_ri_bounds_entry(params_saleh['m_brc_k']['start'], params_saleh['m_brc_k']['end'], params_saleh['m_brc_w']['start'], params_saleh['m_brc_w']['end'], wavelength),
        'ri_res_poa': _create_ri_bounds_entry(params_saleh['m_brc_k']['start'], params_saleh['m_brc_k']['end'], params_saleh['m_brc_w']['start'], params_saleh['m_brc_w']['end'], wavelength),
        'ri_res_soa': _create_ri_bounds_entry(params_saleh['w_brc_k']['start'], params_saleh['w_brc_k']['end'], params_saleh['w_brc_w']['start'], params_saleh['w_brc_w']['end'], wavelength),
        'ri_shp_poa': _create_ri_bounds_entry(params_saleh['w_brc_k']['start'], params_saleh['w_brc_k']['end'], params_saleh['w_brc_w']['start'], params_saleh['w_brc_w']['end'], wavelength),
        'ri_shp_soa': _create_ri_bounds_entry(params_saleh['vw_brc_k']['start'], params_saleh['vw_brc_k']['end'], params_saleh['vw_brc_w']['start'], params_saleh['vw_brc_w']['end'], wavelength),
        'ri_trf_poa': _create_ri_bounds_entry(params_saleh['vw_brc_k']['start'], params_saleh['vw_brc_k']['end'], params_saleh['vw_brc_w']['start'], params_saleh['vw_brc_w']['end'], wavelength),
        'ri_trf_soa': _create_ri_bounds_entry(params_saleh['vw_brc_k']['start'], params_saleh['vw_brc_k']['end'], params_saleh['vw_brc_w']['start'], params_saleh['vw_brc_w']['end'], wavelength),
        'ri_oth_poa': _create_ri_bounds_entry(params_saleh['vw_brc_k']['start'], params_saleh['vw_brc_k']['end'], params_saleh['vw_brc_w']['start'], params_saleh['vw_brc_w']['end'], wavelength),
        'ri_oth_soa': _create_ri_bounds_entry(params_saleh['vw_brc_k']['start'], params_saleh['vw_brc_k']['end'], params_saleh['vw_brc_w']['start'], params_saleh['vw_brc_w']['end'], wavelength),
        }"""

    ri_bounds = {
        'ri_gfs_poa': _create_ri_bounds_entry(dict_of_params[kwargs.get('poa_gfas')]['k']['start'], dict_of_params[kwargs.get('poa_gfas')]['k']['end'],
                                                dict_of_params[kwargs.get('poa_gfas')]['w']['start'], dict_of_params[kwargs.get('poa_gfas')]['w']['end'], wavelength),
                                                
        'ri_gfs_soa': _create_ri_bounds_entry(dict_of_params[kwargs.get('soa_gfas')]['k']['start'], dict_of_params[kwargs.get('soa_gfas')]['k']['end'],
                                                dict_of_params[kwargs.get('soa_gfas')]['w']['start'], dict_of_params[kwargs.get('soa_gfas')]['w']['end'], wavelength),

        'ri_res_poa': _create_ri_bounds_entry(dict_of_params[kwargs.get('poa_res')]['k']['start'], dict_of_params[kwargs.get('poa_res')]['k']['end'],
                                                dict_of_params[kwargs.get('poa_res')]['w']['start'], dict_of_params[kwargs.get('poa_res')]['w']['end'], wavelength),
                                                
        'ri_res_soa': _create_ri_bounds_entry(dict_of_params[kwargs.get('soa_res')]['k']['start'], dict_of_params[kwargs.get('soa_res')]['k']['end'],
                                                dict_of_params[kwargs.get('soa_res')]['w']['start'], dict_of_params[kwargs.get('soa_res')]['w']['end'], wavelength),

        'ri_shp_poa': _create_ri_bounds_entry(dict_of_params[kwargs.get('poa_shp')]['k']['start'], dict_of_params[kwargs.get('poa_shp')]['k']['end'],
                                                dict_of_params[kwargs.get('poa_shp')]['w']['start'], dict_of_params[kwargs.get('poa_shp')]['w']['end'], wavelength),

        'ri_shp_soa': _create_ri_bounds_entry(dict_of_params[kwargs.get('soa_shp')]['k']['start'], dict_of_params[kwargs.get('soa_shp')]['k']['end'],
                                                dict_of_params[kwargs.get('soa_shp')]['w']['start'], dict_of_params[kwargs.get('soa_shp')]['w']['end'], wavelength),

        'ri_trf_poa': _create_ri_bounds_entry(dict_of_params[kwargs.get('poa_trf')]['k']['start'], dict_of_params[kwargs.get('poa_trf')]['k']['end'],
                                                dict_of_params[kwargs.get('poa_trf')]['w']['start'], dict_of_params[kwargs.get('poa_trf')]['w']['end'], wavelength),

        'ri_trf_soa': _create_ri_bounds_entry(dict_of_params[kwargs.get('soa_trf')]['k']['start'], dict_of_params[kwargs.get('soa_trf')]['k']['end'],
                                                dict_of_params[kwargs.get('soa_trf')]['w']['start'], dict_of_params[kwargs.get('soa_trf')]['w']['end'], wavelength),

        'ri_oth_poa': _create_ri_bounds_entry(dict_of_params[kwargs.get('poa_oth')]['k']['start'], dict_of_params[kwargs.get('poa_oth')]['k']['end'],
                                                dict_of_params[kwargs.get('poa_oth')]['w']['start'], dict_of_params[kwargs.get('poa_oth')]['w']['end'], wavelength),

        'ri_oth_soa': _create_ri_bounds_entry(dict_of_params[kwargs.get('soa_oth')]['k']['start'], dict_of_params[kwargs.get('soa_oth')]['k']['end'],
                                                dict_of_params[kwargs.get('soa_oth')]['w']['start'], dict_of_params[kwargs.get('soa_oth')]['w']['end'], wavelength),
    }

    saleh_class = {
        'ri_gfs_poa': names[kwargs.get('poa_gfas')],
        'ri_gfs_soa': names[kwargs.get('soa_gfas')],
        'ri_res_poa': names[kwargs.get('poa_res')],
        'ri_res_soa': names[kwargs.get('soa_res')],
        'ri_shp_poa': names[kwargs.get('poa_shp')],
        'ri_shp_soa': names[kwargs.get('soa_shp')],
        'ri_trf_poa': names[kwargs.get('poa_trf')],
        'ri_trf_soa': names[kwargs.get('soa_trf')],
        'ri_oth_poa': names[kwargs.get('poa_oth')],
        'ri_oth_soa': names[kwargs.get('soa_oth')],
    }

    return ri_bounds, saleh_class


if __name__ == '__main__':

    ri_bounds, tags = get_ri_bounds(370,
                                    poa_gfas='m', soa_gfas='m',
                                    poa_res='m', soa_res='w',
                                    poa_shp='w', soa_shp='vw',
                                    poa_trf='vw', soa_trf='vw',
                                    poa_oth='vw', soa_oth='vw')
    print(ri_bounds, tags)