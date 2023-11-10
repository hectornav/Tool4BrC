"""
Brown Carbon Refractive Index Boundary Calculator

This script calculates the refractive index bounds for various types of aerosols 
using the Saleh methodology. It is designed to provide accurate refractive index 
estimations for different aerosol sources like GFAS (Global Fire Assimilation System), 
residential, shipping, traffic, and others. The script is especially useful in 
atmospheric science research where precise knowledge of aerosol optical properties 
is crucial.

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
    's_brc_k': {'start': 1e-1, 'end': 0.5},
    'vw_brc_w': {'start': 6, 'end': 9},
    'w_brc_w': {'start': 4, 'end': 7},
    'm_brc_w': {'start': 1.5, 'end': 4},
    's_brc_w': {'start': 0.5, 'end': 1.5}
}

def get_ri_bounds(wavelength):
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

    ri_bounds = {
        'ri_gfs_poa': _create_ri_bounds_entry(params_saleh['s_brc_k']['start'], params_saleh['s_brc_k']['end'], params_saleh['s_brc_w']['start'], params_saleh['s_brc_w']['end'], wavelength),
        'ri_gfs_soa': _create_ri_bounds_entry(params_saleh['m_brc_k']['start'], params_saleh['m_brc_k']['end'], params_saleh['m_brc_w']['start'], params_saleh['m_brc_w']['end'], wavelength),
        'ri_res_poa': _create_ri_bounds_entry(params_saleh['m_brc_k']['start'], params_saleh['m_brc_k']['end'], params_saleh['m_brc_w']['start'], params_saleh['m_brc_w']['end'], wavelength),
        'ri_res_soa': _create_ri_bounds_entry(params_saleh['w_brc_k']['start'], params_saleh['w_brc_k']['end'], params_saleh['w_brc_w']['start'], params_saleh['w_brc_w']['end'], wavelength),
        'ri_shp_poa': _create_ri_bounds_entry(params_saleh['w_brc_k']['start'], params_saleh['w_brc_k']['end'], params_saleh['w_brc_w']['start'], params_saleh['w_brc_w']['end'], wavelength),
        'ri_shp_soa': _create_ri_bounds_entry(params_saleh['vw_brc_k']['start'], params_saleh['vw_brc_k']['end'], params_saleh['w_brc_w']['start'], params_saleh['w_brc_w']['end'], wavelength),
        'ri_trf_poa': _create_ri_bounds_entry(params_saleh['vw_brc_k']['start'], params_saleh['vw_brc_k']['end'], params_saleh['w_brc_w']['start'], params_saleh['w_brc_w']['end'], wavelength),
        'ri_trf_soa': _create_ri_bounds_entry(params_saleh['vw_brc_k']['start'], params_saleh['vw_brc_k']['end'], params_saleh['w_brc_w']['start'], params_saleh['w_brc_w']['end'], wavelength),
        'ri_oth_poa': _create_ri_bounds_entry(params_saleh['vw_brc_k']['start'], params_saleh['vw_brc_k']['end'], params_saleh['vw_brc_w']['start'], params_saleh['vw_brc_w']['end'], wavelength),
        'ri_oth_soa': _create_ri_bounds_entry(params_saleh['vw_brc_k']['start'], params_saleh['vw_brc_k']['end'], params_saleh['vw_brc_w']['start'], params_saleh['vw_brc_w']['end'], wavelength),
    }

    return ri_bounds


if __name__ == '__main__':
    ri_bounds = get_ri_bounds(370)
    print(ri_bounds)