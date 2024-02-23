"""
Constants Module

This module defines various constants used in the optimization of refractive index values 
for aerosol absorption. These constants include specific aerosol species, relative 
humidity values, and wavelengths crucial for the model calculations.
"""
import numpy as np

SPECIES = [
    'poagfs', 'soagfs', 'poares', 'soares', 'poashp', 'soashp',
    'poatrf', 'soatrf', 'poaoth', 'soaoth'
]

SPECIESNOSECONDARY = ['gfas', 'resi', 'ship', 'traf', 'othr']

RELATIVE_HUMIDITY = np.array([0.5])
WAVELENGTH = 370

#initial_ri_values = [12e-2, 6e-2, 1e-1, 43e-3, 1e-1, 1e-3, 1e-2, 1e-3, 1e-2, 1e-3]
INITIAL_RI_VALUES = [0.12,0.06,0.10,0.043,0.10,0.001,0.01,0.001,0.01,0.001]
INITIAL_RI_VALUES_NOSECONDARY = [0.12, #GFAS
                                 0.1, #RESI
                                 0.1, #SHIP
                                 0.01, #TRAF
                                 0.01 #OTHR
                                 ]
#give random values

INITIAL_RI_VALUES_OA = [0.12]
SPECIES_OA = ['oa']