"""
Constants Module

This module defines various constants used in the optimization of refractive index values 
for aerosol absorption models. These constants include specific aerosol species, relative 
humidity values, and wavelengths crucial for the model calculations.
"""
import numpy as np

SPECIES = [
    'poagfs', 'soagfs', 'poares', 'soares', 'poashp', 'soashp',
    'poatrf', 'soatrf', 'poaoth', 'soaoth'
]

RELATIVE_HUMIDITY = np.array([0.5])
WAVELENGTH = 370