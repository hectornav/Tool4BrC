"""
Aerosol Optical Properties Calculator

This script calculates the physical and optical properties of various aerosol types,
including sulfate, brown carbon, and other specific categories. It utilizes empirical
formulas and parameters to estimate the refractive index, hygroscopic growth factors,
and other relevant physical properties of aerosols. These calculations are essential
for understanding the interaction of aerosols with light, which is crucial in climate
modeling and atmospheric science studies.

The script defines constants for hygroscopic growth based on established research
and provides functions to calculate physical properties for a range of aerosol types.
It allows for customization of properties for specific aerosol categories, particularly
useful for detailed atmospheric modeling.

Functions:
    calculate_physical_properties: Calculates properties for a range of aerosol types.
    calculate_brown_carbon_properties: Specifically calculates properties for Brown Carbon aerosols.
"""

import numpy as np

# Constants for SULFATE Hygroscopic Growth
# Fitzgerald [1975]; Haywood and Ramaswamy [1998]
SATURATION_LEVELS = np.array([0.80, 0.90, 0.95, 0.99])
F99 = 0.0155 * (SATURATION_LEVELS[3] - 0.97) / (1.02 - SATURATION_LEVELS[3]**1.4)
PHI = np.array([1.058, 1.058, 1.058, 1.058 - F99])
ALPHA = 1.2 * np.exp(0.066 * SATURATION_LEVELS / (PHI - SATURATION_LEVELS))
BETA = np.exp(0.00077 * SATURATION_LEVELS / (1.009 - SATURATION_LEVELS))
ALPHA_VALUES = ALPHA
BETA_VALUES = BETA

# Additional Constants
A30 = 1.0
slope = (ALPHA_VALUES[0] - A30) / (0.80 - 0.30)
A50 = A30 + slope * (0.50 - 0.30)
A70 = A30 + slope * (0.70 - 0.30)
B30, B50, B70 = 1.0, 1.0, 1.0

def calculate_physical_properties(ri_gfas_poa, ri_gfas_soa, ri_resi_poa, ri_resi_soa, ri_shp_poa,
                     ri_shp_soa, ri_trf_poa, ri_trf_soa, ri_othr_poa, ri_othr_soa):
        # rgeo, sdev, rmin, rmax, real, imag, dens, alfa[0%:99%], beta[0%:99%]
        phys = {'WATE':[1.336, 1.0e-08],
                'DUB1':[0.2986, 2.000, 0.100, 0.180, 1.531, 2.5e-03, 2.500, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'DUB2':[0.2986, 2.000, 0.180, 0.300, 1.531, 2.5e-03, 2.500, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'DUB3':[0.2986, 2.000, 0.300, 0.600, 1.531, 2.5e-03, 2.500, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'DUB4':[0.2986, 2.000, 0.600, 1.000, 1.531, 2.5e-03, 2.500, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'DUB5':[0.2986, 2.000, 1.000, 1.800, 1.531, 2.5e-03, 2.650, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'DUB6':[0.2986, 2.000, 1.800, 3.000, 1.531, 2.5e-03, 2.650, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'DUB7':[0.2986, 2.000, 3.000, 6.000, 1.531, 2.5e-03, 2.650, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'DUB8':[0.2986, 2.000, 6.000, 10.00, 1.531, 2.5e-03, 2.650, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'SSB1':[0.1500, 2.800, 0.100, 0.180, 1.557, 1.0e-08, 2.160, [1.0,1.6,1.8,2.0,2.4,2.9,4.8], [1.0 for h in np.arange(7)]],
                'SSB2':[0.1500, 2.800, 0.180, 0.300, 1.557, 1.0e-08, 2.160, [1.0,1.6,1.8,2.0,2.4,2.9,4.8], [1.0 for h in np.arange(7)]],
                'SSB3':[0.1500, 2.800, 0.300, 0.600, 1.557, 1.0e-08, 2.160, [1.0,1.6,1.8,2.0,2.4,2.9,4.8], [1.0 for h in np.arange(7)]],
                'SSB4':[0.1500, 2.800, 0.600, 1.000, 1.557, 1.0e-08, 2.160, [1.0,1.6,1.8,2.0,2.4,2.9,4.8], [1.0 for h in np.arange(7)]],
                'SSB5':[0.1500, 2.800, 1.000, 1.800, 1.557, 1.0e-08, 2.160, [1.0,1.6,1.8,2.0,2.4,2.9,4.8], [1.0 for h in np.arange(7)]],
                'SSB6':[0.1500, 2.800, 1.800, 3.000, 1.557, 1.0e-08, 2.160, [1.0,1.6,1.8,2.0,2.4,2.9,4.8], [1.0 for h in np.arange(7)]],
                'SSB7':[0.1500, 2.800, 3.000, 6.000, 1.557, 1.0e-08, 2.160, [1.0,1.6,1.8,2.0,2.4,2.9,4.8], [1.0 for h in np.arange(7)]],
                'SSB8':[0.1500, 2.800, 6.000, 15.00, 1.557, 1.0e-08, 2.160, [1.0,1.6,1.8,2.0,2.4,2.9,4.8], [1.0 for h in np.arange(7)]],
                'POA0':[0.0212, 2.200, 0.005, 20.00, 1.501, 1.0e-08, 1.800, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'POA1':[0.0212, 2.200, 0.005, 20.00, 1.501, 1.0e-04, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'POW0':[0.0212, 2.200, 0.005, 20.00, 1.501, 1.0e-01, 1.800, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'POW1':[0.0212, 2.200, 0.005, 20.00, 1.501, 1.0e-03, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'SOB1':[0.0212, 2.200, 0.005, 20.00, 1.486, 2.5e-05, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'SOB5':[0.0212, 2.200, 0.005, 20.00, 1.486, 2.5e-05, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'SOB6':[0.0212, 2.200, 0.005, 20.00, 1.486, 2.5e-05, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'SOB7':[0.0212, 2.200, 0.005, 20.00, 1.486, 2.5e-05, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'SBC0':[0.0118, 2.000, 0.005, 20.00, 1.850, 7.1e-01, 1.000, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'SBC1':[0.0118, 2.000, 0.005, 20.00, 1.850, 7.1e-01, 1.000, [1.0,1.0,1.0,1.2,1.4,1.5,1.9], [1.0 for h in np.arange(7)]],
                'SUL1':[0.0695, 2.030, 0.005, 20.00, 1.546, 1.0e-08, 1.700, [A30,A50,A70,A80,A90,A95,A99], [B30,B50,B70,B80,B90,B95,B99]],
                'DUPM10':[0.2986, 2.000, 0.100, 5.00, 1.531, 2.5e-03, 2.650, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'SSPM10':[0.1500, 2.800, 0.100, 5.00, 1.557, 1.0e-08, 2.160, [1.0,1.6,1.8,2.0,2.4,2.9,4.8], [1.0 for h in np.arange(7)]],
                'EPOA':[0.0212, 2.200, 0.005, 20.00, 1.501, 1.0e-08, 1.800, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'OPOA':[0.0212, 2.200, 0.005, 20.00, 1.501, 1.0e-04, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'ASOA':[0.0212, 2.200, 0.005, 20.00, 1.486, 2.5e-05, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'TSOA':[0.0212, 2.200, 0.005, 20.00, 1.486, 2.5e-05, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'ISOA':[0.0212, 2.200, 0.005, 20.00, 1.486, 2.5e-05, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'FSOA':[0.0212, 2.200, 0.005, 20.00, 1.486, 2.5e-05, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'NIPM10':[0.0695, 2.030, 0.005, 20.00, 1.546, 1.0e-08, 1.700, [A30,A50,A70,A80,A90,A95,A99], [B30,B50,B70,B80,B90,B95,B99]],
                'NH1':[0.0695, 2.030, 0.005, 20.00, 1.546, 1.0e-08, 1.700, [A30,A50,A70,A80,A90,A95,A99], [B30,B50,B70,B80,B90,B95,B99]],
                'SUL1':[0.0695, 2.030, 0.005, 20.00, 1.546, 1.0e-08, 1.700, [A30,A50,A70,A80,A90,A95,A99], [B30,B50,B70,B80,B90,B95,B99]],
                'UNSPM10':[0.2986, 2.000, 0.100, 5.00, 1.531, 2.5e-03, 2.650, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'POAGFS':[0.0212, 2.200, 0.005, 20.00, 1.501, ri_gfas_poa, 1.800, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'SOAGFS':[0.0212, 2.200, 0.005, 20.00, 1.486, ri_gfas_soa, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'POARES':[0.0212, 2.200, 0.005, 20.00, 1.501, ri_resi_poa, 1.800, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'SOARES':[0.0212, 2.200, 0.005, 20.00, 1.486, ri_resi_soa, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'POASHP':[0.0212, 2.200, 0.005, 20.00, 1.501, ri_shp_poa, 1.800, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'SOASHP':[0.0212, 2.200, 0.005, 20.00, 1.486, ri_shp_soa, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'POATRF':[0.0212, 2.200, 0.005, 20.00, 1.501, ri_trf_poa, 1.800, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'SOATRF':[0.0212, 2.200, 0.005, 20.00, 1.486, ri_trf_soa, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                'POAOTH':[0.0212, 2.200, 0.005, 20.00, 1.501, ri_othr_poa, 1.800, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'SOAOTH':[0.0212, 2.200, 0.005, 20.00, 1.486, ri_othr_soa, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]],
                
                }
        return phys

def calculate_brown_carbon_properties():
        #crear el BrC
        phys = {
                # rgeo, sdev, rmin, rmax, real, imag, dens, alfa[0%:99%], beta[0%:99%]
                'WATE':[1.336, 1.0e-08],
                'POABRC':[0.0212, 2.200, 0.005, 20.00, 1.501, 0.122, 1.800, [1.0 for h in np.arange(7)],  [1.0 for h in np.arange(7)]],
                'SOABRC':[0.0212, 2.200, 0.005, 20.00, 1.486, 0.06, 1.800, [1.0,1.2,1.4,1.5,1.6,1.8,2.2], [1.0 for h in np.arange(7)]]
        }
        return phys