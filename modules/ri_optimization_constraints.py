"""
Refractive Index Optimization Constraints

This module defines a set of constraints used in the optimization of refractive index (RI) values 
for atmospheric aerosols. These constraints ensure that the optimized RI values adhere to 
expected physical and empirical boundaries. They are crucial for maintaining the physical 
plausibility and accuracy of aerosol models in atmospheric science.

Each function in this module represents a specific constraint and returns a value indicating 
whether the constraint is satisfied. A positive return value signifies that the constraint is met, 
while a negative value indicates a violation of the constraint.

Functions:
- constraint1: Ensures the POA of GFAS is always greater than all others.
- constraint2: Ensures PoA values are strictly greater than SoA values for each aerosol type.
- constraint3: Ensures 'other' and 'traffic' POA and SOA values are less than all other values.
"""
def constraint1(ri_values):
    """
    Ensure that the POA of GFAS is always greater than all others.
    
    The constraint is that ri_gfs_poa should be greater than:
    - ri_gfs_soa, ri_res_poa, ri_res_soa
    - ri_shp_poa, ri_shp_soa
    - ri_trf_poa, ri_trf_soa
    - ri_oth_poa, ri_oth_soa
    
    Parameters:
    ri_values (list): List of refractive index values for different species.
    
    Returns:
    float: The difference between ri_gfs_poa and the maximum of the others.
    """
    return ri_values[0] - max(ri_values[1:])

def constraint2(ri_values):
    '''
    Ensures that the values for Primary Organic Aerosol (POA) are always strictly greater than 
    the values for Secondary Organic Aerosol (SOA). The pairs of values are assumed to be 
    ordered in the sequence of PoA and SOA for each category.
    '''
    # Calculate the differences between POA and SOA for each pair
    poa_soa_differences = [ri_values[i] - ri_values[i + 1] for i in range(0, len(ri_values), 2)]

    # Return the minimum difference, subtract a small epsilon to enforce strict inequality
    return min(poa_soa_differences) - 1e-4

def constraint3(ri_values):
    """
    Esta función ajustada establece una restricción donde se asegura que los valores de 
    POA OTH, SOA OTH, POA TRF y SOA TRF sean los más bajos, pero manteniendo la 
    jerarquía requerida por las otras restricciones.

    Parámetros:
    ri_values (list): Lista de valores de índice de refracción para diferentes especies y fuentes.
    
    Retorna:
    float: El valor mínimo que debe ser mayor o igual a cero para que la restricción se considere cumplida.
    """
    poa_oth = ri_values[8]
    soa_oth = ri_values[9]
    poa_trf = ri_values[6]
    soa_trf = ri_values[7]

    # Encuentra el mínimo de los valores de RI excluyendo OTH y TRF
    min_otros = min([valor for i, valor in enumerate(ri_values) if i not in [6, 7, 8, 9]])

    # Asegurarse de que los valores de OTH y TRF sean los más bajos, respetando PoA > SoA
    return min(min_otros - poa_oth, min_otros - soa_oth, poa_oth - soa_oth, poa_trf - soa_trf)

'''
def constraint2(ri_values):
    
    PoA values are always greater than SoA values:
   ri_gfs_poa > ri_gfs_soa, ri_res_poa > ri_res_soa, ri_shp_poa > ri_shp_soa, ri_trf_poa > ri_trf_soa, ri_oth_poa > ri_oth_soa
    
    return min(ri_values[0] - ri_values[1], ri_values[2] - ri_values[3], ri_values[4] - ri_values[5],
            ri_values[6] - ri_values[7], ri_values[8] - ri_values[9])
'''



#res el segundo que mas absorbe, quitar que el gfas absorbe mas 
#
"""
def constraint2(ri_values):
    '''
    Los valores de SoA siempre son estrictamente mayores que los valores de PoA:
    ri_gfs_soa > ri_gfs_poa, ri_res_soa > ri_res_poa, ri_shp_soa > ri_shp_poa, 
    ri_trf_soa > ri_trf_poa, ri_oth_soa > ri_oth_poa
    '''
    # Calcular las diferencias entre SoA y PoA para cada especie y asegurarse de que sean estrictamente mayores a 0
    soa_poa_diferencias = [ri_values[i + 1] - ri_values[i] for i in range(0, len(ri_values), 2)]
    # Devolver la diferencia mínima - esta debe ser mayor que cero para que la restricción se cumpla
    return min(soa_poa_diferencias) - 1e-6  # Restar un pequeño épsilon para garantizar la desigualdad estricta
"""