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
- constraint3: Ensures 'other' POA and SOA values are less than all other values.
"""

'''
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
'''
"""
def constraint2(ri_values):
    '''
    PoA values are always strictly greater than SoA values:
    ri_gfs_poa > ri_gfs_soa, ri_res_poa > ri_res_soa, ri_shp_poa > ri_shp_soa, 
    ri_trf_poa > ri_trf_soa, ri_oth_poa > ri_oth_soa
    '''
    # Calculate the differences between PoA and SoA for each species and ensure they are strictly greater than 0
    poa_soa_differences = [ri_values[i] - ri_values[i + 1] for i in range(0, len(ri_values), 2)]
    # Return the minimum difference - this must be greater than zero for the constraint to be satisfied
    return min(poa_soa_differences) - 1e-6  # Subtract a small epsilon to enforce strict inequality
"""
'''
def constraint2(ri_values):
    
    PoA values are always greater than SoA values:
   ri_gfs_poa > ri_gfs_soa, ri_res_poa > ri_res_soa, ri_shp_poa > ri_shp_soa, ri_trf_poa > ri_trf_soa, ri_oth_poa > ri_oth_soa
    
    return min(ri_values[0] - ri_values[1], ri_values[2] - ri_values[3], ri_values[4] - ri_values[5],
            ri_values[6] - ri_values[7], ri_values[8] - ri_values[9])
'''

#cambiar a traffic
'''
def constraint3(ri_values):
    """Ensure 'ri_oth_poa' and 'ri_oth_soa' are less than all other values.
    
    Args:
        ri_values (list): List of refractive index values.
    
    Returns:
        float: The minimum difference between the 'other' POA/SOA values 
               and all other values.
    """
    # 'ri_oth_poa' less than all other values
    poa_constraint = max(ri_values[:-2]) - ri_values[-2]
    # 'ri_oth_soa' less than all other values
    soa_constraint = max(ri_values[:-1]) - ri_values[-1]
    
    return min(poa_constraint, soa_constraint)

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

def constraint_poa_gfs_menos_absorbente(ri_values):
    """
    Asegura que el POA de GFAS sea siempre el menos absorbente entre todas las categorías.

    La restricción es que ri_gfs_poa debe ser menor que:
    - ri_gfs_soa, ri_res_poa, ri_res_soa
    - ri_shp_poa, ri_shp_soa
    - ri_trf_poa, ri_trf_soa
    - ri_oth_poa, ri_oth_soa
    
    Parámetros:
    ri_values (list): Lista de valores del índice de refracción para diferentes especies.
    
    Devuelve:
    float: La diferencia entre ri_gfs_poa y el máximo de los valores de las otras categorías.
    """
    # El valor de POA de GFAS
    poa_gfs = ri_values[0]

    # El máximo de los valores restantes
    max_others = max(ri_values[1:])

    # La restricción se satisface si el POA de GFAS es menor que el máximo de los otros
    return poa_gfs - max_others
