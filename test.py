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
