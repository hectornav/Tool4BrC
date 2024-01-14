import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import matplotlib
from matplotlib.lines import Line2D
import calendar as cal
import plotly.figure_factory as ff

from modules import utils as ut
from modules import data_processing as dp


#model_data_mass = dr.get_model_data_4monarch(mass_data='all')
observed_data_opt = pd.read_csv('../absorption/NInventory/obs/absorption/absorption_brc370.csv')
#observed_data_mass = pd.read_csv('../absorption/NInventory/obs/mass_concentration/oa_mass_concentration.csv')

# Get a list of unique stations from the observed data.
unique_stations = observed_data_opt[['station_name']].drop_duplicates()['station_name'].unique()
# Mapping from Spanish month abbreviations to English



#stn = ['Barcelona_PalauReial_DJF']
#mode = 'by_station'
#points = 'complete'

def get_stations():
    path_stns = 'ri_tables/monarch_best_moderately/RI_by_station_season/'
    files = os.listdir(path_stns)
    files = [file for file in files if file.endswith('.csv')]
    #obtain the station names
    stations = [file.split('_')[2] for file in files]
    stations = list(set(stations))
    return stations

def save_df_as_png(model, mode):
    df = ut.calculate_stats(unique_stations, mode=mode, points='complete', model=model)
    df.reset_index(inplace=True)
    #save the dataframe with ff
    fig = ff.create_table(df)
    path = f'ri_stats/{model}/'
    os.makedirs(path, exist_ok=True)
    #mantener indices
    fig.update_layout(
        autosize=False,
        width=900,
        height=300,
        margin=dict(l=0, r=0, b=0, t=0)
    )
    # Ajustar el tamaño de la fuente (opcional)
    for i in range(len(fig.layout.annotations)):
        fig.layout.annotations[i].font.size = 10

    fig.write_image(f"{path}{model}_{mode}.png")

#save_df_as_png('monarch_best_strongly', 'all')

#ut.plot_heatmap_ri_by_stn('monarch_best?weakly', 'gfas_weakly_by_stn')
#ut.plot_line_seasonal_abs(stn, points, mode, model='monarch_best_moderately')
#ut.plot_heatmap_ri_by_stn_season(get_stations(), 'monarch_best_weakly_s', 'weakly_by_stn_season')

"""
for stn in unique_stations:
    ut.plot_line_seasonal_abs(stn, points='complete', mode='by_station', model='monarch_best?weakly')
    ut.plot_bar_seasonal_abs(stn, points='complete', mode='by_station', model='monarch_best?weakly')
    #plot_line_mass(stn, points='complete')
    #plot_bar_mass(stn)

"""

#calculate statistics 
def calcular_y_combinar_estadisticas(unique_stations):
    """
    Calcula y combina estadísticas para casos weakly, moderately y strongly.

    Args:
    unique_stations: Lista de estaciones únicas.

    Returns:
    DataFrame combinado con estadísticas.
    """
    print('Calculando estadísticas')

    # Calcular estadísticas
    weakly = ut.calculate_stats(unique_stations, mode='by_station', points='complete', model='monarch_best_weakly')
    moderately = ut.calculate_stats(unique_stations, mode='by_station', points='complete', model='monarch_best_moderately')
    strongly = ut.calculate_stats(unique_stations, mode='by_station', points='complete', model='monarch_best_strongly')

    # Inicializar un nuevo diccionario para los datos combinados
    datos_combinados = {}

    def asignar_datos(caso_datos, caso_sufijo):
        """
        Extrae y asigna los datos del caso especificado al diccionario de datos combinados.

        Args:
        caso_datos: Datos del caso específico.
        caso_sufijo: Sufijo para identificar el caso en los nombres de las columnas.
        """
        for estacion, datos in caso_datos.items():
            if estacion not in datos_combinados:
                datos_combinados[estacion] = {}
            datos_estacion = datos[list(datos.keys())[0]]  # Obtener los datos de la estación
            for param in ['CORR', 'FB', 'FAC2']:
                datos_combinados[estacion][f'{param}_{caso_sufijo}'] = datos_estacion.get(param)

    # Extraer y asignar datos de cada caso
    asignar_datos(weakly, 'w')
    asignar_datos(moderately, 'm')
    asignar_datos(strongly, 's')

    # Crear y ordenar el DataFrame combinado
    df_combinado = pd.DataFrame.from_dict(datos_combinados, orient='index')
    columnas_ordenadas = [f'{param}_{caso}' for param in ['CORR', 'FB', 'FAC2'] for caso in ['w', 'm', 's']]
    df_combinado = df_combinado[columnas_ordenadas]

    return df_combinado

#uso:
# df_resultado = calcular_y_combinar_estadisticas(lista_de_estaciones_unicas)
# print(df_resultado)
def mostrar_datos_formateados(unique_stations):
    """
    Calcula estadísticas para estaciones únicas, formatea los datos y muestra un DataFrame.

    Args:
    unique_stations: Lista de estaciones únicas.

    Returns:
    Imprime un DataFrame con los datos formateados.
    """
    # Calcular estadísticas para las estaciones únicas
    d = ut.calculate_stats4oa(unique_stations, mode='by_station', points='complete', model='monarch_best_oa')

    # Convertir el diccionario a un formato adecuado para un DataFrame
    formatted_data = [{'Station': location, 'FAC2': details['oa']['FAC2']} for location, details in d.items()]

    # Crear DataFrame
    df = pd.DataFrame(formatted_data)

    # Mostrar el DataFrame
    print(df)

# uso:
#mostrar_datos_formateados(unique_stations)
    
def calcular_moderately_fac2(unique_stations):
    """
    Calcula estadísticas para el caso 'moderately' y extrae solo el valor de 'FAC2'.

    Args:
    unique_stations: Lista de estaciones únicas.

    Returns:
    DataFrame con los valores de 'FAC2' para el caso 'moderately'.
    """
    print('Calculando estadísticas para el caso moderately...')

    # Calcular estadísticas para el caso 'moderately'
    moderately = ut.calculate_stats(unique_stations, mode='by_station', points='complete', model='monarch_best_moderately')

    datos_moderately_fac2 = {estacion: detalles['moderately']['FAC2'] for estacion, detalles in moderately.items()}

    # Crear DataFrame a partir de los datos extraídos
    df_moderately_fac2 = pd.DataFrame(list(datos_moderately_fac2.items()), columns=['Station', 'FAC2'])

    return df_moderately_fac2

# Ejemplo de uso:
df_resultado = calcular_moderately_fac2(unique_stations)
print(df_resultado)