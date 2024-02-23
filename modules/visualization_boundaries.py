"""
Visualization Module

This module provides functionality for visualizing the results of the aerosol absorption 
model optimization. It includes functions to plot and compare optimized refractive index 
values against initial values, enhancing the interpretability of the optimization results.

Functions:
- plot_optimized_vs_initial: Plots a comparison of optimized refractive indices vs. initial values.
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


from . import data_processing as dp


def plot_optimized_vs_initial(stations, observed_data, method, mode, initial_refractive_indices, optimization_result):
    """
    Plot the comparison of optimized refractive indices versus initial values for each station.
    
    Parameters:
    stations (list): List of station names.
    observed_data (DataFrame): The observed data as a DataFrame.
    method (str): The optimization method used.
    mode (str): The mode of optimization.
    initial_refractive_indices (list): The initial guess for refractive indices.
    optimization_result (Result): The result of the optimization process.
    """
    directory = "optimized_vs_initial"
    os.makedirs(directory, exist_ok=True)

    # Remove negative values from observed_data
    observed_data = observed_data[observed_data['AbsBrC370'] > 0]

    for station in stations:
        # Extract observed data with date index
        observed = observed_data[observed_data['station_name'] == station][['AbsBrC370', 'time']]
        observed.set_index('time', inplace=True)
        observed.index = pd.to_datetime(observed.index)

        # Extract initial modeled data (this might need adjustment depending on your model dataframe structure)
        initial_modeled = dp.calculate_abs_modeled(station, initial_refractive_indices)
        
        # Calculate optimized values (this will depend on how you get the optimized data using the optimized refractive indices)
        optimized_modeled = dp.calculate_abs_modeled(station, optimization_result.x)

        plt.figure(figsize=(10, 6))
        # Plot observed and modeled data
        plt.plot(observed.index, observed['AbsBrC370'], label='Observed')
        plt.plot(initial_modeled.index, initial_modeled['AbsBrC370'], label='Initial Modeled', marker='o')
        plt.plot(optimized_modeled.index, optimized_modeled['AbsBrC370'], label='Optimized')

        plt.title(f"Model vs Observed: {station}")
        plt.xlabel('Time')
        plt.ylabel('Absorption Coefficient (Mm$^{-1}$)')
        plt.legend()
        plt.savefig(os.path.join(directory, f"{station}-{method}-{mode}.png"))
        plt.close()

def plot_boundaries(data, saleh_class):
    """
    Plot the boundaries for refractive index values.

    Parameters:
    data (dict): The data to plot.
    saleh_class (dict): The Saleh classification for each refractive index type.
    """
    # Convertir los datos a DataFrame y transponer
    df = pd.DataFrame(data).T

    # Asegurarse de que los datos son numéricos y aplicar formato
    for col in df.columns:
        df[col] = df[col].apply(lambda x: f"{x:.4f}" if isinstance(x, (int, float)) else x)

    # Añadir la clasificación de Saleh como una nueva columna
    df['Saleh Classification'] = df.index.map(saleh_class)

    # Crear una carpeta para guardar las imágenes, si aún no existe
    path = 'ri_boundaries'
    os.makedirs(path, exist_ok=True)

    # Usar Seaborn para establecer el estilo de la figura
    sns.set_theme(style="white")

    # Crear la figura con un tamaño más ajustado
    fig, ax = plt.subplots(figsize=(5.2, 3.1))  # Ajustar el tamaño según sea necesario
    ax.axis('off')
    
    # Crear la tabla con un ajuste en el espaciado
    tabla = ax.table(cellText=df.values, colLabels=df.columns, rowLabels=df.index, 
                     loc='center', cellLoc='center', colColours=["#FFD700"]*len(df.columns),
                     bbox=[0, 0, 1, 1])

    # Ajustar el tamaño de la fuente
    tabla.auto_set_font_size(False)
    tabla.set_fontsize(10)
    tabla.scale(.5, .5)  # Escalar la tabla (ancho, alto)

    # Añadir un título a la figura
    #plt.title("Limites del Índice de Refracción para Diferentes Tipos de Aerosoles", fontsize=14, pad=20)
    #give a filename according to saleh classification
    #join the saleh classification values withthe first letter of each word

    saleh_class = saleh_class.values()
    saleh_class = [x[0] for x in saleh_class]
    saleh_class = ''.join(saleh_class)
    filename = f"ri_boundaries_{saleh_class}.png"
    # Guardar la figura
    plt.savefig(os.path.join(path, filename), bbox_inches='tight', dpi=300)


def plotBoundariesNoSecondary(data, saleh_class):
    """
    Plot the boundaries for refractive index values.

    Parameters:
    data (dict): The data to plot.
    saleh_class (dict): The Saleh classification for each refractive index type.
    """
    # Convertir los datos
    df = pd.DataFrame(data).T

    # Asegurarse de que los datos son numéricos y aplicar formato
    for col in df.columns:
        df[col] = df[col].apply(lambda x: f"{x:.4f}" if isinstance(x, (int, float)) else x)

    # Añadir la clasificación de Saleh como una nueva columna
    df['Saleh Classification'] = df.index.map(saleh_class)

    # Crear una carpeta para guardar las imágenes, si aún no existe
    path = 'ri_boundaries_no_secondary'
    os.makedirs(path, exist_ok=True)

    # Usar Seaborn para establecer el estilo de la figura
    sns.set_theme(style="white")

    # Crear la figura con un tamaño más ajustado
    fig, ax = plt.subplots(figsize=(5.2, 3.1))  # Ajustar el tamaño según sea necesario
    ax.axis('off')
    
    # Crear la tabla con un ajuste en el espaciado
    tabla = ax.table(cellText=df.values, colLabels=df.columns, rowLabels=df.index, 
                     loc='center', cellLoc='center', colColours=["#FFD700"]*len(df.columns),
                     bbox=[0, 0, 1, 1])

    # Ajustar el tamaño de la fuente
    tabla.auto_set_font_size(False)
    tabla.set_fontsize(10)
    tabla.scale(.5, .5)  # Escalar la tabla (ancho, alto)

    # Añadir un título a la figura
    #plt.title("Limites del Índice de Refracción para Diferentes Tipos de Aerosoles", fontsize=14, pad=20)
    #give a filename according to saleh classification
    #join the saleh classification values withthe first letter of each word

    saleh_class = saleh_class.values()
    saleh_class = [x[0] for x in saleh_class]
    saleh_class = ''.join(saleh_class)
    filename = f"ri_boundaries_{saleh_class}.png"
    # Guardar la figura
    plt.savefig(os.path.join(path, filename), bbox_inches='tight', dpi=300)