"""
Utilities Module

This module contains various utility functions that support the main functionalities of 
the aerosol absorption model optimization. These include general-purpose functions for 
printing results, handling data conversions, and other supportive tasks.

Functions:
- print_table: Prints a table of refractive index values post-optimization.
"""
import os
import pandas as pd
import plotly.figure_factory as ff
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import matplotlib
from matplotlib.lines import Line2D
import calendar as cal
import numpy as np
import seaborn as sns
import json 

from . import constants
from . import data_processing as dp
from . import sync_obs_model_data as sync
from . import seasonal_data_grouper as sbys
from . import data_retrieval as dr
from . import conc2abs as ca
from . import statistics as st 

def convention_names_stations():
    dict_of_names = {}
    with open('modules/conv_stat_names.json') as file:
        dict_of_names = json.load(file)
    return dict_of_names


def print_table(optimization_result, mode, method, **kwargs):
    try:
        # Prepare the DataFrame
        ri_values = optimization_result.x
        df_result = pd.DataFrame({
            'RI_name': constants.SPECIES,  # This assumes SPECIES is a predefined list of species names
            'ri_values': ri_values
        })
        df_result['ri_values'] = df_result['ri_values'].round(4)

        # Determine the directory path
        model_name = kwargs.get("model", "default_model")
        directory_path = os.path.join('ri_tables', model_name, f'RI_{mode}')
        os.makedirs(directory_path, exist_ok=True)

        # constantsruct the file path based on the mode
        file_parts = ['RI', method]
        if mode == 'all':
            file_parts.append(kwargs.get("mass_data"))
        elif mode in ['by_category', 'by_season', 'by_station', 'by_station_season']:
            file_parts.extend(kwargs.get(part) for part in mode.split('_')[1:] if kwargs.get(part))
            file_parts.append(kwargs.get("mass_data"))
        
        file_name = '_'.join(filter(None, file_parts))
        file_path = os.path.join(directory_path, file_name)
        
        print(f"Saving table to {file_path}.png and {file_path}.csv")

        # Save the table as an image and CSV
        fig = ff.create_table(df_result)
        fig.update_layout(autosize=False, width=500, height=300)
        fig.write_image(f'{file_path}.png', scale=2)
        df_result.to_csv(f'{file_path}.csv', index=False)

        print("Files saved successfully.")
        return df_result
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def TableNoSecondary(optimization_result, mode, method, **kwargs):
    try:
        # Prepare the DataFrame
        ri_values = optimization_result.x
        df_result = pd.DataFrame({
            'RI_name': constants.SPECIESNOSECONDARY,  # This assumes SPECIES is a predefined list of species names
            'ri_values': ri_values
        })
        df_result['ri_values'] = df_result['ri_values'].round(4)

        # Determine the directory path
        model_name = kwargs.get("model", "default_model")
        directory_path = os.path.join('ri_tables', model_name, f'RI_{mode}', 'no_secondary')
        os.makedirs(directory_path, exist_ok=True)

        # constantsruct the file path based on the mode
        file_parts = ['RI', method]
        if mode == 'all':
            file_parts.append(kwargs.get("mass_data"))
        elif mode in ['by_category', 'by_season', 'by_station', 'by_station_season']:
            file_parts.extend(kwargs.get(part) for part in mode.split('_')[1:] if kwargs.get(part))
            file_parts.append(kwargs.get("mass_data"))
        
        file_name = '_'.join(filter(None, file_parts))
        file_path = os.path.join(directory_path, file_name)
        
        print(f"Saving table to {file_path}.png and {file_path}.csv")

        # Save the table as an image and CSV
        fig = ff.create_table(df_result)
        fig.update_layout(autosize=False, width=500, height=300)
        fig.write_image(f'{file_path}.png', scale=2)
        df_result.to_csv(f'{file_path}.csv', index=False)

        print("Files saved successfully.")
        return df_result
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def print_table4oa(optimization_result, mode, method, **kwargs):
    try:
        # Prepare the DataFrame
        ri_values = optimization_result.x
        df_result = pd.DataFrame({
            'RI_name': constants.SPECIES_OA,  # This assumes SPECIES is a predefined list of species names
            'ri_values': ri_values
        })
        df_result['ri_values'] = df_result['ri_values'].round(4)

        # Determine the directory path
        model_name = kwargs.get("model", "default_model")
        directory_path = os.path.join('ri_tables', model_name, f'RI_{mode}')
        os.makedirs(directory_path, exist_ok=True)

        # constantsruct the file path based on the mode
        file_parts = ['RI', method]
        if mode == 'all':
            file_parts.append(kwargs.get("mass_data"))
        elif mode in ['by_category', 'by_season', 'by_station', 'by_station_season']:
            file_parts.extend(kwargs.get(part) for part in mode.split('_')[1:] if kwargs.get(part))
            file_parts.append(kwargs.get("mass_data"))
        
        file_name = '_'.join(filter(None, file_parts))
        file_path = os.path.join(directory_path, file_name)
        
        print(f"Saving table to {file_path}.png and {file_path}.csv")

        # Save the table as an image and CSV
        fig = ff.create_table(df_result)
        fig.update_layout(autosize=False, width=500, height=300)
        fig.write_image(f'{file_path}.png', scale=2)
        df_result.to_csv(f'{file_path}.csv', index=False)

        print("Files saved successfully.")
        return df_result
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def format_spanish_date_to_english(date_series):
    spanish_to_english_months = {
    'ene': 'Jan', 'feb': 'Feb', 'mar': 'Mar', 'abr': 'Apr', 
    'may': 'May', 'jun': 'Jun', 'jul': 'Jul', 'ago': 'Aug', 
    'sep': 'Sep', 'oct': 'Oct', 'nov': 'Nov', 'dic': 'Dec'
    }
    return date_series.map(lambda x: f"{spanish_to_english_months[x[:3].lower()]}-{x[4:]}")

def set_plot_configuration(i, ax, season, mod, x_mod, ylabel=None, title=None, **kwargs):
    ax.set_title(season, fontweight='bold')
    if ylabel:
        ax.set_ylabel(ylabel)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.grid(True, which='both', axis='both', color='gray', linewidth=0.5, linestyle=':', alpha=0.5)
    ax.tick_params(axis='x', labelsize=10, labelrotation=45, colors='gray')
    ax.tick_params(axis='y', labelsize=10, colors='gray')
    #cambiar frecuencia de ticks x
    if len(mod.index) > 50:
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=15))
        #rotar los tick x
        ax.tick_params(axis='x', rotation=40)
    else:
        if len(mod.index) > 10:
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
            #rotar los tick x
            ax.tick_params(axis='x', rotation=40)
        else:
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
            #rotar los tick x
            ax.tick_params(axis='x', rotation=40)

    #hacer que el eje inicie x en -1
    ax.set_xlim([-1, len(x_mod)+0.5])
    #set title
    ax.set_title(season, fontweight='bold')
    #set legend for plot number 1
    if kwargs.get('legend')=='bars':
        pass
    else:
        if i == 1:
            ax.legend()
        else:
            if ax.get_legend() is not None:
                ax.get_legend().remove()
    
def sync_y_axis_limits(axes):
    max_y = max(ax.get_ylim()[1] for ax in axes)
    for ax in axes:
        ax.set_ylim([0, max_y])

def generate_colors():
    original_cmap = mpl.colormaps['tab10']
    colors = original_cmap(np.linspace(0, 1, 10))
    new_cmap = ListedColormap(colors)
    column_colors = [matplotlib.colors.rgb2hex(new_cmap(k)[:3]) for k in range(new_cmap.N)]
    return column_colors

def save_plot(fig, path, filename):
    os.makedirs(path, exist_ok=True)
    fig.savefig(os.path.join(path, filename), dpi=300, bbox_inches='tight')


#get data for each station
def calculate_absorption(station, best=True, mass_mode='all',SA=False, **kwargs):
    if best and kwargs.get('points')=='best':
        df_mod_stn = pd.read_csv('../absorption/NInventory/mod/4brc/2018/BestModObs/best_' + station + '.csv', index_col=0,\
                          parse_dates=True)
    elif best and kwargs.get('points')=='complete':
        df_mod_stn = pd.read_csv('../absorption/NInventory/mod/4brc/2018/4brc_' + station + '.csv', index_col=0,\
                            parse_dates=True)

    stn_m = ca.get_absorption(df_mod_stn, 
                              WAVELENGTH=constants.WAVELENGTH, 
                              REL_HUM=constants.RELATIVE_HUMIDITY, 
                              method= 'SLSQP', 
                              mode = kwargs.get('mode'), 
                              mass_mode=mass_mode, 
                              SA=SA, 
                              model=kwargs.get('model'),
                              station=station,
                              )
    return stn_m

#get data for each station and case 
def calculate_absorption_cases(station, best=True, mass_mode='all',SA=False, **kwargs):
    if best and kwargs.get('points')=='best':
        df_mod_stn = pd.read_csv('../absorption/NInventory/mod/4brc/2018/BestModObs/best_' + station + '.csv', index_col=0,\
                          parse_dates=True)
    elif best and kwargs.get('points')=='complete':
        df_mod_stn = pd.read_csv('../absorption/NInventory/mod/4brc/2018/4brc_' + station + '.csv', index_col=0,\
                            parse_dates=True)

    stn_m = ca.get_absorption_cases(df_mod_stn, 
                              WAVELENGTH=constants.WAVELENGTH, 
                              REL_HUM=constants.RELATIVE_HUMIDITY, 
                              method= 'SLSQP', 
                              mode = kwargs.get('mode'), 
                              mass_mode=mass_mode, 
                              SA=SA, 
                              model=kwargs.get('model'),
                              station=station,
                              case=kwargs.get('case'),
                              )
    return stn_m

def calculate_absorption4oa(station, best=True, mass_mode='all', **kwargs):
    if best and kwargs.get('points')=='best':
        df_mod_stn = pd.read_csv('../absorption/NInventory/mod/4brc/2018/BestModObs/best_' + station + '.csv', index_col=0,\
                          parse_dates=True)
    elif best and kwargs.get('points')=='complete':
        df_mod_stn = pd.read_csv('../absorption/NInventory/mod/4brc/2018/4brc_' + station + '.csv', index_col=0,\
                            parse_dates=True)
        #sum columns with start wit poa and soa and put in oa column mantain nan values as nan
        df_mod_stn['oa'] = df_mod_stn.filter(regex='^poa|^soa').sum(axis=1, skipna=False)
        #now remove poa and soa columns
        df_mod_stn = df_mod_stn.drop(columns=['poagfs','soagfs','poares','soares','poashp','soashp','poatrf','soatrf','poaoth','soaoth'])

    stn_m = ca.get_absorption4oa(df_mod_stn, 
                              WAVELENGTH=constants.WAVELENGTH, 
                              REL_HUM=constants.RELATIVE_HUMIDITY, 
                              method= 'SLSQP', 
                              mode = kwargs.get('mode'), 
                              mass_mode=mass_mode, 
                              model=kwargs.get('model'),
                              station=station,
                              )
    return stn_m


def plot_line_seasonal_abs(stn, points, mode, model):
    fig, ax = plt.subplots(2, 2, figsize=(10, 5))
    # Adjust subplots
    fig.subplots_adjust(hspace=0.45, wspace=0.11)
    ax = ax.flatten()

    model_abs = calculate_absorption(stn, best=True, mass_mode='best', SA=False, mode=mode, points=points, model=model)
    obser_abs = dr.get_obsabs370(stn, remove_negatives=True)
    mod_apr = calculate_absorption(stn, best=True, mass_mode='best', SA=False, mode='apriori', points=points, model=model)

    for i, season in enumerate(['DJF', 'MAM', 'JJA', 'SON']):
        obs = sbys.season(obser_abs)[season]
        mod = sbys.season(model_abs)[season]
        apr = sbys.season(mod_apr)[season]

        # Check if observation data for the season is empty
        if not obs.empty:
            threshold = 10
            linestyle = '-' if len(obs.index) > threshold else ''
            x_obs = format_spanish_date_to_english(obs.index.strftime('%b-%d'))
            ax[i].plot(x_obs,
                       obs,
                       color='green',
                       label='Obs',
                       linewidth=.8,
                       alpha=0.8,
                       marker='^',
                       linestyle=linestyle,
                       markeredgewidth=.8,
                       ms=3)
        x_mod = format_spanish_date_to_english(mod.index.strftime('%b-%d'))
        ax[i].plot(x_mod, apr, color='blue', label='mod-init', ms=3,
                   linewidth=.8, alpha=0.3, marker='s', linestyle='-', markeredgewidth=.8)

        ax[i].plot(x_mod, mod, color='red', label='mod-opt', 
                   linewidth=.8, alpha=0.5, marker='d', linestyle=':', markeredgewidth=.8, markersize=3)

        set_plot_configuration(i, ax[i], season, mod, x_mod, ylabel='Absorption [M$m^-1$]' if i in [0, 2] else None)

    # Set consistent y-axis limits
    #sync_y_axis_limits(ax)
    #set title for the whole figure
    fig.suptitle(stn, fontweight='bold')
    #make a folder to save plots
    path = 'imgs/seasonal_abs/ts/'
    filename = stn +'_'+mode+'_ts.png'
    save_plot(fig, path, filename)

    return plt.close()


def plot_bar_seasonal_abs(stn, points, mode, model):
    fig, ax = plt.subplots(2, 2, figsize=(25, 9))
    fig.subplots_adjust(hspace=0.25, wspace=0.05)
    ax = ax.flatten()

    sorted_cols = ['poagfs', 'soagfs', 'poares', 'soares', 'poashp', 'soashp',
                   'poatrf', 'soatrf', 'poaoth', 'soaoth']
    model_abs = calculate_absorption(stn, best=True, mass_mode='best', SA=True, mode=mode, points=points, model=model)
    obser_abs = dr.get_obsabs370(stn, remove_negatives=True)

    handles, labels = [], []

    for i, season in enumerate(['DJF', 'MAM', 'JJA', 'SON']):
        obs = sbys.season(obser_abs)[season]
        mod = sbys.season(model_abs)[season][sorted_cols]

        if not obs.empty:
            threshold = 10
            linestyle = '-' if len(obs.index) > threshold else ''
            x_obs = format_spanish_date_to_english(obs.index.strftime('%b-%d'))
            ax[i].plot(
                x_obs, obs, 
                marker='o', color='black', markersize=4, 
                markeredgecolor='black', markeredgewidth=1, 
                linestyle=linestyle, label='Obs', lw=0.5, alpha=0.8
            )

        column_colors = generate_colors()
        x_mod = format_spanish_date_to_english(mod.index.strftime('%b-%d'))

        bottom = pd.Series([0] * mod.shape[0], index=mod.index)

        for j, column in enumerate(mod.columns):
            hatch = '/' if column.startswith('soa') else ''
            bar = ax[i].bar(
                x_mod, mod[column], bottom=bottom, 
                color=column_colors[j % len(column_colors)], 
                label=column, hatch=hatch, edgecolor='black', 
                linewidth=0.1, alpha=0.8
            )
            bottom += mod[column]
            if column not in labels:
                handles.append(bar[0])
                labels.append(column)

        set_plot_configuration(i, ax[i], season, mod, x_mod, ylabel='Absorption [M$m^-1$]' if i in [0, 2] else None, legend='bars')

    #sync_y_axis_limits(ax)
    fig.suptitle(stn, fontweight='bold')
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.77, .97), ncol=5)
    #set a second legend for observed absorption at opposite side of the first legend observatios are black circles
    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10, label='Obs')]
    fig.legend(handles=legend_elements, loc='upper center', ncol=3, bbox_to_anchor=(0.2, .97))
    #make a folder to save plots
    path = 'imgs/seasonal_abs/stacked_ts/'
    filename = stn +'_'+mode+'_stacked_ts.png'
    save_plot(fig, path, filename)

    return plt.close()

def plot_line_mass(stn):
    fig, ax = plt.subplots(2, 2, figsize=(10, 5))
    # adjust subplots
    fig.subplots_adjust(hspace=0.45, wspace=0.05)
    ax = ax.flatten()
    
    obs_mass = dr.get_mass_obs(convention_names_stations()["acr_first"][stn], remove_negatives=True)
    mod_mass = pd.DataFrame(dr.get_mass_mod(convention_names_stations()["acr_first"][stn]).sum(axis=1), columns=['pm2p5oa'])

    for i, season in enumerate(['DJF', 'MAM', 'JJA', 'SON']):
        obs = sbys.season(obs_mass)[season]
        mod = sbys.season(mod_mass)[season]

        #x_mod = format_spanish_date_to_english(mod.index.strftime('%b-%d'))
        x_mod = mod.index.strftime('%b-%d')
        # Plot Model Data
        ax[i].plot(x_mod,
                   mod,
                   color='red',
                   label='mod',
                   linewidth=.8,
                   alpha=0.5,
                   marker='d',
                   linestyle=':',
                   markeredgewidth=.8,
                   markersize=3)

        # Plot Observation Data if Not Empty
        if not obs.empty:
            threshold = 10
            linestyle = '-' if len(obs.index) > threshold else ''
            #x_obs = format_spanish_date_to_english(obs.index.strftime('%b-%d'))
            x_obs = obs.index.strftime('%b-%d')
            ax[i].plot(x_obs,
                       obs,
                       color='green',
                       label='Obs',
                       linewidth=.8,
                       alpha=0.8,
                       marker='^',
                       linestyle=linestyle,
                       markeredgewidth=.8,
                       ms=3)

        # Set plot configuration
        set_plot_configuration(i, ax[i], season, mod, x_mod, ylabel='OA [$\mu$g m$^{-3}$]' if i in [0, 2] else None)

    sync_y_axis_limits(ax)
    #set title for the whole figure
    fig.suptitle(stn, fontweight='bold')
    path = 'imgs/seasonal_mass/ts/'
    filename = stn +'_mass_ts.png'
    save_plot(fig, path, filename)
    return plt.show()

def plot_bar_mass(stn):
    fig, ax = plt.subplots(2, 2, figsize=(25, 9))
    # Adjust subplots
    fig.subplots_adjust(hspace=0.25, wspace=0.04)
    ax = ax.flatten()
    #sorted columns
    sorted_cols = ['poagfs', 'soagfs', 'poares', 'soares', 'poashp', 'soashp',
                       'poatrf', 'soatrf', 'poaoth', 'soaoth']
    model_mass = dr.get_mass_mod(stn)
    obs_mass = dr.get_mass_obs(stn, remove_negatives=True)
    #mod_apr = calculate_absorption(stn, best=True, mass_mode='best', SA=False, mode='apriori', points=points)

    handles = []
    labels = []

    for i, season in enumerate(['DJF', 'MAM', 'JJA', 'SON']):
        obs = sbys.season(obs_mass)[season]
        mod = sbys.season(model_mass)[season]
        mod = mod[sorted_cols]

        #apr = sbys.season(mod_apr)[season]
        if not obs.empty:
            threshold = 10
            if len(obs.index) > threshold:
                linestyle = '-'
            else:
                linestyle = ''
            x_obs = format_spanish_date_to_english(obs.index.strftime('%b-%d'))
            ax[i].plot(
                x_obs, obs, 
                marker='o',  # Estilo de marcador redondo
                color='black',  # Color de relleno del marcador
                markersize=4,  # Tamaño del marcador
                markeredgecolor='black',  # Color del borde del marcador
                markeredgewidth=1,  # Ancho del borde del marcador
                linestyle=linestyle,  # Sin línea conectando los puntos
                label='Obs',
                lw=0.5,
                alpha=0.8
            )
            
        column_colors = generate_colors()

        #plot stacked bar for each season
        x_mod = format_spanish_date_to_english(mod.index.strftime('%b-%d'))
        bottom = pd.Series([0] * mod.shape[0], index=mod.index)
        # Graficar cada columna apilada
        for j, column in enumerate(mod.columns):
            hatch = '/' if column.startswith('soa') else ''
            bar = ax[i].bar(
                x_mod, mod[column], 
                bottom=bottom, 
                color=column_colors[j % len(column_colors)], 
                label=column, 
                hatch=hatch, 
                edgecolor='black',
                linewidth=0.1,
                alpha=0.8
            )
            bottom += mod[column]
            if column not in labels:
                    handles.append(bar[0])  # Añadir el primer manejador de la barra
                    labels.append(column)
        set_plot_configuration(i, ax[i], season, mod, x_mod, ylabel='OA [$\mu$g m$^{-3}$]' if i in [0, 2] else None, legend='bars')

    sync_y_axis_limits(ax)
    fig.suptitle(stn, fontweight='bold')
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.77, .97), ncol=5)
    #set a second legend for observed absorption at opposite side of the first legend observatios are black circles
    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10, label='Obs')]
    fig.legend(handles=legend_elements, loc='upper center', ncol=3, bbox_to_anchor=(0.2, .97))

    path = 'imgs/seasonal_mass/mass_stacked_ts/'
    filename = stn +'_mass_stacked_ts.png'
    save_plot(fig, path, filename)
    
    return plt.close()

def calculate_stats(stns, mode, points, model):
    #list to store stats for each station
    stats_list = []
    for stn in stns:
        model_abs = calculate_absorption(stn, best=True, mass_mode='best', SA=False, mode=mode, points=points, model=model)
        obser_abs = dr.get_obsabs370(stn, remove_negatives=True)
        mod_apr = calculate_absorption(stn, best=True, mass_mode='best', SA=False, mode='apriori', points=points, model=model)

        #make a table with stn as index and columns as stats, add multiindex to add apriori and optimized
        stats = pd.DataFrame(index=[stn], columns=['CORR', 'FB', 'FAC2', 'CORR_i', 'FB_i', 'FAC2_i'])
        #calculate stats for optimized
        stats.loc[stn, 'CORR'] = st.calculate_corr(obser_abs['AbsBrC370'], model_abs['AbsBrC370'])
        stats.loc[stn, 'FB'] = st.calculate_fb(obser_abs['AbsBrC370'], model_abs['AbsBrC370'])
        stats.loc[stn, 'FAC2'] = st.calculate_fac2(obser_abs['AbsBrC370'], model_abs['AbsBrC370'])
        #calculate stats for apriori
        stats.loc[stn, 'CORR_i'] = st.calculate_corr(obser_abs['AbsBrC370'], mod_apr['AbsBrC370'])
        stats.loc[stn, 'FB_i'] = st.calculate_fb(obser_abs['AbsBrC370'], mod_apr['AbsBrC370'])
        stats.loc[stn, 'FAC2_i'] = st.calculate_fac2(obser_abs['AbsBrC370'], mod_apr['AbsBrC370'])
        
        #store stats for each station
        stats_list.append(stats)
    
    #concatenate all stats in one dataframe
    sts = pd.concat([i for i in stats_list], axis=0)
    #add multiindex to stats dataframe
    path = 'stats/'
    os.makedirs(path, exist_ok=True)
    filename = mode + '_stats.csv'
    sts.to_csv(os.path.join(path, filename))    
    return sts

def calculate_stats(stations, mode, points, model):
    #list to store stats for each station
    stats = {}
    for stn in stations:
        model_abs = calculate_absorption(stn, best=True, mass_mode='best', SA=False, mode=mode, points=points, model=model)
        obser_abs = dr.get_obsabs370(stn, remove_negatives=True)
        #mod_apr = calculate_absorption(stn, best=True, mass_mode='best', SA=False, mode='apriori', points=points, model=model)
        if stn not in stats:
            stats[stn] = {}
        if model not in stats[stn]:
            stats[stn][model.split("_")[-1]] = {}
        #calculate stats for optimized add key stn and model
        stats[stn][model.split("_")[-1]]['CORR'] = st.calculate_corr(obser_abs['AbsBrC370'], model_abs['AbsBrC370'])
        stats[stn][model.split("_")[-1]]['FB'] = st.calculate_fb(obser_abs['AbsBrC370'], model_abs['AbsBrC370'])
        stats[stn][model.split("_")[-1]]['FAC2'] = st.calculate_fac2(obser_abs['AbsBrC370'], model_abs['AbsBrC370'])

    #return the dictionary with stats for each station
    return stats

#create a function to get a heatmap of ris for each station
def plot_heatmap_ri_by_stn(model, filename):

    #adding convention acronim names
    import json
    dict_of_names = {}
    with open('modules/conv_stat_names.json') as file:
        dict_of_names = json.load(file)
    
    path = 'ri_tables/' + model + '/RI_by_station/'

    files = os.listdir(path)
    #files = [i for i in files if i.endswith('.csv')]
    files = [i for i in files if 'SLSQP' in i and i.endswith('.csv')]
    df_stns = pd.DataFrame()
    for file in files:
        df = pd.read_csv(path + file)
        df = pd.read_csv(path+file)
        # agregar una columna con el nombre de la estacion en el dataframe actual
        df['station'] = file.replace('.csv', '')  # quitar la extension .csv y agregar la estacion
        indexes = df.RI_name.unique()
        df = df.reset_index()
        #agregar el dataframe al dataframe combinado
        df_stns = pd.concat([df_stns, df], ignore_index=True)
    #set station as columns 
    df_stns = df_stns.pivot_table(index='index', columns='station', values='ri_values')
    stns = [col.split('_')[2] for col in df_stns.columns]
    df_stns.columns = [dict_of_names['stn'][stns[i]] for i in range(len(stns))]
    #sort by column names
    df_stns = df_stns.reindex(sorted(df_stns.columns), axis=1)
    #set index as the ri_names
    df_stns.index = indexes
    fig, ax = plt.subplots(figsize=(18, 8))
    ax.set_title(model.split("_")[-1].title(), fontweight='bold')
    sns.heatmap(df_stns, cmap='coolwarm', annot=True, fmt='.3f', linewidths=0.5, cbar=False)
    path2sv = 'imgs/ri_tables/' + model + '/'
    os.makedirs(path2sv, exist_ok=True)
    save_plot(fig, path2sv, filename + '.png')
    return plt.show()

def convert_table_df(station, model):    
    #read files just for one station
    path_stns = f'ri_tables/{model}/RI_by_station_season/'
    files = os.listdir(path_stns)
    #file = f'RI_SLSQP_{}_{}_best.csv'
    files_csv = [file for file in files if file.endswith('.csv')]
    #match files with station name which is the third element of the file name after spliting by _
    files_stn = [file for file in files_csv if file.split('_')[2] == station]
    # inicializar el dataframe
    df_stns_season = pd.DataFrame()
    # iterar sobre los archivos y agregarlos al dataframe
    for file in files_stn:
        df = pd.read_csv(path_stns+file)
        # agregar una columna con el nombre de la estacion en el dataframe actual
        df['station'] = file.replace('.csv', '')  # quitar la extension .csv y agregar la estacion
        df['season'] = df['station'].str.split('_').str[-2]

        #agregar el dataframe al dataframe combinado
        df_stns_season = pd.concat([df_stns_season, df], ignore_index=True)

    df_stns_season['station']= [col.split('_')[2] for col in df_stns_season['station']]
    #sort by ri_name
    df_stns_season = df_stns_season.sort_values(by=['RI_name'])
    df_pivot = df_stns_season.pivot(index='RI_name', columns=['station', 'season'], values='ri_values')
    #ri_name  always should be in order ['poagfs', 'soagfs', poares', 'soares', 'poashp', 'soashp', 'poatrf', 'soatrf', 'poaoth', 'soaoth']
    df_pivot = df_pivot.reindex(['poagfs', 'soagfs', 'poares', 'soares', 'poashp', 'soashp', 'poatrf', 'soatrf', 'poaoth', 'soaoth'])
    #change order of season values to ['DJF', 'MAM', 'JJA', 'SON']      
    df_pivot = df_pivot.reindex(['DJF', 'MAM', 'JJA', 'SON'], axis=1, level=1)

    return df_pivot

def plot_heatmap_ri_by_stn_season(stations, model, filename):
   
    df_stns = []
    for station in stations:
        try:
            df_stns.append(convert_table_df(station, model))
        except:
            continue


    #concatenate all dataframes in one
    df_stns_seasons = pd.concat(df_stns, axis=1)

    #hacer un htmap
    fig, ax = plt.subplots(figsize=(30, 8))
    sns.heatmap(df_stns_seasons, cmap='coolwarm', annot=True, fmt='.3f', linewidths=0.5, cbar=False)
    path2sv = 'imgs/ri_tables/' + model + '/'
    os.makedirs(path2sv, exist_ok=True)
    save_plot(fig, path2sv, filename + '.png')

    return plt.close()

def plot_line_seasonal_abs_brc(stn, ri_brc_strng, ri_brc_blchd):
    fig, ax = plt.subplots(2, 2, figsize=(10, 5))
    # Adjust subplots
    fig.subplots_adjust(hspace=0.45, wspace=0.11)
    ax = ax.flatten()

    model = dr.get_explct4brcmass(stn, remove_negatives=True)
    model_abs = ca.get_absorption4brc(model, WAVELENGTH=constants.WAVELENGTH,\
                                        REL_HUM=constants.RELATIVE_HUMIDITY,\
                                        ri_brc_strng=ri_brc_strng, ri_brc_blchd=ri_brc_blchd, SA=False)
    obser_abs = dr.get_obsabs370(stn, remove_negatives=True)

    for i, season in enumerate(['DJF', 'MAM', 'JJA', 'SON']):
        obs = sbys.season(obser_abs)[season]
        mod = sbys.season(model_abs)[season]

        # Check if observation data for the season is empty
        if not obs.empty:
            threshold = 10
            linestyle = '-' if len(obs.index) > threshold else ''
            x_obs = format_spanish_date_to_english(obs.index.strftime('%b-%d'))
            ax[i].plot(x_obs,
                       obs,
                       color='green',
                       label='Obs',
                       linewidth=.8,
                       alpha=0.8,
                       marker='^',
                       linestyle=linestyle,
                       markeredgewidth=.8,
                       ms=3)
        x_mod = format_spanish_date_to_english(mod.index.strftime('%b-%d'))

        ax[i].plot(x_mod, mod, color='red', label='mod', 
                   linewidth=.8, alpha=0.5, marker='d', linestyle=':', markeredgewidth=.8, markersize=3)
        ax[i].plot(x_mod, mod*10, color='blue', label='modx10', 
                   linewidth=.8, alpha=0.5, marker='d', linestyle=':', markeredgewidth=.8, markersize=3)

        set_plot_configuration(i, ax[i], season, mod, x_mod, ylabel='Absorption [M$m^-1$]' if i in [0, 2] else None)

    # Set consistent y-axis limits
    #sync_y_axis_limits(ax)
    #set title for the whole figure
    fig.suptitle(stn, fontweight='bold')
    #make a folder to save plots
    path = 'imgs/seasonal_abs_brcexplicit/ts/'
    filename = stn +'_BrC_ts.png'
    save_plot(fig, path, filename)

    return plt.close()

def calculate_stats4oa(stations, mode, points, model):
    #list to store stats for each station
    stats = {}
    for stn in stations:

        model_abs = calculate_absorption4oa(stn, best=True, mass_mode='best', mode=mode, points=points, model=model)
        obser_abs = dr.get_obsabs370(stn, remove_negatives=True)
        #mod_apr = calculate_absorption(stn, best=True, mass_mode='best', SA=False, mode='apriori', points=points, model=model)
        if stn not in stats:
            stats[stn] = {}
        if model not in stats[stn]:
            stats[stn][model.split("_")[-1]] = {}
        #calculate stats for optimized add key stn and model
        #stats[stn][model.split("_")[-1]]['CORR'] = st.calculate_corr(obser_abs['AbsBrC370'], model_abs['AbsBrC370'])
        #stats[stn][model.split("_")[-1]]['FB'] = st.calculate_fb(obser_abs['AbsBrC370'], model_abs['AbsBrC370'])
        stats[stn][model.split("_")[-1]]['FAC2'] = st.calculate_fac2(obser_abs['AbsBrC370'], model_abs['AbsBrC370'])

    #return the dictionary with stats for each station
    return stats