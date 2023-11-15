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

from modules import constants
from modules import data_processing as dp
from modules import sync_obs_model_data as sync
from modules import seasonal_data_grouper as sbys
from modules import data_retrieval as dr
from modules import conc2abs as ca


model_data_mass = dr.get_model_data_4monarch(mass_data='all')
observed_data_opt = pd.read_csv('../absorption/NInventory/obs/absorption/absorption_brc370.csv')
observed_data_mass = pd.read_csv('../absorption/NInventory/obs/mass_concentration/oa_mass_concentration.csv')

# Get a list of unique stations from the observed data.
unique_stations = observed_data_opt[['station_name']].drop_duplicates()['station_name'].unique()
# Mapping from Spanish month abbreviations to English

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
                              model='monarch_best',
                              station=station,
                              )
    return stn_m


def plot_line_seasonal_abs(stn, points, mode):
    fig, ax = plt.subplots(2, 2, figsize=(10, 5), sharey=True)
    # Adjust subplots
    fig.subplots_adjust(hspace=0.45, wspace=0.2)
    ax = ax.flatten()

    model_abs = calculate_absorption(stn, best=True, mass_mode='best', SA=False, mode=mode, points=points)
    obser_abs = dr.get_obsabs370(stn, remove_negatives=True)
    mod_apr = calculate_absorption(stn, best=True, mass_mode='best', SA=False, mode='apriori', points=points)

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
                   linewidth=.8, alpha=0.5, marker='s', linestyle='-', markeredgewidth=.8)

        ax[i].plot(x_mod, mod, color='red', label='mod-opt', 
                   linewidth=.8, alpha=0.5, marker='d', linestyle=':', markeredgewidth=.8, markersize=3)

        set_plot_configuration(i, ax[i], season, mod, x_mod, ylabel='Absorption [M$m^-1$]' if i in [0, 2] else None)

    # Set consistent y-axis limits
    sync_y_axis_limits(ax)
    #set title for the whole figure
    fig.suptitle(stn, fontweight='bold')
    #make a folder to save plots
    path = 'imgs/seasonal_abs/ts/'
    save_plot(fig, path, stn +'_ts.png')

    return plt.close()


def plot_bar_seasonal_abs(stn, points, mode):
    fig, ax = plt.subplots(2, 2, figsize=(25, 9), sharey=True)
    fig.subplots_adjust(hspace=0.25, wspace=0.05)
    ax = ax.flatten()

    sorted_cols = ['poagfs', 'soagfs', 'poares', 'soares', 'poashp', 'soashp',
                   'poatrf', 'soatrf', 'poaoth', 'soaoth']
    model_abs = calculate_absorption(stn, best=True, mass_mode='best', SA=True, mode=mode, points=points)
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

    sync_y_axis_limits(ax)
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

def plot_line_mass(stn, points):
    fig, ax = plt.subplots(2, 2, figsize=(10, 5), sharey=True)
    # adjust subplots
    fig.subplots_adjust(hspace=0.45, wspace=0.05)
    ax = ax.flatten()

    obs_mass = dr.get_mass_obs(stn, remove_negatives=True)
    mod_mass = pd.DataFrame(dr.get_mass_mod(stn).sum(axis=1), columns=['pm2p5oa'])

    for i, season in enumerate(['DJF', 'MAM', 'JJA', 'SON']):
        obs = sbys.season(obs_mass)[season]
        mod = sbys.season(mod_mass)[season]

        x_mod = format_spanish_date_to_english(mod.index.strftime('%b-%d'))
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

        # Set plot configuration
        set_plot_configuration(i, ax[i], season, mod, x_mod, ylabel='OA [$\mu$g m$^{-3}$]' if i in [0, 2] else None)

    sync_y_axis_limits(ax)
    #set title for the whole figure
    fig.suptitle(stn, fontweight='bold')
    path = 'imgs/seasonal_mass/ts/'
    filename = stn +'_mass_ts.png'
    save_plot(fig, path, filename)
    return plt.close()

def plot_bar_mass(stn):
    fig, ax = plt.subplots(2, 2, figsize=(25, 9), sharey=True)
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

for stn in unique_stations:
    plot_line_seasonal_abs(stn, points='complete', mode='by_station')
    plot_bar_seasonal_abs(stn, points='complete', mode='by_station')
    #plot_line_mass(stn, points='complete')
    #plot_bar_mass(stn)
