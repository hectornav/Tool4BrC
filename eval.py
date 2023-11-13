import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

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

model_abs = calculate_absorption('Rigi', best=True, mass_mode='best', SA=False, mode='by_station', points='complete')
obser_abs = dr.get_obsabs370('Rigi')

def plot_line_seasonal(stn, points):
    fig, ax = plt.subplots(2, 2, figsize=(10, 5))
    # Adjust subplots
    fig.subplots_adjust(hspace=0.45, wspace=0.2)
    ax = ax.flatten()

    model_abs = calculate_absorption(stn, best=True, mass_mode='best', SA=False, mode='by_station', points=points)
    obser_abs = dr.get_obsabs370(stn, remove_negatives=True)
    mod_apr = calculate_absorption(stn, best=True, mass_mode='best', SA=False, mode='apriori', points=points)

    for i, season in enumerate(['DJF', 'MAM', 'JJA', 'SON']):
        obs = sbys.season(obser_abs)[season]
        mod = sbys.season(model_abs)[season]
        apr = sbys.season(mod_apr)[season]

        # Check if observation data for the season is empty
        if obs.empty:
            ax[i].text(0.5, 0.5, f'No data for {season}', horizontalalignment='center', verticalalignment='center', transform=ax[i].transAxes)
            ax[i].set_title(season)
            continue

        ax[i].plot(apr.index.strftime('%b-%d'), apr, color='blue', label='mod-init', ms=3,
                   linewidth=.8, alpha=0.5, marker='s', linestyle='-', markeredgewidth=.8)
        ax[i].plot(obs.index.strftime('%b-%d'), obs, color='green', label='Obs', 
                   linewidth=.8, alpha=0.8, marker='^', linestyle='--', markeredgewidth=.8, ms=3)
        ax[i].plot(mod.index.strftime('%b-%d'), mod, color='red', label='mod-opt', 
                   linewidth=.8, alpha=0.5, marker='d', linestyle=':', markeredgewidth=.8, markersize=3)

        ax[i].set_title(season)
        if i in [0, 2]:
            ax[i].set_ylabel('Absorption [M$m^-1$]')

        # Setting legend for one plot and handling x-axis labels
        if i == 1:
            ax[i].legend()
        else:
            if ax[i].get_legend() is not None:
                ax[i].get_legend().remove()

        # Handling x-axis ticks based on data length
        if len(mod.index) > 50:
            ax[i].xaxis.set_major_locator(mdates.DayLocator(interval=15))
        elif len(mod.index) > 10:
            ax[i].xaxis.set_major_locator(mdates.DayLocator(interval=5))
        else:
            ax[i].xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax[i].tick_params(axis='x', rotation=40)
        ax[i].grid(alpha=0.5, linestyle=':')
        #set title
        ax[i].set_title(season, fontweight='bold')
        fig.suptitle(stn, fontweight='bold')

    return plt.show()


plot_line_seasonal('Barcelona_PalauReial', points='complete')