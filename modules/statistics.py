import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import datetime as dt

def preproc_data(obs, mod):
    # Procesar los datos
    # Primero verifica que el índice esté en formato datetime, de lo contrario conviértelo
    if not isinstance(obs.index, pd.DatetimeIndex):
        obs.index = pd.to_datetime(obs.index)
    
    if not isinstance(mod.index, pd.DatetimeIndex):
        mod.index = pd.to_datetime(mod.index)

    # Asegúrate de que estás trabajando con las mismas columnas en caso de DataFrames
    if isinstance(obs, pd.DataFrame) and isinstance(mod, pd.DataFrame):
        # Alinear columnas
        common_columns = obs.columns.intersection(mod.columns)
        obs = obs[common_columns]
        mod = mod[common_columns]

    # Eliminar filas con NaN en cualquier columna para DataFrames
    # o filas con NaN para Series
    obs = obs.dropna()
    mod = mod.dropna()

    # Alinear los datos en obs y mod basándose en sus índices
    obs, mod = obs.align(mod, join='inner')

    # Asegurar de que obs y mod sean Series o arrays 1D antes de devolverlos
    if isinstance(obs, pd.DataFrame):
        obs = obs.squeeze()
    if isinstance(mod, pd.DataFrame):
        mod = mod.squeeze()

    return obs, mod

def preproc_data_collocation(obs, mod, filled=True):
    # Convertir índices a pd.DatetimeIndex si es necesario
    if not isinstance(obs.index, pd.DatetimeIndex):
        obs.index = pd.to_datetime(obs.index)
    if not isinstance(mod.index, pd.DatetimeIndex):
        mod.index = pd.to_datetime(mod.index)

    # Eliminar filas con NaN
    obs = obs.dropna()
    mod = mod.dropna()

    # Alinear los datos en obs y mod basándose en sus índices de fecha/hora
    obs, mod_ = obs.align(mod, join='inner', axis=0)

    if filled:
        #check en que mes obs no tiene datos
        all_months = pd.date_range(start=min(obs.index.min(), mod.index.min()), 
                            end=max(obs.index.max(), mod.index.max()), freq='MS')

        # Identificar los meses vacíos
        meses_vacios = [mes for mes in all_months if mes not in obs.index]
        mod_collocated = pd.DataFrame()

        for mes in meses_vacios:
            if mod_.index.month.isin([mes.month]).any():
                pass
            else:
                df = mod[mod.index.month==mes.month]
                mod_collocated = pd.concat([mod_collocated, df])

        #concatenar el mod_collocated con el mod_ para tener el mod completo
        mod_c = pd.concat([mod_collocated, mod_])
        #sort index
        mod_c = mod_c.sort_index()

        return obs, mod_c
    else:
        return obs, mod_

def calculate_corr(obs, mod):
    """
    Calculate the correlation between the observed and modeled data.

    Parameters
    ----------
    obs : array_like
        The observed data.
    mod : array_like
        The modeled data.

    Returns
    -------
    float
        The correlation between the observed and modeled data.
    """
    obs, mod = preproc_data(obs, mod)
    #check that is with plot
    #fig, ax = plt.subplots()
    #ax.plot(obs.index, obs.values, label='obs')
    #ax.plot(mod.index, mod.values, label='mod')
    #ax.legend()
    #plt.show()
    #calculate correlation
    corr, _ = stats.pearsonr(obs.values, mod.values)
    return np.round(corr, 2)

#calculate normalized mean bias

def calculate_nmb(obs, mod):
    """
    Calculate the normalized mean bias between the observed and modeled data.

    Parameters
    ----------
    obs : array_like
        The observed data.
    mod : array_like
        The modeled data.

    Returns
    -------
    float
        The normalized mean bias between the observed and modeled data.
    """
    obs, mod = preproc_data(obs, mod)
    #calculate nmb
    nmb = (np.sum(mod-obs)/np.sum(obs))*100
    return np.round(nmb, 2)

#calculate normalized mean error
def calculate_nme(obs, mod):
    """
    Calculate the normalized mean error between the observed and modeled data.

    Parameters
    ----------
    obs : array_like
        The observed data.
    mod : array_like
        The modeled data.

    Returns
    -------
    float
        The normalized mean error between the observed and modeled data.
    """
    obs, mod = preproc_data(obs, mod)
    #calculate nme
    nme = (np.sum(np.abs(mod-obs))/np.sum(obs))*100
    return np.round(nme, 2)

#calculate FAC2
def calculate_fac2(obs, mod):
    """
    Calculate the fraction of predictions within a factor of two (FAC2) between the observed and modeled data.

    Parameters
    ----------
    obs : array_like
        The observed data.
    mod : array_like
        The modeled data.

    Returns
    -------
    float
        The fraction of predictions within a factor of two (FAC2) between the observed and modeled data.
    """
    obs, mod = preproc_data(obs, mod)
    #calculate fac2
    frac = mod / obs
    fac2 =  (100.0 / len(frac)) * len(frac[(frac >= 0.5) & (frac <= 2.0)])
    return np.round(fac2, 2)

#calculate (mean fractional bias)

def calculate_fb(obs, mod):
    """
    Calculate the mean fractional bias between the observed and modeled data.

    Parameters
    ----------
    obs : array_like
        The observed data.
    mod : array_like
        The modeled data.

    Returns
    -------
    float
        The mean fractional bias between the observed and modeled data.
    """
    #2*((np.mean(m)-np.mean(o))/(np.mean(m)+np.mean(o)))*100
    obs, mod = preproc_data(obs, mod)
    #calculate mfb
    #mfb = np.mean((mod-obs) / ((mod+obs) / 2.))
    fb = 2*((np.mean(mod)-np.mean(obs))/(np.mean(mod)+np.mean(obs)))*100 #fractional mean bias according to (Herring et. al 2018)
    return np.round(fb, 2)

#calculate (root mean square error)
def calculate_rmse(obs, mod):
    """
    Calculate the root mean square error between the observed and modeled data.

    Parameters
    ----------
    obs : array_like
        The observed data.
    mod : array_like
        The modeled data.

    Returns
    -------
    float
        The root mean square error between the observed and modeled data.
    """
    obs, mod = preproc_data(obs, mod)
    #calculate rmse
    rmse = np.sqrt(np.mean((mod-obs)**2))
    return np.round(rmse, 2)

#calculate (mean absolute error)
def calculate_mae(obs, mod):
    """
    Calculate the mean absolute error between the observed and modeled data.

    Parameters
    ----------
    obs : array_like
        The observed data.
    mod : array_like
        The modeled data.

    Returns
    -------
    float
        The mean absolute error between the observed and modeled data.
    """
    obs, mod = preproc_data(obs, mod)
    #calculate mae
    mae = np.mean(np.abs(mod-obs)) 
    return np.round(mae, 2)