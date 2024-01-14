import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

def preproc_data(obs, mod):
    #process the data
    #first check that index is in datetime format otherwise convert it
    if isinstance(obs.index, pd.DatetimeIndex) == False:
        obs.index = pd.to_datetime(obs.index)
    
    if isinstance(mod.index, pd.DatetimeIndex) == False:
        mod.index = pd.to_datetime(mod.index)

    #select data not nan in obs
    obs = obs[~np.isnan(obs)]
    #select data not nan in mod
    mod = mod[~np.isnan(mod)]
    #select data that is in both obs and mod
    obs, mod = obs.align(mod, join='inner')
    return obs, mod

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