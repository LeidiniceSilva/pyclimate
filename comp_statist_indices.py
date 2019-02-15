# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "05/26/2018"
__description__ = "Statistical indices to see performance of the model"


from netCDF4 import Dataset

import numpy as np
import scipy.stats as st


def check_dims(model, obs):
    
    if not model.ndim == obs.ndim:
        print('Dims are not equals!')
        exit(1)


def filter_nan(model, obs):
    
    data = np.array([model.flatten(), obs.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]
    
    return data[:, 0], data[:, 1]


def compute_corr(model, obs):

    """
    The input arrays must have the same dimentions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Pearson Linear Correlation
    """
    
    # check_dims(model, obs)
    # model, obs = filter_nan(model, obs)
    corr = np.corrcoef(model, obs)[0][1]
    
    return corr


def compute_mae(model, obs):

    """
    The input arrays must have the same dimentions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Mean Absoluty Error
    """

    check_dims(model, obs)
    model, obs = filter_nan(model, obs)
    mae = np.mean(abs(model, obs))
    
    return mae
    

def compute_rmse(model, obs):

    """
    The input arrays must have the same dimentions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Root Mean Square Error
    """

    # check_dims(model, obs)
    # model, obs = filter_nan(model, obs)
    rmse = np.sqrt(((np.array(model) - np.array(obs)) ** 2).mean()) 
    
    return rmse
    
     
def compute_bias(model, obs):

    """
    The input arrays must have the same dimentions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Mean Bias Error
    """

    check_dims(model, obs)
    model, obs = filter_nan(model, obs)
    bias = np.nanmean((np.array(model) - np.array(obs)), axis=0)
    
    return bias


def compute_pbias(model, obs):

    """
    The input arrays must have the same dimentions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Percentage Bias
    """

    # check_dims(model, obs)
    # model, obs = filter_nan(model, obs)
    pbias = 100.0 * sum(np.array(model) - np.array(obs)) / sum(np.array(obs))
    
    return pbias
        
    
def compute_apb(model, obs):

    """
    The input arrays must have the same dimentions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Absolute Percent Bias
    """

    check_dims(model, obs)
    model, obs = filter_nan(model, obs)
    apb = 100.0 * sum(abs(model, obs)) / sum(obs)
    return apb
    
    
def compute_anomaly(model, obs):

    """
    The input arrays must have the same dimentions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Anomaly and Standard Anomaly
    """

    check_dims(model, obs)
    model, obs = filter_nan(model, obs)
    clim_mean = np.nanmean(obs, axis=0)
    clim_std = np.nanstd(obs, axis=0)
    anomaly = model - clim_mean
    standard_anomaly = (model - clim_mean)/clim_std
    
    return anomaly, standard_anomaly
   
    
def compute_fcst_correc(model, obs, fcst):

    """
    The input arrays must have the same dimentions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Forecast Data Correction
    """

    check_dims(model, obs)
    model, obs = filter_nan(model, obs)

    sim = np.sort(model)
    alpha_mod, loc_mod, beta_mod = ss.gamma.fit(sim, loc=0)
    obs = np.sort(obs)
    alpha_obs, loc_obs, beta_obs = ss.gamma.fit(obs, loc=0)

    fcst_fcst_correc = []
    for i in fcst:
        prob = ss.gamma.cdf(i, alpha_mod, scale=beta_mod)
        fcst_correc.append(ss.gamma.ppf(prob, alpha_obs, scale=beta_obs))
        
    return fcst_correc


def compute_efficient_coefficient(model, obs):

    """
    The input arrays must have the same dimentions
    :Param model: Numpy array with model data
    :Param obs: Numpy array with obs data
    :Return: Nash–Sutcliffe Efficient Coefficient
    """

    check_dims(model, obs)
    model, obs = filter_nan(model, obs)
    nash = 1 - sum((model - obs) ** 2) / sum((obs - np.mean(obs)) ** 2)
    
    return nash
