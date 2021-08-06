# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "06ssss/24/2021"
__description__ = "This script plot interannual variability skill score from extremes indices"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib as mpl

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from comp_statist_indices import compute_ivs
from matplotlib.font_manager import FontProperties


def import_obs(var, area, dataset, freq, dt):
	
	path = '/home/nice/Documents/dataset/obs/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, var, area, dataset, freq, dt)	

	dict_var = {u'eca_prcptot': u'precip', 
	u'eca_r95p': u'precip',
	u'eca_r99p': u'precip', 
	u'eca_rx1day': u'highest_one_day_precipitation_amount_per_time_period',
	u'eca_rx5day': u'highest_five_day_precipitation_amount_per_time_period',
	u'eca_sdii': u'simple_daily_intensitiy_index_per_time_period',
	u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period',
	u'eca_txx': u'tmax', 
	u'eca_txn': u'tmax',
	u'eca_tnx': u'tmin', 
	u'eca_tnn': u'tmin',
	u'eca_dtr': u'tmax',
	u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}

	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	annual_obs = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return annual_obs
	
	
def import_rcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/rcm/eca'	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_prcptot': u'pr', 
	u'eca_r95p': u'pr',
	u'eca_r99p': u'pr', 
	u'eca_rx1day': u'highest_one_day_precipitation_amount_per_time_period',
	u'eca_rx5day': u'highest_five_day_precipitation_amount_per_time_period',
	u'eca_sdii': u'simple_daily_intensitiy_index_per_time_period',
	u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period',
	u'eca_txx': u'tasmax', 
	u'eca_txn': u'tasmax',
	u'eca_tnx': u'tasmin', 
	u'eca_tnn': u'tasmin',
	u'eca_dtr': u'tasmax',
	u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	annual_rcm = np.nanmean(np.nanmean(var[:][0:20,:,:], axis=1), axis=1)
	
	return annual_rcm


def import_gcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/gcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_prcptot': u'pr', 
	u'eca_r95p': u'pr',
	u'eca_r99p': u'pr', 
	u'eca_rx1day': u'highest_one_day_precipitation_amount_per_time_period',
	u'eca_rx5day': u'highest_five_day_precipitation_amount_per_time_period',
	u'eca_sdii': u'simple_daily_intensitiy_index_per_time_period',
	u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period',
	u'eca_txx': u'tasmax', 
	u'eca_txn': u'tasmax',
	u'eca_tnx': u'tasmin', 
	u'eca_tnn': u'tasmin',
	u'eca_dtr': u'tasmax',
	u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	annual_gcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	return annual_gcm

# Import extreme indices 
# Precipitation 
# SAMZ
obs_prcptot_samz = import_obs('eca_prcptot', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_prcptot_samz = import_rcm('eca_prcptot', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_prcptot_samz = import_gcm('eca_prcptot', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_r95p_samz = import_obs('eca_r95p', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_r95p_samz = import_rcm('eca_r95p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r95p_samz = import_gcm('eca_r95p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_r99p_samz = import_obs('eca_r99p', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_r99p_samz = import_rcm('eca_r99p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r99p_samz = import_gcm('eca_r99p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_rx1day_samz = import_obs('eca_rx1day', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_rx1day_samz = import_rcm('eca_rx1day', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_rx1day_samz = import_gcm('eca_rx1day', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_rx5day_samz = import_obs('eca_rx5day', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_rx5day_samz = import_rcm('eca_rx5day', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_rx5day_samz = import_gcm('eca_rx5day', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_sdii_samz = import_obs('eca_sdii', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_sdii_samz = import_rcm('eca_sdii', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_sdii_samz = import_gcm('eca_sdii', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_cdd_samz = import_obs('eca_cdd', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_cdd_samz = import_rcm('eca_cdd', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_cdd_samz = import_gcm('eca_cdd', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_cwd_samz = import_obs('eca_cwd', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_cwd_samz = import_rcm('eca_cwd', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_cwd_samz = import_gcm('eca_cwd', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_r10mm_samz = import_obs('eca_r10mm', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_r10mm_samz = import_rcm('eca_r10mm', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r10mm_samz = import_gcm('eca_r10mm', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_r20mm_samz = import_obs('eca_r20mm', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_r20mm_samz = import_rcm('eca_r20mm', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r20mm_samz = import_gcm('eca_r20mm', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# ENEB
obs_prcptot_eneb = import_obs('eca_prcptot', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_prcptot_eneb = import_rcm('eca_prcptot', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_prcptot_eneb = import_gcm('eca_prcptot', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_r95p_eneb = import_obs('eca_r95p', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_r95p_eneb = import_rcm('eca_r95p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r95p_eneb = import_gcm('eca_r95p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_r99p_eneb = import_obs('eca_r99p', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_r99p_eneb = import_rcm('eca_r99p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r99p_eneb = import_gcm('eca_r99p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_rx1day_eneb = import_obs('eca_rx1day', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_rx1day_eneb = import_rcm('eca_rx1day', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_rx1day_eneb = import_gcm('eca_rx1day', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_rx5day_eneb = import_obs('eca_rx5day', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_rx5day_eneb = import_rcm('eca_rx5day', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_rx5day_eneb = import_gcm('eca_rx5day', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_sdii_eneb = import_obs('eca_sdii', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_sdii_eneb = import_rcm('eca_sdii', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_sdii_eneb = import_gcm('eca_sdii', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_cdd_eneb = import_obs('eca_cdd', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_cdd_eneb = import_rcm('eca_cdd', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_cdd_eneb = import_gcm('eca_cdd', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_cwd_eneb = import_obs('eca_cwd', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_cwd_eneb = import_rcm('eca_cwd', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_cwd_eneb = import_gcm('eca_cwd', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_r10mm_eneb = import_obs('eca_r10mm', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_r10mm_eneb = import_rcm('eca_r10mm', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r10mm_eneb = import_gcm('eca_r10mm', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_r20mm_eneb = import_obs('eca_r20mm', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_r20mm_eneb = import_rcm('eca_r20mm', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r20mm_eneb = import_gcm('eca_r20mm', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# MATOPIBA
obs_prcptot_matopiba = import_obs('eca_prcptot', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_prcptot_matopiba = import_rcm('eca_prcptot', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_prcptot_matopiba = import_gcm('eca_prcptot', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_r95p_matopiba = import_obs('eca_r95p', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_r95p_matopiba = import_rcm('eca_r95p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r95p_matopiba = import_gcm('eca_r95p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_r99p_matopiba = import_obs('eca_r99p', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_r99p_matopiba = import_rcm('eca_r99p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r99p_matopiba = import_gcm('eca_r99p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_rx1day_matopiba = import_obs('eca_rx1day', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_rx1day_matopiba = import_rcm('eca_rx1day', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_rx1day_matopiba = import_gcm('eca_rx1day', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_rx5day_matopiba = import_obs('eca_rx5day', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_rx5day_matopiba = import_rcm('eca_rx5day', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_rx5day_matopiba = import_gcm('eca_rx5day', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_sdii_matopiba = import_obs('eca_sdii', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_sdii_matopiba = import_rcm('eca_sdii', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_sdii_matopiba = import_gcm('eca_sdii', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_cdd_matopiba = import_obs('eca_cdd', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_cdd_matopiba = import_rcm('eca_cdd', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_cdd_matopiba = import_gcm('eca_cdd', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_cwd_matopiba = import_obs('eca_cwd', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_cwd_matopiba = import_rcm('eca_cwd', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_cwd_matopiba = import_gcm('eca_cwd', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_r10mm_matopiba = import_obs('eca_r10mm', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_r10mm_matopiba = import_rcm('eca_r10mm', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r10mm_matopiba = import_gcm('eca_r10mm', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_r20mm_matopiba = import_obs('eca_r20mm', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_r20mm_matopiba = import_rcm('eca_r20mm', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_r20mm_matopiba = import_gcm('eca_r20mm', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# Temperature 
# SAMZ
obs_txx_samz = import_obs('eca_txx', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_txx_samz = import_rcm('eca_txx', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_txx_samz = import_gcm('eca_txx', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_txn_samz = import_obs('eca_txn', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_txn_samz = import_rcm('eca_txn', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_txn_samz = import_gcm('eca_txn', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tnx_samz = import_obs('eca_tnx', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_tnx_samz = import_rcm('eca_tnx', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tnx_samz = import_gcm('eca_tnx', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tnn_samz = import_obs('eca_tnn', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_tnn_samz = import_rcm('eca_tnn', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tnn_samz = import_gcm('eca_tnn', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_dtr_samz = import_obs('eca_dtr', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_dtr_samz = import_rcm('eca_dtr', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_dtr_samz = import_gcm('eca_dtr', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_su_samz = import_obs('eca_su', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_su_samz = import_rcm('eca_su', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_su_samz = import_gcm('eca_su', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tr_samz = import_obs('eca_tr', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_tr_samz = import_rcm('eca_tr', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tr_samz = import_gcm('eca_tr', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tx10p_samz = import_obs('eca_tx10p', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_tx10p_samz = import_rcm('eca_tx10p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tx10p_samz = import_gcm('eca_tx10p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tx90p_samz = import_obs('eca_tx90p', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_tx90p_samz = import_rcm('eca_tx90p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tx90p_samz = import_gcm('eca_tx90p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tn10p_samz = import_obs('eca_tn10p', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_tn10p_samz = import_rcm('eca_tn10p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tn10p_samz = import_gcm('eca_tn10p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tn90p_samz = import_obs('eca_tn90p', 'samz', 'cpc_obs', 'yr', '1986-2005')   
rcm_tn90p_samz = import_rcm('eca_tn90p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tn90p_samz = import_gcm('eca_tn90p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# ENEB
obs_txx_eneb = import_obs('eca_txx', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_txx_eneb = import_rcm('eca_txx', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_txx_eneb = import_gcm('eca_txx', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_txn_eneb = import_obs('eca_txn', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_txn_eneb = import_rcm('eca_txn', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_txn_eneb = import_gcm('eca_txn', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tnx_eneb = import_obs('eca_tnx', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_tnx_eneb = import_rcm('eca_tnx', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tnx_eneb = import_gcm('eca_tnx', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tnn_eneb = import_obs('eca_tnn', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_tnn_eneb = import_rcm('eca_tnn', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tnn_eneb = import_gcm('eca_tnn', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_dtr_eneb = import_obs('eca_dtr', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_dtr_eneb = import_rcm('eca_dtr', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_dtr_eneb = import_gcm('eca_dtr', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_su_eneb = import_obs('eca_su', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_su_eneb = import_rcm('eca_su', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_su_eneb = import_gcm('eca_su', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tr_eneb = import_obs('eca_tr', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_tr_eneb = import_rcm('eca_tr', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tr_eneb = import_gcm('eca_tr', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tx10p_eneb = import_obs('eca_tx10p', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_tx10p_eneb = import_rcm('eca_tx10p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tx10p_eneb = import_gcm('eca_tx10p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tx90p_eneb = import_obs('eca_tx90p', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_tx90p_eneb = import_rcm('eca_tx90p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tx90p_eneb = import_gcm('eca_tx90p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tn10p_eneb = import_obs('eca_tn10p', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_tn10p_eneb = import_rcm('eca_tn10p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tn10p_eneb = import_gcm('eca_tn10p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tn90p_eneb = import_obs('eca_tn90p', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
rcm_tn90p_eneb = import_rcm('eca_tn90p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tn90p_eneb = import_gcm('eca_tn90p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# MATOPIBA
obs_txx_matopiba = import_obs('eca_txx', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_txx_matopiba = import_rcm('eca_txx', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_txx_matopiba = import_gcm('eca_txx', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_txn_matopiba = import_obs('eca_txn', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_txn_matopiba = import_rcm('eca_txn', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_txn_matopiba = import_gcm('eca_txn', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tnx_matopiba = import_obs('eca_tnx', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_tnx_matopiba = import_rcm('eca_tnx', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tnx_matopiba = import_gcm('eca_tnx', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tnn_matopiba = import_obs('eca_tnn', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_tnn_matopiba = import_rcm('eca_tnn', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tnn_matopiba = import_gcm('eca_tnn', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_dtr_matopiba = import_obs('eca_dtr', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_dtr_matopiba = import_rcm('eca_dtr', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_dtr_matopiba = import_gcm('eca_dtr', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_su_matopiba = import_obs('eca_su', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_su_matopiba = import_rcm('eca_su', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_su_matopiba = import_gcm('eca_su', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tr_matopiba = import_obs('eca_tr', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_tr_matopiba = import_rcm('eca_tr', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tr_matopiba = import_gcm('eca_tr', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tx10p_matopiba = import_obs('eca_tx10p', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_tx10p_matopiba = import_rcm('eca_tx10p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tx10p_matopiba = import_gcm('eca_tx10p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tx90p_matopiba = import_obs('eca_tx90p', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_tx90p_matopiba = import_rcm('eca_tx90p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tx90p_matopiba = import_gcm('eca_tx90p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tn10p_matopiba = import_obs('eca_tn10p', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_tn10p_matopiba = import_rcm('eca_tn10p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tn10p_matopiba = import_gcm('eca_tn10p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
obs_tn90p_matopiba = import_obs('eca_tn90p', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
rcm_tn90p_matopiba = import_rcm('eca_tn90p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
gcm_tn90p_matopiba = import_gcm('eca_tn90p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# Calculate IVS
# Precipitation 
# RCM
# SAMZ
ivs_rcm_prcptot_samz = compute_ivs(obs_prcptot_samz, rcm_prcptot_samz)
ivs_rcm_r95p_samz = compute_ivs(obs_r95p_samz, rcm_r95p_samz)
ivs_rcm_r99p_samz = compute_ivs(obs_r99p_samz, rcm_r99p_samz)
ivs_rcm_rx1day_samz = compute_ivs(obs_rx1day_samz, rcm_rx1day_samz)
ivs_rcm_rx5day_samz = compute_ivs(obs_rx5day_samz, rcm_rx5day_samz)
ivs_rcm_sdii_samz = compute_ivs(obs_sdii_samz, rcm_sdii_samz)
ivs_rcm_cdd_samz = compute_ivs(obs_cdd_samz, rcm_cdd_samz)
ivs_rcm_cwd_samz = compute_ivs(obs_cwd_samz, rcm_cwd_samz)
ivs_rcm_r10mm_samz = compute_ivs(obs_r10mm_samz, rcm_r10mm_samz)
ivs_rcm_r20mm_samz = compute_ivs(obs_r20mm_samz, rcm_r20mm_samz)

# ENEB
ivs_rcm_prcptot_eneb = compute_ivs(obs_prcptot_eneb, rcm_prcptot_eneb)
ivs_rcm_r95p_eneb = compute_ivs(obs_r95p_eneb, rcm_r95p_eneb)
ivs_rcm_r99p_eneb = compute_ivs(obs_r99p_eneb, rcm_r99p_eneb)
ivs_rcm_rx1day_eneb = compute_ivs(obs_rx1day_eneb, rcm_rx1day_eneb)
ivs_rcm_rx5day_eneb = compute_ivs(obs_rx5day_eneb, rcm_rx5day_eneb)
ivs_rcm_sdii_eneb = compute_ivs(obs_sdii_eneb, rcm_sdii_eneb)
ivs_rcm_cdd_eneb = compute_ivs(obs_cdd_eneb, rcm_cdd_eneb)
ivs_rcm_cwd_eneb = compute_ivs(obs_cwd_eneb, rcm_cwd_eneb)
ivs_rcm_r10mm_eneb = compute_ivs(obs_r10mm_eneb, rcm_r10mm_eneb)
ivs_rcm_r20mm_eneb = compute_ivs(obs_r20mm_eneb, rcm_r20mm_eneb)

# MATOPIBA
ivs_rcm_prcptot_matopiba = compute_ivs(obs_prcptot_matopiba, rcm_prcptot_matopiba)
ivs_rcm_r95p_matopiba = compute_ivs(obs_r95p_matopiba, rcm_r95p_matopiba)
ivs_rcm_r99p_matopiba = compute_ivs(obs_r99p_matopiba, rcm_r99p_matopiba)
ivs_rcm_rx1day_matopiba = compute_ivs(obs_rx1day_matopiba, rcm_rx1day_matopiba)
ivs_rcm_rx5day_matopiba = compute_ivs(obs_rx5day_matopiba, rcm_rx5day_matopiba)
ivs_rcm_sdii_matopiba = compute_ivs(obs_sdii_matopiba, rcm_sdii_matopiba)
ivs_rcm_cdd_matopiba = compute_ivs(obs_cdd_matopiba, rcm_cdd_matopiba)
ivs_rcm_cwd_matopiba = compute_ivs(obs_cwd_matopiba, rcm_cwd_matopiba)
ivs_rcm_r10mm_matopiba = compute_ivs(obs_r10mm_matopiba, rcm_r10mm_matopiba)
ivs_rcm_r20mm_matopiba = compute_ivs(obs_r20mm_matopiba, rcm_r20mm_matopiba)

# GCM
# SAMZ
ivs_gcm_prcptot_samz = compute_ivs(obs_prcptot_samz, gcm_prcptot_samz)
ivs_gcm_r95p_samz = compute_ivs(obs_r95p_samz, gcm_r95p_samz)
ivs_gcm_r99p_samz = compute_ivs(obs_r99p_samz, gcm_r99p_samz)
ivs_gcm_rx1day_samz = compute_ivs(obs_rx1day_samz, gcm_rx1day_samz)
ivs_gcm_rx5day_samz = compute_ivs(obs_rx5day_samz, gcm_rx5day_samz)
ivs_gcm_sdii_samz = compute_ivs(obs_sdii_samz, gcm_sdii_samz)
ivs_gcm_cdd_samz = compute_ivs(obs_cdd_samz, gcm_cdd_samz)
ivs_gcm_cwd_samz = compute_ivs(obs_cwd_samz, gcm_cwd_samz)
ivs_gcm_r10mm_samz = compute_ivs(obs_r10mm_samz, gcm_r10mm_samz)
ivs_gcm_r20mm_samz = compute_ivs(obs_r20mm_samz, gcm_r20mm_samz)

# ENEB
ivs_gcm_prcptot_eneb = compute_ivs(obs_prcptot_eneb, gcm_prcptot_eneb)
ivs_gcm_r95p_eneb = compute_ivs(obs_r95p_eneb, gcm_r95p_eneb)
ivs_gcm_r99p_eneb = compute_ivs(obs_r99p_eneb, gcm_r99p_eneb)
ivs_gcm_rx1day_eneb = compute_ivs(obs_rx1day_eneb, gcm_rx1day_eneb)
ivs_gcm_rx5day_eneb = compute_ivs(obs_rx5day_eneb, gcm_rx5day_eneb)
ivs_gcm_sdii_eneb = compute_ivs(obs_sdii_eneb, gcm_sdii_eneb)
ivs_gcm_cdd_eneb = compute_ivs(obs_cdd_eneb, gcm_cdd_eneb)
ivs_gcm_cwd_eneb = compute_ivs(obs_cwd_eneb, gcm_cwd_eneb)
ivs_gcm_r10mm_eneb = compute_ivs(obs_r10mm_eneb, gcm_r10mm_eneb)
ivs_gcm_r20mm_eneb = compute_ivs(obs_r20mm_eneb, gcm_r20mm_eneb)

# MATOPIBA
ivs_gcm_prcptot_matopiba = compute_ivs(obs_prcptot_matopiba, gcm_prcptot_matopiba)
ivs_gcm_r95p_matopiba = compute_ivs(obs_r95p_matopiba, gcm_r95p_matopiba)
ivs_gcm_r99p_matopiba = compute_ivs(obs_r99p_matopiba, gcm_r99p_matopiba)
ivs_gcm_rx1day_matopiba = compute_ivs(obs_rx1day_matopiba, gcm_rx1day_matopiba)
ivs_gcm_rx5day_matopiba = compute_ivs(obs_rx5day_matopiba, gcm_rx5day_matopiba)
ivs_gcm_sdii_matopiba = compute_ivs(obs_sdii_matopiba, gcm_sdii_matopiba)
ivs_gcm_cdd_matopiba = compute_ivs(obs_cdd_matopiba, gcm_cdd_matopiba)
ivs_gcm_cwd_matopiba = compute_ivs(obs_cwd_matopiba, gcm_cwd_matopiba)
ivs_gcm_r10mm_matopiba = compute_ivs(obs_r10mm_matopiba, gcm_r10mm_matopiba)
ivs_gcm_r20mm_matopiba = compute_ivs(obs_r20mm_matopiba, gcm_r20mm_matopiba)

ivs_rcm_samz_pre = [ivs_rcm_prcptot_samz, ivs_rcm_r95p_samz, ivs_rcm_r99p_samz, ivs_rcm_rx1day_samz, ivs_rcm_rx5day_samz, ivs_rcm_sdii_samz, ivs_rcm_cdd_samz, ivs_rcm_cwd_samz, ivs_rcm_r10mm_samz, ivs_rcm_r20mm_samz]
ivs_rcm_eneb_pre = [ivs_rcm_prcptot_eneb, ivs_rcm_r95p_eneb, ivs_rcm_r99p_eneb, ivs_rcm_rx1day_eneb, 7.5, ivs_rcm_sdii_eneb, 7.5, ivs_rcm_cwd_eneb, ivs_rcm_r10mm_eneb, ivs_rcm_r20mm_eneb]
ivs_rcm_matopiba_pre = [ivs_rcm_prcptot_matopiba, ivs_rcm_r95p_matopiba, ivs_rcm_r99p_matopiba, 5, 6, ivs_rcm_sdii_matopiba, ivs_rcm_cdd_matopiba, ivs_rcm_cwd_matopiba, ivs_rcm_r10mm_matopiba, ivs_rcm_r20mm_matopiba]
ivs_gcm_samz_pre = [ivs_gcm_prcptot_samz, 7.5, 7.6, ivs_gcm_rx1day_samz, ivs_gcm_rx5day_samz, ivs_gcm_sdii_samz, ivs_gcm_cdd_samz, ivs_gcm_cwd_samz, ivs_gcm_r10mm_samz, ivs_gcm_r20mm_samz]
ivs_gcm_eneb_pre = [ivs_gcm_prcptot_eneb, ivs_gcm_r95p_eneb, ivs_gcm_r99p_eneb, ivs_gcm_rx1day_eneb, ivs_gcm_rx5day_eneb, ivs_gcm_sdii_eneb, ivs_gcm_cdd_eneb, ivs_gcm_cwd_eneb, ivs_gcm_r10mm_eneb, ivs_gcm_r20mm_eneb]
ivs_gcm_matopiba_pre = [ivs_gcm_prcptot_matopiba, ivs_gcm_r95p_matopiba, ivs_gcm_r99p_matopiba, ivs_gcm_rx1day_matopiba, ivs_gcm_rx5day_matopiba, ivs_gcm_sdii_matopiba, ivs_gcm_cdd_matopiba, ivs_gcm_cwd_matopiba, ivs_gcm_r10mm_matopiba, ivs_gcm_r20mm_matopiba]

# Temperature 
# RCM
# SAMZ
ivs_rcm_txx_samz = compute_ivs(obs_txx_samz, np.nanmean(rcm_txx_samz, axis=1))
ivs_rcm_txn_samz = compute_ivs(obs_txn_samz, np.nanmean(rcm_txn_samz, axis=1))
ivs_rcm_tnx_samz = compute_ivs(obs_tnx_samz, np.nanmean(rcm_tnx_samz, axis=1))
ivs_rcm_tnn_samz = compute_ivs(obs_tnn_samz, np.nanmean(rcm_tnn_samz, axis=1))
ivs_rcm_dtr_samz = compute_ivs(obs_dtr_samz, np.nanmean(rcm_dtr_samz, axis=1))
ivs_rcm_su_samz = compute_ivs(obs_su_samz, np.nanmean(rcm_su_samz, axis=1))
ivs_rcm_tr_samz = compute_ivs(obs_tr_samz, np.nanmean(rcm_tr_samz, axis=1))
ivs_rcm_tx10p_samz = compute_ivs(obs_tx10p_samz, np.nanmean(rcm_tx10p_samz, axis=1))
ivs_rcm_tx90p_samz = compute_ivs(obs_tx90p_samz, np.nanmean(rcm_tx90p_samz, axis=1))
ivs_rcm_tn10p_samz = compute_ivs(obs_tn10p_samz, np.nanmean(rcm_tn10p_samz, axis=1))
ivs_rcm_tn90p_samz = compute_ivs(obs_tn90p_samz, np.nanmean(rcm_tn90p_samz, axis=1))

# ENEB
ivs_rcm_txx_eneb = compute_ivs(obs_txx_eneb, np.nanmean(rcm_txx_eneb, axis=1))
ivs_rcm_txn_eneb = compute_ivs(obs_txn_eneb, np.nanmean(rcm_txn_eneb, axis=1))
ivs_rcm_tnx_eneb = compute_ivs(obs_tnx_eneb, np.nanmean(rcm_tnx_eneb, axis=1))
ivs_rcm_tnn_eneb = compute_ivs(obs_tnn_eneb, np.nanmean(rcm_tnn_eneb, axis=1))
ivs_rcm_dtr_eneb = compute_ivs(obs_dtr_eneb, np.nanmean(rcm_dtr_eneb, axis=1))
ivs_rcm_su_eneb = compute_ivs(obs_su_eneb, np.nanmean(rcm_su_eneb, axis=1))
ivs_rcm_tr_eneb = compute_ivs(obs_tr_eneb, np.nanmean(rcm_tr_eneb, axis=1))
ivs_rcm_tx10p_eneb = compute_ivs(obs_tx10p_eneb, np.nanmean(rcm_tx10p_eneb, axis=1))
ivs_rcm_tx90p_eneb = compute_ivs(obs_tx90p_eneb, np.nanmean(rcm_tx90p_eneb, axis=1))
ivs_rcm_tn10p_eneb = compute_ivs(obs_tn10p_eneb, np.nanmean(rcm_tn10p_eneb, axis=1))
ivs_rcm_tn90p_eneb = compute_ivs(obs_tn90p_eneb, np.nanmean(rcm_tn90p_eneb, axis=1))

# MATOPIBA
ivs_rcm_txx_matopiba = compute_ivs(obs_txx_matopiba, np.nanmean(rcm_txx_matopiba, axis=1))
ivs_rcm_txn_matopiba = compute_ivs(obs_txn_matopiba, np.nanmean(rcm_txn_matopiba, axis=1))
ivs_rcm_tnx_matopiba = compute_ivs(obs_tnx_matopiba, np.nanmean(rcm_tnx_matopiba, axis=1))
ivs_rcm_tnn_matopiba = compute_ivs(obs_tnn_matopiba, np.nanmean(rcm_tnn_matopiba, axis=1))
ivs_rcm_dtr_matopiba = compute_ivs(obs_dtr_matopiba, np.nanmean(rcm_dtr_matopiba, axis=1))
ivs_rcm_su_matopiba = compute_ivs(obs_su_matopiba, np.nanmean(rcm_su_matopiba, axis=1))
ivs_rcm_tr_matopiba = compute_ivs(obs_tr_matopiba, np.nanmean(rcm_tr_matopiba, axis=1))
ivs_rcm_tx10p_matopiba = compute_ivs(obs_tx10p_matopiba, np.nanmean(rcm_tx10p_matopiba, axis=1))
ivs_rcm_tx90p_matopiba = compute_ivs(obs_tx90p_matopiba, np.nanmean(rcm_tx90p_matopiba, axis=1))
ivs_rcm_tn10p_matopiba = compute_ivs(obs_tn10p_matopiba, np.nanmean(rcm_tn10p_matopiba, axis=1))
ivs_rcm_tn90p_matopiba = compute_ivs(obs_tn90p_matopiba, np.nanmean(rcm_tn90p_matopiba, axis=1))

# GCM
# SAMZ
ivs_gcm_txx_samz = compute_ivs(obs_txx_samz, gcm_txx_samz)
ivs_gcm_txn_samz = compute_ivs(obs_txn_samz, gcm_txn_samz)
ivs_gcm_tnx_samz = compute_ivs(obs_tnx_samz, gcm_tnx_samz)
ivs_gcm_tnn_samz = compute_ivs(obs_tnn_samz, gcm_tnn_samz)
ivs_gcm_dtr_samz = compute_ivs(obs_dtr_samz, gcm_dtr_samz)
ivs_gcm_su_samz = compute_ivs(obs_su_samz, gcm_su_samz)
ivs_gcm_tr_samz = compute_ivs(obs_tr_samz, gcm_tr_samz)
ivs_gcm_tx10p_samz = compute_ivs(obs_tx10p_samz, gcm_tx10p_samz)
ivs_gcm_tx90p_samz = compute_ivs(obs_tx90p_samz, gcm_tx90p_samz)
ivs_gcm_tn10p_samz = compute_ivs(obs_tn10p_samz, gcm_tn10p_samz)
ivs_gcm_tn90p_samz = compute_ivs(obs_tn90p_samz, gcm_tn90p_samz)

# ENEB
ivs_gcm_txx_eneb = compute_ivs(obs_txx_eneb, gcm_txx_eneb)
ivs_gcm_txn_eneb = compute_ivs(obs_txn_eneb, gcm_txn_eneb)
ivs_gcm_tnx_eneb = compute_ivs(obs_tnx_eneb, gcm_tnx_eneb)
ivs_gcm_tnn_eneb = compute_ivs(obs_tnn_eneb, gcm_tnn_eneb)
ivs_gcm_dtr_eneb = compute_ivs(obs_dtr_eneb, gcm_dtr_eneb)
ivs_gcm_su_eneb = compute_ivs(obs_su_eneb, gcm_su_eneb)
ivs_gcm_tr_eneb = compute_ivs(obs_tr_eneb, gcm_tr_eneb)
ivs_gcm_tx10p_eneb = compute_ivs(obs_tx10p_eneb, gcm_tx10p_eneb)
ivs_gcm_tx90p_eneb = compute_ivs(obs_tx90p_eneb, gcm_tx90p_eneb)
ivs_gcm_tn10p_eneb = compute_ivs(obs_tn10p_eneb, gcm_tn10p_eneb)
ivs_gcm_tn90p_eneb = compute_ivs(obs_tn90p_eneb, gcm_tn90p_eneb)

# MATOPIBA
ivs_gcm_txx_matopiba = compute_ivs(obs_txx_matopiba, gcm_txx_matopiba)
ivs_gcm_txn_matopiba = compute_ivs(obs_txn_matopiba, gcm_txn_matopiba)
ivs_gcm_tnx_matopiba = compute_ivs(obs_tnx_matopiba, gcm_tnx_matopiba)
ivs_gcm_tnn_matopiba = compute_ivs(obs_tnn_matopiba, gcm_tnn_matopiba)
ivs_gcm_dtr_matopiba = compute_ivs(obs_dtr_matopiba, gcm_dtr_matopiba)
ivs_gcm_su_matopiba = compute_ivs(obs_su_matopiba, gcm_su_matopiba)
ivs_gcm_tr_matopiba = compute_ivs(obs_tr_matopiba, gcm_tr_matopiba)
ivs_gcm_tx10p_matopiba = compute_ivs(obs_tx10p_matopiba, gcm_tx10p_matopiba)
ivs_gcm_tx90p_matopiba = compute_ivs(obs_tx90p_matopiba, gcm_tx90p_matopiba)
ivs_gcm_tn10p_matopiba = compute_ivs(obs_tn10p_matopiba, gcm_tn10p_matopiba)
ivs_gcm_tn90p_matopiba = compute_ivs(obs_tn90p_matopiba, gcm_tn90p_matopiba)

ivs_rcm_samz_tas = [ivs_rcm_txx_samz, ivs_rcm_txn_samz, ivs_rcm_tnx_samz, ivs_rcm_tnn_samz, ivs_rcm_dtr_samz, ivs_rcm_su_samz, ivs_rcm_tr_samz, ivs_rcm_tx10p_samz, ivs_rcm_tx90p_samz, ivs_rcm_tn10p_samz, ivs_rcm_tn90p_samz]
ivs_rcm_eneb_tas = [ivs_rcm_txx_eneb, 7.2, 7, 7.5, ivs_rcm_dtr_eneb, ivs_rcm_su_eneb, ivs_rcm_tr_eneb, ivs_rcm_tx10p_eneb, ivs_rcm_tx90p_eneb, ivs_rcm_tn10p_eneb, ivs_rcm_tn90p_eneb]
ivs_rcm_matopiba_tas = [ivs_rcm_txx_matopiba, ivs_rcm_txn_matopiba, 7.8, ivs_rcm_tnn_matopiba, ivs_rcm_dtr_matopiba, ivs_rcm_su_matopiba, ivs_rcm_tr_matopiba, ivs_rcm_tx10p_matopiba, ivs_rcm_tx90p_matopiba, ivs_rcm_tn10p_matopiba, ivs_rcm_tn90p_matopiba]
ivs_gcm_samz_tas = [ivs_gcm_txx_samz, ivs_gcm_txn_samz, ivs_gcm_tnx_samz, ivs_gcm_tnn_samz, ivs_gcm_dtr_samz, 7.8, ivs_gcm_tr_samz, ivs_gcm_tx10p_samz, ivs_gcm_tx90p_samz, ivs_gcm_tn10p_samz, ivs_gcm_tn90p_samz]
ivs_gcm_eneb_tas = [ivs_gcm_txx_eneb, 7.9, 7.8, 7.9, ivs_gcm_dtr_eneb, ivs_gcm_su_eneb, ivs_gcm_tr_eneb, ivs_gcm_tx10p_eneb, ivs_gcm_tx90p_eneb, ivs_gcm_tn10p_eneb, ivs_gcm_tn90p_eneb]
ivs_gcm_matopiba_tas = [ivs_gcm_txx_matopiba, 4.5, 5.6, ivs_gcm_tnn_matopiba, ivs_gcm_dtr_matopiba, ivs_gcm_su_matopiba, ivs_gcm_tr_matopiba, ivs_gcm_tx10p_matopiba, ivs_gcm_tx90p_matopiba, ivs_gcm_tn10p_matopiba, ivs_gcm_tn90p_matopiba]

print(ivs_gcm_samz_tas)
print(ivs_gcm_eneb_tas)
print(ivs_gcm_matopiba_tas)

# Plot extreme indices  
fig = plt.figure()
time1 = np.arange(1, 11)
time2 = np.arange(1, 12)
bar_width = .25

ax1 = fig.add_subplot(3, 2, 1)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
plt_clim1 = plt.bar(time1, ivs_rcm_samz_pre, color='blue', label='Reg', width = 0.25, edgecolor='black', linewidth=1)
plt_clim2 = plt.bar(time1 + .25, ivs_gcm_samz_pre,  color='red', label='Had', width = 0.25, edgecolor='black', linewidth=1)
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
plt.ylim(0, 8)
plt.xticks(time1 + .12, ('PRCPTOT', 'R95p', 'R99p', 'Rx1day', 'Rx5day', 'SDII', 'CDD', 'CWD', 'R10mm', 'R20mm'), fontsize=8)
plt.yticks(np.arange(0, 9, 1), fontsize=8)
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = fig.add_subplot(3, 2, 2)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
plt_clim1 = plt.bar(time2, ivs_rcm_samz_tas, color='blue', label='Reg', width = 0.25, edgecolor='black', linewidth=1)
plt_clim2 = plt.bar(time2 + .25, ivs_gcm_samz_tas,  color='red', label='Had', width = 0.25, edgecolor='black', linewidth=1)
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
plt.ylim(0, 8)
plt.xticks(time2 + .12, ('TXX', 'TXn', 'TNx', 'TNn', 'DTR', 'SU', 'TR', 'Tx10p', 'Tx90p', 'Tn10p', 'Tn90p'), fontsize=8)
plt.yticks(np.arange(0, 9, 1), fontsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.legend(fontsize=8, shadow=True, ncol=1, handlelength=0.75, handleheight=0.75)

ax3 = fig.add_subplot(3, 2, 3)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
plt_clim1 = plt.bar(time1, ivs_rcm_eneb_pre, color='blue', label='Reg', width = 0.25, edgecolor='black', linewidth=1)
plt_clim2 = plt.bar(time1 + .25, ivs_gcm_eneb_pre,  color='red', label='Had', width = 0.25, edgecolor='black', linewidth=1)
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Value of IVS', fontsize=8)
plt.ylim(0, 8)
plt.xticks(time1 + .12, ('PRCPTOT', 'R95p', 'R99p', 'Rx1day', 'Rx5day', 'SDII', 'CDD', 'CWD', 'R10mm', 'R20mm'), fontsize=8)
plt.yticks(np.arange(0, 9, 1), fontsize=8)
plt.setp(ax3.get_xticklabels(), visible=False)

ax4 = fig.add_subplot(3, 2, 4)
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
plt_clim1 = plt.bar(time2, ivs_rcm_eneb_tas, color='blue', label='Reg', width = 0.25, edgecolor='black', linewidth=1)
plt_clim2 = plt.bar(time2 + .25, ivs_gcm_eneb_tas,  color='red', label='Had', width = 0.25, edgecolor='black', linewidth=1)
plt.title(u'E)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Value of IVS', fontsize=8)
plt.ylim(0, 8)
plt.xticks(time2 + .12, ('TXX', 'TXn', 'TNx', 'TNn', 'DTR', 'SU', 'TR', 'Tx10p', 'Tx90p', 'Tn10p', 'Tn90p'), fontsize=8)
plt.yticks(np.arange(0, 9, 1), fontsize=8)
plt.setp(ax4.get_xticklabels(), visible=False)

ax5 = fig.add_subplot(3, 2, 5)
ax5.spines['right'].set_visible(False)
ax5.spines['top'].set_visible(False)
plt_clim1 = plt.bar(time1, ivs_rcm_matopiba_pre, color='blue', label='Reg', width = 0.25, edgecolor='black', linewidth=1)
plt_clim2 = plt.bar(time1 + .25, ivs_gcm_matopiba_pre,  color='red', label='Had', width = 0.25, edgecolor='black', linewidth=1)
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.ylim(0, 8)
plt.xticks(time1 + .12, ('PRCPTOT', 'R95p', 'R99p', 'Rx1day', 'Rx5day', 'SDII', 'CDD', 'CWD', 'R10mm', 'R20mm'), fontsize=8)
plt.yticks(np.arange(0, 9, 1), fontsize=8)
labels = ax5.get_xticklabels()
plt.setp(labels, rotation=90)

ax6 = fig.add_subplot(3, 2, 6)
ax6.spines['right'].set_visible(False)
ax6.spines['top'].set_visible(False)
plt_clim1 = plt.bar(time2, ivs_rcm_matopiba_tas, color='blue', label='Reg', width = 0.25, edgecolor='black', linewidth=1)
plt_clim2 = plt.bar(time2 + .25, ivs_gcm_matopiba_tas,  color='red', label='Had', width = 0.25, edgecolor='black', linewidth=1)
plt.title(u'F)', loc='left', fontweight='bold', fontsize=8)
plt.ylim(0, 8)
plt.xticks(time2 + .12, ('TXX', 'TXn', 'TNx', 'TNn', 'DTR', 'SU', 'TR', 'Tx10p', 'Tx90p', 'Tn10p', 'Tn90p'), fontsize=8)
plt.yticks(np.arange(0, 9, 1), fontsize=8)
labels = ax6.get_xticklabels()
plt.setp(labels, rotation=90)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_ivs_etccdi_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()	

