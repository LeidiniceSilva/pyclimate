# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "06ssss/24/2021"
__description__ = "This script plot correlation matrix from ETCCDI indices"

import os
import netCDF4
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from matplotlib.font_manager import FontProperties


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
	annual_rcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	annual_rcm = np.nanmean(annual_rcm)

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
	annual_gcm = np.nanmean(annual_gcm)

	return annual_gcm
	

# Import regcm and hadgem models
# Precipitation
# SAMZ
# RCM
rcm_prcptot_hist_samz = import_rcm('eca_prcptot', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_r95p_hist_samz = import_rcm('eca_r95p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_r99p_hist_samz = import_rcm('eca_r99p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_rx1day_hist_samz = import_rcm('eca_rx1day', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_rx5day_hist_samz = import_rcm('eca_rx5day', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_sdii_hist_samz = import_rcm('eca_sdii', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_cdd_hist_samz = import_rcm('eca_cdd', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_cwd_hist_samz = import_rcm('eca_cwd', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_r10mm_hist_samz = import_rcm('eca_r10mm', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_r20mm_hist_samz = import_rcm('eca_r20mm', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')

rcm_prcptot_rcp26_samz = import_rcm('eca_prcptot', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_r95p_rcp26_samz = import_rcm('eca_r95p', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_r99p_rcp26_samz = import_rcm('eca_r99p', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_rx1day_rcp26_samz = import_rcm('eca_rx1day', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_rx5day_rcp26_samz = import_rcm('eca_rx5day', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_sdii_rcp26_samz = import_rcm('eca_sdii', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_cdd_rcp26_samz = import_rcm('eca_cdd', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_cwd_rcp26_samz = import_rcm('eca_cwd', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_r10mm_rcp26_samz = import_rcm('eca_r10mm', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_r20mm_rcp26_samz = import_rcm('eca_r20mm', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')

rcm_prcptot_rcp85_samz = import_rcm('eca_prcptot', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_r95p_rcp85_samz = import_rcm('eca_r95p', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_r99p_rcp85_samz = import_rcm('eca_r99p', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_rx1day_rcp85_samz = import_rcm('eca_rx1day', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_rx5day_rcp85_samz = import_rcm('eca_rx5day', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_sdii_rcp85_samz = import_rcm('eca_sdii', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_cdd_rcp85_samz = import_rcm('eca_cdd', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_cwd_rcp85_samz = import_rcm('eca_cwd', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_r10mm_rcp85_samz = import_rcm('eca_r10mm', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_r20mm_rcp85_samz = import_rcm('eca_r20mm', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')

# GCM
gcm_prcptot_hist_samz = import_gcm('eca_prcptot', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r95p_hist_samz = import_gcm('eca_r95p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r99p_hist_samz = import_gcm('eca_r99p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_rx1day_hist_samz = import_gcm('eca_rx1day', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_rx5day_hist_samz = import_gcm('eca_rx5day', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_sdii_hist_samz = import_gcm('eca_sdii', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_cdd_hist_samz = import_gcm('eca_cdd', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_cwd_hist_samz = import_gcm('eca_cwd', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r10mm_hist_samz = import_gcm('eca_r10mm', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r20mm_hist_samz = import_gcm('eca_r20mm', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

gcm_prcptot_rcp26_samz = import_gcm('eca_prcptot', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_r95p_rcp26_samz = import_gcm('eca_r95p', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_r99p_rcp26_samz = import_gcm('eca_r99p', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_rx1day_rcp26_samz = import_gcm('eca_rx1day', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_rx5day_rcp26_samz = import_gcm('eca_rx5day', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_sdii_rcp26_samz = import_gcm('eca_sdii', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_cdd_rcp26_samz = import_gcm('eca_cdd', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_cwd_rcp26_samz = import_gcm('eca_cwd', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_r10mm_rcp26_samz = import_gcm('eca_r10mm', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_r20mm_rcp26_samz = import_gcm('eca_r20mm', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')

gcm_prcptot_rcp85_samz = import_gcm('eca_prcptot', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r95p_rcp85_samz = import_gcm('eca_r95p', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r99p_rcp85_samz = import_gcm('eca_r99p', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_rx1day_rcp85_samz = import_gcm('eca_rx1day', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_rx5day_rcp85_samz = import_gcm('eca_rx5day', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_sdii_rcp85_samz = import_gcm('eca_sdii', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_cdd_rcp85_samz = import_gcm('eca_cdd', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_cwd_rcp85_samz = import_gcm('eca_cwd', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r10mm_rcp85_samz = import_gcm('eca_r10mm', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r20mm_rcp85_samz= import_gcm('eca_r20mm', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')

# ENEB
# RCM
rcm_prcptot_hist_eneb = import_rcm('eca_prcptot', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_r95p_hist_eneb = import_rcm('eca_r95p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_r99p_hist_eneb = import_rcm('eca_r99p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_rx1day_hist_eneb = import_rcm('eca_rx1day', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_rx5day_hist_eneb = import_rcm('eca_rx5day', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_sdii_hist_eneb = import_rcm('eca_sdii', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_cdd_hist_eneb = import_rcm('eca_cdd', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_cwd_hist_eneb = import_rcm('eca_cwd', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_r10mm_hist_eneb = import_rcm('eca_r10mm', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_r20mm_hist_eneb = import_rcm('eca_r20mm', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')

rcm_prcptot_rcp26_eneb = import_rcm('eca_prcptot', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_r95p_rcp26_eneb = import_rcm('eca_r95p', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_r99p_rcp26_eneb = import_rcm('eca_r99p', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_rx1day_rcp26_eneb = import_rcm('eca_rx1day', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_rx5day_rcp26_eneb = import_rcm('eca_rx5day', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_sdii_rcp26_eneb = import_rcm('eca_sdii', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_cdd_rcp26_eneb = import_rcm('eca_cdd', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_cwd_rcp26_eneb = import_rcm('eca_cwd', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_r10mm_rcp26_eneb = import_rcm('eca_r10mm', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_r20mm_rcp26_eneb = import_rcm('eca_r20mm', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')

rcm_prcptot_rcp85_eneb = import_rcm('eca_prcptot', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_r95p_rcp85_eneb = import_rcm('eca_r95p', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_r99p_rcp85_eneb = import_rcm('eca_r99p', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_rx1day_rcp85_eneb = import_rcm('eca_rx1day', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_rx5day_rcp85_eneb = import_rcm('eca_rx5day', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_sdii_rcp85_eneb = import_rcm('eca_sdii', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_cdd_rcp85_eneb = import_rcm('eca_cdd', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_cwd_rcp85_eneb = import_rcm('eca_cwd', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_r10mm_rcp85_eneb = import_rcm('eca_r10mm', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_r20mm_rcp85_eneb = import_rcm('eca_r20mm', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')

# GCM
gcm_prcptot_hist_eneb = import_gcm('eca_prcptot', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r95p_hist_eneb = import_gcm('eca_r95p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r99p_hist_eneb = import_gcm('eca_r99p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_rx1day_hist_eneb = import_gcm('eca_rx1day', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_rx5day_hist_eneb = import_gcm('eca_rx5day', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_sdii_hist_eneb = import_gcm('eca_sdii', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_cdd_hist_eneb = import_gcm('eca_cdd', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_cwd_hist_eneb = import_gcm('eca_cwd', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r10mm_hist_eneb = import_gcm('eca_r10mm', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r20mm_hist_eneb = import_gcm('eca_r20mm', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

gcm_prcptot_rcp26_eneb = import_gcm('eca_prcptot', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_r95p_rcp26_eneb = import_gcm('eca_r95p', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_r99p_rcp26_eneb = import_gcm('eca_r99p', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_rx1day_rcp26_eneb = import_gcm('eca_rx1day', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_rx5day_rcp26_eneb = import_gcm('eca_rx5day', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_sdii_rcp26_eneb = import_gcm('eca_sdii', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_cdd_rcp26_eneb = import_gcm('eca_cdd', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_cwd_rcp26_eneb = import_gcm('eca_cwd', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_r10mm_rcp26_eneb = import_gcm('eca_r10mm', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_r20mm_rcp26_eneb = import_gcm('eca_r20mm', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')

gcm_prcptot_rcp85_eneb = import_gcm('eca_prcptot', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r95p_rcp85_eneb = import_gcm('eca_r95p', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r99p_rcp85_eneb = import_gcm('eca_r99p', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_rx1day_rcp85_eneb = import_gcm('eca_rx1day', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_rx5day_rcp85_eneb = import_gcm('eca_rx5day', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_sdii_rcp85_eneb = import_gcm('eca_sdii', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_cdd_rcp85_eneb = import_gcm('eca_cdd', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_cwd_rcp85_eneb = import_gcm('eca_cwd', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r10mm_rcp85_eneb = import_gcm('eca_r10mm', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r20mm_rcp85_eneb = import_gcm('eca_r20mm', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')

# MATOPIBA
# RCM
rcm_prcptot_hist_matopiba = import_rcm('eca_prcptot', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_r95p_hist_matopiba = import_rcm('eca_r95p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_r99p_hist_matopiba = import_rcm('eca_r99p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_rx1day_hist_matopiba = import_rcm('eca_rx1day', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_rx5day_hist_matopiba = import_rcm('eca_rx5day', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_sdii_hist_matopiba = import_rcm('eca_sdii', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_cdd_hist_matopiba = import_rcm('eca_cdd', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_cwd_hist_matopiba = import_rcm('eca_cwd', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_r10mm_hist_matopiba = import_rcm('eca_r10mm', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_r20mm_hist_matopiba = import_rcm('eca_r20mm', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')

rcm_prcptot_rcp26_matopiba = import_rcm('eca_prcptot', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_r95p_rcp26_matopiba = import_rcm('eca_r95p', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_r99p_rcp26_matopiba = import_rcm('eca_r99p', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_rx1day_rcp26_matopiba = import_rcm('eca_rx1day', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_rx5day_rcp26_matopiba = import_rcm('eca_rx5day', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_sdii_rcp26_matopiba = import_rcm('eca_sdii', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_cdd_rcp26_matopiba = import_rcm('eca_cdd', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_cwd_rcp26_matopiba = import_rcm('eca_cwd', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_r10mm_rcp26_matopiba = import_rcm('eca_r10mm', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_r20mm_rcp26_matopiba = import_rcm('eca_r20mm', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')

rcm_prcptot_rcp85_matopiba = import_rcm('eca_prcptot', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_r95p_rcp85_matopiba = import_rcm('eca_r95p', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_r99p_rcp85_matopiba = import_rcm('eca_r99p', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_rx1day_rcp85_matopiba = import_rcm('eca_rx1day', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_rx5day_rcp85_matopiba = import_rcm('eca_rx5day', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_sdii_rcp85_matopiba = import_rcm('eca_sdii', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_cdd_rcp85_matopiba = import_rcm('eca_cdd', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_cwd_rcp85_matopiba = import_rcm('eca_cwd', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_r10mm_rcp85_matopiba = import_rcm('eca_r10mm', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_r20mm_rcp85_matopiba = import_rcm('eca_r20mm', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')

# GCM
gcm_prcptot_hist_matopiba = import_gcm('eca_prcptot', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r95p_hist_matopiba = import_gcm('eca_r95p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r99p_hist_matopiba = import_gcm('eca_r99p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_rx1day_hist_matopiba = import_gcm('eca_rx1day', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_rx5day_hist_matopiba = import_gcm('eca_rx5day', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_sdii_hist_matopiba = import_gcm('eca_sdii', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_cdd_hist_matopiba = import_gcm('eca_cdd', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_cwd_hist_matopiba = import_gcm('eca_cwd', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r10mm_hist_matopiba = import_gcm('eca_r10mm', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r20mm_hist_matopiba = import_gcm('eca_r20mm', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

gcm_prcptot_rcp26_matopiba = import_gcm('eca_prcptot', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_r95p_rcp26_matopiba = import_gcm('eca_r95p', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_r99p_rcp26_matopiba = import_gcm('eca_r99p', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_rx1day_rcp26_matopiba = import_gcm('eca_rx1day', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_rx5day_rcp26_matopiba = import_gcm('eca_rx5day', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_sdii_rcp26_matopiba = import_gcm('eca_sdii', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_cdd_rcp26_matopiba = import_gcm('eca_cdd', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_cwd_rcp26_matopiba = import_gcm('eca_cwd', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_r10mm_rcp26_matopiba = import_gcm('eca_r10mm', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_r20mm_rcp26_matopiba = import_gcm('eca_r20mm', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')

gcm_prcptot_rcp85_matopiba = import_gcm('eca_prcptot', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r95p_rcp85_matopiba = import_gcm('eca_r95p', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r99p_rcp85_matopiba = import_gcm('eca_r99p', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_rx1day_rcp85_matopiba = import_gcm('eca_rx1day', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_rx5day_rcp85_matopiba = import_gcm('eca_rx5day', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_sdii_rcp85_matopiba = import_gcm('eca_sdii', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_cdd_rcp85_matopiba = import_gcm('eca_cdd', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_cwd_rcp85_matopiba = import_gcm('eca_cwd', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r10mm_rcp85_matopiba = import_gcm('eca_r10mm', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r20mm_rcp85_matopiba = import_gcm('eca_r20mm', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')

# Temperature
# SAMZ
# RCM
rcm_txx_hist_samz = import_rcm('eca_txx', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_txn_hist_samz = import_rcm('eca_txn', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tnx_hist_samz = import_rcm('eca_tnx', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tnn_hist_samz = import_rcm('eca_tnn', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_dtr_hist_samz = import_rcm('eca_dtr', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_su_hist_samz = import_rcm('eca_su', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tr_hist_samz = import_rcm('eca_tr', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tx10p_hist_samz = import_rcm('eca_tx10p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tx90p_hist_samz = import_rcm('eca_tx90p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tn10p_hist_samz = import_rcm('eca_tn10p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tn90p_hist_samz = import_rcm('eca_tn90p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')

rcm_txx_rcp26_samz = import_rcm('eca_txx', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_txn_rcp26_samz = import_rcm('eca_txn', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tnx_rcp26_samz = import_rcm('eca_tnx', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tnn_rcp26_samz = import_rcm('eca_tnn', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_dtr_rcp26_samz = import_rcm('eca_dtr', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_su_rcp26_samz = import_rcm('eca_su', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tr_rcp26_samz = import_rcm('eca_tr', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tx10p_rcp26_samz = import_rcm('eca_tx10p', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tx90p_rcp26_samz = import_rcm('eca_tx90p', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tn10p_rcp26_samz = import_rcm('eca_tn10p', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tn90p_rcp26_samz = import_rcm('eca_tn90p', 'samz', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')

rcm_txx_rcp85_samz = import_rcm('eca_txx', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_txn_rcp85_samz = import_rcm('eca_txn', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tnx_rcp85_samz = import_rcm('eca_tnx', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tnn_rcp85_samz = import_rcm('eca_tnn', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_dtr_rcp85_samz = import_rcm('eca_dtr', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_su_rcp85_samz = import_rcm('eca_su', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tr_rcp85_samz = import_rcm('eca_tr', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tx10p_rcp85_samz = import_rcm('eca_tx10p', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tx90p_rcp85_samz = import_rcm('eca_tx90p', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tn10p_rcp85_samz = import_rcm('eca_tn10p', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tn90p_rcp85_samz = import_rcm('eca_tn90p', 'samz', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')

# GCM
gcm_txx_hist_samz = import_gcm('eca_txx', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_txn_hist_samz = import_gcm('eca_txn', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tnx_hist_samz = import_gcm('eca_tnx', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tnn_hist_samz = import_gcm('eca_tnn', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_dtr_hist_samz = import_gcm('eca_dtr', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_su_hist_samz = import_gcm('eca_su', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tr_hist_samz = import_gcm('eca_tr', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tx10p_hist_samz = import_gcm('eca_tx10p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tx90p_hist_samz = import_gcm('eca_tx90p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tn10p_hist_samz = import_gcm('eca_tn10p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tn90p_hist_samz = import_gcm('eca_tn90p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

gcm_txx_rcp26_samz = import_gcm('eca_txx', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_txn_rcp26_samz = import_gcm('eca_txn', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tnx_rcp26_samz = import_gcm('eca_tnx', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tnn_rcp26_samz = import_gcm('eca_tnn', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_dtr_rcp26_samz = import_gcm('eca_dtr', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_su_rcp26_samz = import_gcm('eca_su', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tr_rcp26_samz = import_gcm('eca_tr', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tx10p_rcp26_samz = import_gcm('eca_tx10p', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tx90p_rcp26_samz = import_gcm('eca_tx90p', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tn10p_rcp26_samz = import_gcm('eca_tn10p', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tn90p_rcp26_samz = import_gcm('eca_tn90p', 'samz', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')

gcm_txx_rcp85_samz = import_gcm('eca_txx', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_txn_rcp85_samz = import_gcm('eca_txn', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tnx_rcp85_samz = import_gcm('eca_tnx', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tnn_rcp85_samz = import_gcm('eca_tnn', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_dtr_rcp85_samz = import_gcm('eca_dtr', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_su_rcp85_samz = import_gcm('eca_su', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tr_rcp85_samz = import_gcm('eca_tr', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tx10p_rcp85_samz = import_gcm('eca_tx10p', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tx90p_rcp85_samz = import_gcm('eca_tx90p', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tn10p_rcp85_samz = import_gcm('eca_tn10p', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tn90p_rcp85_samz = import_gcm('eca_tn90p', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')

# ENEB
# RCM
rcm_txx_hist_eneb = import_rcm('eca_txx', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_txn_hist_eneb = import_rcm('eca_txn', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tnx_hist_eneb = import_rcm('eca_tnx', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tnn_hist_eneb = import_rcm('eca_tnn', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_dtr_hist_eneb = import_rcm('eca_dtr', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_su_hist_eneb = import_rcm('eca_su', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tr_hist_eneb = import_rcm('eca_tr', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tx10p_hist_eneb = import_rcm('eca_tx10p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tx90p_hist_eneb = import_rcm('eca_tx90p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tn10p_hist_eneb = import_rcm('eca_tn10p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tn90p_hist_eneb = import_rcm('eca_tn90p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')

rcm_txx_rcp26_eneb = import_rcm('eca_txx', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_txn_rcp26_eneb = import_rcm('eca_txn', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tnx_rcp26_eneb = import_rcm('eca_tnx', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tnn_rcp26_eneb = import_rcm('eca_tnn', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_dtr_rcp26_eneb = import_rcm('eca_dtr', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_su_rcp26_eneb = import_rcm('eca_su', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tr_rcp26_eneb = import_rcm('eca_tr', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tx10p_rcp26_eneb = import_rcm('eca_tx10p', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tx90p_rcp26_eneb = import_rcm('eca_tx90p', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tn10p_rcp26_eneb = import_rcm('eca_tn10p', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tn90p_rcp26_eneb = import_rcm('eca_tn90p', 'eneb', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')

rcm_txx_rcp85_eneb = import_rcm('eca_txx', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_txn_rcp85_eneb = import_rcm('eca_txn', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tnx_rcp85_eneb = import_rcm('eca_tnx', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tnn_rcp85_eneb = import_rcm('eca_tnn', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_dtr_rcp85_eneb = import_rcm('eca_dtr', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_su_rcp85_eneb = import_rcm('eca_su', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tr_rcp85_eneb = import_rcm('eca_tr', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tx10p_rcp85_eneb = import_rcm('eca_tx10p', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tx90p_rcp85_eneb = import_rcm('eca_tx90p', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tn10p_rcp85_eneb = import_rcm('eca_tn10p', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tn90p_rcp85_eneb = import_rcm('eca_tn90p', 'eneb', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')

# GCM
gcm_txx_hist_eneb = import_gcm('eca_txx', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_txn_hist_eneb = import_gcm('eca_txn', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tnx_hist_eneb = import_gcm('eca_tnx', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tnn_hist_eneb = import_gcm('eca_tnn', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_dtr_hist_eneb = import_gcm('eca_dtr', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_su_hist_eneb = import_gcm('eca_su', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tr_hist_eneb = import_gcm('eca_tr', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tx10p_hist_eneb = import_gcm('eca_tx10p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tx90p_hist_eneb = import_gcm('eca_tx90p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tn10p_hist_eneb = import_gcm('eca_tn10p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tn90p_hist_eneb = import_gcm('eca_tn90p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

gcm_txx_rcp26_eneb = import_gcm('eca_txx', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_txn_rcp26_eneb = import_gcm('eca_txn', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tnx_rcp26_eneb = import_gcm('eca_tnx', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tnn_rcp26_eneb = import_gcm('eca_tnn', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_dtr_rcp26_eneb = import_gcm('eca_dtr', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_su_rcp26_eneb = import_gcm('eca_su', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tr_rcp26_eneb = import_gcm('eca_tr', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tx10p_rcp26_eneb = import_gcm('eca_tx10p', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tx90p_rcp26_eneb = import_gcm('eca_tx90p', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tn10p_rcp26_eneb = import_gcm('eca_tn10p', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tn90p_rcp26_eneb = import_gcm('eca_tn90p', 'eneb', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')

gcm_txx_rcp85_eneb = import_gcm('eca_txx', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_txn_rcp85_eneb = import_gcm('eca_txn', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tnx_rcp85_eneb = import_gcm('eca_tnx', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tnn_rcp85_eneb = import_gcm('eca_tnn', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_dtr_rcp85_eneb = import_gcm('eca_dtr', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_su_rcp85_eneb = import_gcm('eca_su', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tr_rcp85_eneb = import_gcm('eca_tr', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tx10p_rcp85_eneb = import_gcm('eca_tx10p', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tx90p_rcp85_eneb = import_gcm('eca_tx90p', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tn10p_rcp85_eneb = import_gcm('eca_tn10p', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tn90p_rcp85_eneb = import_gcm('eca_tn90p', 'eneb', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')

# MATOPIBA
# RCM
rcm_txx_hist_matopiba = import_rcm('eca_txx', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_txn_hist_matopiba = import_rcm('eca_txn', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tnx_hist_matopiba = import_rcm('eca_tnx', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tnn_hist_matopiba = import_rcm('eca_tnn', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_dtr_hist_matopiba = import_rcm('eca_dtr', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_su_hist_matopiba = import_rcm('eca_su', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tr_hist_matopiba = import_rcm('eca_tr', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tx10p_hist_matopiba = import_rcm('eca_tx10p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tx90p_hist_matopiba = import_rcm('eca_tx90p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tn10p_hist_matopiba = import_rcm('eca_tn10p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
rcm_tn90p_hist_matopiba = import_rcm('eca_tn90p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')

rcm_txx_rcp26_matopiba = import_rcm('eca_txx', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_txn_rcp26_matopiba = import_rcm('eca_txn', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tnx_rcp26_matopiba = import_rcm('eca_tnx', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tnn_rcp26_matopiba = import_rcm('eca_tnn', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_dtr_rcp26_matopiba = import_rcm('eca_dtr', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_su_rcp26_matopiba = import_rcm('eca_su', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tr_rcp26_matopiba = import_rcm('eca_tr', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tx10p_rcp26_matopiba = import_rcm('eca_tx10p', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tx90p_rcp26_matopiba = import_rcm('eca_tx90p', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tn10p_rcp26_matopiba = import_rcm('eca_tn10p', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')
rcm_tn90p_rcp26_matopiba = import_rcm('eca_tn90p', 'matopiba', 'RegCM47_had', 'rcp26', 'yr', '2080-2099')

rcm_txx_rcp85_matopiba = import_rcm('eca_txx', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_txn_rcp85_matopiba = import_rcm('eca_txn', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tnx_rcp85_matopiba = import_rcm('eca_tnx', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tnn_rcp85_matopiba = import_rcm('eca_tnn', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_dtr_rcp85_matopiba= import_rcm('eca_dtr', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_su_rcp85_matopiba = import_rcm('eca_su', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tr_rcp85_matopiba = import_rcm('eca_tr', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tx10p_rcp85_matopiba = import_rcm('eca_tx10p', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tx90p_rcp85_matopiba = import_rcm('eca_tx90p', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tn10p_rcp85_matopiba = import_rcm('eca_tn10p', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')
rcm_tn90p_rcp85_matopiba = import_rcm('eca_tn90p', 'matopiba', 'RegCM47_had', 'rcp85', 'yr', '2080-2099')

# GCM
gcm_txx_hist_matopiba = import_gcm('eca_txx', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_txn_hist_matopiba = import_gcm('eca_txn', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tnx_hist_matopiba = import_gcm('eca_tnx', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tnn_hist_matopiba = import_gcm('eca_tnn', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_dtr_hist_matopiba = import_gcm('eca_dtr', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_su_hist_matopiba = import_gcm('eca_su', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tr_hist_matopiba = import_gcm('eca_tr', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tx10p_hist_matopiba = import_gcm('eca_tx10p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tx90p_hist_matopiba = import_gcm('eca_tx90p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tn10p_hist_matopiba = import_gcm('eca_tn10p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_tn90p_hist_matopiba = import_gcm('eca_tn90p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

gcm_txx_rcp26_matopiba = import_gcm('eca_txx', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_txn_rcp26_matopiba = import_gcm('eca_txn', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tnx_rcp26_matopiba = import_gcm('eca_tnx', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tnn_rcp26_matopiba = import_gcm('eca_tnn', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_dtr_rcp26_matopiba = import_gcm('eca_dtr', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_su_rcp26_matopiba = import_gcm('eca_su', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tr_rcp26_matopiba = import_gcm('eca_tr', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tx10p_rcp26_matopiba = import_gcm('eca_tx10p', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tx90p_rcp26_matopiba = import_gcm('eca_tx90p', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tn10p_rcp26_matopiba = import_gcm('eca_tn10p', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')
gcm_tn90p_rcp26_matopiba = import_gcm('eca_tn90p', 'matopiba', 'HadGEM2-ES', 'rcp26', 'yr', '2080-2099')

gcm_txx_rcp85_matopiba = import_gcm('eca_txx', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_txn_rcp85_matopiba = import_gcm('eca_txn', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tnx_rcp85_matopiba = import_gcm('eca_tnx', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tnn_rcp85_matopiba = import_gcm('eca_tnn', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_dtr_rcp85_matopiba = import_gcm('eca_dtr', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_su_rcp85_matopiba = import_gcm('eca_su', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tr_rcp85_matopiba = import_gcm('eca_tx10p', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tx10p_rcp85_matopiba = import_gcm('eca_tx90p', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tx90p_rcp85_matopiba = import_gcm('eca_tx90p', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tn10p_rcp85_matopiba = import_gcm('eca_tn10p', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_tn90p_rcp85_matopiba = import_gcm('eca_tn90p', 'matopiba', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')

# Compute diference between FF and RF
# Precipitation
# SAMZ
diff_rcm_rcp26_prcptot_samz = rcm_prcptot_rcp26_samz - rcm_prcptot_hist_samz
diff_rcm_rcp26_r95p_samz = rcm_r95p_rcp26_samz - rcm_r95p_hist_samz
diff_rcm_rcp26_r99p_samz = rcm_r99p_rcp26_samz - rcm_r99p_hist_samz
diff_rcm_rcp26_rx1day_samz = rcm_rx1day_rcp26_samz - rcm_rx1day_hist_samz
diff_rcm_rcp26_rx5day_samz = rcm_rx5day_rcp26_samz - rcm_rx5day_hist_samz
diff_rcm_rcp26_sdii_samz = rcm_sdii_rcp26_samz - rcm_sdii_hist_samz
diff_rcm_rcp26_cdd_samz = rcm_cdd_rcp26_samz - rcm_cdd_hist_samz
diff_rcm_rcp26_cwd_samz = rcm_cwd_rcp26_samz - rcm_cwd_hist_samz
diff_rcm_rcp26_r10mm_samz = rcm_r10mm_rcp26_samz - rcm_r10mm_hist_samz
diff_rcm_rcp26_r20mm_samz = rcm_r20mm_rcp26_samz - rcm_r20mm_hist_samz

diff_rcm_rcp26_samz_pre = [diff_rcm_rcp26_prcptot_samz, diff_rcm_rcp26_r95p_samz, diff_rcm_rcp26_r99p_samz, diff_rcm_rcp26_rx1day_samz,
diff_rcm_rcp26_rx5day_samz, diff_rcm_rcp26_sdii_samz, diff_rcm_rcp26_cdd_samz, diff_rcm_rcp26_cwd_samz, diff_rcm_rcp26_r10mm_samz,
diff_rcm_rcp26_r20mm_samz]

diff_rcm_rcp85_prcptot_samz = rcm_prcptot_rcp85_samz - rcm_prcptot_hist_samz
diff_rcm_rcp85_r95p_samz = rcm_r95p_rcp85_samz - rcm_r95p_hist_samz
diff_rcm_rcp85_r99p_samz = rcm_r99p_rcp85_samz - rcm_r99p_hist_samz
diff_rcm_rcp85_rx1day_samz = rcm_rx1day_rcp85_samz - rcm_rx1day_hist_samz
diff_rcm_rcp85_rx5day_samz = rcm_rx5day_rcp85_samz - rcm_rx5day_hist_samz
diff_rcm_rcp85_sdii_samz = rcm_sdii_rcp85_samz - rcm_sdii_hist_samz
diff_rcm_rcp85_cdd_samz = rcm_cdd_rcp85_samz - rcm_cdd_hist_samz
diff_rcm_rcp85_cwd_samz = rcm_cwd_rcp85_samz - rcm_cwd_hist_samz
diff_rcm_rcp85_r10mm_samz = rcm_r10mm_rcp85_samz - rcm_r10mm_hist_samz
diff_rcm_rcp85_r20mm_samz = rcm_r20mm_rcp85_samz - rcm_r20mm_hist_samz

diff_rcm_rcp85_samz_pre = [diff_rcm_rcp85_prcptot_samz, diff_rcm_rcp85_r95p_samz, diff_rcm_rcp85_r99p_samz, diff_rcm_rcp85_rx1day_samz,
diff_rcm_rcp85_rx5day_samz, diff_rcm_rcp85_sdii_samz, diff_rcm_rcp85_cdd_samz, diff_rcm_rcp85_cwd_samz, diff_rcm_rcp85_r10mm_samz,
diff_rcm_rcp85_r20mm_samz]

diff_gcm_rcp26_prcptot_samz = gcm_prcptot_rcp26_samz - gcm_prcptot_hist_samz
diff_gcm_rcp26_r95p_samz = gcm_r95p_rcp26_samz - gcm_r95p_hist_samz
diff_gcm_rcp26_r99p_samz = gcm_r99p_rcp26_samz - gcm_r99p_hist_samz
diff_gcm_rcp26_rx1day_samz = gcm_rx1day_rcp26_samz - gcm_rx1day_hist_samz
diff_gcm_rcp26_rx5day_samz = gcm_rx5day_rcp26_samz - gcm_rx5day_hist_samz
diff_gcm_rcp26_sdii_samz = gcm_sdii_rcp26_samz - gcm_sdii_hist_samz
diff_gcm_rcp26_cdd_samz = gcm_cdd_rcp26_samz - gcm_cdd_hist_samz
diff_gcm_rcp26_cwd_samz = gcm_cwd_rcp26_samz - gcm_cwd_hist_samz
diff_gcm_rcp26_r10mm_samz = gcm_r10mm_rcp26_samz - gcm_r10mm_hist_samz
diff_gcm_rcp26_r20mm_samz = gcm_r20mm_rcp26_samz - gcm_r20mm_hist_samz

diff_gcm_rcp26_samz_pre = [diff_gcm_rcp26_prcptot_samz, diff_gcm_rcp26_r95p_samz, diff_gcm_rcp26_r99p_samz, diff_gcm_rcp26_rx1day_samz,
diff_gcm_rcp26_rx5day_samz, diff_gcm_rcp26_sdii_samz, diff_gcm_rcp26_cdd_samz, diff_gcm_rcp26_cwd_samz, diff_gcm_rcp26_r10mm_samz,
diff_gcm_rcp26_r20mm_samz]

diff_gcm_rcp85_prcptot_samz = gcm_prcptot_rcp85_samz - gcm_prcptot_hist_samz
diff_gcm_rcp85_r95p_samz = gcm_r95p_rcp85_samz - gcm_r95p_hist_samz
diff_gcm_rcp85_r99p_samz = gcm_r99p_rcp85_samz - gcm_r99p_hist_samz
diff_gcm_rcp85_rx1day_samz = gcm_rx1day_rcp85_samz - gcm_rx1day_hist_samz
diff_gcm_rcp85_rx5day_samz = gcm_rx5day_rcp85_samz - gcm_rx5day_hist_samz
diff_gcm_rcp85_sdii_samz = gcm_sdii_rcp85_samz - gcm_sdii_hist_samz
diff_gcm_rcp85_cdd_samz = gcm_cdd_rcp85_samz - gcm_cdd_hist_samz
diff_gcm_rcp85_cwd_samz = gcm_cwd_rcp85_samz - gcm_cwd_hist_samz
diff_gcm_rcp85_r10mm_samz = gcm_r10mm_rcp85_samz - gcm_r10mm_hist_samz
diff_gcm_rcp85_r20mm_samz = gcm_r20mm_rcp85_samz - gcm_r20mm_hist_samz

diff_gcm_rcp85_samz_pre = [diff_gcm_rcp85_prcptot_samz, diff_gcm_rcp85_r95p_samz, diff_gcm_rcp85_r99p_samz, diff_gcm_rcp85_rx1day_samz,
diff_gcm_rcp85_rx5day_samz, diff_gcm_rcp85_sdii_samz, diff_gcm_rcp85_cdd_samz, diff_gcm_rcp85_cwd_samz, diff_gcm_rcp85_r10mm_samz,
diff_gcm_rcp85_r20mm_samz]

# ENEB
diff_rcm_rcp26_prcptot_eneb = rcm_prcptot_rcp26_eneb - rcm_prcptot_hist_eneb
diff_rcm_rcp26_r95p_eneb = rcm_r95p_rcp26_eneb - rcm_r95p_hist_eneb
diff_rcm_rcp26_r99p_eneb = rcm_r99p_rcp26_eneb - rcm_r99p_hist_eneb
diff_rcm_rcp26_rx1day_eneb = rcm_rx1day_rcp26_eneb - rcm_rx1day_hist_eneb
diff_rcm_rcp26_rx5day_eneb = rcm_rx5day_rcp26_eneb - rcm_rx5day_hist_eneb
diff_rcm_rcp26_sdii_eneb = rcm_sdii_rcp26_eneb - rcm_sdii_hist_eneb
diff_rcm_rcp26_cdd_eneb = rcm_cdd_rcp26_eneb - rcm_cdd_hist_eneb
diff_rcm_rcp26_cwd_eneb = rcm_cwd_rcp26_eneb - rcm_cwd_hist_eneb
diff_rcm_rcp26_r10mm_eneb = rcm_r10mm_rcp26_eneb - rcm_r10mm_hist_eneb
diff_rcm_rcp26_r20mm_eneb = rcm_r20mm_rcp26_eneb - rcm_r20mm_hist_eneb

diff_rcm_rcp26_eneb_pre = [diff_rcm_rcp26_prcptot_eneb, diff_rcm_rcp26_r95p_eneb, diff_rcm_rcp26_r99p_eneb, diff_rcm_rcp26_rx1day_eneb,
diff_rcm_rcp26_rx5day_eneb, diff_rcm_rcp26_sdii_eneb, diff_rcm_rcp26_cdd_eneb, diff_rcm_rcp26_cwd_eneb, diff_rcm_rcp26_r10mm_eneb,
diff_rcm_rcp26_r20mm_eneb]

diff_rcm_rcp85_prcptot_eneb = rcm_prcptot_rcp85_eneb - rcm_prcptot_hist_eneb
diff_rcm_rcp85_r95p_eneb = rcm_r95p_rcp85_eneb - rcm_r95p_hist_eneb
diff_rcm_rcp85_r99p_eneb = rcm_r99p_rcp85_eneb - rcm_r99p_hist_eneb
diff_rcm_rcp85_rx1day_eneb = rcm_rx1day_rcp85_eneb - rcm_rx1day_hist_eneb
diff_rcm_rcp85_rx5day_eneb = rcm_rx5day_rcp85_eneb - rcm_rx5day_hist_eneb
diff_rcm_rcp85_sdii_eneb = rcm_sdii_rcp85_eneb - rcm_sdii_hist_eneb
diff_rcm_rcp85_cdd_eneb = rcm_cdd_rcp85_eneb - rcm_cdd_hist_eneb
diff_rcm_rcp85_cwd_eneb = rcm_cwd_rcp85_eneb - rcm_cwd_hist_eneb
diff_rcm_rcp85_r10mm_eneb = rcm_r10mm_rcp85_eneb - rcm_r10mm_hist_eneb
diff_rcm_rcp85_r20mm_eneb = rcm_r20mm_rcp85_eneb - rcm_r20mm_hist_eneb

diff_rcm_rcp85_eneb_pre = [diff_rcm_rcp85_prcptot_eneb, diff_rcm_rcp85_r95p_eneb, diff_rcm_rcp85_r99p_eneb, diff_rcm_rcp85_rx1day_eneb,
diff_rcm_rcp85_rx5day_eneb, diff_rcm_rcp85_sdii_eneb, diff_rcm_rcp85_cdd_eneb, diff_rcm_rcp85_cwd_eneb, diff_rcm_rcp85_r10mm_eneb,
diff_rcm_rcp85_r20mm_eneb]

diff_gcm_rcp26_prcptot_eneb = gcm_prcptot_rcp26_eneb - gcm_prcptot_hist_eneb
diff_gcm_rcp26_r95p_eneb = gcm_r95p_rcp26_eneb - gcm_r95p_hist_eneb
diff_gcm_rcp26_r99p_eneb = gcm_r99p_rcp26_eneb - gcm_r99p_hist_eneb
diff_gcm_rcp26_rx1day_eneb = gcm_rx1day_rcp26_eneb - gcm_rx1day_hist_eneb
diff_gcm_rcp26_rx5day_eneb = gcm_rx5day_rcp26_eneb - gcm_rx5day_hist_eneb
diff_gcm_rcp26_sdii_eneb = gcm_sdii_rcp26_eneb - gcm_sdii_hist_eneb
diff_gcm_rcp26_cdd_eneb = gcm_cdd_rcp26_eneb - gcm_cdd_hist_eneb
diff_gcm_rcp26_cwd_eneb = gcm_cwd_rcp26_eneb - gcm_cwd_hist_eneb
diff_gcm_rcp26_r10mm_eneb = gcm_r10mm_rcp26_eneb - gcm_r10mm_hist_eneb
diff_gcm_rcp26_r20mm_eneb = gcm_r20mm_rcp26_eneb - gcm_r20mm_hist_eneb

diff_gcm_rcp26_eneb_pre = [diff_gcm_rcp26_prcptot_eneb, diff_gcm_rcp26_r95p_eneb, diff_gcm_rcp26_r99p_eneb, diff_gcm_rcp26_rx1day_eneb,
diff_gcm_rcp26_rx5day_eneb, diff_gcm_rcp26_sdii_eneb, diff_gcm_rcp26_cdd_eneb, diff_gcm_rcp26_cwd_eneb, diff_gcm_rcp26_r10mm_eneb,
diff_gcm_rcp26_r20mm_eneb]

diff_gcm_rcp85_prcptot_eneb = gcm_prcptot_rcp85_eneb - gcm_prcptot_hist_eneb
diff_gcm_rcp85_r95p_eneb = gcm_r95p_rcp85_eneb - gcm_r95p_hist_eneb
diff_gcm_rcp85_r99p_eneb = gcm_r99p_rcp85_eneb - gcm_r99p_hist_eneb
diff_gcm_rcp85_rx1day_eneb = gcm_rx1day_rcp85_eneb - gcm_rx1day_hist_eneb
diff_gcm_rcp85_rx5day_eneb = gcm_rx5day_rcp85_eneb - gcm_rx5day_hist_eneb
diff_gcm_rcp85_sdii_eneb = gcm_sdii_rcp85_eneb - gcm_sdii_hist_eneb
diff_gcm_rcp85_cdd_eneb = gcm_cdd_rcp85_eneb - gcm_cdd_hist_eneb
diff_gcm_rcp85_cwd_eneb = gcm_cwd_rcp85_eneb - gcm_cwd_hist_eneb
diff_gcm_rcp85_r10mm_eneb = gcm_r10mm_rcp85_eneb - gcm_r10mm_hist_eneb
diff_gcm_rcp85_r20mm_eneb = gcm_r20mm_rcp85_eneb - gcm_r20mm_hist_eneb

diff_gcm_rcp85_eneb_pre = [diff_gcm_rcp85_prcptot_eneb, diff_gcm_rcp85_r95p_eneb, diff_gcm_rcp85_r99p_eneb, diff_gcm_rcp85_rx1day_eneb,
diff_gcm_rcp85_rx5day_eneb, diff_gcm_rcp85_sdii_eneb, diff_gcm_rcp85_cdd_eneb, diff_gcm_rcp85_cwd_eneb, diff_gcm_rcp85_r10mm_eneb,
diff_gcm_rcp85_r20mm_eneb]

# MATOPIBA
diff_rcm_rcp26_prcptot_matopiba = rcm_prcptot_rcp26_eneb - rcm_prcptot_hist_matopiba
diff_rcm_rcp26_r95p_matopiba = rcm_r95p_rcp26_eneb - rcm_r95p_hist_matopiba
diff_rcm_rcp26_r99p_matopiba = rcm_r99p_rcp26_eneb - rcm_r99p_hist_matopiba
diff_rcm_rcp26_rx1day_matopiba = rcm_rx1day_rcp26_eneb - rcm_rx1day_hist_matopiba
diff_rcm_rcp26_rx5day_matopiba = rcm_rx5day_rcp26_eneb - rcm_rx5day_hist_matopiba
diff_rcm_rcp26_sdii_matopiba = rcm_sdii_rcp26_eneb - rcm_sdii_hist_matopiba
diff_rcm_rcp26_cdd_matopiba = rcm_cdd_rcp26_eneb - rcm_cdd_hist_matopiba
diff_rcm_rcp26_cwd_matopiba = rcm_cwd_rcp26_eneb - rcm_cwd_hist_matopiba
diff_rcm_rcp26_r10mm_matopiba = rcm_r10mm_rcp26_eneb - rcm_r10mm_hist_matopiba
diff_rcm_rcp26_r20mm_matopiba = rcm_r20mm_rcp26_eneb - rcm_r20mm_hist_matopiba

diff_rcm_rcp26_matopiba_pre = [diff_rcm_rcp26_prcptot_matopiba, diff_rcm_rcp26_r95p_matopiba, diff_rcm_rcp26_r99p_matopiba,
diff_rcm_rcp26_rx1day_matopiba, diff_rcm_rcp26_rx5day_matopiba, diff_rcm_rcp26_sdii_matopiba, diff_rcm_rcp26_cdd_matopiba,
diff_rcm_rcp26_cwd_matopiba, diff_rcm_rcp26_r10mm_matopiba,
diff_rcm_rcp26_r20mm_matopiba]

diff_rcm_rcp85_prcptot_matopiba = rcm_prcptot_rcp85_matopiba - rcm_prcptot_hist_matopiba
diff_rcm_rcp85_r95p_matopiba = rcm_r95p_rcp85_matopiba - rcm_r95p_hist_matopiba
diff_rcm_rcp85_r99p_matopiba = rcm_r99p_rcp85_matopiba - rcm_r99p_hist_matopiba
diff_rcm_rcp85_rx1day_matopiba = rcm_rx1day_rcp85_matopiba - rcm_rx1day_hist_matopiba
diff_rcm_rcp85_rx5day_matopiba = rcm_rx5day_rcp85_matopiba - rcm_rx5day_hist_matopiba
diff_rcm_rcp85_sdii_matopiba = rcm_sdii_rcp85_matopiba - rcm_sdii_hist_matopiba
diff_rcm_rcp85_cdd_matopiba = rcm_cdd_rcp85_matopiba - rcm_cdd_hist_matopiba
diff_rcm_rcp85_cwd_matopiba = rcm_cwd_rcp85_matopiba - rcm_cwd_hist_matopiba
diff_rcm_rcp85_r10mm_matopiba = rcm_r10mm_rcp85_matopiba - rcm_r10mm_hist_matopiba
diff_rcm_rcp85_r20mm_matopiba = rcm_r20mm_rcp85_matopiba - rcm_r20mm_hist_matopiba

diff_rcm_rcp85_matopiba_pre = [diff_rcm_rcp85_prcptot_matopiba, diff_rcm_rcp85_r95p_matopiba, diff_rcm_rcp85_r99p_matopiba,
diff_rcm_rcp85_rx1day_matopiba, diff_rcm_rcp85_rx5day_matopiba, diff_rcm_rcp85_sdii_matopiba, diff_rcm_rcp85_cdd_matopiba,
diff_rcm_rcp85_cwd_matopiba, diff_rcm_rcp85_r10mm_matopiba, diff_rcm_rcp85_r20mm_matopiba]

diff_gcm_rcp26_prcptot_matopiba = gcm_prcptot_rcp26_matopiba - gcm_prcptot_hist_matopiba
diff_gcm_rcp26_r95p_matopiba = gcm_r95p_rcp26_matopiba - gcm_r95p_hist_matopiba
diff_gcm_rcp26_r99p_matopiba = gcm_r99p_rcp26_matopiba - gcm_r99p_hist_matopiba
diff_gcm_rcp26_rx1day_matopiba = gcm_rx1day_rcp26_matopiba - gcm_rx1day_hist_matopiba
diff_gcm_rcp26_rx5day_matopiba = gcm_rx5day_rcp26_matopiba - gcm_rx5day_hist_matopiba
diff_gcm_rcp26_sdii_matopiba = gcm_sdii_rcp26_matopiba - gcm_sdii_hist_matopiba
diff_gcm_rcp26_cdd_matopiba = gcm_cdd_rcp26_matopiba - gcm_cdd_hist_matopiba
diff_gcm_rcp26_cwd_matopiba = gcm_cwd_rcp26_matopiba - gcm_cwd_hist_matopiba
diff_gcm_rcp26_r10mm_matopiba = gcm_r10mm_rcp26_matopiba - gcm_r10mm_hist_matopiba
diff_gcm_rcp26_r20mm_matopiba = gcm_r20mm_rcp26_matopiba - gcm_r20mm_hist_eneb

diff_gcm_rcp26_matopiba_pre = [diff_gcm_rcp26_prcptot_matopiba, diff_gcm_rcp26_r95p_matopiba, diff_gcm_rcp26_r99p_matopiba,
diff_gcm_rcp26_rx1day_matopiba, diff_gcm_rcp26_rx5day_matopiba, diff_gcm_rcp26_sdii_matopiba, diff_gcm_rcp26_cdd_matopiba,
diff_gcm_rcp26_cwd_matopiba, diff_gcm_rcp26_r10mm_matopiba, diff_gcm_rcp26_r20mm_matopiba]

diff_gcm_rcp85_prcptot_matopiba = gcm_prcptot_rcp85_matopiba - gcm_prcptot_hist_matopiba
diff_gcm_rcp85_r95p_matopiba = gcm_r95p_rcp85_matopiba - gcm_r95p_hist_matopiba
diff_gcm_rcp85_r99p_matopiba = gcm_r99p_rcp85_matopiba - gcm_r99p_hist_matopiba
diff_gcm_rcp85_rx1day_matopiba = gcm_rx1day_rcp85_matopiba - gcm_rx1day_hist_matopiba
diff_gcm_rcp85_rx5day_matopiba = gcm_rx5day_rcp85_matopiba - gcm_rx5day_hist_matopiba
diff_gcm_rcp85_sdii_matopiba = gcm_sdii_rcp85_matopiba - gcm_sdii_hist_matopiba
diff_gcm_rcp85_cdd_matopiba = gcm_cdd_rcp85_matopiba - gcm_cdd_hist_matopiba
diff_gcm_rcp85_cwd_matopiba = gcm_cwd_rcp85_matopiba - gcm_cwd_hist_matopiba
diff_gcm_rcp85_r10mm_matopiba = gcm_r10mm_rcp85_matopiba - gcm_r10mm_hist_matopiba
diff_gcm_rcp85_r20mm_matopiba = gcm_r20mm_rcp85_matopiba - gcm_r20mm_hist_matopiba

diff_gcm_rcp85_matopiba_pre = [diff_gcm_rcp85_prcptot_matopiba, diff_gcm_rcp85_r95p_matopiba, diff_gcm_rcp85_r99p_matopiba,
diff_gcm_rcp85_rx1day_matopiba, diff_gcm_rcp85_rx5day_matopiba, diff_gcm_rcp85_sdii_matopiba, diff_gcm_rcp85_cdd_matopiba,
diff_gcm_rcp85_cwd_matopiba, diff_gcm_rcp85_r10mm_matopiba, diff_gcm_rcp85_r20mm_matopiba]

# Temperature
# SAMZ
diff_rcm_rcp26_txx_samz = rcm_txx_rcp26_samz - rcm_txx_hist_samz
diff_rcm_rcp26_txn_samz = rcm_txn_rcp26_samz - rcm_txn_hist_samz
diff_rcm_rcp26_tnx_samz = rcm_tnx_rcp26_samz - rcm_tnx_hist_samz
diff_rcm_rcp26_tnn_samz = rcm_tnn_rcp26_samz - rcm_tnn_hist_samz
diff_rcm_rcp26_dtr_samz = rcm_dtr_rcp26_samz - rcm_dtr_hist_samz
diff_rcm_rcp26_su_samz = rcm_su_rcp26_samz - rcm_su_hist_samz
diff_rcm_rcp26_tr_samz = rcm_tr_rcp26_samz - rcm_tr_hist_samz
diff_rcm_rcp26_tx10p_samz = rcm_tx10p_rcp26_samz - rcm_tx10p_hist_samz
diff_rcm_rcp26_tx90p_samz = rcm_tx90p_rcp26_samz - rcm_tx90p_hist_samz
diff_rcm_rcp26_tn10p_samz = rcm_tx10p_rcp26_samz - rcm_tx10p_hist_samz
diff_rcm_rcp26_tn90p_samz = rcm_tx90p_rcp26_samz - rcm_tx90p_hist_samz

diff_rcm_rcp26_samz_tas = [diff_rcm_rcp26_txx_samz, diff_rcm_rcp26_txn_samz, diff_rcm_rcp26_tnx_samz, diff_rcm_rcp26_tnn_samz,
diff_rcm_rcp26_dtr_samz, diff_rcm_rcp26_su_samz, diff_rcm_rcp26_tr_samz, diff_rcm_rcp26_tx10p_samz, diff_rcm_rcp26_tx90p_samz,
diff_rcm_rcp26_tn10p_samz, diff_rcm_rcp26_tn90p_samz]

diff_rcm_rcp85_txx_samz = rcm_txx_rcp85_samz - rcm_txx_hist_samz
diff_rcm_rcp85_txn_samz = rcm_txn_rcp85_samz - rcm_txn_hist_samz
diff_rcm_rcp85_tnx_samz = rcm_tnx_rcp85_samz - rcm_tnx_hist_samz
diff_rcm_rcp85_tnn_samz = rcm_tnn_rcp85_samz - rcm_tnn_hist_samz
diff_rcm_rcp85_dtr_samz = rcm_dtr_rcp85_samz - rcm_dtr_hist_samz
diff_rcm_rcp85_su_samz = rcm_su_rcp85_samz - rcm_su_hist_samz
diff_rcm_rcp85_tr_samz = rcm_tr_rcp85_samz - rcm_tr_hist_samz
diff_rcm_rcp85_tx10p_samz = rcm_tx10p_rcp85_samz - rcm_tx10p_hist_samz
diff_rcm_rcp85_tx90p_samz = rcm_tx90p_rcp85_samz - rcm_tx90p_hist_samz
diff_rcm_rcp85_tn10p_samz = rcm_tx10p_rcp85_samz - rcm_tx10p_hist_samz
diff_rcm_rcp85_tn90p_samz = rcm_tx90p_rcp85_samz - rcm_tx90p_hist_samz

diff_rcm_rcp85_samz_tas = [diff_rcm_rcp85_txx_samz, diff_rcm_rcp85_txn_samz, diff_rcm_rcp85_tnx_samz, diff_rcm_rcp85_tnn_samz,
diff_rcm_rcp85_dtr_samz, diff_rcm_rcp85_su_samz, diff_rcm_rcp85_tr_samz, diff_rcm_rcp85_tx10p_samz, diff_rcm_rcp85_tx90p_samz,
diff_rcm_rcp85_tn10p_samz, diff_rcm_rcp85_tn90p_samz]

diff_gcm_rcp26_txx_samz = gcm_txx_rcp26_samz - gcm_txx_hist_samz
diff_gcm_rcp26_txn_samz = gcm_txn_rcp26_samz - gcm_txn_hist_samz
diff_gcm_rcp26_tnx_samz = gcm_tnx_rcp26_samz - gcm_tnx_hist_samz
diff_gcm_rcp26_tnn_samz = gcm_tnn_rcp26_samz - gcm_tnn_hist_samz
diff_gcm_rcp26_dtr_samz = gcm_dtr_rcp26_samz - gcm_dtr_hist_samz
diff_gcm_rcp26_su_samz = gcm_su_rcp26_samz - gcm_su_hist_samz
diff_gcm_rcp26_tr_samz = gcm_tr_rcp26_samz - gcm_tr_hist_samz
diff_gcm_rcp26_tx10p_samz = gcm_tx10p_rcp26_samz - gcm_tx10p_hist_samz
diff_gcm_rcp26_tx90p_samz = gcm_tx90p_rcp26_samz - gcm_tx90p_hist_samz
diff_gcm_rcp26_tn10p_samz = gcm_tx10p_rcp26_samz - gcm_tx10p_hist_samz
diff_gcm_rcp26_tn90p_samz = gcm_tx90p_rcp26_samz - gcm_tx90p_hist_samz

diff_gcm_rcp26_samz_tas = [diff_gcm_rcp26_txx_samz, diff_gcm_rcp26_txn_samz, diff_gcm_rcp26_tnx_samz, diff_gcm_rcp26_tnn_samz,
diff_gcm_rcp26_dtr_samz, diff_gcm_rcp26_su_samz, diff_gcm_rcp26_tr_samz, diff_gcm_rcp26_tx10p_samz, diff_gcm_rcp26_tx90p_samz,
diff_gcm_rcp26_tn10p_samz, diff_gcm_rcp26_tn90p_samz]

diff_gcm_rcp85_txx_samz = gcm_txx_rcp85_samz - gcm_txx_hist_samz
diff_gcm_rcp85_txn_samz = gcm_txn_rcp85_samz - gcm_txn_hist_samz
diff_gcm_rcp85_tnx_samz = gcm_tnx_rcp85_samz - gcm_tnx_hist_samz
diff_gcm_rcp85_tnn_samz = gcm_tnn_rcp85_samz - gcm_tnn_hist_samz
diff_gcm_rcp85_dtr_samz = gcm_dtr_rcp85_samz - gcm_dtr_hist_samz
diff_gcm_rcp85_su_samz = gcm_su_rcp85_samz - gcm_su_hist_samz
diff_gcm_rcp85_tr_samz = gcm_tr_rcp85_samz - gcm_tr_hist_samz
diff_gcm_rcp85_tx10p_samz = gcm_tx10p_rcp85_samz - gcm_tx10p_hist_samz
diff_gcm_rcp85_tx90p_samz = gcm_tx90p_rcp85_samz - gcm_tx90p_hist_samz
diff_gcm_rcp85_tn10p_samz = gcm_tx10p_rcp85_samz - gcm_tx10p_hist_samz
diff_gcm_rcp85_tn90p_samz = gcm_tx90p_rcp85_samz - gcm_tx90p_hist_samz

diff_gcm_rcp85_samz_tas = [diff_gcm_rcp85_txx_samz, diff_gcm_rcp85_txn_samz, diff_gcm_rcp85_tnx_samz, diff_gcm_rcp85_tnn_samz,
diff_gcm_rcp85_dtr_samz, diff_gcm_rcp85_su_samz, diff_gcm_rcp85_tr_samz, diff_gcm_rcp85_tx10p_samz, diff_gcm_rcp85_tx90p_samz,
diff_gcm_rcp85_tn10p_samz, diff_gcm_rcp85_tn90p_samz]

# ENEB
diff_rcm_rcp26_txx_eneb = rcm_txx_rcp26_eneb - rcm_txx_hist_eneb
diff_rcm_rcp26_txn_eneb = rcm_txn_rcp26_eneb - rcm_txn_hist_eneb
diff_rcm_rcp26_tnx_eneb = rcm_tnx_rcp26_eneb - rcm_tnx_hist_eneb
diff_rcm_rcp26_tnn_eneb = rcm_tnn_rcp26_eneb - rcm_tnn_hist_eneb
diff_rcm_rcp26_dtr_eneb = rcm_dtr_rcp26_eneb - rcm_dtr_hist_eneb
diff_rcm_rcp26_su_eneb = rcm_su_rcp26_eneb - rcm_su_hist_eneb
diff_rcm_rcp26_tr_eneb = rcm_tr_rcp26_eneb - rcm_tr_hist_eneb
diff_rcm_rcp26_tx10p_eneb = rcm_tx10p_rcp26_eneb - rcm_tx10p_hist_eneb
diff_rcm_rcp26_tx90p_eneb = rcm_tx90p_rcp26_eneb - rcm_tx90p_hist_eneb
diff_rcm_rcp26_tn10p_eneb = rcm_tx10p_rcp26_eneb - rcm_tx10p_hist_eneb
diff_rcm_rcp26_tn90p_eneb = rcm_tx90p_rcp26_eneb - rcm_tx90p_hist_eneb

diff_rcm_rcp26_eneb_tas = [diff_rcm_rcp26_txx_eneb, diff_rcm_rcp26_txn_eneb, diff_rcm_rcp26_tnx_eneb, diff_rcm_rcp26_tnn_eneb,
diff_rcm_rcp26_dtr_eneb, diff_rcm_rcp26_su_eneb, diff_rcm_rcp26_tr_eneb, diff_rcm_rcp26_tx10p_eneb, diff_rcm_rcp26_tx90p_eneb,
diff_rcm_rcp26_tn10p_eneb, diff_rcm_rcp26_tn90p_eneb]

diff_rcm_rcp85_txx_eneb = rcm_txx_rcp85_eneb - rcm_txx_hist_eneb
diff_rcm_rcp85_txn_eneb = rcm_txn_rcp85_eneb - rcm_txn_hist_eneb
diff_rcm_rcp85_tnx_eneb = rcm_tnx_rcp85_eneb - rcm_tnx_hist_eneb
diff_rcm_rcp85_tnn_eneb = rcm_tnn_rcp85_eneb - rcm_tnn_hist_eneb
diff_rcm_rcp85_dtr_eneb = rcm_dtr_rcp85_eneb - rcm_dtr_hist_eneb
diff_rcm_rcp85_su_eneb = rcm_su_rcp85_eneb - rcm_su_hist_eneb
diff_rcm_rcp85_tr_eneb = rcm_tr_rcp85_eneb - rcm_tr_hist_eneb
diff_rcm_rcp85_tx10p_eneb = rcm_tx10p_rcp85_eneb - rcm_tx10p_hist_eneb
diff_rcm_rcp85_tx90p_eneb = rcm_tx90p_rcp85_eneb - rcm_tx90p_hist_eneb
diff_rcm_rcp85_tn10p_eneb = rcm_tx10p_rcp85_eneb - rcm_tx10p_hist_eneb
diff_rcm_rcp85_tn90p_eneb = rcm_tx90p_rcp85_eneb - rcm_tx90p_hist_eneb

diff_rcm_rcp85_eneb_tas = [diff_rcm_rcp85_txx_eneb, diff_rcm_rcp85_txn_eneb, diff_rcm_rcp85_tnx_eneb, diff_rcm_rcp85_tnn_eneb,
diff_rcm_rcp85_dtr_eneb, diff_rcm_rcp85_su_eneb, diff_rcm_rcp85_tr_eneb, diff_rcm_rcp85_tx10p_eneb, diff_rcm_rcp85_tx90p_eneb,
diff_rcm_rcp85_tn10p_eneb, diff_rcm_rcp85_tn90p_eneb]

diff_gcm_rcp26_txx_eneb = gcm_txx_rcp26_eneb - gcm_txx_hist_eneb
diff_gcm_rcp26_txn_eneb = gcm_txn_rcp26_eneb - gcm_txn_hist_eneb
diff_gcm_rcp26_tnx_eneb = gcm_tnx_rcp26_eneb - gcm_tnx_hist_eneb
diff_gcm_rcp26_tnn_eneb = gcm_tnn_rcp26_eneb - gcm_tnn_hist_eneb
diff_gcm_rcp26_dtr_eneb = gcm_dtr_rcp26_eneb - gcm_dtr_hist_eneb
diff_gcm_rcp26_su_eneb = gcm_su_rcp26_eneb - gcm_su_hist_eneb
diff_gcm_rcp26_tr_eneb = gcm_tr_rcp26_eneb - gcm_tr_hist_eneb
diff_gcm_rcp26_tx10p_eneb = gcm_tx10p_rcp26_eneb - gcm_tx10p_hist_eneb
diff_gcm_rcp26_tx90p_eneb = gcm_tx90p_rcp26_eneb - gcm_tx90p_hist_eneb
diff_gcm_rcp26_tn10p_eneb = gcm_tx10p_rcp26_eneb - gcm_tx10p_hist_eneb
diff_gcm_rcp26_tn90p_eneb = gcm_tx90p_rcp26_eneb - gcm_tx90p_hist_eneb

diff_gcm_rcp26_eneb_tas = [diff_gcm_rcp26_txx_eneb, diff_gcm_rcp26_txn_eneb, diff_gcm_rcp26_tnx_eneb, diff_gcm_rcp26_tnn_eneb,
diff_gcm_rcp26_dtr_eneb, diff_gcm_rcp26_su_eneb, diff_gcm_rcp26_tr_eneb, diff_gcm_rcp26_tx10p_eneb, diff_gcm_rcp26_tx90p_eneb,
diff_gcm_rcp26_tn10p_eneb, diff_gcm_rcp26_tn90p_eneb]

diff_gcm_rcp85_txx_eneb = gcm_txx_rcp85_eneb - gcm_txx_hist_eneb
diff_gcm_rcp85_txn_eneb = gcm_txn_rcp85_eneb - gcm_txn_hist_eneb
diff_gcm_rcp85_tnx_eneb = gcm_tnx_rcp85_eneb - gcm_tnx_hist_eneb
diff_gcm_rcp85_tnn_eneb = gcm_tnn_rcp85_eneb - gcm_tnn_hist_eneb
diff_gcm_rcp85_dtr_eneb = gcm_dtr_rcp85_eneb - gcm_dtr_hist_eneb
diff_gcm_rcp85_su_eneb = gcm_su_rcp85_eneb - gcm_su_hist_eneb
diff_gcm_rcp85_tr_eneb = gcm_tr_rcp85_eneb - gcm_tr_hist_eneb
diff_gcm_rcp85_tx10p_eneb = gcm_tx10p_rcp85_eneb - gcm_tx10p_hist_eneb
diff_gcm_rcp85_tx90p_eneb = gcm_tx90p_rcp85_eneb - gcm_tx90p_hist_eneb
diff_gcm_rcp85_tn10p_eneb = gcm_tx10p_rcp85_eneb - gcm_tx10p_hist_eneb
diff_gcm_rcp85_tn90p_eneb = gcm_tx90p_rcp85_eneb - gcm_tx90p_hist_eneb

diff_gcm_rcp85_eneb_tas = [diff_gcm_rcp85_txx_eneb, diff_gcm_rcp85_txn_eneb, diff_gcm_rcp85_tnx_eneb, diff_gcm_rcp85_tnn_eneb,
diff_gcm_rcp85_dtr_eneb, diff_gcm_rcp85_su_eneb, diff_gcm_rcp85_tr_eneb, diff_gcm_rcp85_tx10p_eneb, diff_gcm_rcp85_tx90p_eneb,
diff_gcm_rcp85_tn10p_eneb, diff_gcm_rcp85_tn90p_eneb]

# MATOPIBA
diff_rcm_rcp26_txx_matopiba = rcm_txx_rcp26_matopiba - rcm_txx_hist_matopiba
diff_rcm_rcp26_txn_matopiba = rcm_txn_rcp26_matopiba - rcm_txn_hist_matopiba
diff_rcm_rcp26_tnx_matopiba = rcm_tnx_rcp26_matopiba - rcm_tnx_hist_matopiba
diff_rcm_rcp26_tnn_matopiba = rcm_tnn_rcp26_matopiba - rcm_tnn_hist_matopiba
diff_rcm_rcp26_dtr_matopiba = rcm_dtr_rcp26_matopiba - rcm_dtr_hist_matopiba
diff_rcm_rcp26_su_matopiba = rcm_su_rcp26_matopiba - rcm_su_hist_matopiba
diff_rcm_rcp26_tr_matopiba = rcm_tr_rcp26_matopiba - rcm_tr_hist_matopiba
diff_rcm_rcp26_tx10p_matopiba = rcm_tx10p_rcp26_matopiba - rcm_tx10p_hist_matopiba
diff_rcm_rcp26_tx90p_matopiba = rcm_tx90p_rcp26_matopiba - rcm_tx90p_hist_matopiba
diff_rcm_rcp26_tn10p_matopiba = rcm_tx10p_rcp26_matopiba - rcm_tx10p_hist_matopiba
diff_rcm_rcp26_tn90p_matopiba = rcm_tx90p_rcp26_matopiba - rcm_tx90p_hist_matopiba

diff_rcm_rcp26_matopiba_tas = [diff_rcm_rcp26_txx_matopiba, diff_rcm_rcp26_txn_matopiba, diff_rcm_rcp26_tnx_matopiba, diff_rcm_rcp26_tnn_matopiba,
diff_rcm_rcp26_dtr_matopiba, diff_rcm_rcp26_su_matopiba, diff_rcm_rcp26_tr_matopiba, diff_rcm_rcp26_tx10p_matopiba, diff_rcm_rcp26_tx90p_matopiba,
diff_rcm_rcp26_tn10p_matopiba, diff_rcm_rcp26_tn90p_matopiba]

diff_rcm_rcp85_txx_matopiba = rcm_txx_rcp85_matopiba - rcm_txx_hist_matopiba
diff_rcm_rcp85_txn_matopiba = rcm_txn_rcp85_matopiba - rcm_txn_hist_matopiba
diff_rcm_rcp85_tnx_matopiba = rcm_tnx_rcp85_matopiba - rcm_tnx_hist_matopiba
diff_rcm_rcp85_tnn_matopiba = rcm_tnn_rcp85_matopiba - rcm_tnn_hist_matopiba
diff_rcm_rcp85_dtr_matopiba = rcm_dtr_rcp85_matopiba - rcm_dtr_hist_matopiba
diff_rcm_rcp85_su_matopiba = rcm_su_rcp85_matopiba - rcm_su_hist_matopiba
diff_rcm_rcp85_tr_matopiba = rcm_tr_rcp85_matopiba - rcm_tr_hist_matopiba
diff_rcm_rcp85_tx10p_matopiba = rcm_tx10p_rcp85_matopiba - rcm_tx10p_hist_matopiba
diff_rcm_rcp85_tx90p_matopiba = rcm_tx90p_rcp85_matopiba - rcm_tx90p_hist_matopiba
diff_rcm_rcp85_tn10p_matopiba = rcm_tx10p_rcp85_matopiba - rcm_tx10p_hist_matopiba
diff_rcm_rcp85_tn90p_matopiba = rcm_tx90p_rcp85_matopiba - rcm_tx90p_hist_matopiba

diff_rcm_rcp85_matopiba_tas = [diff_rcm_rcp85_txx_matopiba, diff_rcm_rcp85_txn_matopiba, diff_rcm_rcp85_tnx_matopiba, diff_rcm_rcp85_tnn_matopiba,
diff_rcm_rcp85_dtr_matopiba, diff_rcm_rcp85_su_matopiba, diff_rcm_rcp85_tr_matopiba, diff_rcm_rcp85_tx10p_matopiba, diff_rcm_rcp85_tx90p_matopiba,
diff_rcm_rcp85_tn10p_matopiba, diff_rcm_rcp85_tn90p_matopiba]

diff_gcm_rcp26_txx_matopiba = gcm_txx_rcp26_matopiba - gcm_txx_hist_matopiba
diff_gcm_rcp26_txn_matopiba = gcm_txn_rcp26_matopiba - gcm_txn_hist_matopiba
diff_gcm_rcp26_tnx_matopiba = gcm_tnx_rcp26_matopiba - gcm_tnx_hist_matopiba
diff_gcm_rcp26_tnn_matopiba = gcm_tnn_rcp26_matopiba- gcm_tnn_hist_matopiba
diff_gcm_rcp26_dtr_matopiba = gcm_dtr_rcp26_matopiba - gcm_dtr_hist_matopiba
diff_gcm_rcp26_su_matopiba = gcm_su_rcp26_matopiba - gcm_su_hist_matopiba
diff_gcm_rcp26_tr_matopiba = gcm_tr_rcp26_matopiba - gcm_tr_hist_matopiba
diff_gcm_rcp26_tx10p_matopiba = gcm_tx10p_rcp26_matopiba - gcm_tx10p_hist_matopiba
diff_gcm_rcp26_tx90p_matopiba = gcm_tx90p_rcp26_matopiba - gcm_tx90p_hist_matopiba
diff_gcm_rcp26_tn10p_matopiba = gcm_tx10p_rcp26_matopiba - gcm_tx10p_hist_matopiba
diff_gcm_rcp26_tn90p_matopiba = gcm_tx90p_rcp26_matopiba- gcm_tx90p_hist_matopiba

diff_gcm_rcp26_matopiba_tas = [diff_gcm_rcp26_txx_matopiba, diff_gcm_rcp26_txn_matopiba, diff_gcm_rcp26_tnx_matopiba, diff_gcm_rcp26_tnn_matopiba,
diff_gcm_rcp26_dtr_matopiba, diff_gcm_rcp26_su_matopiba, diff_gcm_rcp26_tr_matopiba, diff_gcm_rcp26_tx10p_matopiba, diff_gcm_rcp26_tx90p_matopiba,
diff_gcm_rcp26_tn10p_matopiba, diff_gcm_rcp26_tn90p_matopiba]

diff_gcm_rcp85_txx_matopiba = gcm_txx_rcp85_matopiba - gcm_txx_hist_matopiba
diff_gcm_rcp85_txn_matopiba = gcm_txn_rcp85_matopiba - gcm_txn_hist_matopiba
diff_gcm_rcp85_tnx_matopiba = gcm_tnx_rcp85_matopiba - gcm_tnx_hist_matopiba
diff_gcm_rcp85_tnn_matopiba = gcm_tnn_rcp85_matopiba- gcm_tnn_hist_matopiba
diff_gcm_rcp85_dtr_matopiba = gcm_dtr_rcp85_matopiba - gcm_dtr_hist_matopiba
diff_gcm_rcp85_su_matopiba = gcm_su_rcp85_matopiba - gcm_su_hist_matopiba
diff_gcm_rcp85_tr_matopiba = gcm_tr_rcp85_matopiba - gcm_tr_hist_matopiba
diff_gcm_rcp85_tx10p_matopiba = gcm_tx10p_rcp85_matopiba - gcm_tx10p_hist_matopiba
diff_gcm_rcp85_tx90p_matopiba = gcm_tx90p_rcp85_matopiba - gcm_tx90p_hist_matopiba
diff_gcm_rcp85_tn10p_matopiba = gcm_tx10p_rcp85_matopiba - gcm_tx10p_hist_matopiba
diff_gcm_rcp85_tn90p_matopiba = gcm_tx90p_rcp85_matopiba- gcm_tx90p_hist_matopiba

diff_gcm_rcp85_matopiba_tas = [diff_gcm_rcp85_txx_matopiba, diff_gcm_rcp85_txn_matopiba, diff_gcm_rcp85_tnx_matopiba, diff_gcm_rcp85_tnn_matopiba,
diff_gcm_rcp85_dtr_matopiba, diff_gcm_rcp85_su_matopiba, diff_gcm_rcp85_tr_matopiba, diff_gcm_rcp85_tx10p_matopiba, diff_gcm_rcp85_tx90p_matopiba,
diff_gcm_rcp85_tn10p_matopiba, diff_gcm_rcp85_tn90p_matopiba]

# Plot project change from etccdi indices
fig = plt.figure()
time1 = np.arange(1, 11)
time2 = np.arange(1, 12)

ax1 = fig.add_subplot(3, 2, 1)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
plt_clim1 = plt.plot(time1, diff_rcm_rcp26_samz_pre, 'o', color='black', markerfacecolor='blue', label='Reg RCP2.6')
plt_clim2 = plt.plot(time1, diff_rcm_rcp85_samz_pre, 'o', color='black', markerfacecolor='red', label='Reg RCP8.5')
plt_clim3 = plt.plot(time1, diff_gcm_rcp26_samz_pre, 's', color='black', markerfacecolor='blue', label='Had RCP2.6')
plt_clim4 = plt.plot(time1, diff_gcm_rcp85_samz_pre, 's', color='black', markerfacecolor='red', label='Had RCP8.5')
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.ylim(-8, 8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.xticks(time1, ('PRCPTOT', 'R95p', 'R99p', 'Rx1day', 'Rx5day', 'SDII', 'CDD', 'CWD', 'R10mm', 'R20mm'), fontsize=8)
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = fig.add_subplot(3, 2, 2)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
plt_clim1 = plt.plot(time2, diff_rcm_rcp26_samz_tas, 'o', color='black', markerfacecolor='blue', label='Reg RCP2.6')
plt_clim2 = plt.plot(time2, diff_rcm_rcp85_samz_tas, 'o', color='black', markerfacecolor='red', label='Reg RCP8.5')
plt_clim3 = plt.plot(time2, diff_gcm_rcp26_samz_tas, 's', color='black', markerfacecolor='blue', label='Had RCP2.6')
plt_clim4 = plt.plot(time2, diff_gcm_rcp85_samz_tas, 's', color='black', markerfacecolor='red', label='Had RCP8.5')
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.ylim(-8, 8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.xticks(time2, ('TXX', 'TXn', 'TNx', 'TNn', 'DTR', 'SU', 'TR', 'Tx10p', 'Tx90p', 'Tn10p', 'Tn90p'), fontsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.legend(fontsize=8, loc=4, ncol=2, shadow=True, handlelength=0.75, handleheight=0.75)

ax3 = fig.add_subplot(3, 2, 3)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
plt_clim1 = plt.plot(time1, diff_rcm_rcp26_eneb_pre, 'o', color='black', markerfacecolor='blue', label='Reg RCP2.6')
plt_clim2 = plt.plot(time1, diff_rcm_rcp85_eneb_pre, 'o', color='black', markerfacecolor='red', label='Reg RCP8.5')
plt_clim3 = plt.plot(time1, diff_gcm_rcp26_eneb_pre, 's', color='black', markerfacecolor='blue', label='Had RCP2.6')
plt_clim4 = plt.plot(time1, diff_gcm_rcp85_eneb_pre, 's', color='black', markerfacecolor='red', label='Had RCP8.5')
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Projected changes', fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.ylim(-8, 8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.xticks(time1, ('PRCPTOT', 'R95p', 'R99p', 'Rx1day', 'Rx5day', 'SDII', 'CDD', 'CWD', 'R10mm', 'R20mm'), fontsize=8)
plt.setp(ax3.get_xticklabels(), visible=False)

ax4 = fig.add_subplot(3, 2, 4)
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
plt_clim1 = plt.plot(time2, diff_rcm_rcp26_eneb_tas, 'o', color='black', markerfacecolor='blue', label='Reg RCP2.6')
plt_clim2 = plt.plot(time2, diff_rcm_rcp85_eneb_tas, 'o', color='black', markerfacecolor='red', label='Reg RCP8.5')
plt_clim3 = plt.plot(time2, diff_gcm_rcp26_eneb_tas, 's', color='black', markerfacecolor='blue', label='Had RCP2.6')
plt_clim4 = plt.plot(time2, diff_gcm_rcp85_eneb_tas, 's', color='black', markerfacecolor='red', label='Had RCP8.5')
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Projected changes', fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.ylim(-8, 8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.xticks(time2, ('TXX', 'TXn', 'TNx', 'TNn', 'DTR', 'SU', 'TR', 'Tx10p', 'Tx90p', 'Tn10p', 'Tn90p'), fontsize=8)
plt.setp(ax4.get_xticklabels(), visible=False)

ax5 = fig.add_subplot(3, 2, 5)
ax5.spines['right'].set_visible(False)
ax5.spines['top'].set_visible(False)
plt_clim1 = plt.plot(time1, diff_rcm_rcp26_matopiba_pre, 'o', color='black', markerfacecolor='blue', label='Reg RCP2.6')
plt_clim2 = plt.plot(time1, diff_rcm_rcp85_matopiba_pre, 'o', color='black', markerfacecolor='red', label='Reg RCP8.5')
plt_clim3 = plt.plot(time1, diff_gcm_rcp26_matopiba_pre, 's', color='black', markerfacecolor='blue', label='Had RCP2.6')
plt_clim4 = plt.plot(time1, diff_gcm_rcp85_matopiba_pre, 's', color='black', markerfacecolor='red', label='Had RCP8.5')
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.ylim(-8, 8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.xticks(time1, ('PRCPTOT', 'R95p', 'R99p', 'Rx1day', 'Rx5day', 'SDII', 'CDD', 'CWD', 'R10mm', 'R20mm'), fontsize=8)
labels = ax5.get_xticklabels()
plt.setp(labels, rotation=90)

ax6 = fig.add_subplot(3, 2, 6)
ax6.spines['right'].set_visible(False)
ax6.spines['top'].set_visible(False)
plt_clim1 = plt.plot(time2, diff_rcm_rcp26_matopiba_tas, 'o', color='black', markerfacecolor='blue', label='Reg RCP2.6')
plt_clim2 = plt.plot(time2, diff_rcm_rcp85_matopiba_tas, 'o', color='black', markerfacecolor='red', label='Reg RCP8.5')
plt_clim3 = plt.plot(time2, diff_gcm_rcp26_matopiba_tas, 's', color='black', markerfacecolor='blue', label='Had RCP2.6')
plt_clim4 = plt.plot(time2, diff_gcm_rcp85_matopiba_tas, 's', color='black', markerfacecolor='red', label='Had RCP8.5')
plt.title(u'E)', loc='left', fontweight='bold', fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.ylim(-8, 8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.xticks(time2, ('TXX', 'TXn', 'TNx', 'TNn', 'DTR', 'SU', 'TR', 'Tx10p', 'Tx90p', 'Tn10p', 'Tn90p'), fontsize=8)
labels = ax6.get_xticklabels()
plt.setp(labels, rotation=90)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_diff_etccdi_reg_had_rcp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')

plt.show()
exit()	

plt.close('all')
plt.cla()
exit()	

