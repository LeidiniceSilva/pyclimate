# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "06ssss/24/2021"
__description__ = "This script plot correlation matrix from ETCCDI indices"

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
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}
	
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
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	annual_gcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	return annual_gcm
	

# Import regcm exp and cru databases 
# SAMZ
# Historical period
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

# RCP2.6
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

# RCP8.5
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

gcm_prcptot_rcp85_samz = import_gcm('eca_prcptot', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r95p_rcp85_samz = import_gcm('eca_r95p', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r99p_rcp85_samz = import_gcm('eca_r99p', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_rx1day_rcp85_samz = import_gcm('eca_rx1day', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_rx5day_rcp85_samz = import_gcm('eca_rx5day', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_sdii_rcp85_samz = import_gcm('eca_sdii', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_cdd_rcp85_samz = import_gcm('eca_cdd', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_cwd_rcp85_samz = import_gcm('eca_cwd', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r10mm_rcp85_samz = import_gcm('eca_r10mm', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')
gcm_r20mm_rcp85 _samz= import_gcm('eca_r20mm', 'samz', 'HadGEM2-ES', 'rcp85', 'yr', '2080-2099')

# ENEB
# Historical Period
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

# RCP2.6
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

# RCP8.5
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
# Historical Period
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

gcm_prcptot_hist_matopiba = import_gcm('eca_prcptot', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r95p_hist_matopiba= import_gcm('eca_r95p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r99p_hist_matopiba = import_gcm('eca_r99p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_rx1day_hist_matopiba = import_gcm('eca_rx1day', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_rx5day_hist_matopiba = import_gcm('eca_rx5day', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_sdii_hist_matopiba = import_gcm('eca_sdii', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_cdd_hist_matopiba = import_gcm('eca_cdd', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_cwd_hist_matopiba = import_gcm('eca_cwd', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r10mm_hist_matopiba = import_gcm('eca_r10mm', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
gcm_r20mm_hist_matopiba = import_gcm('eca_r20mm', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# RCP2.6
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

# RCP8.5
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

# Plot correlation matrix 
fig = plt.figure()
time1 = np.arange(1, 11)
time2 = np.arange(1, 12)
bar_width = .25

ax1 = fig.add_subplot(3, 2, 1)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
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
ax2.spines['bottom'].set_visible(False)
plt_clim1 = plt.bar(time2, ivs_rcm_samz_tas, color='blue', label='Reg', width = 0.25, edgecolor='black', linewidth=1)
plt_clim2 = plt.bar(time2 + .25, ivs_gcm_samz_tas,  color='red', label='Had', width = 0.25, edgecolor='black', linewidth=1)
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Value of IVS', fontsize=8)
plt.ylim(0, 8)
plt.xticks(time2 + .12, ('TXX', 'TXn', 'TNx', 'TNn', 'DTR', 'SU', 'TR', 'Tx10p', 'Tx90p', 'Tn10p', 'Tn90p'), fontsize=8)
plt.yticks(np.arange(0, 9, 1), fontsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.legend(fontsize=8, shadow=True, ncol=1, handlelength=0.75, handleheight=0.75)

ax3 = fig.add_subplot(3, 2, 3)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
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
ax4.spines['bottom'].set_visible(False)
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
ax5.spines['bottom'].set_visible(False)
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
ax6.spines['bottom'].set_visible(False)
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

plt.close('all')
plt.cla()
exit()	

