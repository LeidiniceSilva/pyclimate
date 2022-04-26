# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "06/24/2021"
__description__ = "This script plot change projected from extremes indices"

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
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}
	
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
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	annual_gcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	annual_gcm = np.nanmean(annual_gcm)

	return annual_gcm
	

# Import extreme indices 
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

rcm_rcp26_prcptot = [diff_rcm_rcp26_prcptot_samz, diff_rcm_rcp26_prcptot_eneb, diff_rcm_rcp26_prcptot_matopiba]
rcm_rcp85_prcptot = [diff_rcm_rcp85_prcptot_samz, diff_rcm_rcp85_prcptot_eneb, diff_rcm_rcp85_prcptot_matopiba]
gcm_rcp26_prcptot = [diff_gcm_rcp26_prcptot_samz, diff_gcm_rcp26_prcptot_eneb, diff_gcm_rcp26_prcptot_matopiba]
gcm_rcp85_prcptot = [diff_gcm_rcp85_prcptot_samz, diff_gcm_rcp85_prcptot_eneb, diff_gcm_rcp85_prcptot_matopiba]

rcm_rcp26_r95p = [diff_rcm_rcp26_r95p_samz, diff_rcm_rcp26_r95p_eneb, diff_rcm_rcp26_r95p_matopiba]
rcm_rcp85_r95p = [diff_rcm_rcp85_r95p_samz, diff_rcm_rcp85_r95p_eneb, diff_rcm_rcp85_r95p_matopiba]
gcm_rcp26_r95p = [diff_gcm_rcp26_r95p_samz, diff_gcm_rcp26_r95p_eneb, diff_gcm_rcp26_r95p_matopiba]
gcm_rcp85_r95p = [diff_gcm_rcp85_r95p_samz, diff_gcm_rcp85_r95p_eneb, diff_gcm_rcp85_r95p_matopiba]

rcm_rcp26_r99p = [diff_rcm_rcp26_r95p_samz, diff_rcm_rcp26_r99p_eneb, diff_rcm_rcp26_r99p_matopiba]
rcm_rcp85_r99p = [diff_rcm_rcp85_r95p_samz, diff_rcm_rcp85_r99p_eneb, diff_rcm_rcp85_r99p_matopiba]
gcm_rcp26_r99p = [diff_gcm_rcp26_r95p_samz, diff_gcm_rcp85_r99p_eneb, diff_gcm_rcp85_r99p_matopiba]
gcm_rcp85_r99p = [diff_gcm_rcp85_r95p_samz, diff_gcm_rcp85_r99p_eneb, diff_gcm_rcp85_r99p_matopiba]

rcm_rcp26_rx1day = [diff_rcm_rcp26_rx1day_samz, diff_rcm_rcp26_rx1day_eneb, diff_rcm_rcp26_rx1day_matopiba]
rcm_rcp85_rx1day = [diff_rcm_rcp85_rx1day_samz, diff_rcm_rcp85_rx1day_eneb, diff_rcm_rcp85_rx1day_matopiba]
gcm_rcp26_rx1day = [diff_gcm_rcp26_rx1day_samz, diff_gcm_rcp26_rx1day_eneb, diff_gcm_rcp26_rx1day_matopiba]
gcm_rcp85_rx1day = [diff_gcm_rcp85_rx1day_samz, diff_gcm_rcp85_rx1day_eneb, diff_gcm_rcp85_rx1day_matopiba]

rcm_rcp26_rx5day = [diff_rcm_rcp26_rx5day_samz, diff_rcm_rcp26_rx5day_eneb, diff_rcm_rcp26_rx5day_matopiba]
rcm_rcp85_rx5day = [diff_rcm_rcp85_rx5day_samz, diff_rcm_rcp85_rx5day_eneb, diff_rcm_rcp85_rx5day_matopiba]
gcm_rcp26_rx5day = [diff_gcm_rcp26_rx5day_samz, diff_gcm_rcp26_rx5day_eneb, diff_gcm_rcp26_rx5day_matopiba]
gcm_rcp85_rx5day = [diff_gcm_rcp85_rx5day_samz, diff_gcm_rcp85_rx5day_eneb, diff_gcm_rcp85_rx5day_matopiba]

rcm_rcp26_sdii = [diff_rcm_rcp26_sdii_samz, diff_rcm_rcp26_sdii_eneb, diff_rcm_rcp26_sdii_matopiba]
rcm_rcp85_sdii = [diff_rcm_rcp85_sdii_samz, diff_rcm_rcp85_sdii_eneb, diff_rcm_rcp85_sdii_matopiba]
gcm_rcp26_sdii = [diff_gcm_rcp26_sdii_samz, diff_gcm_rcp26_sdii_eneb, diff_gcm_rcp26_sdii_matopiba]
gcm_rcp85_sdii = [diff_gcm_rcp85_sdii_samz, diff_gcm_rcp85_sdii_eneb, diff_gcm_rcp85_sdii_matopiba]

rcm_rcp26_cdd = [diff_rcm_rcp26_cdd_samz, diff_rcm_rcp26_cdd_eneb, diff_rcm_rcp26_cdd_matopiba]
rcm_rcp85_cdd = [diff_rcm_rcp85_cdd_samz, diff_rcm_rcp85_cdd_eneb, diff_rcm_rcp85_cdd_matopiba]
gcm_rcp26_cdd = [diff_gcm_rcp26_cdd_samz, diff_gcm_rcp26_cdd_eneb, diff_gcm_rcp26_cdd_matopiba]
gcm_rcp85_cdd = [diff_gcm_rcp85_cdd_samz, diff_gcm_rcp85_cdd_eneb, diff_gcm_rcp85_cdd_matopiba]

rcm_rcp26_cwd = [diff_rcm_rcp26_cwd_samz, diff_rcm_rcp26_cwd_eneb, diff_rcm_rcp26_cwd_matopiba]
rcm_rcp85_cwd = [diff_rcm_rcp85_cwd_samz, diff_rcm_rcp85_cwd_eneb, diff_rcm_rcp85_cwd_matopiba]
gcm_rcp26_cwd = [diff_gcm_rcp26_cwd_samz, diff_gcm_rcp26_cwd_eneb, diff_gcm_rcp26_cwd_matopiba]
gcm_rcp85_cwd = [diff_gcm_rcp85_cwd_samz, diff_gcm_rcp85_cwd_eneb, diff_gcm_rcp85_cwd_matopiba]

rcm_rcp26_r10mm = [diff_rcm_rcp26_r10mm_samz, diff_rcm_rcp26_r10mm_eneb, diff_rcm_rcp26_r10mm_matopiba]
rcm_rcp85_r10mm = [diff_rcm_rcp85_r10mm_samz, diff_rcm_rcp85_r10mm_eneb, diff_rcm_rcp85_r10mm_matopiba]
gcm_rcp26_r10mm = [diff_gcm_rcp26_r10mm_samz, diff_gcm_rcp26_r10mm_eneb, diff_gcm_rcp26_r10mm_matopiba]
gcm_rcp85_r10mm = [diff_gcm_rcp85_r10mm_samz, diff_gcm_rcp85_r10mm_eneb, diff_gcm_rcp85_r10mm_matopiba]

rcm_rcp26_r20mm = [diff_rcm_rcp26_r20mm_samz, diff_rcm_rcp26_r20mm_eneb, diff_rcm_rcp26_r20mm_matopiba]
rcm_rcp85_r20mm = [diff_rcm_rcp85_r20mm_samz, diff_rcm_rcp85_r20mm_eneb, diff_rcm_rcp85_r20mm_matopiba]
gcm_rcp26_r20mm = [diff_gcm_rcp26_r20mm_samz, diff_gcm_rcp26_r20mm_eneb, diff_gcm_rcp26_r20mm_matopiba]
gcm_rcp85_r20mm = [diff_gcm_rcp85_r20mm_samz, diff_gcm_rcp85_r20mm_eneb, diff_gcm_rcp85_r20mm_matopiba]

# Plot extreme indices 
fig = plt.figure()
time = np.arange(1, 4)

ax1 = fig.add_subplot(4, 3, 1)
plt1 = plt.plot(time, rcm_rcp26_prcptot, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_prcptot, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP2.6')
plt3 = plt.plot(time, gcm_rcp26_prcptot, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_prcptot, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP2.6')
plt.title(u'A) PRCPTOT (mm)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-1000, 1000)
plt.yticks(np.arange(-1000, 1500, 500), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax2 = fig.add_subplot(4, 3, 2)
plt1 = plt.plot(time, rcm_rcp26_r95p, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_r95p, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_r95p, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_r95p, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'B) R95p (mm)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-10, 10)
plt.yticks(np.arange(-10, 15, 5), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax3 = fig.add_subplot(4, 3, 3)
plt1 = plt.plot(time, rcm_rcp26_r99p, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_r99p, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP2.6')
plt3 = plt.plot(time, gcm_rcp26_r99p, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP8.5')
plt4 = plt.plot(time, gcm_rcp85_r99p, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'C) R99p (mm)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-10, 10)
plt.yticks(np.arange(-10, 15, 5), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax3.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax4 = fig.add_subplot(4, 3, 4)
plt1 = plt.plot(time, rcm_rcp26_rx1day, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_rx1day, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_rx1day, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_rx1day, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'D) Rx1day (mm)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-10, 10)
plt.yticks(np.arange(-10, 15, 5), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax4.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax5 = fig.add_subplot(4, 3, 5)
plt1 = plt.plot(time, rcm_rcp26_rx5day, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_rx5day, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_rx5day, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_rx5day, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'E) Rx5day (mm)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-20, 20)
plt.yticks(np.arange(-20, 30, 10), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax5.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax6 = fig.add_subplot(4, 3, 6)
plt1 = plt.plot(time, rcm_rcp26_sdii, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_sdii, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_sdii, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_sdii, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'F) SDII (mm d⁻¹)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-2, 2)
plt.yticks(np.arange(-2, 3, 1), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax6.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax7 = fig.add_subplot(4, 3, 7)
plt1 = plt.plot(time, rcm_rcp26_cdd, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_cdd, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_cdd, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_cdd, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'G) CDD (dias)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-200, 200)
plt.yticks(np.arange(-200, 300, 100), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax7.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax8 = fig.add_subplot(4, 3, 8)
plt1 = plt.plot(time, rcm_rcp26_cwd, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_cwd, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_cwd, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_cwd, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'H) CWD (dias)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-50, 50)
plt.yticks(np.arange(-50, 75, 25), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

ax9 = fig.add_subplot(4, 3, 9)
plt1 = plt.plot(time, rcm_rcp26_r10mm, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_r10mm, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_r10mm, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_r10mm, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'I) R10mm (dias)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-50, 50)
plt.yticks(np.arange(-50, 75, 25), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

ax10 = fig.add_subplot(4, 3, 10)
plt1 = plt.plot(time, rcm_rcp26_r20mm, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_r20mm, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_r20mm, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_r20mm, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'J) R20mm (dias)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-10, 10)
plt.yticks(np.arange(-10, 15, 5), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
legend(bbox_to_anchor=(1.05, 1.), loc=2, ncol=1, fontsize=7, frameon=False)
plt.grid(True, which='major', linestyle='--')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_diff_etccdi_reg_had_rcp_pre.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()	
