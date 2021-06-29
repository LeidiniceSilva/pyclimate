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
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	annual_rcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
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
	annual_gcm = np.nanmean(np.nanmean(var[:][20:,:,:], axis=1), axis=1)

	return annual_gcm
	

# Import regcm and hadgem models
# SAMZ
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

rcm_rcp26_samz = {'PRCPTOT': rcm_prcptot_rcp26_samz,
'R95p': rcm_r95p_rcp26_samz,
'R99p': rcm_r99p_rcp26_samz,
'Rx1day': rcm_rx1day_rcp26_samz,
'Rx5day': rcm_rx5day_rcp26_samz,
'SDII': rcm_sdii_rcp26_samz,
'CDD': rcm_cdd_rcp26_samz,
'CWD': rcm_cwd_rcp26_samz, 
'R10mm': rcm_r10mm_rcp26_samz,
'R20mm': rcm_r20mm_rcp26_samz}

rcm_rcp85_samz = {'PRCPTOT': rcm_prcptot_rcp85_samz,
'R95p': rcm_r95p_rcp85_samz,
'R99p': rcm_r99p_rcp85_samz,
'Rx1day': rcm_rx1day_rcp85_samz,
'Rx5day': rcm_rx5day_rcp85_samz,
'SDII': rcm_sdii_rcp85_samz,
'CDD': rcm_cdd_rcp85_samz,
'CWD': rcm_cwd_rcp85_samz, 
'R10mm': rcm_r10mm_rcp85_samz,
'R20mm': rcm_r20mm_rcp85_samz}

gcm_rcp26_samz = {'PRCPTOT': gcm_prcptot_rcp26_samz,
'R95p': gcm_r95p_rcp26_samz,
'R99p': gcm_r99p_rcp26_samz,
'Rx1day': gcm_rx1day_rcp26_samz,
'Rx5day': gcm_rx5day_rcp26_samz,
'SDII': gcm_sdii_rcp26_samz,
'CDD': gcm_cdd_rcp26_samz,
'CWD': gcm_cwd_rcp26_samz, 
'R10mm': gcm_r10mm_rcp26_samz,
'R20mm': gcm_r20mm_rcp26_samz}

gcm_rcp85_samz = {'PRCPTOT': gcm_prcptot_rcp85_samz,
'R95p': gcm_r95p_rcp85_samz,
'R99p': gcm_r99p_rcp85_samz,
'Rx1day': gcm_rx1day_rcp85_samz,
'Rx5day': gcm_rx5day_rcp85_samz,
'SDII': gcm_sdii_rcp85_samz,
'CDD': gcm_cdd_rcp85_samz,
'CWD': gcm_cwd_rcp85_samz, 
'R10mm': gcm_r10mm_rcp85_samz,
'R20mm': gcm_r20mm_rcp85_samz}

rcm_rcp26_eneb = {'PRCPTOT': rcm_prcptot_rcp26_eneb,
'R95p': rcm_r95p_rcp26_eneb,
'R99p': rcm_r99p_rcp26_eneb,
'Rx1day': rcm_rx1day_rcp26_eneb,
'Rx5day': rcm_rx5day_rcp26_eneb,
'SDII': rcm_sdii_rcp26_eneb,
'CDD': rcm_cdd_rcp26_eneb,
'CWD': rcm_cwd_rcp26_eneb, 
'R10mm': rcm_r10mm_rcp26_eneb,
'R20mm': rcm_r20mm_rcp26_eneb}

rcm_rcp85_eneb = {'PRCPTOT': rcm_prcptot_rcp85_eneb,
'R95p': rcm_r95p_rcp85_eneb,
'R99p': rcm_r99p_rcp85_eneb,
'Rx1day': rcm_rx1day_rcp85_eneb,
'Rx5day': rcm_rx5day_rcp85_eneb,
'SDII': rcm_sdii_rcp85_eneb,
'CDD': rcm_cdd_rcp85_eneb,
'CWD': rcm_cwd_rcp85_eneb, 
'R10mm': rcm_r10mm_rcp85_eneb,
'R20mm': rcm_r20mm_rcp85_eneb}

gcm_rcp26_eneb = {'PRCPTOT': gcm_prcptot_rcp26_eneb,
'R95p': gcm_r95p_rcp26_eneb,
'R99p': gcm_r99p_rcp26_eneb,
'Rx1day': gcm_rx1day_rcp26_eneb,
'Rx5day': gcm_rx5day_rcp26_eneb,
'SDII': gcm_sdii_rcp26_eneb,
'CDD': gcm_cdd_rcp26_eneb,
'CWD': gcm_cwd_rcp26_eneb, 
'R10mm': gcm_r10mm_rcp26_eneb,
'R20mm': gcm_r20mm_rcp26_eneb}

gcm_rcp85_eneb = {'PRCPTOT': gcm_prcptot_rcp85_eneb,
'R95p': gcm_r95p_rcp85_eneb,
'R99p': gcm_r99p_rcp85_eneb,
'Rx1day': gcm_rx1day_rcp85_eneb,
'Rx5day': gcm_rx5day_rcp85_eneb,
'SDII': gcm_sdii_rcp85_eneb,
'CDD': gcm_cdd_rcp85_eneb,
'CWD': gcm_cwd_rcp85_eneb, 
'R10mm': gcm_r10mm_rcp85_eneb,
'R20mm': gcm_r20mm_rcp85_eneb}

rcm_rcp26_matopiba = {'PRCPTOT': rcm_prcptot_rcp26_matopiba,
'R95p': rcm_r95p_rcp26_matopiba,
'R99p': rcm_r99p_rcp26_matopiba,
'Rx1day': rcm_rx1day_rcp26_matopiba,
'Rx5day': rcm_rx5day_rcp26_matopiba,
'SDII': rcm_sdii_rcp26_matopiba,
'CDD': rcm_cdd_rcp26_matopiba,
'CWD': rcm_cwd_rcp26_matopiba, 
'R10mm': rcm_r10mm_rcp26_matopiba,
'R20mm': rcm_r20mm_rcp26_matopiba}

rcm_rcp85_matopiba = {'PRCPTOT': rcm_prcptot_rcp85_matopiba,
'R95p': rcm_r95p_rcp85_matopiba,
'R99p': rcm_r99p_rcp85_matopiba,
'Rx1day': rcm_rx1day_rcp85_matopiba,
'Rx5day': rcm_rx5day_rcp85_matopiba,
'SDII': rcm_sdii_rcp85_matopiba,
'CDD': rcm_cdd_rcp85_matopiba,
'CWD': rcm_cwd_rcp85_matopiba, 
'R10mm': rcm_r10mm_rcp85_matopiba,
'R20mm': rcm_r20mm_rcp85_matopiba}

gcm_rcp26_matopiba = {'PRCPTOT': gcm_prcptot_rcp26_matopiba,
'R95p': gcm_r95p_rcp26_matopiba,
'R99p': gcm_r99p_rcp26_matopiba,
'Rx1day': gcm_rx1day_rcp26_matopiba,
'Rx5day': gcm_rx5day_rcp26_matopiba,
'SDII': gcm_sdii_rcp26_matopiba,
'CDD': gcm_cdd_rcp26_matopiba,
'CWD': gcm_cwd_rcp26_matopiba, 
'R10mm': gcm_r10mm_rcp26_matopiba,
'R20mm': gcm_r20mm_rcp26_matopiba}

gcm_rcp85_matopiba = {'PRCPTOT': gcm_prcptot_rcp85_matopiba,
'R95p': gcm_r95p_rcp85_matopiba,
'R99p': gcm_r99p_rcp85_matopiba,
'Rx1day': gcm_rx1day_rcp85_matopiba,
'Rx5day': gcm_rx5day_rcp85_matopiba,
'SDII': gcm_sdii_rcp85_matopiba,
'CDD': gcm_cdd_rcp85_matopiba,
'CWD': gcm_cwd_rcp85_matopiba, 
'R10mm': gcm_r10mm_rcp85_matopiba,
'R20mm': gcm_r20mm_rcp85_matopiba}

rcm_rcp26_samz = pd.DataFrame(rcm_rcp26_samz,columns=['PRCPTOT','R95p','R99p','Rx1day','Rx5day','SDII','CDD','CWD','R10mm','R20mm'])
rcm_rcp85_samz = pd.DataFrame(rcm_rcp85_samz,columns=['PRCPTOT','R95p','R99p','Rx1day','Rx5day','SDII','CDD','CWD','R10mm','R20mm'])
#~ gcm_rcp26_samz = pd.DataFrame(gcm_rcp26_samz,columns=['PRCPTOT','R95p','R99p','Rx1day','Rx5day','SDII','CDD','CWD','R10mm','R20mm'])
gcm_rcp85_samz = pd.DataFrame(gcm_rcp85_samz,columns=['PRCPTOT','R95p','R99p','Rx1day','Rx5day','SDII','CDD','CWD','R10mm','R20mm'])

rcm_rcp26_eneb = pd.DataFrame(rcm_rcp26_eneb,columns=['PRCPTOT','R95p','R99p','Rx1day','Rx5day','SDII','CDD','CWD','R10mm','R20mm'])
rcm_rcp85_eneb = pd.DataFrame(rcm_rcp85_eneb,columns=['PRCPTOT','R95p','R99p','Rx1day','Rx5day','SDII','CDD','CWD','R10mm','R20mm'])
#~ gcm_rcp26_eneb = pd.DataFrame(gcm_rcp26_eneb,columns=['PRCPTOT','R95p','R99p','Rx1day','Rx5day','SDII','CDD','CWD','R10mm','R20mm'])
gcm_rcp85_eneb = pd.DataFrame(gcm_rcp85_eneb,columns=['PRCPTOT','R95p','R99p','Rx1day','Rx5day','SDII','CDD','CWD','R10mm','R20mm'])

rcm_rcp26_matopiba = pd.DataFrame(rcm_rcp26_matopiba,columns=['PRCPTOT','R95p','R99p','Rx1day','Rx5day','SDII','CDD','CWD','R10mm','R20mm'])
rcm_rcp85_matopiba = pd.DataFrame(rcm_rcp85_matopiba,columns=['PRCPTOT','R95p','R99p','Rx1day','Rx5day','SDII','CDD','CWD','R10mm','R20mm'])
#~ gcm_rcp26_matopiba = pd.DataFrame(gcm_rcp26_matopiba,columns=['PRCPTOT','R95p','R99p','Rx1day','Rx5day','SDII','CDD','CWD','R10mm','R20mm'])
gcm_rcp85_matopiba = pd.DataFrame(gcm_rcp85_matopiba,columns=['PRCPTOT','R95p','R99p','Rx1day','Rx5day','SDII','CDD','CWD','R10mm','R20mm'])

rcm_rcp26_samz_corr = rcm_rcp26_samz.corr()
rcm_rcp85_samz_corr = rcm_rcp85_samz.corr()
#~ gcm_rcp26_samz_corr = gcm_rcp26_samz.corr()
gcm_rcp85_samz_corr = gcm_rcp85_samz.corr()

rcm_rcp26_eneb_corr = rcm_rcp26_eneb.corr()
rcm_rcp85_eneb_corr = rcm_rcp85_eneb.corr()
#~ gcm_rcp26_eneb_corr = gcm_rcp26_eneb.corr()
gcm_rcp85_eneb_corr = gcm_rcp85_eneb.corr()

rcm_rcp26_matopiba_corr = rcm_rcp26_matopiba.corr()
rcm_rcp85_matopiba_corr = rcm_rcp85_matopiba.corr()
#~ gcm_rcp26_matopiba_corr = gcm_rcp26_matopiba.corr()
gcm_rcp85_matopiba_corr = gcm_rcp85_matopiba.corr()

# Plot correlation matrix 
fig = plt.figure(figsize=(9, 10))

ax1 = fig.add_subplot(4, 3, 1)
mask = np.zeros_like(rcm_rcp26_samz_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp26_samz_corr, cmap='BrBG', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5}, linewidths=.6, ax=ax1)
heatmap.set_title('A)', fontdict={'fontsize':8}, loc='left', fontweight='bold')
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = fig.add_subplot(4, 3, 2)
mask = np.zeros_like(rcm_rcp26_eneb_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp26_eneb_corr, cmap='BrBG', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5}, linewidths=.6, ax=ax2)
heatmap.set_title('B)', fontdict={'fontsize':8}, loc='left', fontweight='bold')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)

ax3 = fig.add_subplot(4, 3, 3)
mask = np.zeros_like(rcm_rcp26_matopiba_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp26_matopiba_corr, cmap='BrBG', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5}, linewidths=.6, ax=ax3)
heatmap.set_title('C)', fontdict={'fontsize':8}, loc='left', fontweight='bold')
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)

ax4 = fig.add_subplot(4, 3, 4)
mask = np.zeros_like(rcm_rcp85_samz_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp85_samz_corr, cmap='BrBG', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5}, linewidths=.6, ax=ax4)
heatmap.set_title('D)', fontdict={'fontsize':8}, loc='left', fontweight='bold')
plt.setp(ax4.get_xticklabels(), visible=False)

ax5 = fig.add_subplot(4, 3, 5)
mask = np.zeros_like(rcm_rcp85_eneb_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp85_eneb_corr, cmap='BrBG', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5}, linewidths=.6, ax=ax5)
heatmap.set_title('E)', fontdict={'fontsize':8}, loc='left', fontweight='bold')
plt.setp(ax5.get_xticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)

ax6 = fig.add_subplot(4, 3, 6)
mask = np.zeros_like(rcm_rcp85_matopiba_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp85_matopiba_corr, cmap='BrBG', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5}, linewidths=.6, ax=ax6)
heatmap.set_title('F)', fontdict={'fontsize':8}, loc='left', fontweight='bold')
plt.setp(ax6.get_xticklabels(), visible=False)
plt.setp(ax6.get_yticklabels(), visible=False)

ax7 = fig.add_subplot(4, 3, 7)
mask = np.zeros_like(rcm_rcp26_eneb_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp26_eneb_corr, cmap='BrBG', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5}, linewidths=.6, ax=ax7)
heatmap.set_title('G)', fontdict={'fontsize':8}, loc='left', fontweight='bold')
plt.setp(ax7.get_xticklabels(), visible=False)

ax8 = fig.add_subplot(4, 3, 8)
mask = np.zeros_like(rcm_rcp26_eneb_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp26_eneb_corr, cmap='BrBG', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5}, linewidths=.6, ax=ax8)
heatmap.set_title('H)', fontdict={'fontsize':8}, loc='left', fontweight='bold')
plt.setp(ax8.get_xticklabels(), visible=False)
plt.setp(ax8.get_yticklabels(), visible=False)

ax9 = fig.add_subplot(4, 3, 9)
mask = np.zeros_like(rcm_rcp26_matopiba_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp26_matopiba_corr, cmap='BrBG', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5}, linewidths=.6, ax=ax9)
heatmap.set_title('I)', fontdict={'fontsize':8}, loc='left', fontweight='bold')
plt.setp(ax9.get_xticklabels(), visible=False)
plt.setp(ax9.get_yticklabels(), visible=False)

ax10 = fig.add_subplot(4, 3, 10)
mask = np.zeros_like(gcm_rcp85_samz_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(gcm_rcp85_samz_corr, cmap='BrBG', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5}, linewidths=.6, ax=ax10)
heatmap.set_title('J)', fontdict={'fontsize':8}, loc='left', fontweight='bold')

ax11 = fig.add_subplot(4, 3, 11)
mask = np.zeros_like(gcm_rcp85_eneb_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(gcm_rcp85_eneb_corr, cmap='BrBG', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5}, linewidths=.6, ax=ax11)
heatmap.set_title('K)', fontdict={'fontsize':8}, loc='left', fontweight='bold')
plt.setp(ax11.get_yticklabels(), visible=False)

ax12 = fig.add_subplot(4, 3, 12)
mask = np.zeros_like(gcm_rcp85_matopiba_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(gcm_rcp85_matopiba_corr, cmap='BrBG', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5}, linewidths=.6, ax=ax12)
heatmap.set_title('L)', fontdict={'fontsize':8}, loc='left', fontweight='bold')
plt.setp(ax12.get_yticklabels(), visible=False)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_matrix_corr_etccdi_pre_reg_had_rcp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')

plt.close('all')
plt.cla()
exit()	

