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

	dict_var = {u'eca_txx': u'tasmax', 
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

	dict_var = {u'eca_txx': u'tasmax', 
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
	

# Import extreme indices 
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
diff_rcm_rcp26_tn10p_samz = rcm_tn10p_rcp26_samz - rcm_tn10p_hist_samz
diff_rcm_rcp26_tn90p_samz = rcm_tn90p_rcp26_samz - rcm_tn90p_hist_samz

diff_rcm_rcp85_txx_samz = rcm_txx_rcp85_samz - rcm_txx_hist_samz
diff_rcm_rcp85_txn_samz = rcm_txn_rcp85_samz - rcm_txn_hist_samz
diff_rcm_rcp85_tnx_samz = rcm_tnx_rcp85_samz - rcm_tnx_hist_samz
diff_rcm_rcp85_tnn_samz = rcm_tnn_rcp85_samz - rcm_tnn_hist_samz
diff_rcm_rcp85_dtr_samz = rcm_dtr_rcp85_samz - rcm_dtr_hist_samz
diff_rcm_rcp85_su_samz = rcm_su_rcp85_samz - rcm_su_hist_samz
diff_rcm_rcp85_tr_samz = rcm_tr_rcp85_samz - rcm_tr_hist_samz
diff_rcm_rcp85_tx10p_samz = rcm_tx10p_rcp85_samz - rcm_tx10p_hist_samz 
diff_rcm_rcp85_tx90p_samz = rcm_tx90p_rcp85_samz - rcm_tx90p_hist_samz 
diff_rcm_rcp85_tn10p_samz = rcm_tn10p_rcp85_samz - rcm_tn10p_hist_samz
diff_rcm_rcp85_tn90p_samz = rcm_tn90p_rcp85_samz - rcm_tn90p_hist_samz

diff_gcm_rcp26_txx_samz = gcm_txx_rcp26_samz - gcm_txx_hist_samz
diff_gcm_rcp26_txn_samz = gcm_txn_rcp26_samz - gcm_txn_hist_samz
diff_gcm_rcp26_tnx_samz = gcm_tnx_rcp26_samz - gcm_tnx_hist_samz
diff_gcm_rcp26_tnn_samz = gcm_tnn_rcp26_samz - gcm_tnn_hist_samz
diff_gcm_rcp26_dtr_samz = gcm_dtr_rcp26_samz - gcm_dtr_hist_samz
diff_gcm_rcp26_su_samz = gcm_su_rcp26_samz - gcm_su_hist_samz
diff_gcm_rcp26_tr_samz = gcm_tr_rcp26_samz - gcm_tr_hist_samz
diff_gcm_rcp26_tx10p_samz = gcm_tx10p_rcp26_samz - gcm_tx10p_hist_samz
diff_gcm_rcp26_tx90p_samz = gcm_tx90p_rcp26_samz - gcm_tx90p_hist_samz
diff_gcm_rcp26_tn10p_samz = gcm_tn10p_rcp26_samz - gcm_tn10p_hist_samz 
diff_gcm_rcp26_tn90p_samz = gcm_tn90p_rcp26_samz - gcm_tn90p_hist_samz

diff_gcm_rcp85_txx_samz = gcm_txx_rcp85_samz - gcm_txx_hist_samz
diff_gcm_rcp85_txn_samz = gcm_txn_rcp85_samz - gcm_txn_hist_samz
diff_gcm_rcp85_tnx_samz = gcm_tnx_rcp85_samz - gcm_tnx_hist_samz
diff_gcm_rcp85_tnn_samz = gcm_tnn_rcp85_samz - gcm_tnn_hist_samz
diff_gcm_rcp85_dtr_samz = gcm_dtr_rcp85_samz - gcm_dtr_hist_samz
diff_gcm_rcp85_su_samz = gcm_su_rcp85_samz - gcm_su_hist_samz
diff_gcm_rcp85_tr_samz = gcm_tr_rcp85_samz - gcm_tr_hist_samz
diff_gcm_rcp85_tx10p_samz = gcm_tx10p_rcp85_samz - gcm_tx10p_hist_samz
diff_gcm_rcp85_tx90p_samz = gcm_tx90p_rcp85_samz - gcm_tx90p_hist_samz 
diff_gcm_rcp85_tn10p_samz = gcm_tn10p_rcp85_samz - gcm_tn10p_hist_samz 
diff_gcm_rcp85_tn90p_samz = gcm_tn90p_rcp85_samz - gcm_tn90p_hist_samz 

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
diff_rcm_rcp26_tn10p_eneb = rcm_tn10p_rcp26_eneb - rcm_tn10p_hist_eneb
diff_rcm_rcp26_tn90p_eneb = rcm_tn90p_rcp26_eneb - rcm_tn90p_hist_eneb

diff_rcm_rcp85_txx_eneb = rcm_txx_rcp85_eneb - rcm_txx_hist_eneb
diff_rcm_rcp85_txn_eneb = rcm_txn_rcp85_eneb - rcm_txn_hist_eneb
diff_rcm_rcp85_tnx_eneb = rcm_tnx_rcp85_eneb - rcm_tnx_hist_eneb
diff_rcm_rcp85_tnn_eneb = rcm_tnn_rcp85_eneb - rcm_tnn_hist_eneb
diff_rcm_rcp85_dtr_eneb = rcm_dtr_rcp85_eneb - rcm_dtr_hist_eneb
diff_rcm_rcp85_su_eneb = rcm_su_rcp85_eneb - rcm_su_hist_eneb
diff_rcm_rcp85_tr_eneb = rcm_tr_rcp85_eneb - rcm_tr_hist_eneb
diff_rcm_rcp85_tx10p_eneb = rcm_tx10p_rcp85_eneb - rcm_tx10p_hist_eneb
diff_rcm_rcp85_tx90p_eneb = rcm_tx90p_rcp85_eneb - rcm_tx90p_hist_eneb 
diff_rcm_rcp85_tn10p_eneb = rcm_tn10p_rcp85_eneb - rcm_tn10p_hist_eneb
diff_rcm_rcp85_tn90p_eneb = rcm_tn90p_rcp85_eneb - rcm_tn90p_hist_eneb

diff_gcm_rcp26_txx_eneb = gcm_txx_rcp26_eneb - gcm_txx_hist_eneb
diff_gcm_rcp26_txn_eneb = gcm_txn_rcp26_eneb - gcm_txn_hist_eneb
diff_gcm_rcp26_tnx_eneb = gcm_tnx_rcp26_eneb - gcm_tnx_hist_eneb
diff_gcm_rcp26_tnn_eneb = gcm_tnn_rcp26_eneb - gcm_tnn_hist_eneb
diff_gcm_rcp26_dtr_eneb = gcm_dtr_rcp26_eneb - gcm_dtr_hist_eneb
diff_gcm_rcp26_su_eneb = gcm_su_rcp26_eneb - gcm_su_hist_eneb
diff_gcm_rcp26_tr_eneb = gcm_tr_rcp26_eneb - gcm_tr_hist_eneb
diff_gcm_rcp26_tx10p_eneb = gcm_tx10p_rcp26_eneb - gcm_tx10p_hist_eneb
diff_gcm_rcp26_tx90p_eneb = gcm_tx90p_rcp26_eneb - gcm_tx90p_hist_eneb 
diff_gcm_rcp26_tn10p_eneb = gcm_tn10p_rcp26_eneb - gcm_tn10p_hist_eneb
diff_gcm_rcp26_tn90p_eneb = gcm_tn90p_rcp26_eneb - gcm_tn90p_hist_eneb 

diff_gcm_rcp85_txx_eneb = gcm_txx_rcp85_eneb - gcm_txx_hist_eneb
diff_gcm_rcp85_txn_eneb = gcm_txn_rcp85_eneb - gcm_txn_hist_eneb
diff_gcm_rcp85_tnx_eneb = gcm_tnx_rcp85_eneb - gcm_tnx_hist_eneb
diff_gcm_rcp85_tnn_eneb = gcm_tnn_rcp85_eneb - gcm_tnn_hist_eneb
diff_gcm_rcp85_dtr_eneb = gcm_dtr_rcp85_eneb - gcm_dtr_hist_eneb
diff_gcm_rcp85_su_eneb = gcm_su_rcp85_eneb - gcm_su_hist_eneb
diff_gcm_rcp85_tr_eneb = gcm_tr_rcp85_eneb - gcm_tr_hist_eneb
diff_gcm_rcp85_tx10p_eneb = gcm_tx10p_rcp85_eneb - gcm_tx10p_hist_eneb 
diff_gcm_rcp85_tx90p_eneb = gcm_tx90p_rcp85_eneb - gcm_tx90p_hist_eneb
diff_gcm_rcp85_tn10p_eneb = gcm_tn10p_rcp85_eneb - gcm_tn10p_hist_eneb 
diff_gcm_rcp85_tn90p_eneb = gcm_tn90p_rcp85_eneb - gcm_tn90p_hist_eneb 

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
diff_rcm_rcp26_tn10p_matopiba = rcm_tn10p_rcp26_matopiba - rcm_tn10p_hist_matopiba
diff_rcm_rcp26_tn90p_matopiba = rcm_tn90p_rcp26_matopiba - rcm_tn90p_hist_matopiba

diff_rcm_rcp85_txx_matopiba = rcm_txx_rcp85_matopiba - rcm_txx_hist_matopiba
diff_rcm_rcp85_txn_matopiba = rcm_txn_rcp85_matopiba - rcm_txn_hist_matopiba
diff_rcm_rcp85_tnx_matopiba = rcm_tnx_rcp85_matopiba - rcm_tnx_hist_matopiba
diff_rcm_rcp85_tnn_matopiba = rcm_tnn_rcp85_matopiba - rcm_tnn_hist_matopiba
diff_rcm_rcp85_dtr_matopiba = rcm_dtr_rcp85_matopiba - rcm_dtr_hist_matopiba
diff_rcm_rcp85_su_matopiba = rcm_su_rcp85_matopiba - rcm_su_hist_matopiba
diff_rcm_rcp85_tr_matopiba = rcm_tr_rcp85_matopiba - rcm_tr_hist_matopiba
diff_rcm_rcp85_tx10p_matopiba = rcm_tx10p_rcp85_matopiba - rcm_tx10p_hist_matopiba 
diff_rcm_rcp85_tx90p_matopiba = rcm_tx90p_rcp85_matopiba - rcm_tx90p_hist_matopiba
diff_rcm_rcp85_tn10p_matopiba = rcm_tn10p_rcp85_matopiba - rcm_tn10p_hist_matopiba
diff_rcm_rcp85_tn90p_matopiba = rcm_tn90p_rcp85_matopiba - rcm_tn90p_hist_matopiba

diff_gcm_rcp26_txx_matopiba = gcm_txx_rcp26_matopiba - gcm_txx_hist_matopiba
diff_gcm_rcp26_txn_matopiba = gcm_txn_rcp26_matopiba - gcm_txn_hist_matopiba
diff_gcm_rcp26_tnx_matopiba = gcm_tnx_rcp26_matopiba - gcm_tnx_hist_matopiba
diff_gcm_rcp26_tnn_matopiba = gcm_tnn_rcp26_matopiba- gcm_tnn_hist_matopiba
diff_gcm_rcp26_dtr_matopiba = gcm_dtr_rcp26_matopiba - gcm_dtr_hist_matopiba
diff_gcm_rcp26_su_matopiba = gcm_su_rcp26_matopiba - gcm_su_hist_matopiba
diff_gcm_rcp26_tr_matopiba = gcm_tr_rcp26_matopiba - gcm_tr_hist_matopiba
diff_gcm_rcp26_tx10p_matopiba = gcm_tx10p_rcp26_matopiba - gcm_tx10p_hist_matopiba
diff_gcm_rcp26_tx90p_matopiba = gcm_tx90p_rcp26_matopiba - gcm_tx90p_hist_matopiba 
diff_gcm_rcp26_tn10p_matopiba = gcm_tn10p_rcp26_matopiba - gcm_tn10p_hist_matopiba
diff_gcm_rcp26_tn90p_matopiba = gcm_tn90p_rcp26_matopiba - gcm_tn90p_hist_matopiba 

diff_gcm_rcp85_txx_matopiba = gcm_txx_rcp85_matopiba - gcm_txx_hist_matopiba
diff_gcm_rcp85_txn_matopiba = gcm_txn_rcp85_matopiba - gcm_txn_hist_matopiba
diff_gcm_rcp85_tnx_matopiba = gcm_tnx_rcp85_matopiba - gcm_tnx_hist_matopiba
diff_gcm_rcp85_tnn_matopiba = gcm_tnn_rcp85_matopiba- gcm_tnn_hist_matopiba
diff_gcm_rcp85_dtr_matopiba = gcm_dtr_rcp85_matopiba - gcm_dtr_hist_matopiba
diff_gcm_rcp85_su_matopiba = gcm_su_rcp85_matopiba - gcm_su_hist_matopiba
diff_gcm_rcp85_tr_matopiba = gcm_tr_rcp85_matopiba - gcm_tr_hist_matopiba
diff_gcm_rcp85_tx10p_matopiba = gcm_tx10p_rcp85_matopiba - gcm_tx10p_hist_matopiba
diff_gcm_rcp85_tx90p_matopiba = gcm_tx90p_rcp85_matopiba - gcm_tx90p_hist_matopiba
diff_gcm_rcp85_tn10p_matopiba = gcm_tn10p_rcp85_matopiba - gcm_tn10p_hist_matopiba
diff_gcm_rcp85_tn90p_matopiba = gcm_tn90p_rcp85_matopiba - gcm_tn90p_hist_matopiba 

rcm_rcp26_txx = [diff_rcm_rcp26_txx_samz, diff_rcm_rcp26_txx_eneb, diff_rcm_rcp26_txx_matopiba]
rcm_rcp85_txx = [diff_rcm_rcp85_txx_samz, diff_rcm_rcp85_txx_eneb, diff_rcm_rcp85_txx_matopiba]
gcm_rcp26_txx = [diff_gcm_rcp26_txx_samz, diff_gcm_rcp26_txx_eneb, diff_gcm_rcp26_txx_matopiba]
gcm_rcp85_txx = [diff_gcm_rcp85_txx_samz, diff_gcm_rcp85_txx_eneb, diff_gcm_rcp85_txx_matopiba]

rcm_rcp26_txn = [diff_rcm_rcp26_txn_samz, diff_rcm_rcp26_txn_eneb, diff_rcm_rcp26_txn_matopiba]
rcm_rcp85_txn = [diff_rcm_rcp85_txn_samz, diff_rcm_rcp85_txn_eneb, diff_rcm_rcp85_txn_matopiba]
gcm_rcp26_txn = [diff_gcm_rcp26_txn_samz, diff_gcm_rcp26_txn_eneb, diff_gcm_rcp26_txn_matopiba]
gcm_rcp85_txn = [diff_gcm_rcp85_txn_samz, diff_gcm_rcp85_txn_eneb, diff_gcm_rcp85_txn_matopiba]

rcm_rcp26_tnx = [diff_rcm_rcp26_tnx_samz, diff_rcm_rcp26_tnx_eneb, diff_rcm_rcp26_tnx_matopiba]
rcm_rcp85_tnx = [diff_rcm_rcp85_tnx_samz, diff_rcm_rcp85_tnx_eneb, diff_rcm_rcp85_tnx_matopiba]
gcm_rcp26_tnx = [diff_gcm_rcp26_tnx_samz, diff_gcm_rcp26_tnx_eneb, diff_gcm_rcp26_tnx_matopiba]
gcm_rcp85_tnx = [diff_gcm_rcp85_tnx_samz, diff_gcm_rcp85_tnx_eneb, diff_gcm_rcp85_tnx_matopiba]

rcm_rcp26_tnn = [diff_rcm_rcp26_tnn_samz, diff_rcm_rcp26_tnn_eneb, diff_rcm_rcp26_tnn_matopiba]
rcm_rcp85_tnn = [diff_rcm_rcp85_tnn_samz, diff_rcm_rcp85_tnn_eneb, diff_rcm_rcp85_tnn_matopiba]
gcm_rcp26_tnn = [diff_gcm_rcp26_tnn_samz, diff_gcm_rcp26_tnn_eneb, diff_gcm_rcp26_tnn_matopiba]
gcm_rcp85_tnn = [diff_gcm_rcp85_tnn_samz, diff_gcm_rcp85_tnn_eneb, diff_gcm_rcp85_tnn_matopiba]

rcm_rcp26_dtr = [diff_rcm_rcp26_dtr_samz, diff_rcm_rcp26_dtr_eneb, diff_rcm_rcp26_dtr_matopiba]
rcm_rcp85_dtr = [diff_rcm_rcp85_dtr_samz, diff_rcm_rcp85_dtr_eneb, diff_rcm_rcp85_dtr_matopiba]
gcm_rcp26_dtr = [diff_gcm_rcp26_dtr_samz, diff_gcm_rcp26_dtr_eneb, diff_gcm_rcp26_dtr_matopiba]
gcm_rcp85_dtr = [diff_gcm_rcp85_dtr_samz, diff_gcm_rcp85_dtr_eneb, diff_gcm_rcp85_dtr_matopiba]

rcm_rcp26_su = [diff_rcm_rcp26_su_samz, diff_rcm_rcp26_su_eneb, diff_rcm_rcp26_su_matopiba]
rcm_rcp85_su = [diff_rcm_rcp85_su_samz, diff_rcm_rcp85_su_eneb, diff_rcm_rcp85_su_matopiba]
gcm_rcp26_su = [diff_gcm_rcp26_su_samz, diff_gcm_rcp26_su_eneb, diff_gcm_rcp26_su_matopiba]
gcm_rcp85_su = [diff_gcm_rcp85_su_samz, diff_gcm_rcp85_su_eneb, diff_gcm_rcp85_su_matopiba]

rcm_rcp26_tr = [diff_rcm_rcp26_tr_samz, diff_rcm_rcp26_tr_eneb, diff_rcm_rcp26_tr_matopiba]
rcm_rcp85_tr = [diff_rcm_rcp85_tr_samz, diff_rcm_rcp85_tr_eneb, diff_rcm_rcp85_tr_matopiba]
gcm_rcp26_tr = [diff_gcm_rcp26_tr_samz, diff_gcm_rcp26_tr_eneb, diff_gcm_rcp26_tr_matopiba]
gcm_rcp85_tr = [diff_gcm_rcp85_tr_samz, diff_gcm_rcp85_tr_eneb, diff_gcm_rcp85_tr_matopiba]

rcm_rcp26_tx10p = [diff_rcm_rcp26_tx10p_samz, diff_rcm_rcp26_tx10p_eneb, diff_rcm_rcp26_tx10p_matopiba]
rcm_rcp85_tx10p = [diff_rcm_rcp85_tx10p_samz, diff_rcm_rcp85_tx10p_eneb, diff_rcm_rcp85_tx10p_matopiba]
gcm_rcp26_tx10p = [diff_gcm_rcp26_tx10p_samz, diff_gcm_rcp26_tx10p_eneb, diff_gcm_rcp26_tx10p_matopiba]
gcm_rcp85_tx10p = [diff_gcm_rcp85_tx10p_samz, diff_gcm_rcp85_tx10p_eneb, diff_gcm_rcp85_tx10p_matopiba]

rcm_rcp26_tx90p = [diff_rcm_rcp26_tx90p_samz, diff_rcm_rcp26_tx90p_eneb, diff_rcm_rcp26_tx90p_matopiba]
rcm_rcp85_tx90p = [diff_rcm_rcp85_tx90p_samz, diff_rcm_rcp85_tx90p_eneb, diff_rcm_rcp85_tx90p_matopiba]
gcm_rcp26_tx90p = [diff_gcm_rcp26_tx90p_samz, diff_gcm_rcp26_tx90p_eneb, diff_gcm_rcp26_tx90p_matopiba]
gcm_rcp85_tx90p = [diff_gcm_rcp85_tx90p_samz, diff_gcm_rcp85_tx90p_eneb, diff_gcm_rcp85_tx90p_matopiba]

rcm_rcp26_tn10p = [diff_rcm_rcp26_tn10p_samz, diff_rcm_rcp26_tn10p_eneb, diff_rcm_rcp26_tn10p_matopiba]
rcm_rcp85_tn10p = [diff_rcm_rcp85_tn10p_samz, diff_rcm_rcp85_tn10p_eneb, diff_rcm_rcp85_tn10p_matopiba]
gcm_rcp26_tn10p = [diff_gcm_rcp26_tn10p_samz, diff_gcm_rcp26_tn10p_eneb, diff_gcm_rcp26_tn10p_matopiba]
gcm_rcp85_tn10p = [diff_gcm_rcp85_tn10p_samz, diff_gcm_rcp85_tn10p_eneb, diff_gcm_rcp85_tn10p_matopiba]

rcm_rcp26_tn90p = [diff_rcm_rcp26_tn90p_samz, diff_rcm_rcp26_tn90p_eneb, diff_rcm_rcp26_tn90p_matopiba]
rcm_rcp85_tn90p = [diff_rcm_rcp85_tn90p_samz, diff_rcm_rcp85_tn90p_eneb, diff_rcm_rcp85_tn90p_matopiba]
gcm_rcp26_tn90p = [diff_gcm_rcp26_tn90p_samz, diff_gcm_rcp26_tn90p_eneb, diff_gcm_rcp26_tn90p_matopiba]
gcm_rcp85_tn90p = [diff_gcm_rcp85_tn90p_samz, diff_gcm_rcp85_tn90p_eneb, diff_gcm_rcp85_tn90p_matopiba]

# Plot extreme indices 
fig = plt.figure()
time = np.arange(1, 4)

ax1 = fig.add_subplot(4, 3, 1)
plt1 = plt.plot(time, rcm_rcp26_txx, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_txx, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP2.6')
plt3 = plt.plot(time, gcm_rcp26_txx, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_txx, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP2.6')
plt.title(u'A) TXx (°C)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-10, 10)
plt.yticks(np.arange(-10, 15, 5), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax2 = fig.add_subplot(4, 3, 2)
plt1 = plt.plot(time, rcm_rcp26_txn, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_txn, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_txn, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_txn, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'B) TXn (°C)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-10, 10)
plt.yticks(np.arange(-10, 15, 5), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax3 = fig.add_subplot(4, 3, 3)
plt1 = plt.plot(time, rcm_rcp26_tnx, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_tnx, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP2.6')
plt3 = plt.plot(time, gcm_rcp26_tnx, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP8.5')
plt4 = plt.plot(time, gcm_rcp85_tnx, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'C) TNx (°C)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-10, 10)
plt.yticks(np.arange(-10, 15, 5), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax3.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax4 = fig.add_subplot(4, 3, 4)
plt1 = plt.plot(time, rcm_rcp26_tnn, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_tnn, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_tnn, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_tnn, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'D) TNn (°C)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-10, 10)
plt.yticks(np.arange(-10, 15, 5), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax4.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax5 = fig.add_subplot(4, 3, 5)
plt1 = plt.plot(time, rcm_rcp26_dtr, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_dtr, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_dtr, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_dtr, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'E) DTR (°C)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-2, 2)
plt.yticks(np.arange(-2, 3, 1), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax5.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax6 = fig.add_subplot(4, 3, 6)
plt1 = plt.plot(time, rcm_rcp26_su, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_su, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_su, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_su, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'F) SU (days)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-20, 20)
plt.yticks(np.arange(-20, 30, 10), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax6.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax7 = fig.add_subplot(4, 3, 7)
plt1 = plt.plot(time, rcm_rcp26_tr, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_tr, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_tr, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_tr, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'G) TR (days)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-100, 100)
plt.yticks(np.arange(-100, 150, 50), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax7.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax8 = fig.add_subplot(4, 3, 8)
plt1 = plt.plot(time, rcm_rcp26_tx10p, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_tx10p, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_tx10p, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_tx10p, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'H) TX10p (%)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-50, 50)
plt.yticks(np.arange(-50, 75, 25), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax8.get_xticklabels(), visible=False)
plt.tick_params(bottom = False)
plt.grid(True, which='major', linestyle='--')

ax9 = fig.add_subplot(4, 3, 9)
plt1 = plt.plot(time, rcm_rcp26_tx90p, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_tx90p, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_tx90p, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_tx90p, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'I) TX90p (%)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-100, 100)
plt.yticks(np.arange(-100, 150, 50), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

ax10 = fig.add_subplot(4, 3, 10)
plt1 = plt.plot(time, rcm_rcp26_tn10p, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_tn10p, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_tn10p, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_tn10p, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'J)  TN10p (%)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-50, 50)
plt.yticks(np.arange(-50, 75, 25), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

ax10 = fig.add_subplot(4, 3, 11)
plt1 = plt.plot(time, rcm_rcp26_tn90p, 'o', color='black', markerfacecolor='blue', label='RegCM4.7 RCP2.6')
plt2 = plt.plot(time, rcm_rcp85_tn90p, 'o', color='black', markerfacecolor='red',  label='RegCM4.7 RCP8.5')
plt3 = plt.plot(time, gcm_rcp26_tn90p, 's', color='black', markerfacecolor='blue', label='HadGEM2-ES RCP2.6')
plt4 = plt.plot(time, gcm_rcp85_tn90p, 's', color='black', markerfacecolor='red',  label='HadGEM2-ES RCP8.5')
plt.title(u'K)  TN90p (%)', loc='left', fontweight='bold', fontsize=7)
plt.ylim(-100, 100)
plt.yticks(np.arange(-100, 150, 50), fontsize=7)
plt.xticks(time, ('SAMZ', 'ENEB', 'MATOPIBA'), fontsize=7)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
legend(bbox_to_anchor=(1.05, 1.), loc=2, ncol=1, fontsize=7, frameon=False)
plt.grid(True, which='major', linestyle='--')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_diff_etccdi_reg_had_rcp_tas.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()	
