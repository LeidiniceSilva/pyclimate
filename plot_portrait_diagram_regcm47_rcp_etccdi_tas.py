# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "06/24/2021"
__description__ = "This script plot portrait diagram from extremes indices to rcp"

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
	annual_rcm = np.nanmean(np.nanmean(var[:][0:20,:,:], axis=1), axis=1)
	annual_rcm = np.nanmean(annual_rcm, axis=1)
	
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
	annual_gcm = np.nanmean(np.nanmean(var[:][0:20,:,:], axis=1), axis=1)

	return annual_gcm
	

# Import extreme indices 
# SAMZ
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

rcm_rcp26_samz = {'TXx': rcm_txx_rcp26_samz,
'TXn': rcm_txn_rcp26_samz,
'TNx': rcm_tnx_rcp26_samz,
'TNn': rcm_tnn_rcp26_samz,
'DTR': rcm_dtr_rcp26_samz,
'SU': rcm_su_rcp26_samz,
'TR': rcm_tr_rcp26_samz,
'TX10p': rcm_tx10p_rcp26_samz, 
'TX90p': rcm_tx90p_rcp26_samz,
'TN10p': rcm_tn10p_rcp26_samz,
'TN90p': rcm_tn90p_rcp26_samz}

rcm_rcp85_samz = {'TXx': rcm_txx_rcp85_samz,
'TXn': rcm_txn_rcp85_samz,
'TNx': rcm_tnx_rcp85_samz,
'TNn': rcm_tnn_rcp85_samz,
'DTR': rcm_dtr_rcp85_samz,
'SU': rcm_su_rcp85_samz,
'TR': rcm_tr_rcp85_samz,
'TX10p': rcm_tx10p_rcp85_samz, 
'TX90p': rcm_tx90p_rcp85_samz,
'TN10p': rcm_tn10p_rcp85_samz,
'TN90p': rcm_tn90p_rcp85_samz}

gcm_rcp26_samz = {'TXx': gcm_txx_rcp26_samz,
'TXn': gcm_txn_rcp26_samz,
'TNx': gcm_tnx_rcp26_samz,
'TNn': gcm_tnn_rcp26_samz,
'DTR': gcm_dtr_rcp26_samz,
'SU': gcm_su_rcp26_samz,
'TR': gcm_tr_rcp26_samz,
'TX10p': gcm_tx10p_rcp26_samz, 
'TX90p': gcm_tx90p_rcp26_samz,
'TN10p': gcm_tn10p_rcp26_samz,
'TN90p': gcm_tn90p_rcp26_samz}

gcm_rcp85_samz = {'TXx': gcm_txx_rcp85_samz,
'TXn': gcm_txn_rcp85_samz,
'TNx': gcm_tnx_rcp85_samz,
'TNn': gcm_tnn_rcp85_samz,
'DTR': gcm_dtr_rcp85_samz,
'SU': gcm_su_rcp85_samz,
'TR': gcm_tr_rcp85_samz,
'TX10p': gcm_tx10p_rcp85_samz, 
'TX90p': gcm_tx90p_rcp85_samz,
'TN10p': gcm_tn10p_rcp85_samz,
'TN90p': gcm_tn90p_rcp85_samz}

rcm_rcp26_eneb = {'TXx': rcm_txx_rcp26_eneb,
'TXn': rcm_txn_rcp26_eneb,
'TNx': rcm_tnx_rcp26_eneb,
'TNn': rcm_tnn_rcp26_eneb,
'DTR': rcm_dtr_rcp26_eneb,
'SU': rcm_su_rcp26_eneb,
'TR': rcm_tr_rcp26_eneb,
'TX10p': rcm_tx10p_rcp26_eneb, 
'TX90p': rcm_tx90p_rcp26_eneb,
'TN10p': rcm_tn10p_rcp26_eneb,
'TN90p': rcm_tn90p_rcp26_eneb}

rcm_rcp85_eneb = {'TXx': rcm_txx_rcp85_eneb,
'TXn': rcm_txn_rcp85_eneb,
'TNx': rcm_tnx_rcp85_eneb,
'TNn': rcm_tnn_rcp85_eneb,
'DTR': rcm_dtr_rcp85_eneb,
'SU': rcm_su_rcp85_eneb,
'TR': rcm_tr_rcp85_eneb,
'TX10p': rcm_tx10p_rcp85_eneb, 
'TX90p': rcm_tx90p_rcp85_eneb,
'TN10p': rcm_tn10p_rcp85_eneb,
'TN90p': rcm_tn90p_rcp85_eneb}

gcm_rcp26_eneb = {'TXx': gcm_txx_rcp26_eneb,
'TXn': gcm_txn_rcp26_eneb,
'TNx': gcm_tnx_rcp26_eneb,
'TNn': gcm_tnn_rcp26_eneb,
'DTR': gcm_dtr_rcp26_eneb,
'SU': gcm_su_rcp26_eneb,
'TR': gcm_tr_rcp26_eneb,
'TX10p': gcm_tx10p_rcp26_eneb, 
'TX90p': gcm_tx90p_rcp26_eneb,
'TN10p': gcm_tn10p_rcp26_eneb,
'TN90p': gcm_tn90p_rcp26_eneb}

gcm_rcp85_eneb = {'TXx': gcm_txx_rcp85_eneb,
'TXn': gcm_txn_rcp85_eneb,
'TNx': gcm_tnx_rcp85_eneb,
'TNn': gcm_tnn_rcp85_eneb,
'DTR': gcm_dtr_rcp85_eneb,
'SU': gcm_su_rcp85_eneb,
'TR': gcm_tr_rcp85_eneb,
'TX10p': gcm_tx10p_rcp85_eneb, 
'TX90p': gcm_tx90p_rcp85_eneb,
'TN10p': gcm_tn10p_rcp85_eneb,
'TN90p': gcm_tn90p_rcp85_eneb}

rcm_rcp26_matopiba = {'TXx': rcm_txx_rcp26_matopiba,
'TXn': rcm_txn_rcp26_matopiba,
'TNx': rcm_tnx_rcp26_matopiba,
'TNn': rcm_tnn_rcp26_matopiba,
'DTR': rcm_dtr_rcp26_matopiba,
'SU': rcm_su_rcp26_matopiba,
'TR': rcm_tr_rcp26_matopiba,
'TX10p': rcm_tx10p_rcp26_matopiba, 
'TX90p': rcm_tx90p_rcp26_matopiba,
'TN10p': rcm_tn10p_rcp26_matopiba,
'TN90p': rcm_tn90p_rcp26_matopiba}

rcm_rcp85_matopiba = {'TXx': rcm_txx_rcp85_matopiba,
'TXn': rcm_txn_rcp85_matopiba,
'TNx': rcm_tnx_rcp85_matopiba,
'TNn': rcm_tnn_rcp85_matopiba,
'DTR': rcm_dtr_rcp85_matopiba,
'SU': rcm_su_rcp85_matopiba,
'TR': rcm_tr_rcp85_matopiba,
'TX10p': rcm_tx10p_rcp85_matopiba, 
'TX90p': rcm_tx90p_rcp85_matopiba,
'TN10p': rcm_tn10p_rcp85_matopiba,
'TN90p': rcm_tn90p_rcp85_matopiba}

gcm_rcp26_matopiba = {'TXx': gcm_txx_rcp26_matopiba,
'TXn': gcm_txn_rcp26_matopiba,
'TNx': gcm_tnx_rcp26_matopiba,
'TNn': gcm_tnn_rcp26_matopiba,
'DTR': gcm_dtr_rcp26_matopiba,
'SU': gcm_su_rcp26_matopiba,
'TR': gcm_tr_rcp26_matopiba,
'TX10p': gcm_tx10p_rcp26_matopiba, 
'TX90p': gcm_tx90p_rcp26_matopiba,
'TN10p': gcm_tn10p_rcp26_matopiba,
'TN90p': gcm_tn90p_rcp26_matopiba}

gcm_rcp85_matopiba = {'TXx': gcm_txx_rcp85_matopiba,
'TXn': gcm_txn_rcp85_matopiba,
'TNx': gcm_tnx_rcp85_matopiba,
'TNn': gcm_tnn_rcp85_matopiba,
'DTR': gcm_dtr_rcp85_matopiba,
'SU': gcm_su_rcp85_matopiba,
'TR': gcm_tr_rcp85_matopiba,
'TX10p': gcm_tx10p_rcp85_matopiba, 
'TX90p': gcm_tx90p_rcp85_matopiba,
'TN10p': gcm_tn10p_rcp85_matopiba,
'TN90p': gcm_tn90p_rcp85_matopiba}

rcm_rcp26_samz = pd.DataFrame(rcm_rcp26_samz,columns=['TXx','TXn','TNx','TNn','DTR','SU','TR','TX10p','TX90p','TN10p','TN90p'])
rcm_rcp85_samz = pd.DataFrame(rcm_rcp85_samz,columns=['TXx','TXn','TNx','TNn','DTR','SU','TR','TX10p','TX90p','TN10p','TN90p'])
gcm_rcp26_samz = pd.DataFrame(gcm_rcp26_samz,columns=['TXx','TXn','TNx','TNn','DTR','SU','TR','TX10p','TX90p','TN10p','TN90p'])
gcm_rcp85_samz = pd.DataFrame(gcm_rcp85_samz,columns=['TXx','TXn','TNx','TNn','DTR','SU','TR','TX10p','TX90p','TN10p','TN90p'])

rcm_rcp26_eneb = pd.DataFrame(rcm_rcp26_eneb,columns=['TXx','TXn','TNx','TNn','DTR','SU','TR','TX10p','TX90p','TN10p','TN90p'])
rcm_rcp85_eneb = pd.DataFrame(rcm_rcp85_eneb,columns=['TXx','TXn','TNx','TNn','DTR','SU','TR','TX10p','TX90p','TN10p','TN90p'])
gcm_rcp26_eneb = pd.DataFrame(gcm_rcp26_eneb,columns=['TXx','TXn','TNx','TNn','DTR','SU','TR','TX10p','TX90p','TN10p','TN90p'])
gcm_rcp85_eneb = pd.DataFrame(gcm_rcp85_eneb,columns=['TXx','TXn','TNx','TNn','DTR','SU','TR','TX10p','TX90p','TN10p','TN90p'])

rcm_rcp26_matopiba = pd.DataFrame(rcm_rcp26_matopiba,columns=['TXx','TXn','TNx','TNn','DTR','SU','TR','TX10p','TX90p','TN10p','TN90p'])
rcm_rcp85_matopiba = pd.DataFrame(rcm_rcp85_matopiba,columns=['TXx','TXn','TNx','TNn','DTR','SU','TR','TX10p','TX90p','TN10p','TN90p'])
gcm_rcp26_matopiba = pd.DataFrame(gcm_rcp26_matopiba,columns=['TXx','TXn','TNx','TNn','DTR','SU','TR','TX10p','TX90p','TN10p','TN90p'])
gcm_rcp85_matopiba = pd.DataFrame(gcm_rcp85_matopiba,columns=['TXx','TXn','TNx','TNn','DTR','SU','TR','TX10p','TX90p','TN10p','TN90p'])

rcm_rcp26_samz_corr = rcm_rcp26_samz.corr()
rcm_rcp85_samz_corr = rcm_rcp85_samz.corr()
gcm_rcp26_samz_corr = gcm_rcp26_samz.corr()
gcm_rcp85_samz_corr = gcm_rcp85_samz.corr()

rcm_rcp26_eneb_corr = rcm_rcp26_eneb.corr()
rcm_rcp85_eneb_corr = rcm_rcp85_eneb.corr()
gcm_rcp26_eneb_corr = gcm_rcp26_eneb.corr()
gcm_rcp85_eneb_corr = gcm_rcp85_eneb.corr()

rcm_rcp26_matopiba_corr = rcm_rcp26_matopiba.corr()
rcm_rcp85_matopiba_corr = rcm_rcp85_matopiba.corr()
gcm_rcp26_matopiba_corr = gcm_rcp26_matopiba.corr()
gcm_rcp85_matopiba_corr = gcm_rcp85_matopiba.corr()

# Plot extreme indices  
fig = plt.figure(figsize=(9, 10))

ax1 = fig.add_subplot(4, 3, 1)
mask = np.zeros_like(rcm_rcp26_samz_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp26_samz_corr, cmap='bwr', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5}, linewidths=.6, ax=ax1)
heatmap.set_title('A)', fontdict={'fontsize':9}, loc='left', fontweight='bold')
heatmap.set_ylabel('RegCM4.7 RCP2.6', fontdict={'fontsize':9}, fontweight='bold')
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = fig.add_subplot(4, 3, 2)
mask = np.zeros_like(rcm_rcp26_eneb_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp26_eneb_corr, cmap='bwr', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5},linewidths=.6, ax=ax2)
heatmap.set_title('B)', fontdict={'fontsize':9}, loc='left', fontweight='bold')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)

ax3 = fig.add_subplot(4, 3, 3)
mask = np.zeros_like(rcm_rcp26_matopiba_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp26_matopiba_corr, cmap='bwr', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5},linewidths=.6, ax=ax3)
heatmap.set_title('C)', fontdict={'fontsize':9}, loc='left', fontweight='bold')
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)

ax4 = fig.add_subplot(4, 3, 4)
mask = np.zeros_like(rcm_rcp85_samz_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp85_samz_corr, cmap='bwr', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5},linewidths=.6, ax=ax4)
heatmap.set_title('D)', fontdict={'fontsize':9}, loc='left', fontweight='bold')
heatmap.set_ylabel('RegCM4.7 RCP8.5', fontdict={'fontsize':9}, fontweight='bold')
plt.setp(ax4.get_xticklabels(), visible=False)

ax5 = fig.add_subplot(4, 3, 5)
mask = np.zeros_like(rcm_rcp85_eneb_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp85_eneb_corr, cmap='bwr', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5},linewidths=.6, ax=ax5)
heatmap.set_title('E)', fontdict={'fontsize':9}, loc='left', fontweight='bold')
plt.setp(ax5.get_xticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)

ax6 = fig.add_subplot(4, 3, 6)
mask = np.zeros_like(rcm_rcp85_matopiba_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(rcm_rcp85_matopiba_corr, cmap='bwr', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5},linewidths=.6, ax=ax6)
heatmap.set_title('F)', fontdict={'fontsize':9}, loc='left', fontweight='bold')
plt.setp(ax6.get_xticklabels(), visible=False)
plt.setp(ax6.get_yticklabels(), visible=False)

ax7 = fig.add_subplot(4, 3, 7)
mask = np.zeros_like(gcm_rcp26_samz_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(gcm_rcp26_samz_corr, cmap='bwr', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5},linewidths=.6, ax=ax7)
heatmap.set_title('G)', fontdict={'fontsize':9}, loc='left', fontweight='bold')
heatmap.set_ylabel('HadGEM2.6 RCP2.6', fontdict={'fontsize':9}, fontweight='bold')
plt.setp(ax7.get_xticklabels(), visible=False)

ax8 = fig.add_subplot(4, 3, 8)
mask = np.zeros_like(gcm_rcp26_eneb_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(gcm_rcp26_eneb_corr, cmap='bwr', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5},linewidths=.6, ax=ax8)
heatmap.set_title('H)', fontdict={'fontsize':9}, loc='left', fontweight='bold')
plt.setp(ax8.get_xticklabels(), visible=False)
plt.setp(ax8.get_yticklabels(), visible=False)

ax9 = fig.add_subplot(4, 3, 9)
mask = np.zeros_like(gcm_rcp26_matopiba_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(gcm_rcp26_matopiba_corr, cmap='bwr', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5},linewidths=.6, ax=ax9)
heatmap.set_title('I)', fontdict={'fontsize':9}, loc='left', fontweight='bold')
plt.setp(ax9.get_xticklabels(), visible=False)
plt.setp(ax9.get_yticklabels(), visible=False)

ax10 = fig.add_subplot(4, 3, 10)
mask = np.zeros_like(gcm_rcp85_samz_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(gcm_rcp85_samz_corr, cmap='bwr', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5},linewidths=.6, ax=ax10)
heatmap.set_title('J)', fontdict={'fontsize':9}, loc='left', fontweight='bold')
heatmap.set_ylabel('HadGEM2-ES RCP8.5', fontdict={'fontsize':9}, fontweight='bold')

ax11 = fig.add_subplot(4, 3, 11)
mask = np.zeros_like(gcm_rcp85_eneb_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(gcm_rcp85_eneb_corr, cmap='bwr', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5},linewidths=.6, ax=ax11)
heatmap.set_title('K)', fontdict={'fontsize':9}, loc='left', fontweight='bold')
plt.setp(ax11.get_yticklabels(), visible=False)

ax12 = fig.add_subplot(4, 3, 12)
mask = np.zeros_like(gcm_rcp85_matopiba_corr, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
heatmap = sns.heatmap(gcm_rcp85_matopiba_corr, cmap='bwr', vmin=-1, vmax=1, center=0, mask=mask, annot=True, fmt='.1f', annot_kws={"size":6.5},linewidths=.6, ax=ax12)
heatmap.set_title('L)', fontdict={'fontsize':9}, loc='left', fontweight='bold')
plt.setp(ax12.get_yticklabels(), visible=False)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_matrix_corr_etccdi_tas_reg_had_rcp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()	



