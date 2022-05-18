# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/08/2019"
__description__ = "This script compute statistical indices from extremes indices"

import os
import netCDF4
import numpy as np
from comp_statist_indices import compute_corr, compute_rmse


def import_obs(var, area, dataset, freq, dt):
	
	path = '/home/nice/Documents/dataset/obs/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, var, area, dataset, freq, dt)	

	dict_var = {u'eca_txx': u'tmax', 
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
	annual_obs = np.nanmean(np.nanmean(var[:][0:20,:,:], axis=1), axis=1)
	
	return annual_obs
	
	
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
	
# Compute correlation 
r_samz_rcm_txx = compute_corr(obs_txx_samz, np.nanmean(rcm_txx_samz, axis=1))
r_samz_rcm_txn = compute_corr(obs_txn_samz, np.nanmean(rcm_txn_samz, axis=1))
r_samz_rcm_tnx = compute_corr(obs_tnx_samz, np.nanmean(rcm_tnx_samz, axis=1))
r_samz_rcm_tnn = compute_corr(obs_tnn_samz, np.nanmean(rcm_tnn_samz, axis=1))
r_samz_rcm_dtr = compute_corr(obs_dtr_samz, np.nanmean(rcm_dtr_samz, axis=1))
r_samz_rcm_su = compute_corr(obs_su_samz, np.nanmean(rcm_su_samz, axis=1))
r_samz_rcm_tr = compute_corr(obs_tr_samz, np.nanmean(rcm_tr_samz, axis=1))
r_samz_rcm_tx10p = compute_corr(obs_tx10p_samz, np.nanmean(rcm_tx10p_samz, axis=1))
r_samz_rcm_tx90p = compute_corr(obs_tx90p_samz, np.nanmean(rcm_tx90p_samz, axis=1))
r_samz_rcm_tn10p = compute_corr(obs_tn10p_samz, np.nanmean(rcm_tn10p_samz, axis=1))
r_samz_rcm_tn90p = compute_corr(obs_tn90p_samz, np.nanmean(rcm_tn90p_samz, axis=1))

r_samz_gcm_txx = compute_corr(obs_txx_samz, gcm_txx_samz)
r_samz_gcm_txn = compute_corr(obs_txn_samz, gcm_txn_samz)
r_samz_gcm_tnx = compute_corr(obs_tnx_samz, gcm_tnx_samz)
r_samz_gcm_tnn = compute_corr(obs_tnn_samz, gcm_tnn_samz)
r_samz_gcm_dtr = compute_corr(obs_dtr_samz, gcm_dtr_samz)
r_samz_gcm_su = compute_corr(obs_su_samz, gcm_su_samz)
r_samz_gcm_tr = compute_corr(obs_tr_samz, gcm_tr_samz)
r_samz_gcm_tx10p = compute_corr(obs_tx10p_samz, gcm_tx10p_samz)
r_samz_gcm_tx90p = compute_corr(obs_tx90p_samz, gcm_tx90p_samz)
r_samz_gcm_tn10p = compute_corr(obs_tn10p_samz, gcm_tn10p_samz)
r_samz_gcm_tn90p = compute_corr(obs_tn90p_samz, gcm_tn90p_samz)

print("r_samz_rcm")
print(r_samz_rcm_txx, r_samz_rcm_txn, r_samz_rcm_tnx, r_samz_rcm_tnn, r_samz_rcm_dtr, r_samz_rcm_su,
r_samz_rcm_tr, r_samz_rcm_tx10p, r_samz_rcm_tx90p, r_samz_rcm_tn10p, r_samz_rcm_tn90p)
print("r_samz_gcm")
print(r_samz_gcm_txx, r_samz_gcm_txn, r_samz_gcm_tnx, r_samz_gcm_tnn, r_samz_gcm_dtr, r_samz_gcm_su,
r_samz_gcm_tr, r_samz_gcm_tx10p, r_samz_gcm_tx90p, r_samz_gcm_tn10p, r_samz_gcm_tn90p)

r_eneb_rcm_txx = compute_corr(obs_txx_eneb, np.nanmean(rcm_txx_eneb, axis=1))
r_eneb_rcm_txn = compute_corr(obs_txn_eneb, np.nanmean(rcm_txn_eneb, axis=1))
r_eneb_rcm_tnx = compute_corr(obs_tnx_eneb, np.nanmean(rcm_tnx_eneb, axis=1))
r_eneb_rcm_tnn = compute_corr(obs_tnn_eneb, np.nanmean(rcm_tnn_eneb, axis=1))
r_eneb_rcm_dtr = compute_corr(obs_dtr_eneb, np.nanmean(rcm_dtr_eneb, axis=1))
r_eneb_rcm_su = compute_corr(obs_su_eneb, np.nanmean(rcm_su_eneb, axis=1))
r_eneb_rcm_tr = compute_corr(obs_tr_eneb, np.nanmean(rcm_tr_eneb, axis=1))
r_eneb_rcm_tx10p = compute_corr(obs_tx10p_eneb, np.nanmean(rcm_tx10p_eneb, axis=1))
r_eneb_rcm_tx90p = compute_corr(obs_tx90p_eneb, np.nanmean(rcm_tx90p_eneb, axis=1))
r_eneb_rcm_tn10p = compute_corr(obs_tn10p_eneb, np.nanmean(rcm_tn10p_eneb, axis=1))
r_eneb_rcm_tn90p = compute_corr(obs_tn90p_eneb, np.nanmean(rcm_tn90p_eneb, axis=1))

r_eneb_gcm_txx = compute_corr(obs_txx_eneb, gcm_txx_eneb)
r_eneb_gcm_txn = compute_corr(obs_txn_eneb, gcm_txn_eneb)
r_eneb_gcm_tnx = compute_corr(obs_tnx_eneb, gcm_tnx_eneb)
r_eneb_gcm_tnn = compute_corr(obs_tnn_eneb, gcm_tnn_eneb)
r_eneb_gcm_dtr = compute_corr(obs_dtr_eneb, gcm_dtr_eneb)
r_eneb_gcm_su = compute_corr(obs_su_eneb, gcm_su_eneb)
r_eneb_gcm_tr = compute_corr(obs_tr_eneb, gcm_tr_eneb)
r_eneb_gcm_tx10p = compute_corr(obs_tx10p_eneb, gcm_tx10p_eneb)
r_eneb_gcm_tx90p = compute_corr(obs_tx90p_eneb, gcm_tx90p_eneb)
r_eneb_gcm_tn10p = compute_corr(obs_tn10p_eneb, gcm_tn10p_eneb)
r_eneb_gcm_tn90p = compute_corr(obs_tn90p_eneb, gcm_tn90p_eneb)

print("r_eneb_rcm")
print(r_eneb_rcm_txx, r_eneb_rcm_txn, r_eneb_rcm_tnx, r_eneb_rcm_tnn, r_eneb_rcm_dtr, r_eneb_rcm_su,
r_eneb_rcm_tr, r_eneb_rcm_tx10p, r_eneb_rcm_tx90p, r_eneb_rcm_tn10p, r_eneb_rcm_tn90p)
print("r_eneb_gcm")
print(r_eneb_gcm_txx, r_eneb_gcm_txn, r_eneb_gcm_tnx, r_eneb_gcm_tnn, r_eneb_gcm_dtr, r_eneb_gcm_su,
r_eneb_gcm_tr, r_eneb_gcm_tx10p, r_eneb_gcm_tx90p, r_eneb_gcm_tn10p, r_eneb_gcm_tn90p)

r_matopiba_rcm_txx = compute_corr(obs_txx_matopiba, np.nanmean(rcm_txx_matopiba, axis=1))
r_matopiba_rcm_txn = compute_corr(obs_txn_matopiba, np.nanmean(rcm_txn_matopiba, axis=1))
r_matopiba_rcm_tnx = compute_corr(obs_tnx_matopiba, np.nanmean(rcm_tnx_matopiba, axis=1))
r_matopiba_rcm_tnn = compute_corr(obs_tnn_matopiba, np.nanmean(rcm_tnn_matopiba, axis=1))
r_matopiba_rcm_dtr = compute_corr(obs_dtr_matopiba, np.nanmean(rcm_dtr_matopiba, axis=1))
r_matopiba_rcm_su = compute_corr(obs_su_matopiba, np.nanmean(rcm_su_matopiba, axis=1))
r_matopiba_rcm_tr = compute_corr(obs_tr_matopiba, np.nanmean(rcm_tr_matopiba, axis=1))
r_matopiba_rcm_tx10p = compute_corr(obs_tx10p_matopiba, np.nanmean(rcm_tx10p_matopiba, axis=1))
r_matopiba_rcm_tx90p = compute_corr(obs_tx90p_matopiba, np.nanmean(rcm_tx90p_matopiba, axis=1))
r_matopiba_rcm_tn10p = compute_corr(obs_tn10p_matopiba, np.nanmean(rcm_tn10p_matopiba, axis=1))
r_matopiba_rcm_tn90p = compute_corr(obs_tn90p_matopiba, np.nanmean(rcm_tn90p_matopiba, axis=1))

r_matopiba_gcm_txx = compute_corr(obs_txx_matopiba, gcm_txx_matopiba)
r_matopiba_gcm_txn = compute_corr(obs_txn_matopiba, gcm_txn_matopiba)
r_matopiba_gcm_tnx = compute_corr(obs_tnx_matopiba, gcm_tnx_matopiba)
r_matopiba_gcm_tnn = compute_corr(obs_tnn_matopiba, gcm_tnn_matopiba)
r_matopiba_gcm_dtr = compute_corr(obs_dtr_matopiba, gcm_dtr_matopiba)
r_matopiba_gcm_su = compute_corr(obs_su_matopiba, gcm_su_matopiba)
r_matopiba_gcm_tr = compute_corr(obs_tr_matopiba, gcm_tr_matopiba)
r_matopiba_gcm_tx10p = compute_corr(obs_tx10p_matopiba, gcm_tx10p_matopiba)
r_matopiba_gcm_tx90p = compute_corr(obs_tx90p_matopiba, gcm_tx90p_matopiba)
r_matopiba_gcm_tn10p = compute_corr(obs_tn10p_matopiba, gcm_tn10p_matopiba)
r_matopiba_gcm_tn90p = compute_corr(obs_tn90p_matopiba, gcm_tn90p_matopiba)

print("r_matopiba_rcm")
print(r_matopiba_rcm_txx, r_matopiba_rcm_txn, r_matopiba_rcm_tnx, r_matopiba_rcm_tnn, r_matopiba_rcm_dtr, r_matopiba_rcm_su,
r_matopiba_rcm_tr, r_matopiba_rcm_tx10p, r_matopiba_rcm_tx90p, r_matopiba_rcm_tn10p, r_matopiba_rcm_tn90p)
print("r_matopiba_gcm")
print(r_matopiba_gcm_txx, r_matopiba_gcm_txn, r_matopiba_gcm_tnx, r_matopiba_gcm_tnn, r_matopiba_gcm_dtr, r_matopiba_gcm_su,
r_matopiba_gcm_tr, r_matopiba_gcm_tx10p, r_matopiba_gcm_tx90p, r_matopiba_gcm_tn10p, r_matopiba_gcm_tn90p)
print()

# Compute root mean square error 
rmse_samz_rcm_txx = compute_rmse(np.nanmean(rcm_txx_samz, axis=1), obs_txx_samz)
rmse_samz_rcm_txn = compute_rmse(np.nanmean(rcm_txn_samz, axis=1), obs_txn_samz)
rmse_samz_rcm_tnx = compute_rmse(np.nanmean(rcm_tnx_samz, axis=1), obs_tnx_samz)
rmse_samz_rcm_tnn = compute_rmse(np.nanmean(rcm_tnn_samz, axis=1), obs_tnn_samz)
rmse_samz_rcm_dtr = compute_rmse(np.nanmean(rcm_dtr_samz, axis=1), obs_dtr_samz)
rmse_samz_rcm_su = compute_rmse(np.nanmean(rcm_su_samz, axis=1), obs_su_samz)
rmse_samz_rcm_tr = compute_rmse(np.nanmean(rcm_tr_samz, axis=1), obs_tr_samz)
rmse_samz_rcm_tx10p = compute_rmse(np.nanmean(rcm_tx10p_samz, axis=1), obs_tx10p_samz)
rmse_samz_rcm_tx90p = compute_rmse(np.nanmean(rcm_tx90p_samz, axis=1), obs_tx90p_samz)
rmse_samz_rcm_tn10p = compute_rmse(np.nanmean(rcm_tn10p_samz, axis=1), obs_tn10p_samz)
rmse_samz_rcm_tn90p = compute_rmse(np.nanmean(rcm_tn90p_samz, axis=1), obs_tn90p_samz)

rmse_samz_gcm_txx = compute_rmse(gcm_txx_samz, obs_txx_samz)
rmse_samz_gcm_txn = compute_rmse(gcm_txn_samz, obs_txn_samz)
rmse_samz_gcm_tnx = compute_rmse(gcm_tnx_samz, obs_tnx_samz)
rmse_samz_gcm_tnn = compute_rmse(gcm_tnn_samz, obs_tnn_samz)
rmse_samz_gcm_dtr = compute_rmse(gcm_dtr_samz, obs_dtr_samz)
rmse_samz_gcm_su = compute_rmse(gcm_su_samz, obs_su_samz)
rmse_samz_gcm_tr = compute_rmse(gcm_tr_samz, obs_tr_samz)
rmse_samz_gcm_tx10p = compute_rmse(gcm_tx10p_samz, obs_tx10p_samz)
rmse_samz_gcm_tx90p = compute_rmse(gcm_tx90p_samz, obs_tx90p_samz)
rmse_samz_gcm_tn10p = compute_rmse(gcm_tn10p_samz, obs_tn10p_samz)
rmse_samz_gcm_tn90p = compute_rmse(gcm_tn90p_samz, obs_tn90p_samz)

print("rmse_samz_rcm")
print(rmse_samz_rcm_txx, rmse_samz_rcm_txn, rmse_samz_rcm_tnx, rmse_samz_rcm_tnn, rmse_samz_rcm_dtr, rmse_samz_rcm_su,
rmse_samz_rcm_tr, rmse_samz_rcm_tx10p, rmse_samz_rcm_tx90p, rmse_samz_rcm_tn10p, rmse_samz_rcm_tn90p)
print("rmse_samz_gcm")
print(rmse_samz_gcm_txx, rmse_samz_gcm_txn, rmse_samz_gcm_tnx, rmse_samz_gcm_tnn, rmse_samz_gcm_dtr, rmse_samz_gcm_su,
rmse_samz_gcm_tr, rmse_samz_gcm_tx10p, rmse_samz_gcm_tx90p, rmse_samz_gcm_tn10p, rmse_samz_gcm_tn90p)

rmse_eneb_rcm_txx = compute_rmse(np.nanmean(rcm_txx_eneb, axis=1), obs_txx_eneb)
rmse_eneb_rcm_txn = compute_rmse(np.nanmean(rcm_txn_eneb, axis=1), obs_txn_eneb)
rmse_eneb_rcm_tnx = compute_rmse(np.nanmean(rcm_tnx_eneb, axis=1), obs_tnx_eneb)
rmse_eneb_rcm_tnn = compute_rmse(np.nanmean(rcm_tnn_eneb, axis=1), obs_tnn_eneb)
rmse_eneb_rcm_dtr = compute_rmse(np.nanmean(rcm_dtr_eneb, axis=1), obs_dtr_eneb)
rmse_eneb_rcm_su = compute_rmse(np.nanmean(rcm_su_eneb, axis=1), obs_su_eneb)
rmse_eneb_rcm_tr = compute_rmse(np.nanmean(rcm_tr_eneb, axis=1), obs_tr_eneb)
rmse_eneb_rcm_tx10p = compute_rmse(np.nanmean(rcm_tx10p_eneb, axis=1), obs_tx10p_eneb)
rmse_eneb_rcm_tx90p = compute_rmse(np.nanmean(rcm_tx90p_eneb, axis=1), obs_tx90p_eneb)
rmse_eneb_rcm_tn10p = compute_rmse(np.nanmean(rcm_tn10p_eneb, axis=1), obs_tn10p_eneb)
rmse_eneb_rcm_tn90p = compute_rmse(np.nanmean(rcm_tn90p_eneb, axis=1), obs_tn90p_eneb)

rmse_eneb_gcm_txx = compute_rmse(gcm_txx_eneb, obs_txx_eneb)
rmse_eneb_gcm_txn = compute_rmse(gcm_txn_eneb, obs_txn_eneb)
rmse_eneb_gcm_tnx = compute_rmse(gcm_tnx_eneb, obs_tnx_eneb)
rmse_eneb_gcm_tnn = compute_rmse(gcm_tnn_eneb, obs_tnn_eneb)
rmse_eneb_gcm_dtr = compute_rmse(gcm_dtr_eneb, obs_dtr_eneb)
rmse_eneb_gcm_su = compute_rmse(gcm_su_eneb, obs_su_eneb)
rmse_eneb_gcm_tr = compute_rmse(gcm_tr_eneb, obs_tr_eneb)
rmse_eneb_gcm_tx10p = compute_rmse(gcm_tx10p_eneb, obs_tx10p_eneb)
rmse_eneb_gcm_tx90p = compute_rmse(gcm_tx90p_eneb, obs_tx90p_eneb)
rmse_eneb_gcm_tn10p = compute_rmse(gcm_tn10p_eneb, obs_tn10p_eneb)
rmse_eneb_gcm_tn90p = compute_rmse(gcm_tn90p_eneb, obs_tn90p_eneb)

print("rmse_eneb_rcm")
print(rmse_eneb_rcm_txx, rmse_eneb_rcm_txn, rmse_eneb_rcm_tnx, rmse_eneb_rcm_tnn, rmse_eneb_rcm_dtr, rmse_eneb_rcm_su,
rmse_eneb_rcm_tr, rmse_eneb_rcm_tx10p, rmse_eneb_rcm_tx90p, rmse_eneb_rcm_tn10p, rmse_eneb_rcm_tn90p)
print("rmse_eneb_gcm")
print(rmse_eneb_gcm_txx, rmse_eneb_gcm_txn, rmse_eneb_gcm_tnx, rmse_eneb_gcm_tnn, rmse_eneb_gcm_dtr, rmse_eneb_gcm_su,
rmse_eneb_gcm_tr, rmse_eneb_gcm_tx10p, rmse_eneb_gcm_tx90p, rmse_eneb_gcm_tn10p, rmse_eneb_gcm_tn90p)

rmse_matopiba_rcm_txx = compute_rmse(np.nanmean(rcm_txx_matopiba, axis=1), obs_txx_matopiba)
rmse_matopiba_rcm_txn = compute_rmse(np.nanmean(rcm_txn_matopiba, axis=1), obs_txn_matopiba)
rmse_matopiba_rcm_tnx = compute_rmse(np.nanmean(rcm_tnx_matopiba, axis=1), obs_tnx_matopiba)
rmse_matopiba_rcm_tnn = compute_rmse(np.nanmean(rcm_tnn_matopiba, axis=1), obs_tnn_matopiba)
rmse_matopiba_rcm_dtr = compute_rmse(np.nanmean(rcm_dtr_matopiba, axis=1), obs_dtr_matopiba)
rmse_matopiba_rcm_su = compute_rmse(np.nanmean(rcm_su_matopiba, axis=1), obs_su_matopiba)
rmse_matopiba_rcm_tr = compute_rmse(np.nanmean(rcm_tr_matopiba, axis=1), obs_tr_matopiba)
rmse_matopiba_rcm_tx10p = compute_rmse(np.nanmean(rcm_tx10p_matopiba, axis=1), obs_tx10p_matopiba)
rmse_matopiba_rcm_tx90p = compute_rmse(np.nanmean(rcm_tx90p_matopiba, axis=1), obs_tx90p_matopiba)
rmse_matopiba_rcm_tn10p = compute_rmse(np.nanmean(rcm_tn10p_matopiba, axis=1), obs_tn10p_matopiba)
rmse_matopiba_rcm_tn90p = compute_rmse(np.nanmean(rcm_tn90p_matopiba, axis=1), obs_tn90p_matopiba)

rmse_matopiba_gcm_txx = compute_rmse(gcm_txx_matopiba, obs_txx_matopiba)
rmse_matopiba_gcm_txn = compute_rmse(gcm_txn_matopiba, obs_txn_matopiba)
rmse_matopiba_gcm_tnx = compute_rmse(gcm_tnx_matopiba, obs_tnx_matopiba)
rmse_matopiba_gcm_tnn = compute_rmse(gcm_tnn_matopiba, obs_tnn_matopiba)
rmse_matopiba_gcm_dtr = compute_rmse(gcm_dtr_matopiba, obs_dtr_matopiba)
rmse_matopiba_gcm_su = compute_rmse(gcm_su_matopiba, obs_su_matopiba)
rmse_matopiba_gcm_tr = compute_rmse(gcm_tr_matopiba, obs_tr_matopiba)
rmse_matopiba_gcm_tx10p = compute_rmse(gcm_tx10p_matopiba, obs_tx10p_matopiba)
rmse_matopiba_gcm_tx90p = compute_rmse(gcm_tx90p_matopiba, obs_tx90p_matopiba)
rmse_matopiba_gcm_tn10p = compute_rmse(gcm_tn10p_matopiba, obs_tn10p_matopiba)
rmse_matopiba_gcm_tn90p = compute_rmse(gcm_tn90p_matopiba, obs_tn90p_matopiba)

print("rmse_matopiba_rcm")
print(rmse_matopiba_rcm_txx, rmse_matopiba_rcm_txn, rmse_matopiba_rcm_tnx, rmse_matopiba_rcm_tnn, rmse_matopiba_rcm_dtr, rmse_matopiba_rcm_su,
rmse_matopiba_rcm_tr, rmse_matopiba_rcm_tx10p, rmse_matopiba_rcm_tx90p, rmse_matopiba_rcm_tn10p, rmse_matopiba_rcm_tn90p)
print("rmse_matopiba_gcm")
print(rmse_matopiba_gcm_txx, rmse_matopiba_gcm_txn, rmse_matopiba_gcm_tnx, rmse_matopiba_gcm_tnn, rmse_matopiba_gcm_dtr, rmse_matopiba_gcm_su,
rmse_matopiba_gcm_tr, rmse_matopiba_gcm_tx10p, rmse_matopiba_gcm_tx90p, rmse_matopiba_gcm_tn10p, rmse_matopiba_gcm_tn90p)

