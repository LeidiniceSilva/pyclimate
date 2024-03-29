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

	dict_var = {u'eca_prcptot': u'precip', 
	u'eca_r95p': u'precip',
	u'eca_r99p': u'precip', 
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
	obs  = np.nanmean(np.nanmean(var[:][0:20,:,:], axis=1), axis=1)
	
	return obs
	
	
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
	rcm  = np.nanmean(np.nanmean(var[:][0:20,:,:], axis=1), axis=1)

	return rcm


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
	gcm  = np.nanmean(np.nanmean(var[:][0:20,:,:], axis=1), axis=1)

	return gcm


# Import extreme indices 
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
	
# Compute correlation 
r_samz_rcm_prcptot = compute_corr(obs_prcptot_samz, rcm_prcptot_samz)
r_samz_rcm_r95p = compute_corr(obs_r95p_samz, rcm_r95p_samz)
r_samz_rcm_r99p = compute_corr(obs_r99p_samz, rcm_r99p_samz)
r_samz_rcm_rx1day = compute_corr(obs_rx1day_samz, rcm_rx1day_samz)
r_samz_rcm_rx5day = compute_corr(obs_rx5day_samz, rcm_rx5day_samz)
r_samz_rcm_sdii = compute_corr(obs_sdii_samz, rcm_sdii_samz)
r_samz_rcm_cdd = compute_corr(obs_cdd_samz, rcm_cdd_samz)
r_samz_rcm_cwd = compute_corr(obs_cwd_samz, rcm_cwd_samz)
r_samz_rcm_r10mm = compute_corr(obs_r10mm_samz, rcm_r10mm_samz)
r_samz_rcm_r20mm = compute_corr(obs_r20mm_samz, rcm_r20mm_samz)

r_samz_gcm_prcptot = compute_corr(obs_prcptot_samz, gcm_prcptot_samz)
r_samz_gcm_r95p = compute_corr(obs_r95p_samz, gcm_r95p_samz)
r_samz_gcm_r99p = compute_corr(obs_r99p_samz, gcm_r99p_samz)
r_samz_gcm_rx1day = compute_corr(obs_rx1day_samz, gcm_rx1day_samz)
r_samz_gcm_rx5day = compute_corr(obs_rx5day_samz, gcm_rx5day_samz)
r_samz_gcm_sdii = compute_corr(obs_sdii_samz, gcm_sdii_samz)
r_samz_gcm_cdd = compute_corr(obs_cdd_samz, gcm_cdd_samz)
r_samz_gcm_cwd = compute_corr(obs_cwd_samz, gcm_cwd_samz)
r_samz_gcm_r10mm = compute_corr(obs_r10mm_samz, gcm_r10mm_samz)
r_samz_gcm_r20mm = compute_corr(obs_r20mm_samz, gcm_r20mm_samz)

print("r_samz_rcm")
print(r_samz_rcm_prcptot, r_samz_rcm_r95p, r_samz_rcm_r99p, r_samz_rcm_rx1day, r_samz_rcm_rx5day, r_samz_rcm_sdii,
r_samz_rcm_cdd, r_samz_rcm_cwd, r_samz_rcm_r10mm, r_samz_rcm_r20mm)
print("r_samz_gcm")
print(r_samz_gcm_prcptot, r_samz_gcm_r95p, r_samz_gcm_r99p, r_samz_gcm_rx1day, r_samz_gcm_rx5day, r_samz_gcm_sdii,
r_samz_gcm_cdd, r_samz_gcm_cwd, r_samz_gcm_r10mm, r_samz_gcm_r20mm)

r_eneb_rcm_prcptot = compute_corr(obs_prcptot_eneb, rcm_prcptot_eneb)
r_eneb_rcm_r95p = compute_corr(obs_r95p_eneb, rcm_r95p_eneb)
r_eneb_rcm_r99p = compute_corr(obs_r99p_eneb, rcm_r99p_eneb)
r_eneb_rcm_rx1day = compute_corr(obs_rx1day_eneb, rcm_rx1day_eneb)
r_eneb_rcm_rx5day = compute_corr(obs_rx5day_eneb, rcm_rx5day_eneb)
r_eneb_rcm_sdii = compute_corr(obs_sdii_eneb, rcm_sdii_eneb)
r_eneb_rcm_cdd = compute_corr(obs_cdd_eneb, rcm_cdd_eneb)
r_eneb_rcm_cwd = compute_corr(obs_cwd_eneb, rcm_cwd_eneb)
r_eneb_rcm_r10mm = compute_corr(obs_r10mm_eneb, rcm_r10mm_eneb)
r_eneb_rcm_r20mm = compute_corr(obs_r20mm_eneb, rcm_r20mm_eneb)

r_eneb_gcm_prcptot = compute_corr(obs_prcptot_eneb, gcm_prcptot_eneb)
r_eneb_gcm_r95p = compute_corr(obs_r95p_eneb, gcm_r95p_eneb)
r_eneb_gcm_r99p = compute_corr(obs_r99p_eneb, gcm_r99p_eneb)
r_eneb_gcm_rx1day = compute_corr(obs_rx1day_eneb, gcm_rx1day_eneb)
r_eneb_gcm_rx5day = compute_corr(obs_rx5day_eneb, gcm_rx5day_eneb)
r_eneb_gcm_sdii = compute_corr(obs_sdii_eneb, gcm_sdii_eneb)
r_eneb_gcm_cdd = compute_corr(obs_cdd_eneb, gcm_cdd_eneb)
r_eneb_gcm_cwd = compute_corr(obs_cwd_eneb, gcm_cwd_eneb)
r_eneb_gcm_r10mm = compute_corr(obs_r10mm_eneb, gcm_r10mm_eneb)
r_eneb_gcm_r20mm = compute_corr(obs_r20mm_eneb, gcm_r20mm_eneb)

print("r_eneb_rcm")
print(r_eneb_rcm_prcptot, r_eneb_rcm_r95p, r_eneb_rcm_r99p, r_eneb_rcm_rx1day, r_eneb_rcm_rx5day, r_eneb_rcm_sdii,
r_eneb_rcm_cdd, r_samz_rcm_cwd, r_eneb_rcm_r10mm, r_eneb_rcm_r20mm)

print("r_eneb_gcm")
print(r_eneb_gcm_prcptot, r_eneb_gcm_r95p, r_eneb_gcm_r99p, r_eneb_gcm_rx1day, r_eneb_gcm_rx5day, r_eneb_gcm_sdii,
r_eneb_gcm_cdd, r_eneb_gcm_cwd, r_eneb_gcm_r10mm, r_eneb_gcm_r20mm)

r_matopiba_rcm_prcptot = compute_corr(obs_prcptot_matopiba, rcm_prcptot_matopiba)
r_matopiba_rcm_r95p = compute_corr(obs_r95p_matopiba, rcm_r95p_matopiba)
r_matopiba_rcm_r99p = compute_corr(obs_r99p_matopiba, rcm_r99p_matopiba)
r_matopiba_rcm_rx1day = compute_corr(obs_rx1day_matopiba, rcm_rx1day_matopiba)
r_matopiba_rcm_rx5day = compute_corr(obs_rx5day_matopiba, rcm_rx5day_matopiba)
r_matopiba_rcm_sdii = compute_corr(obs_sdii_matopiba, rcm_sdii_matopiba)
r_matopiba_rcm_cdd = compute_corr(obs_cdd_matopiba, rcm_cdd_matopiba)
r_matopiba_rcm_cwd = compute_corr(obs_cwd_matopiba, rcm_cwd_matopiba)
r_matopiba_rcm_r10mm = compute_corr(obs_r10mm_matopiba, rcm_r10mm_matopiba)
r_matopiba_rcm_r20mm = compute_corr(obs_r20mm_matopiba, rcm_r20mm_matopiba)

r_matopiba_gcm_prcptot = compute_corr(obs_prcptot_matopiba, gcm_prcptot_matopiba)
r_matopiba_gcm_r95p = compute_corr(obs_r95p_matopiba, gcm_r95p_matopiba)
r_matopiba_gcm_r99p = compute_corr(obs_r99p_matopiba, gcm_r99p_matopiba)
r_matopiba_gcm_rx1day = compute_corr(obs_rx1day_matopiba, gcm_rx1day_matopiba)
r_matopiba_gcm_rx5day = compute_corr(obs_rx5day_matopiba, gcm_rx5day_matopiba)
r_matopiba_gcm_sdii = compute_corr(obs_sdii_matopiba, gcm_sdii_matopiba)
r_matopiba_gcm_cdd = compute_corr(obs_cdd_matopiba, gcm_cdd_matopiba)
r_matopiba_gcm_cwd = compute_corr(obs_cwd_matopiba, gcm_cwd_matopiba)
r_matopiba_gcm_r10mm = compute_corr(obs_r10mm_matopiba, gcm_r10mm_matopiba)
r_matopiba_gcm_r20mm = compute_corr(obs_r20mm_matopiba, gcm_r20mm_matopiba)

print("r_matopiba_rcm")
print(r_matopiba_rcm_prcptot, r_matopiba_rcm_r95p, r_matopiba_rcm_r99p, r_matopiba_rcm_rx1day, r_matopiba_rcm_rx5day, r_matopiba_rcm_sdii,
r_matopiba_rcm_cdd, r_matopiba_rcm_cwd, r_matopiba_rcm_r10mm, r_matopiba_rcm_r20mm)
print("r_matopiba_gcm")
print(r_matopiba_gcm_prcptot, r_matopiba_gcm_r95p, r_matopiba_gcm_r99p, r_matopiba_gcm_rx1day, r_matopiba_gcm_rx5day, r_matopiba_gcm_sdii,
r_matopiba_gcm_cdd, r_matopiba_gcm_cwd, r_matopiba_gcm_r10mm, r_matopiba_gcm_r20mm)
print()

# Compute root mean square error 
rmse_samz_rcm_prcptot = compute_rmse(rcm_prcptot_samz, obs_prcptot_samz)
rmse_samz_rcm_r95p = compute_rmse(rcm_r95p_samz, obs_r95p_samz)
rmse_samz_rcm_r99p = compute_rmse(rcm_r99p_samz, obs_r99p_samz)
rmse_samz_rcm_rx1day = compute_rmse(rcm_rx1day_samz, obs_rx1day_samz)
rmse_samz_rcm_rx5day = compute_rmse(rcm_rx5day_samz, obs_rx5day_samz)
rmse_samz_rcm_sdii =compute_rmse(rcm_sdii_samz, obs_sdii_samz)
rmse_samz_rcm_cdd = compute_rmse(rcm_cdd_samz, obs_cdd_samz)
rmse_samz_rcm_cwd = compute_rmse(rcm_cwd_samz, obs_cwd_samz)
rmse_samz_rcm_r10mm = compute_rmse(rcm_r10mm_samz, obs_r10mm_samz)
rmse_samz_rcm_r20mm = compute_rmse(rcm_r10mm_samz, obs_r20mm_samz)

rmse_samz_gcm_prcptot = compute_rmse(gcm_prcptot_samz, obs_prcptot_samz)
rmse_samz_gcm_r95p = compute_rmse(gcm_r95p_samz, obs_r95p_samz)
rmse_samz_gcm_r99p = compute_rmse(gcm_r99p_samz, obs_r99p_samz)
rmse_samz_gcm_rx1day = compute_rmse(gcm_rx1day_samz, obs_rx1day_samz)
rmse_samz_gcm_rx5day = compute_rmse(gcm_rx5day_samz, obs_rx5day_samz)
rmse_samz_gcm_sdii = compute_rmse(gcm_sdii_samz, obs_sdii_samz)
rmse_samz_gcm_cdd = compute_rmse(gcm_cdd_samz, obs_cdd_samz)
rmse_samz_gcm_cwd = compute_rmse(gcm_cwd_samz, obs_cwd_samz)
rmse_samz_gcm_r10mm = compute_rmse(gcm_r10mm_samz, obs_r10mm_samz)
rmse_samz_gcm_r20mm = compute_rmse(gcm_r20mm_samz, obs_r20mm_samz)

print("rmse_samz_rcm")
print(rmse_samz_rcm_prcptot, rmse_samz_rcm_r95p, rmse_samz_rcm_r99p, rmse_samz_rcm_rx1day, rmse_samz_rcm_rx5day, rmse_samz_rcm_sdii,
rmse_samz_rcm_cdd, rmse_samz_rcm_cwd, rmse_samz_rcm_r10mm, rmse_samz_rcm_r20mm)
print("rmse_samz_gcm")
print(rmse_samz_gcm_prcptot, rmse_samz_gcm_r95p, rmse_samz_gcm_r99p, rmse_samz_gcm_rx1day, rmse_samz_gcm_rx5day, rmse_samz_gcm_sdii,
rmse_samz_gcm_cdd, rmse_samz_gcm_cwd, rmse_samz_gcm_r10mm, rmse_samz_gcm_r20mm)

rmse_eneb_rcm_prcptot = compute_rmse(rcm_prcptot_eneb, obs_prcptot_eneb)
rmse_eneb_rcm_r95p = compute_rmse(rcm_r95p_eneb, obs_r95p_eneb)
rmse_eneb_rcm_r99p = compute_rmse(rcm_r99p_eneb, obs_r99p_eneb)
rmse_eneb_rcm_rx1day = compute_rmse(rcm_rx1day_eneb, obs_rx1day_eneb)
rmse_eneb_rcm_rx5day = compute_rmse(rcm_rx5day_eneb, obs_rx5day_eneb)
rmse_eneb_rcm_sdii =compute_rmse(rcm_sdii_eneb, obs_sdii_eneb)
rmse_eneb_rcm_cdd = compute_rmse(rcm_cdd_eneb, obs_cdd_eneb)
rmse_eneb_rcm_cwd = compute_rmse(rcm_cwd_eneb, obs_cwd_eneb)
rmse_eneb_rcm_r10mm = compute_rmse(rcm_r10mm_eneb, obs_r10mm_eneb)
rmse_eneb_rcm_r20mm = compute_rmse(rcm_r10mm_eneb, obs_r20mm_eneb)

rmse_eneb_gcm_prcptot = compute_rmse(gcm_prcptot_samz, obs_prcptot_samz)
rmse_eneb_gcm_r95p = compute_rmse(gcm_r95p_samz, obs_r95p_samz)
rmse_eneb_gcm_r99p = compute_rmse(gcm_r99p_samz, obs_r99p_samz)
rmse_eneb_gcm_rx1day = compute_rmse(gcm_rx1day_samz, obs_rx1day_samz)
rmse_eneb_gcm_rx5day = compute_rmse(gcm_rx5day_samz, obs_rx5day_samz)
rmse_eneb_gcm_sdii = compute_rmse(gcm_sdii_samz, obs_sdii_samz)
rmse_eneb_gcm_cdd = compute_rmse(gcm_cdd_samz, obs_cdd_samz)
rmse_eneb_gcm_cwd = compute_rmse(gcm_cwd_samz, obs_cwd_samz)
rmse_eneb_gcm_r10mm = compute_rmse(gcm_r10mm_samz, obs_r10mm_samz)
rmse_eneb_gcm_r20mm = compute_rmse(gcm_r20mm_samz, obs_r20mm_samz)

print("rmse_eneb_rcm")
print(rmse_eneb_rcm_prcptot, rmse_eneb_rcm_r95p, rmse_eneb_rcm_r99p, rmse_eneb_rcm_rx1day, rmse_eneb_rcm_rx5day, rmse_eneb_rcm_sdii,
rmse_eneb_rcm_cdd, rmse_eneb_rcm_cwd, rmse_eneb_rcm_r10mm, rmse_eneb_rcm_r20mm)
print("rmse_eneb_gcm")
print(rmse_eneb_gcm_prcptot, rmse_eneb_gcm_r95p, rmse_eneb_gcm_r99p, rmse_eneb_gcm_rx1day, rmse_eneb_gcm_rx5day, rmse_eneb_gcm_sdii,
rmse_eneb_gcm_cdd, rmse_eneb_gcm_cwd, rmse_eneb_gcm_r10mm, rmse_eneb_gcm_r20mm)

rmse_matopiba_rcm_prcptot = compute_rmse(rcm_prcptot_matopiba, obs_prcptot_matopiba)
rmse_matopiba_rcm_r95p = compute_rmse(rcm_r95p_matopiba, obs_r95p_matopiba)
rmse_matopiba_rcm_r99p = compute_rmse(rcm_r99p_matopiba, obs_r99p_matopiba)
rmse_matopiba_rcm_rx1day = compute_rmse(rcm_rx1day_matopiba, obs_rx1day_matopiba)
rmse_matopiba_rcm_rx5day = compute_rmse(rcm_rx5day_matopiba, obs_rx5day_matopiba)
rmse_matopiba_rcm_sdii =compute_rmse(rcm_sdii_matopiba, obs_sdii_matopiba)
rmse_matopiba_rcm_cdd = compute_rmse(rcm_cdd_matopiba, obs_cdd_matopiba)
rmse_matopiba_rcm_cwd = compute_rmse(rcm_cwd_matopiba, obs_cwd_matopiba)
rmse_matopiba_rcm_r10mm = compute_rmse(rcm_r10mm_matopiba, obs_r10mm_matopiba)
rmse_matopiba_rcm_r20mm = compute_rmse(rcm_r10mm_matopiba, obs_r20mm_matopiba)

rmse_matopiba_gcm_prcptot = compute_rmse(gcm_prcptot_matopiba, obs_prcptot_matopiba)
rmse_matopiba_gcm_r95p = compute_rmse(gcm_r95p_matopiba, obs_r95p_matopiba)
rmse_matopiba_gcm_r99p = compute_rmse(gcm_r99p_matopiba, obs_r99p_matopiba)
rmse_matopiba_gcm_rx1day = compute_rmse(gcm_rx1day_matopiba, obs_rx1day_matopiba)
rmse_matopiba_gcm_rx5day = compute_rmse(gcm_rx5day_matopiba, obs_rx5day_matopiba)
rmse_matopiba_gcm_sdii = compute_rmse(gcm_sdii_matopiba, obs_sdii_matopiba)
rmse_matopiba_gcm_cdd = compute_rmse(gcm_cdd_matopiba, obs_cdd_matopiba)
rmse_matopiba_gcm_cwd = compute_rmse(gcm_cwd_matopiba, obs_cwd_matopiba)
rmse_matopiba_gcm_r10mm = compute_rmse(gcm_r10mm_matopiba, obs_r10mm_matopiba)
rmse_matopiba_gcm_r20mm = compute_rmse(gcm_r20mm_matopiba, obs_r20mm_matopiba)

print("rmse_matopiba_rcm")
print(rmse_matopiba_rcm_prcptot, rmse_matopiba_rcm_r95p, rmse_matopiba_rcm_r99p, rmse_matopiba_rcm_rx1day, rmse_matopiba_rcm_rx5day, rmse_matopiba_rcm_sdii,
rmse_samz_rcm_cdd, rmse_samz_rcm_cwd, rmse_samz_rcm_r10mm, rmse_samz_rcm_r20mm)
print("rmse_matopiba_gcm")
print(rmse_matopiba_gcm_prcptot, rmse_matopiba_gcm_r95p, rmse_matopiba_gcm_r99p, rmse_matopiba_gcm_rx1day, rmse_matopiba_gcm_rx5day, rmse_matopiba_gcm_sdii,
rmse_matopiba_gcm_cdd, rmse_matopiba_gcm_cwd, rmse_matopiba_gcm_r10mm, rmse_matopiba_gcm_r20mm)


