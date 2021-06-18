# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "02/15/2019"
__description__ = "This script plot Taylor Diagram seasonal from CMIP5 models end OBS basedata"

import os
import netCDF4
import numpy as np
import texttable as tt

from comp_statist_indices import compute_bias, compute_rmse, compute_corr, compute_effic_coeffic


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/rcm'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_obs = value[2:240:3,:,:]
	djf_obs = np.nanmean(np.nanmean(season_obs[3:80:4], axis=1), axis=1)
	mam_obs = np.nanmean(np.nanmean(season_obs[0:80:4], axis=1), axis=1)
	jja_obs = np.nanmean(np.nanmean(season_obs[1:80:4], axis=1), axis=1)
	son_obs = np.nanmean(np.nanmean(season_obs[2:80:4], axis=1), axis=1)
	annual_obs = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1)

	return djf_obs, mam_obs, jja_obs, son_obs, annual_obs
	

def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_rcm = value[2:240:3,:,:]
	djf_rcm = np.nanmean(np.nanmean(season_rcm[3:80:4], axis=1), axis=1)
	mam_rcm = np.nanmean(np.nanmean(season_rcm[0:80:4], axis=1), axis=1)
	jja_rcm = np.nanmean(np.nanmean(season_rcm[1:80:4], axis=1), axis=1)
	son_rcm = np.nanmean(np.nanmean(season_rcm[2:80:4], axis=1), axis=1)
	annual_rcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1)

	return djf_rcm, mam_rcm, jja_rcm, son_rcm, annual_rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_gcm = value[2:240:3,:,:]
	djf_gcm = np.nanmean(np.nanmean(season_gcm[3:80:4], axis=1), axis=1)
	mam_gcm = np.nanmean(np.nanmean(season_gcm[0:80:4], axis=1), axis=1)
	jja_gcm = np.nanmean(np.nanmean(season_gcm[1:80:4], axis=1), axis=1)
	son_gcm = np.nanmean(np.nanmean(season_gcm[2:80:4], axis=1), axis=1)
	annual_gcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1)
	
	return djf_gcm, mam_gcm, jja_gcm, son_gcm, annual_gcm


# Import regcm exp and cru databases 
pre_djf_cru_samz, pre_mam_cru_samz, pre_jja_cru_samz, pre_son_cru_samz, pre_annual_cru_samz = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')	   
pre_djf_cru_eneb, pre_mam_cru_eneb, pre_jja_cru_eneb, pre_son_cru_eneb, pre_annual_cru_eneb = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')	   
pre_djf_cru_matopiba, pre_mam_cru_matopiba, pre_jja_cru_matopiba, pre_son_cru_matopiba, pre_annual_cru_matopiba = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')	   

pre_djf_rcm_samz, pre_mam_rcm_samz, pre_jja_rcm_samz, pre_son_rcm_samz, pre_annual_rcm_samz = import_rcm('pre', 'amz_neb', 'hist', '1986-2005')	   
pre_djf_rcm_eneb, pre_mam_rcm_eneb, pre_jja_rcm_eneb, pre_son_rcm_eneb, pre_annual_rcm_eneb = import_rcm('pre', 'amz_neb', 'hist', '1986-2005')	   
pre_djf_rcm_matopiba, pre_mam_rcm_matopiba, pre_jja_rcm_matopiba, pre_son_rcm_matopiba, pre_annual_rcm_matopiba = import_rcm('pre', 'amz_neb', 'hist', '1986-2005')	   

pre_djf_gcm_samz, pre_mam_gcm_samz, pre_jja_gcm_samz, pre_son_gcm_samz, pre_annual_gcm_samz = import_gcm('pre', 'amz_neb', 'hist', '1986-2005')	   
pre_djf_gcm_eneb, pre_mam_gcm_eneb, pre_jja_gcm_eneb, pre_son_gcm_eneb, pre_annual_gcm_eneb = import_gcm('pre', 'amz_neb', 'hist', '1986-2005')	   
pre_djf_gcm_matopiba, pre_mam_gcm_matopiba, pre_jja_gcm_matopiba, pre_son_gcm_matopiba, pre_annual_gcm_matopiba = import_gcm('pre', 'amz_neb', 'hist', '1986-2005')	   

tmp_djf_cru_samz, tmp_mam_cru_samz, tmp_jja_cru_samz, tmp_son_cru_samz, tmp_annual_cru_samz = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')	   
tmp_djf_cru_eneb, tmp_mam_cru_eneb, tmp_jja_cru_eneb, tmp_son_cru_eneb, tmp_annual_cru_eneb = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')	   
tmp_djf_cru_matopiba, tmp_mam_cru_matopiba, tmp_jja_cru_matopiba, tmp_son_cru_matopiba, tmp_annual_cru_matopiba = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')	   

tas_djf_rcm_samz, tas_mam_rcm_samz, tas_jja_rcm_samz, tas_son_rcm_samz, tas_annual_rcm_samz = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')	   
tas_djf_rcm_eneb, tas_mam_rcm_eneb, tas_jja_rcm_eneb, tas_son_rcm_eneb, tas_annual_rcm_eneb = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')	   
tas_djf_rcm_matopiba, tas_mam_rcm_matopiba, tas_jja_rcm_matopiba, tas_son_rcm_matopiba, tas_annual_rcm_matopiba = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')	   

tas_djf_gcm_samz, tas_mam_gcm_samz, tas_jja_gcm_samz, tas_son_gcm_samz, tas_annual_gcm_samz = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')	   
tas_djf_gcm_eneb, tas_mam_gcm_eneb, tas_jja_gcm_eneb, tas_son_gcm_eneb, tas_annual_gcm_eneb = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')	   
tas_djf_gcm_matopiba, tas_mam_gcm_matopiba, tas_jja_gcm_matopiba, tas_son_gcm_matopiba, tas_annual_gcm_matopiba = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')	   

# Compute statiscts indices from models
# Precipitação
# CRU
mean_pre_djf_cru_samz = np.nanmean(pre_djf_cru_samz)
mean_pre_mam_cru_samz = np.nanmean(pre_mam_cru_samz)
mean_pre_jja_cru_samz = np.nanmean(pre_jja_cru_samz)
mean_pre_son_cru_samz = np.nanmean(pre_son_cru_samz)
mean_pre_annual_cru_samz = np.nanmean(pre_annual_cru_samz)

mean_pre_djf_cru_eneb = np.nanmean(pre_djf_cru_eneb)
mean_pre_mam_cru_eneb = np.nanmean(pre_mam_cru_eneb)
mean_pre_jja_cru_eneb = np.nanmean(pre_jja_cru_eneb)
mean_pre_son_cru_eneb = np.nanmean(pre_son_cru_eneb)
mean_pre_annual_cru_eneb = np.nanmean(pre_annual_cru_eneb)

mean_pre_djf_cru_matopiba = np.nanmean(pre_djf_cru_matopiba)
mean_pre_mam_cru_matopiba = np.nanmean(pre_mam_cru_matopiba)
mean_pre_jja_cru_matopiba = np.nanmean(pre_jja_cru_matopiba)
mean_pre_son_cru_matopiba = np.nanmean(pre_son_cru_matopiba)
mean_pre_annual_cru_matopiba = np.nanmean(pre_annual_cru_matopiba)

# RCM
mean_pre_djf_rcm_samz = np.nanmean(pre_djf_rcm_samz)
mean_pre_mam_rcm_samz = np.nanmean(pre_mam_rcm_samz)
mean_pre_jja_rcm_samz = np.nanmean(pre_jja_rcm_samz)
mean_pre_son_rcm_samz = np.nanmean(pre_son_rcm_samz)
mean_pre_annual_rcm_samz = np.nanmean(pre_annual_rcm_samz)

mean_pre_djf_rcm_eneb = np.nanmean(pre_djf_rcm_eneb)
mean_pre_mam_rcm_eneb = np.nanmean(pre_mam_rcm_eneb)
mean_pre_jja_rcm_eneb = np.nanmean(pre_jja_rcm_eneb)
mean_pre_son_rcm_eneb = np.nanmean(pre_son_rcm_eneb)
mean_pre_annual_rcm_eneb = np.nanmean(pre_annual_rcm_eneb)

mean_pre_djf_rcm_matopiba = np.nanmean(pre_djf_rcm_matopiba)
mean_pre_mam_rcm_matopiba = np.nanmean(pre_mam_rcm_matopiba)
mean_pre_jja_rcm_matopiba = np.nanmean(pre_jja_rcm_matopiba)
mean_pre_son_rcm_matopiba = np.nanmean(pre_son_rcm_matopiba)
mean_pre_annual_rcm_matopiba = np.nanmean(pre_annual_rcm_matopiba)

# GCM
mean_pre_djf_gcm_samz = np.nanmean(pre_djf_gcm_samz)
mean_pre_mam_gcm_samz = np.nanmean(pre_mam_gcm_samz)
mean_pre_jja_gcm_samz = np.nanmean(pre_jja_gcm_samz)
mean_pre_son_gcm_samz = np.nanmean(pre_son_gcm_samz)
mean_pre_annual_gcm_samz = np.nanmean(pre_annual_gcm_samz)

mean_pre_djf_gcm_eneb = np.nanmean(pre_djf_gcm_eneb)
mean_pre_mam_gcm_eneb = np.nanmean(pre_mam_gcm_eneb)
mean_pre_jja_gcm_eneb = np.nanmean(pre_jja_gcm_eneb)
mean_pre_son_gcm_eneb = np.nanmean(pre_son_gcm_eneb)
mean_pre_annual_gcm_eneb = np.nanmean(pre_annual_gcm_eneb)

mean_pre_djf_gcm_matopiba = np.nanmean(pre_djf_gcm_matopiba)
mean_pre_mam_gcm_matopiba = np.nanmean(pre_mam_gcm_matopiba)
mean_pre_jja_gcm_matopiba = np.nanmean(pre_jja_gcm_matopiba)
mean_pre_son_gcm_matopiba = np.nanmean(pre_son_gcm_matopiba)
mean_pre_annual_gcm_matopiba = np.nanmean(pre_annual_gcm_matopiba)

# temperatura
# CRU
mean_tas_djf_cru_samz = np.nanmean(tas_djf_cru_samz)
mean_tas_mam_cru_samz = np.nanmean(tas_mam_cru_samz)
mean_tas_jja_cru_samz = np.nanmean(tas_jja_cru_samz)
mean_tas_son_cru_samz = np.nanmean(tas_son_cru_samz)
mean_tas_annual_cru_samz = np.nanmean(tas_annual_cru_samz)

mean_tas_djf_cru_eneb = np.nanmean(tas_djf_cru_eneb)
mean_tas_mam_cru_eneb = np.nanmean(tas_mam_cru_eneb)
mean_tas_jja_cru_eneb = np.nanmean(tas_jja_cru_eneb)
mean_tas_son_cru_eneb = np.nanmean(tas_son_cru_eneb)
mean_tas_annual_cru_eneb = np.nanmean(tas_annual_cru_eneb)

mean_tas_djf_cru_matopiba = np.nanmean(tas_djf_cru_matopiba)
mean_tas_mam_cru_matopiba = np.nanmean(tas_mam_cru_matopiba)
mean_tas_jja_cru_matopiba = np.nanmean(tas_jja_cru_matopiba)
mean_tas_son_cru_matopiba = np.nanmean(tas_son_cru_matopiba)
mean_tas_annual_cru_matopiba = np.nanmean(tas_annual_cru_matopiba)

# RCM
mean_tas_djf_rcm_samz = np.nanmean(tas_djf_rcm_samz)
mean_tas_mam_rcm_samz = np.nanmean(tas_mam_rcm_samz)
mean_tas_jja_rcm_samz = np.nanmean(tas_jja_rcm_samz)
mean_tas_son_rcm_samz = np.nanmean(tas_son_rcm_samz)
mean_tas_annual_rcm_samz = np.nanmean(tas_annual_rcm_samz)

mean_tas_djf_rcm_eneb = np.nanmean(tas_djf_rcm_eneb)
mean_tas_mam_rcm_eneb = np.nanmean(tas_mam_rcm_eneb)
mean_tas_jja_rcm_eneb = np.nanmean(tas_jja_rcm_eneb)
mean_tas_son_rcm_eneb = np.nanmean(tas_son_rcm_eneb)
mean_tas_annual_rcm_eneb = np.nanmean(tas_annual_rcm_eneb)

mean_tas_djf_rcm_matopiba = np.nanmean(tas_djf_rcm_matopiba)
mean_tas_mam_rcm_matopiba = np.nanmean(tas_mam_rcm_matopiba)
mean_tas_jja_rcm_matopiba = np.nanmean(tas_jja_rcm_matopiba)
mean_tas_son_rcm_matopiba = np.nanmean(tas_son_rcm_matopiba)
mean_tas_annual_rcm_matopiba = np.nanmean(tas_annual_rcm_matopiba)

# GCM
mean_tas_djf_gcm_samz = np.nanmean(tas_djf_gcm_samz)
mean_tas_mam_gcm_samz = np.nanmean(tas_mam_gcm_samz)
mean_tas_jja_gcm_samz = np.nanmean(tas_jja_gcm_samz)
mean_tas_son_gcm_samz = np.nanmean(tas_son_gcm_samz)
mean_tas_annual_gcm_samz = np.nanmean(tas_annual_gcm_samz)

mean_tas_djf_gcm_eneb = np.nanmean(tas_djf_gcm_eneb)
mean_tas_mam_gcm_eneb = np.nanmean(tas_mam_gcm_eneb)
mean_tas_jja_gcm_eneb = np.nanmean(tas_jja_gcm_eneb)
mean_tas_son_gcm_eneb = np.nanmean(tas_son_gcm_eneb)
mean_tas_annual_gcm_eneb = np.nanmean(tas_annual_gcm_eneb)

mean_tas_djf_gcm_matopiba = np.nanmean(tas_djf_gcm_matopiba)
mean_tas_mam_gcm_matopiba = np.nanmean(tas_mam_gcm_matopiba)
mean_tas_jja_gcm_matopiba = np.nanmean(tas_jja_gcm_matopiba)
mean_tas_son_gcm_matopiba = np.nanmean(tas_son_gcm_matopiba)
mean_tas_annual_gcm_matopiba = np.nanmean(tas_annual_gcm_matopiba)

# Precipitation
# RCM
# SAMZ
bias_pre_djf_rcm_samz = compute_bias(pre_djf_rcm_samz, pre_djf_cru_samz)
bias_pre_mam_rcm_samz = compute_bias(pre_mam_rcm_samz, pre_mam_cru_samz)
bias_pre_jja_rcm_samz = compute_bias(pre_jja_rcm_samz, pre_jja_cru_samz)
bias_pre_son_rcm_samz = compute_bias(pre_son_rcm_samz, pre_son_cru_samz)
bias_pre_annual_rcm_samz = compute_bias(pre_annual_rcm_samz, pre_annual_cru_samz)

rmse_pre_djf_rcm_samz = compute_rmse(pre_djf_rcm_samz, pre_djf_cru_samz)
rmse_pre_mam_rcm_samz = compute_rmse(pre_mam_rcm_samz, pre_mam_cru_samz)
rmse_pre_jja_rcm_samz = compute_rmse(pre_jja_rcm_samz, pre_jja_cru_samz)
rmse_pre_son_rcm_samz = compute_rmse(pre_son_rcm_samz, pre_son_cru_samz)
rmse_pre_annual_rcm_samz = compute_rmse(pre_annual_rcm_samz, pre_annual_cru_samz)

r_pre_djf_rcm_samz    = round((np.corrcoef(np.array(pre_djf_cru_samz), np.array(pre_djf_rcm_samz)))[0][1], 3)
r_pre_mam_rcm_samz    = round((np.corrcoef(np.array(pre_mam_cru_samz), np.array(pre_mam_rcm_samz)))[0][1], 3)
r_pre_jja_rcm_samz    = round((np.corrcoef(np.array(pre_jja_cru_samz), np.array(pre_jja_rcm_samz)))[0][1], 3)
r_pre_son_rcm_samz    = round((np.corrcoef(np.array(pre_son_cru_samz), np.array(pre_son_rcm_samz)))[0][1], 3)
r_pre_annual_rcm_samz = round((np.corrcoef(np.array(pre_annual_cru_samz), np.array(pre_annual_rcm_samz)))[0][1], 3)

nse_pre_djf_rcm_samz = compute_compute_effic_coefficrmse(pre_djf_rcm_samz, pre_djf_cru_samz)
nse_pre_mam_rcm_samz = computecompute_effic_coeffic_rmse(pre_mam_rcm_samz, pre_mam_cru_samz)
nse_pre_jja_rcm_samz = compute_compute_effic_coefficrmse(pre_jja_rcm_samz, pre_jja_cru_samz)
nse_pre_son_rcm_samz = compute_compute_effic_coefficrmse(pre_son_rcm_samz, pre_son_cru_samz)
nse_pre_annual_rcm_samz = compute_effic_coeffic(pre_annual_rcm_samz, pre_annual_cru_samz)

# ENEB
bias_pre_djf_rcm_eneb = compute_bias(pre_djf_rcm_samz, pre_djf_cru_samz)
bias_pre_mam_rcm_eneb = compute_bias(pre_mam_rcm_samz, pre_mam_cru_samz)
bias_pre_jja_rcm_eneb = compute_bias(pre_jja_rcm_samz, pre_jja_cru_samz)
bias_pre_son_rcm_eneb = compute_bias(pre_son_rcm_samz, pre_son_cru_samz)
bias_pre_annual_rcm_eneb = compute_bias(pre_annual_rcm_samz, pre_annual_cru_samz)

rmse_pre_djf_rcm_eneb = compute_rmse(pre_djf_rcm_eneb, pre_djf_cru_eneb)
rmse_pre_mam_rcm_eneb = compute_rmse(pre_mam_rcm_eneb, pre_mam_cru_eneb)
rmse_pre_jja_rcm_eneb = compute_rmse(pre_jja_rcm_eneb, pre_jja_cru_eneb)
rmse_pre_son_rcm_eneb = compute_rmse(pre_son_rcm_eneb, pre_son_cru_eneb)
rmse_pre_annual_rcm_eneb = compute_rmse(pre_annual_rcm_eneb, pre_annual_cru_eneb)

r_pre_djf_rcm_eneb    = round((np.corrcoef(np.array(pre_djf_cru_eneb), np.array(pre_djf_rcm_eneb)))[0][1], 3)
r_pre_mam_rcm_eneb    = round((np.corrcoef(np.array(pre_mam_cru_eneb), np.array(pre_mam_rcm_eneb)))[0][1], 3)
r_pre_jja_rcm_eneb = round((np.corrcoef(np.array(pre_jja_cru_eneb), np.array(pre_jja_rcm_eneb)))[0][1], 3)
r_pre_son_rcm_eneb    = round((np.corrcoef(np.array(pre_son_cru_eneb), np.array(pre_son_rcm_eneb)))[0][1], 3)
r_pre_annual_rcm_eneb = round((np.corrcoef(np.array(pre_annual_cru_eneb), np.array(pre_annual_rcm_eneb)))[0][1], 3)

nse_pre_djf_rcm_eneb = compute_compute_effic_coefficrmse(pre_djf_rcm_eneb, pre_djf_cru_eneb)
nse_pre_mam_rcm_eneb = computecompute_effic_coeffic_rmse(pre_mam_rcm_eneb, pre_mam_cru_eneb)
nse_pre_jja_rcm_eneb = compute_compute_effic_coefficrmse(pre_jja_rcm_eneb, pre_jja_cru_eneb)
nse_pre_son_rcm_eneb = compute_compute_effic_coefficrmse(pre_son_rcm_eneb, pre_son_cru_eneb)
nse_pre_annual_rcm_eneb = compute_effic_coeffic(pre_annual_rcm_eneb, pre_annual_cru_eneb)
	
# MATOPIBA
bias_pre_djf_rcm_matopiba = compute_bias(pre_djf_rcm_matopiba, pre_djf_cru_matopibaz)
bias_pre_mam_rcm_matopiba = compute_bias(pre_mam_rcm_matopiba, pre_mam_cru_matopiba)
bias_pre_jja_rcm_matopiba = compute_bias(pre_jja_rcm_matopiba, pre_jja_cru_matopiba)
bias_pre_son_rcm_matopiba = compute_bias(pre_son_rcm_matopiba, pre_son_cru_matopiba)
bias_pre_annual_rcm_matopiba = compute_bias(pre_annual_rcm_matopiba, pre_annual_cru_matopiba

rmse_pre_djf_rcm_matopiba = compute_rmse(pre_djf_rcm_matopiba, pre_djf_cru_matopiba)
rmse_pre_mam_rcm_matopiba = compute_rmse(pre_mam_rcm_matopiba, pre_mam_cru_matopiba)
rmse_pre_jja_rcm_matopiba = compute_rmse(pre_jja_rcm_matopiba, pre_jja_cru_matopiba)
rmse_pre_son_rcm_matopiba = compute_rmse(pre_son_rcm_matopiba, pre_son_cru_matopiba)
rmse_pre_annual_rcm_matopiba = compute_rmse(pre_annual_rcm_matopiba, pre_annual_cru_matopiba)

r_pre_djf_rcm_matopiba    = round((np.corrcoef(np.array(pre_djf_cru_matopiba), np.array(pre_djf_rcm_matopiba)))[0][1], 3)
r_pre_mam_rcm_matopiba    = round((np.corrcoef(np.array(pre_mam_cru_matopiba), np.array(pre_mam_rcm_matopiba)))[0][1], 3)
r_pre_jja_rcm_matopiba    = round((np.corrcoef(np.array(pre_jja_cru_matopiba), np.array(pre_jja_rcm_matopiba)))[0][1], 3)
r_pre_son_rcm_matopiba    = round((np.corrcoef(np.array(pre_son_cru_matopiba), np.array(pre_son_rcm_matopiba)))[0][1], 3)
r_pre_annual_rcm_matopiba = round((np.corrcoef(np.array(pre_annual_cru_matopiba), np.array(pre_annual_rcm_matopiba)))[0][1], 3)

nse_pre_djf_rcm_matopiba = compute_compute_effic_coefficrmse(pre_djf_rcm_matopiba, pre_djf_cru_matopiba)
nse_pre_mam_rcm_matopiba = computecompute_effic_coeffic_rmse(pre_mam_rcm_matopiba, pre_mam_cru_matopiba)
nse_pre_jja_rcm_matopiba = compute_compute_effic_coefficrmse(pre_jja_rcm_matopiba, pre_jja_cru_matopiba)
nse_pre_son_rcm_matopiba = compute_compute_effic_coefficrmse(pre_son_rcm_matopiba, pre_son_cru_matopiba)
nse_pre_annual_rcm_matopiba = compute_effic_coeffic(pre_annual_rcm_matopiba, pre_annual_cru_matopiba)
		







