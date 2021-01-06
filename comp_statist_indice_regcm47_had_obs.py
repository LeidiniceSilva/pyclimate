# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot portrait from Reg and Had models end obs database"

import os
import netCDF4
import matplotlib
import numpy as np
import seaborn as sb 
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from comp_statist_indices import compute_corr, compute_mae, compute_rmse, compute_index_agreement


def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_rcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1) 
	sea_rcm = np.nanmean(np.nanmean(value[0:240:3,:,:], axis=1), axis=1)
	djf_rcm = sea_rcm[0:80:4]
	mam_rcm = sea_rcm[1:80:4]
	jja_rcm = sea_rcm[2:80:4]
	son_rcm = sea_rcm[3:80:4]

	return annual_rcm, djf_rcm, mam_rcm, jja_rcm, son_rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_gcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1) 
	sea_gcm = np.nanmean(np.nanmean(value[0:240:3,:,:], axis=1), axis=1)
	djf_gcm = sea_gcm[0:80:4]
	mam_gcm = sea_gcm[1:80:4]
	jja_gcm = sea_gcm[2:80:4]
	son_gcm = sea_gcm[3:80:4]

	return annual_gcm, djf_gcm, mam_gcm, jja_gcm, son_gcm

	
def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_obs = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1) 
	sea_obs = np.nanmean(np.nanmean(value[0:240:3,:,:], axis=1), axis=1)
	djf_obs = sea_obs[0:80:4]
	mam_obs = sea_obs[1:80:4]
	jja_obs = sea_obs[2:80:4]
	son_obs = sea_obs[3:80:4]

	return annual_obs, djf_obs, mam_obs, jja_obs, son_obs
	
               
# Import regcm exps model end obs database climatology
annual_rcm_pre_samz, djf_rcm_pre_samz, mam_rcm_pre_samz, jja_rcm_pre_samz, son_rcm_pre_samz = import_rcm('pr', 'samz', 'hist', '1986-2005')
annual_gcm_pre_samz, djf_gcm_pre_samz, mam_gcm_pre_samz, jja_gcm_pre_samz, son_gcm_pre_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')
annual_obs_pre_samz, djf_obs_pre_samz, mam_obs_pre_samz, jja_obs_pre_samz, son_obs_pre_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')

annual_rcm_pre_eneb, djf_rcm_pre_eneb, mam_rcm_pre_eneb, jja_rcm_pre_eneb, son_rcm_pre_eneb = import_rcm('pr', 'eneb', 'hist', '1986-2005')
annual_gcm_pre_eneb, djf_gcm_pre_eneb, mam_gcm_pre_eneb, jja_gcm_pre_eneb, son_gcm_pre_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')
annual_obs_pre_eneb, djf_obs_pre_eneb, mam_obs_pre_eneb, jja_obs_pre_eneb, son_obs_pre_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')

annual_rcm_pre_matopiba, djf_rcm_pre_matopiba, mam_rcm_pre_matopiba, jja_rcm_pre_matopiba, son_rcm_pre_matopiba = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
annual_gcm_pre_matopiba, djf_gcm_pre_matopiba, mam_gcm_pre_matopiba, jja_gcm_pre_matopiba, son_gcm_pre_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')
annual_obs_pre_matopiba, djf_obs_pre_matopiba, mam_obs_pre_matopiba, jja_obs_pre_matopiba, son_obs_pre_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')

annual_rcm_tas_samz, djf_rcm_tas_samz, mam_rcm_tas_samz, jja_rcm_tas_samz, son_rcm_tas_samz = import_rcm('tas', 'samz', 'hist', '1986-2005')
annual_gcm_tas_samz, djf_gcm_tas_samz, mam_gcm_tas_samz, jja_gcm_tas_samz, son_gcm_tas_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')
annual_obs_tas_samz, djf_obs_tas_samz, mam_obs_tas_samz, jja_obs_tas_samz, son_obs_tas_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')

annual_rcm_tas_eneb, djf_rcm_tas_eneb, mam_rcm_tas_eneb, jja_rcm_tas_eneb, son_rcm_tas_eneb = import_rcm('tas', 'eneb', 'hist', '1986-2005')
annual_gcm_tas_eneb, djf_gcm_tas_eneb, mam_gcm_tas_eneb, jja_gcm_tas_eneb, son_gcm_tas_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')
annual_obs_tas_eneb, djf_obs_tas_eneb, mam_obs_tas_eneb, jja_obs_tas_eneb, son_obs_tas_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')

annual_rcm_tas_matopiba, djf_rcm_tas_matopiba, mam_rcm_tas_matopiba, jja_rcm_tas_matopiba, son_rcm_tas_matopiba = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
annual_gcm_tas_matopiba, djf_gcm_tas_matopiba, mam_gcm_tas_matopiba, jja_gcm_tas_matopiba, son_gcm_tas_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')
annual_obs_tas_matopiba, djf_obs_tas_matopiba, mam_obs_tas_matopiba, jja_obs_tas_matopiba, son_obs_tas_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')

# Compute statistical indices 
# regcm - pre (samz, eneb, matopiba)
corr_rcm_pre_samz = compute_corr(annual_rcm_pre_samz, annual_obs_pre_samz)
mae_rcm_pre_samz = compute_mae(annual_rcm_pre_samz, annual_obs_pre_samz)
rmse_rcm_pre_samz = compute_rmse(annual_rcm_pre_samz, annual_obs_pre_samz)
icw_rcm_pre_samz = compute_index_agreement(annual_rcm_pre_samz, annual_obs_pre_samz)

corr_rcm_pre_eneb = compute_corr(annual_rcm_pre_eneb, annual_obs_pre_eneb)
mae_rcm_pre_eneb = compute_mae(annual_rcm_pre_eneb, annual_obs_pre_eneb)
rmse_rcm_pre_eneb = compute_rmse(annual_rcm_pre_eneb, annual_obs_pre_eneb)
icw_rcm_pre_eneb = compute_index_agreement(annual_rcm_pre_eneb, annual_obs_pre_eneb)

corr_rcm_pre_matopiba = compute_corr(annual_rcm_pre_matopiba, annual_obs_pre_matopiba)
mae_rcm_pre_matopiba = compute_mae(annual_rcm_pre_matopiba, annual_obs_pre_matopiba)
rmse_rcm_pre_matopiba = compute_rmse(annual_rcm_pre_matopiba, annual_obs_pre_matopiba)
icw_rcm_pre_matopiba = compute_index_agreement(annual_rcm_pre_matopiba, annual_obs_pre_matopiba)

# hadgem - pre (samz, eneb, matopiba)
corr_gcm_pre_samz = compute_corr(annual_gcm_pre_samz, annual_obs_pre_samz)
mae_gcm_pre_samz = compute_mae(annual_gcm_pre_samz, annual_obs_pre_samz)
rmse_gcm_pre_samz = compute_rmse(annual_gcm_pre_samz, annual_obs_pre_samz)
icw_gcm_pre_samz = compute_index_agreement(annual_gcm_pre_samz, annual_obs_pre_samz)

corr_gcm_pre_eneb = compute_corr(annual_gcm_pre_eneb, annual_obs_pre_eneb)
mae_gcm_pre_eneb = compute_mae(annual_gcm_pre_eneb, annual_obs_pre_eneb)
rmse_gcm_pre_eneb = compute_rmse(annual_gcm_pre_eneb, annual_obs_pre_eneb)
icw_gcm_pre_eneb = compute_index_agreement(annual_gcm_pre_eneb, annual_obs_pre_eneb)

corr_gcm_pre_matopiba = compute_corr(annual_gcm_pre_matopiba, annual_obs_pre_matopiba)
mae_gcm_pre_matopiba = compute_mae(annual_gcm_pre_matopiba, annual_obs_pre_matopiba)
rmse_gcm_pre_matopiba = compute_rmse(annual_gcm_pre_matopiba, annual_obs_pre_matopiba)
icw_gcm_pre_matopiba = compute_index_agreement(annual_gcm_pre_matopiba, annual_obs_pre_matopiba)

# regcm - tas (samz, eneb, matopiba)

#~ corr_rcm_tas_samz = compute_corr(annual_rcm_tas_samz, annual_obs_tas_samz)
#~ mae_rcm_tas_samz = compute_mae(annual_rcm_tas_samz, annual_obs_tas_samz)
#~ rmse_rcm_tas_samz = compute_rmse(annual_rcm_tas_samz, annual_obs_tas_samz)
#~ icw_rcm_tas_samz = compute_index_agreement(annual_rcm_tas_samz, annual_obs_tas_samz)

#~ corr_rcm_tas_eneb = compute_corr(annual_rcm_tas_eneb, annual_obs_tas_eneb)
#~ mae_rcm_tas_eneb = compute_mae(annual_rcm_tas_eneb, annual_obs_tas_eneb)
#~ rmse_rcm_tas_eneb = compute_rmse(annual_rcm_tas_eneb, annual_obs_tas_eneb)
#~ icw_rcm_tas_eneb = compute_index_agreement(annual_rcm_tas_eneb, annual_obs_tas_eneb)

#~ corr_rcm_tas_matopiba = compute_corr(annual_rcm_tas_matopiba, annual_obs_tas_matopiba)
#~ mae_rcm_tas_matopiba = compute_mae(annual_rcm_tas_matopiba, annual_obs_tas_matopiba)
#~ rmse_rcm_tas_matopiba = compute_rmse(annual_rcm_tas_matopiba, annual_obs_tas_matopiba)
#~ icw_rcm_tas_matopiba = compute_index_agreement(annual_rcm_tas_matopiba, annual_obs_tas_matopiba)

# hadgem - tas (samz, eneb, matopiba)
corr_gcm_tas_samz = compute_corr(annual_gcm_tas_samz, annual_obs_tas_samz)
mae_gcm_tas_samz = compute_mae(annual_gcm_tas_samz, annual_obs_tas_samz)
rmse_gcm_tas_samz = compute_rmse(annual_gcm_tas_samz, annual_obs_tas_samz)
icw_gcm_tas_samz = compute_index_agreement(annual_gcm_tas_samz, annual_obs_tas_samz)

corr_gcm_tas_eneb = compute_corr(annual_gcm_tas_eneb, annual_obs_tas_eneb)
mae_gcm_tas_eneb = compute_mae(annual_gcm_tas_eneb, annual_obs_tas_eneb)
rmse_gcm_tas_eneb = compute_rmse(annual_gcm_tas_eneb, annual_obs_tas_eneb)
icw_gcm_tas_eneb = compute_index_agreement(annual_gcm_tas_eneb, annual_obs_tas_eneb)

corr_gcm_tas_matopiba = compute_corr(annual_gcm_tas_matopiba, annual_obs_tas_matopiba)
mae_gcm_tas_matopiba = compute_mae(annual_gcm_tas_matopiba, annual_obs_tas_matopiba)
rmse_gcm_tas_matopiba = compute_rmse(annual_gcm_tas_matopiba, annual_obs_tas_matopiba)
icw_gcm_tas_matopiba = compute_index_agreement(annual_gcm_tas_matopiba, annual_obs_tas_matopiba)

# Print results from compute statistical indices
print('regcm --> pre --> samz')
print(corr_rcm_pre_samz)
print(mae_rcm_pre_samz)
print(rmse_rcm_pre_samz)
print(icw_rcm_pre_samz)

print('regcm --> pre --> eneb')
print(corr_rcm_pre_eneb)
print(mae_rcm_pre_eneb)
print(rmse_rcm_pre_eneb)
print(icw_rcm_pre_eneb)

print('regcm --> pre --> matopiba')
print(corr_rcm_pre_matopiba)
print(mae_rcm_pre_matopiba)
print(rmse_rcm_pre_matopiba)
print(icw_rcm_pre_matopiba)

print('hadgem --> pre --> samz')
print(corr_gcm_pre_samz)
print(mae_gcm_pre_samz)
print(rmse_gcm_pre_samz)
print(icw_gcm_pre_samz)

print('hadgem --> pre --> eneb')
print(corr_gcm_pre_eneb)
print(mae_gcm_pre_eneb)
print(rmse_gcm_pre_eneb)
print(icw_gcm_pre_eneb)

print('hadgem --> pre --> matopiba')
print(corr_gcm_pre_matopiba)
print(mae_gcm_pre_matopiba)
print(rmse_gcm_pre_matopiba)
print(icw_gcm_pre_matopiba)

#~ print('regcm --> tas --> samz')
#~ print(corr_rcm_tas_samz)
#~ print(mae_rcm_tas_samz)
#~ print(rmse_rcm_tas_samz)
#~ print(icw_rcm_tas_samz)

#~ print('regcm --> tas --> eneb')
#~ print(corr_rcm_tas_eneb)
#~ print(mae_rcm_tas_eneb)
#~ print(rmse_rcm_tas_eneb)
#~ print(icw_rcm_tas_eneb)

#~ print('regcm --> tas --> matopiba')
#~ print(corr_rcm_tas_matopiba)
#~ print(mae_rcm_tas_matopiba)
#~ print(rmse_rcm_tas_matopiba)
#~ print(icw_rcm_tas_matopiba)

print('hadgem --> tas --> samz')
print(corr_gcm_tas_samz)
print(mae_gcm_tas_samz)
print(rmse_gcm_tas_samz)
print(icw_gcm_tas_samz)

print('hadgem --> tas --> eneb')
print(corr_gcm_tas_eneb)
print(mae_gcm_tas_eneb)
print(rmse_gcm_tas_eneb)
print(icw_gcm_tas_eneb)

print('hadgem --> tas --> matopiba')
print(corr_gcm_tas_matopiba)
print(mae_gcm_tas_matopiba)
print(rmse_gcm_tas_matopiba)
print(icw_gcm_tas_matopiba)

exit()
