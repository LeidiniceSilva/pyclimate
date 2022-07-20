# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot portrait diagram from regcm47 and hadgem models and obs database"

import os
import netCDF4
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from matplotlib.font_manager import FontProperties
from comp_statist_indices import compute_ioa


def import_obs(var, area, dataset, period, dt):
	
	path = '/home/nice/Documents/dataset/obs/reg_exp1'
	arq  = '{0}/{1}_{2}_{3}_obs_{4}_{5}_lonlat.nc'.format(path, var, area, dataset, period, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	mean_obs = np.nanmean(np.nanmean(value, axis=1), axis=1) 
	
	return mean_obs
	
	
def import_rcm(var, area, exp, period, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp1/hist'
	arq  = '{0}/{1}_{2}_reg_had_{3}_{4}_{5}_lonlat_seamask.nc'.format(path, var, area, exp, period, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	mean_rcm = np.nanmean(np.nanmean(value, axis=1), axis=1) 

	return mean_rcm
	

def import_gcm(var, area, exp, period, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp1/hist'
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_{4}_{5}_lonlat_seamask.nc'.format(path, var, area, exp, period, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	mean_gcm = np.nanmean(np.nanmean(value, axis=1), axis=1) 

	return mean_gcm

	              
# Import models and obs database 
# Precipitation
djf_obs_pre_samz = import_obs('pre', 'samz', 'cru_ts4.04', 'djf', '1986-2005')	   
mam_obs_pre_samz = import_obs('pre', 'samz', 'cru_ts4.04', 'mam', '1986-2005')	   
jja_obs_pre_samz = import_obs('pre', 'samz', 'cru_ts4.04', 'jja', '1986-2005')	   
son_obs_pre_samz = import_obs('pre', 'samz', 'cru_ts4.04', 'son', '1986-2005')	   
ann_obs_pre_samz = import_obs('pre', 'samz', 'cru_ts4.04', 'ann', '1986-2005')	   

djf_obs_pre_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', 'djf', '1986-2005')	   
mam_obs_pre_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', 'mam', '1986-2005')	   
jja_obs_pre_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', 'jja', '1986-2005')	   
son_obs_pre_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', 'son', '1986-2005')	   
ann_obs_pre_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', 'ann', '1986-2005')	   

djf_obs_pre_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', 'djf', '1986-2005')	   
mam_obs_pre_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', 'mam', '1986-2005')	   
jja_obs_pre_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', 'jja', '1986-2005')	   
son_obs_pre_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', 'son', '1986-2005')	   
ann_obs_pre_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', 'ann', '1986-2005')	   

djf_rcm_pre_samz = import_rcm('pr', 'samz', 'hist', 'djf', '1986-2005')	   
mam_rcm_pre_samz = import_rcm('pr', 'samz', 'hist', 'mam', '1986-2005')	   
jja_rcm_pre_samz = import_rcm('pr', 'samz', 'hist', 'jja', '1986-2005')	   
son_rcm_pre_samz = import_rcm('pr', 'samz', 'hist', 'son', '1986-2005')	   
ann_rcm_pre_samz = import_rcm('pr', 'samz', 'hist', 'ann', '1986-2005')	  

djf_rcm_pre_eneb = import_rcm('pr', 'eneb', 'hist', 'djf', '1986-2005')	   
mam_rcm_pre_eneb = import_rcm('pr', 'eneb', 'hist', 'mam', '1986-2005')	   
jja_rcm_pre_eneb = import_rcm('pr', 'eneb', 'hist', 'jja', '1986-2005')	   
son_rcm_pre_eneb = import_rcm('pr', 'eneb', 'hist', 'son', '1986-2005')	   
ann_rcm_pre_eneb = import_rcm('pr', 'eneb', 'hist', 'ann', '1986-2005')	  

djf_rcm_pre_matopiba = import_rcm('pr', 'matopiba', 'hist', 'djf', '1986-2005')	   
mam_rcm_pre_matopiba = import_rcm('pr', 'matopiba', 'hist', 'mam', '1986-2005')	   
jja_rcm_pre_matopiba = import_rcm('pr', 'matopiba', 'hist', 'jja', '1986-2005')	   
son_rcm_pre_matopiba = import_rcm('pr', 'matopiba', 'hist', 'son', '1986-2005')	   
ann_rcm_pre_matopiba = import_rcm('pr', 'matopiba', 'hist', 'ann', '1986-2005')	  

djf_gcm_pre_samz = import_gcm('pr', 'samz', 'hist', 'djf', '1986-2005')	   
mam_gcm_pre_samz = import_gcm('pr', 'samz', 'hist', 'mam', '1986-2005')	   
jja_gcm_pre_samz = import_gcm('pr', 'samz', 'hist', 'jja', '1986-2005')	   
son_gcm_pre_samz = import_gcm('pr', 'samz', 'hist', 'son', '1986-2005')	   
ann_gcm_pre_samz = import_gcm('pr', 'samz', 'hist', 'ann', '1986-2005')	  

djf_gcm_pre_eneb = import_gcm('pr', 'eneb', 'hist', 'djf', '1986-2005')	   
mam_gcm_pre_eneb = import_gcm('pr', 'eneb', 'hist', 'mam', '1986-2005')	   
jja_gcm_pre_eneb = import_gcm('pr', 'eneb', 'hist', 'jja', '1986-2005')	   
son_gcm_pre_eneb = import_gcm('pr', 'eneb', 'hist', 'son', '1986-2005')	   
ann_gcm_pre_eneb = import_gcm('pr', 'eneb', 'hist', 'ann', '1986-2005')	  

djf_gcm_pre_matopiba = import_gcm('pr', 'matopiba', 'hist', 'djf', '1986-2005')	   
mam_gcm_pre_matopiba = import_gcm('pr', 'matopiba', 'hist', 'mam', '1986-2005')	   
jja_gcm_pre_matopiba = import_gcm('pr', 'matopiba', 'hist', 'jja', '1986-2005')	   
son_gcm_pre_matopiba = import_gcm('pr', 'matopiba', 'hist', 'son', '1986-2005')	   
ann_gcm_pre_matopiba = import_gcm('pr', 'matopiba', 'hist', 'ann', '1986-2005')	  

# Temperature
djf_obs_tas_samz = import_obs('tmp', 'samz', 'cru_ts4.04', 'djf', '1986-2005')	   
mam_obs_tas_samz = import_obs('tmp', 'samz', 'cru_ts4.04', 'mam', '1986-2005')	   
jja_obs_tas_samz = import_obs('tmp', 'samz', 'cru_ts4.04', 'jja', '1986-2005')	   
son_obs_tas_samz = import_obs('tmp', 'samz', 'cru_ts4.04', 'son', '1986-2005')	   
ann_obs_tas_samz = import_obs('tmp', 'samz', 'cru_ts4.04', 'ann', '1986-2005')	   

djf_obs_tas_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', 'djf', '1986-2005')	   
mam_obs_tas_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', 'mam', '1986-2005')	   
jja_obs_tas_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', 'jja', '1986-2005')	   
son_obs_tas_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', 'son', '1986-2005')	   
ann_obs_tas_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', 'ann', '1986-2005')	   

djf_obs_tas_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', 'djf', '1986-2005')	   
mam_obs_tas_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', 'mam', '1986-2005')	   
jja_obs_tas_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', 'jja', '1986-2005')	   
son_obs_tas_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', 'son', '1986-2005')	   
ann_obs_tas_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', 'ann', '1986-2005')	   

djf_rcm_tas_samz = import_rcm('tas', 'samz', 'hist', 'djf', '1986-2005')	   
mam_rcm_tas_samz = import_rcm('tas', 'samz', 'hist', 'mam', '1986-2005')	   
jja_rcm_tas_samz = import_rcm('tas', 'samz', 'hist', 'jja', '1986-2005')	   
son_rcm_tas_samz = import_rcm('tas', 'samz', 'hist', 'son', '1986-2005')	   
ann_rcm_tas_samz = import_rcm('tas', 'samz', 'hist', 'ann', '1986-2005')	  

djf_rcm_tas_eneb = import_rcm('tas', 'eneb', 'hist', 'djf', '1986-2005')	   
mam_rcm_tas_eneb = import_rcm('tas', 'eneb', 'hist', 'mam', '1986-2005')	   
jja_rcm_tas_eneb = import_rcm('tas', 'eneb', 'hist', 'jja', '1986-2005')	   
son_rcm_tas_eneb = import_rcm('tas', 'eneb', 'hist', 'son', '1986-2005')	   
ann_rcm_tas_eneb = import_rcm('tas', 'eneb', 'hist', 'ann', '1986-2005')	  

djf_rcm_tas_matopiba = import_rcm('tas', 'matopiba', 'hist', 'djf', '1986-2005')	   
mam_rcm_tas_matopiba = import_rcm('tas', 'matopiba', 'hist', 'mam', '1986-2005')	   
jja_rcm_tas_matopiba = import_rcm('tas', 'matopiba', 'hist', 'jja', '1986-2005')	   
son_rcm_tas_matopiba = import_rcm('tas', 'matopiba', 'hist', 'son', '1986-2005')	   
ann_rcm_tas_matopiba = import_rcm('tas', 'matopiba', 'hist', 'ann', '1986-2005')	  

djf_gcm_tas_samz = import_gcm('tas', 'samz', 'hist', 'djf', '1986-2005')	   
mam_gcm_tas_samz = import_gcm('tas', 'samz', 'hist', 'mam', '1986-2005')	   
jja_gcm_tas_samz = import_gcm('tas', 'samz', 'hist', 'jja', '1986-2005')	   
son_gcm_tas_samz = import_gcm('tas', 'samz', 'hist', 'son', '1986-2005')	   
ann_gcm_tas_samz = import_gcm('tas', 'samz', 'hist', 'ann', '1986-2005')	  

djf_gcm_tas_eneb = import_gcm('tas', 'eneb', 'hist', 'djf', '1986-2005')	   
mam_gcm_tas_eneb = import_gcm('tas', 'eneb', 'hist', 'mam', '1986-2005')	   
jja_gcm_tas_eneb = import_gcm('tas', 'eneb', 'hist', 'jja', '1986-2005')	   
son_gcm_tas_eneb = import_gcm('tas', 'eneb', 'hist', 'son', '1986-2005')	   
ann_gcm_tas_eneb = import_gcm('tas', 'eneb', 'hist', 'ann', '1986-2005')	  

djf_gcm_tas_matopiba = import_gcm('tas', 'matopiba', 'hist', 'djf', '1986-2005')	   
mam_gcm_tas_matopiba = import_gcm('tas', 'matopiba', 'hist', 'mam', '1986-2005')	   
jja_gcm_tas_matopiba = import_gcm('tas', 'matopiba', 'hist', 'jja', '1986-2005')	   
son_gcm_tas_matopiba = import_gcm('tas', 'matopiba', 'hist', 'son', '1986-2005')	   
ann_gcm_tas_matopiba = import_gcm('tas', 'matopiba', 'hist', 'ann', '1986-2005')	  

# Compute index of agreement
# Precipitation
# RegCM
nse_djf_rcm_pre_samz = compute_ioa(djf_rcm_pre_samz, djf_obs_pre_samz)
nse_mam_rcm_pre_samz = compute_ioa(mam_rcm_pre_samz, mam_obs_pre_samz)
nse_jja_rcm_pre_samz = compute_ioa(jja_rcm_pre_samz, jja_obs_pre_samz)
nse_son_rcm_pre_samz = compute_ioa(son_rcm_pre_samz, son_obs_pre_samz)
nse_ann_rcm_pre_samz = compute_ioa(ann_rcm_pre_samz, ann_obs_pre_samz)
nse_djf_rcm_pre_eneb = compute_ioa(djf_rcm_pre_eneb, djf_obs_pre_eneb)
nse_mam_rcm_pre_eneb = compute_ioa(mam_rcm_pre_eneb, mam_obs_pre_eneb)
nse_jja_rcm_pre_eneb = compute_ioa(jja_rcm_pre_eneb, jja_obs_pre_eneb)
nse_son_rcm_pre_eneb = compute_ioa(son_rcm_pre_eneb, son_obs_pre_eneb)
nse_ann_rcm_pre_eneb = compute_ioa(ann_rcm_pre_eneb, ann_obs_pre_eneb)
nse_djf_rcm_pre_matopiba = compute_ioa(djf_rcm_pre_matopiba, djf_obs_pre_matopiba)
nse_mam_rcm_pre_matopiba = compute_ioa(mam_rcm_pre_matopiba, mam_obs_pre_matopiba)
nse_jja_rcm_pre_matopiba = compute_ioa(jja_rcm_pre_matopiba, jja_obs_pre_matopiba)
nse_son_rcm_pre_matopiba = compute_ioa(son_rcm_pre_matopiba, son_obs_pre_matopiba)
nse_ann_rcm_pre_matopiba = compute_ioa(ann_rcm_pre_matopiba, ann_obs_pre_matopiba)

# HadGEM2-ES
nse_djf_gcm_pre_samz = compute_ioa(djf_gcm_pre_samz, djf_obs_pre_samz)
nse_mam_gcm_pre_samz = compute_ioa(mam_gcm_pre_samz, mam_obs_pre_samz)
nse_jja_gcm_pre_samz = compute_ioa(jja_gcm_pre_samz, jja_obs_pre_samz)
nse_son_gcm_pre_samz = compute_ioa(son_gcm_pre_samz, son_obs_pre_samz)
nse_ann_gcm_pre_samz = compute_ioa(ann_gcm_pre_samz, ann_obs_pre_samz)
nse_djf_gcm_pre_eneb = compute_ioa(djf_gcm_pre_eneb, djf_obs_pre_eneb)
nse_mam_gcm_pre_eneb = compute_ioa(mam_gcm_pre_eneb, mam_obs_pre_eneb)
nse_jja_gcm_pre_eneb = compute_ioa(jja_gcm_pre_eneb, jja_obs_pre_eneb)
nse_son_gcm_pre_eneb = compute_ioa(son_gcm_pre_eneb, son_obs_pre_eneb)
nse_ann_gcm_pre_eneb = compute_ioa(ann_gcm_pre_eneb, ann_obs_pre_eneb)
nse_djf_gcm_pre_matopiba = compute_ioa(djf_gcm_pre_matopiba, djf_obs_pre_matopiba)
nse_mam_gcm_pre_matopiba = compute_ioa(mam_gcm_pre_matopiba, mam_obs_pre_matopiba)
nse_jja_gcm_pre_matopiba = compute_ioa(jja_gcm_pre_matopiba, jja_obs_pre_matopiba)
nse_son_gcm_pre_matopiba = compute_ioa(son_gcm_pre_matopiba, son_obs_pre_matopiba)
nse_ann_gcm_pre_matopiba = compute_ioa(ann_gcm_pre_matopiba, ann_obs_pre_matopiba)

# Temperature
# RegCM
nse_djf_rcm_tas_samz = compute_ioa(np.nanmean(djf_rcm_tas_samz, axis=1), djf_obs_tas_samz)
nse_mam_rcm_tas_samz = compute_ioa(np.nanmean(mam_rcm_tas_samz, axis=1), mam_obs_tas_samz)
nse_jja_rcm_tas_samz = compute_ioa(np.nanmean(jja_rcm_tas_samz, axis=1), jja_obs_tas_samz)
nse_son_rcm_tas_samz = compute_ioa(np.nanmean(son_rcm_tas_samz, axis=1), son_obs_tas_samz)
nse_ann_rcm_tas_samz = compute_ioa(np.nanmean(ann_rcm_tas_samz, axis=1), ann_obs_tas_samz)
nse_djf_rcm_tas_eneb = compute_ioa(np.nanmean(djf_rcm_tas_eneb, axis=1), djf_obs_tas_eneb)
nse_mam_rcm_tas_eneb = compute_ioa(np.nanmean(mam_rcm_tas_eneb, axis=1), mam_obs_tas_eneb)
nse_jja_rcm_tas_eneb = compute_ioa(np.nanmean(jja_rcm_tas_eneb, axis=1), jja_obs_tas_eneb)
nse_son_rcm_tas_eneb = compute_ioa(np.nanmean(son_rcm_tas_eneb, axis=1), son_obs_tas_eneb)
nse_ann_rcm_tas_eneb = compute_ioa(np.nanmean(ann_rcm_tas_eneb, axis=1), ann_obs_tas_eneb)
nse_djf_rcm_tas_matopiba = compute_ioa(np.nanmean(djf_rcm_tas_matopiba, axis=1), djf_obs_tas_matopiba)
nse_mam_rcm_tas_matopiba = compute_ioa(np.nanmean(mam_rcm_tas_matopiba, axis=1), mam_obs_tas_matopiba)
nse_jja_rcm_tas_matopiba = compute_ioa(np.nanmean(jja_rcm_tas_matopiba, axis=1), jja_obs_tas_matopiba)
nse_son_rcm_tas_matopiba = compute_ioa(np.nanmean(son_rcm_tas_matopiba, axis=1), son_obs_tas_matopiba)
nse_ann_rcm_tas_matopiba = compute_ioa(np.nanmean(ann_rcm_tas_matopiba, axis=1), ann_obs_tas_matopiba)

# HadGEM2-ES
nse_djf_gcm_tas_samz = compute_ioa(djf_gcm_tas_samz, djf_obs_tas_samz)
nse_mam_gcm_tas_samz = compute_ioa(mam_gcm_tas_samz, mam_obs_tas_samz)
nse_jja_gcm_tas_samz = compute_ioa(jja_gcm_tas_samz, jja_obs_tas_samz)
nse_son_gcm_tas_samz = compute_ioa(son_gcm_tas_samz, son_obs_tas_samz)
nse_ann_gcm_tas_samz = compute_ioa(ann_gcm_tas_samz, ann_obs_tas_samz)
nse_djf_gcm_tas_eneb = compute_ioa(djf_gcm_tas_eneb, djf_obs_tas_eneb)
nse_mam_gcm_tas_eneb = compute_ioa(mam_gcm_tas_eneb, mam_obs_tas_eneb)
nse_jja_gcm_tas_eneb = compute_ioa(jja_gcm_tas_eneb, jja_obs_tas_eneb)
nse_son_gcm_tas_eneb = compute_ioa(son_gcm_tas_eneb, son_obs_tas_eneb)
nse_ann_gcm_tas_eneb = compute_ioa(ann_gcm_tas_eneb, ann_obs_tas_eneb)
nse_djf_gcm_tas_matopiba = compute_ioa(djf_gcm_tas_matopiba, djf_obs_tas_matopiba)
nse_mam_gcm_tas_matopiba = compute_ioa(mam_gcm_tas_matopiba, mam_obs_tas_matopiba)
nse_jja_gcm_tas_matopiba = compute_ioa(jja_gcm_tas_matopiba, jja_obs_tas_matopiba)
nse_son_gcm_tas_matopiba = compute_ioa(son_gcm_tas_matopiba, son_obs_tas_matopiba)
nse_ann_gcm_tas_matopiba = compute_ioa(ann_gcm_tas_matopiba, ann_obs_tas_matopiba)

nse_rcm_pre = np.array([[nse_djf_rcm_pre_matopiba, nse_mam_rcm_pre_matopiba, nse_jja_rcm_pre_matopiba, nse_son_rcm_pre_matopiba, nse_ann_rcm_pre_matopiba],
[nse_djf_rcm_pre_eneb, nse_mam_rcm_pre_eneb, nse_jja_rcm_pre_eneb, nse_son_rcm_pre_eneb, nse_ann_rcm_pre_eneb], 
[nse_djf_rcm_pre_samz, nse_mam_rcm_pre_samz, nse_jja_rcm_pre_samz, nse_son_rcm_pre_samz, nse_ann_rcm_pre_samz]])

nse_gcm_pre = np.array([[nse_djf_gcm_pre_matopiba, nse_mam_gcm_pre_matopiba, nse_jja_gcm_pre_matopiba, nse_son_gcm_pre_matopiba, nse_ann_gcm_pre_matopiba],
[nse_djf_gcm_pre_eneb, nse_mam_gcm_pre_eneb, nse_jja_gcm_pre_eneb, nse_son_gcm_pre_eneb, nse_ann_gcm_pre_eneb], 
[nse_djf_gcm_pre_samz, nse_mam_gcm_pre_samz, nse_jja_gcm_pre_samz, nse_son_gcm_pre_samz, nse_ann_gcm_pre_samz]])

nse_rcm_tas = np.array([[nse_djf_rcm_tas_matopiba, nse_mam_rcm_tas_matopiba, nse_jja_rcm_tas_matopiba, nse_son_rcm_tas_matopiba, nse_ann_rcm_tas_matopiba],
[nse_djf_rcm_tas_eneb, nse_mam_rcm_tas_eneb, nse_jja_rcm_tas_eneb, nse_son_rcm_tas_eneb, nse_ann_rcm_tas_eneb], 
[nse_djf_rcm_tas_samz, nse_mam_rcm_tas_samz, nse_jja_rcm_tas_samz, nse_son_rcm_tas_samz, nse_ann_rcm_tas_samz]])

nse_gcm_tas = np.array([[nse_djf_gcm_tas_matopiba, nse_mam_gcm_tas_matopiba, nse_jja_gcm_tas_matopiba, nse_son_gcm_tas_matopiba, nse_ann_gcm_tas_matopiba],
[nse_djf_gcm_tas_eneb, nse_mam_gcm_tas_eneb, nse_jja_gcm_tas_eneb, nse_son_gcm_tas_eneb, nse_ann_gcm_tas_eneb], 
[nse_djf_gcm_tas_samz, nse_mam_gcm_tas_samz, nse_jja_gcm_tas_samz, nse_son_gcm_tas_samz, nse_ann_gcm_tas_samz]])

print(nse_rcm_pre)
print(nse_gcm_pre)

# Plot models and obs database 
fig, axes = plt.subplots(nrows=2, ncols=2, constrained_layout=True)
norm = colors.BoundaryNorm(boundaries=np.arange(0, 1, 0.1), ncolors=256)

xlabels = [u'DJF', u'MAM', u'JJA', u'SON', u'ANN']
ylabels = [u'MATOPIBA', u'ENEB', u'SAMZ']

# First column heatmaps with same colormap
pcm1 = axes[0, 0].pcolormesh(nse_rcm_pre, edgecolors ='k', linewidths = 1.5, norm=norm, cmap='Greys')
axes[0, 0].set_title(u'A) RegCM4.7', loc='left', fontweight='bold', fontsize=8)
axes[0, 0].set_xticks(np.arange(nse_rcm_pre.shape[1]) + 0.5)
axes[0, 0].set_yticks(np.arange(nse_rcm_pre.shape[0]) + 0.5)
axes[0, 0].set_yticklabels(ylabels, fontsize=8)
plt.setp(axes[0, 0].get_xticklabels(), visible=False)
axes[0, 0].tick_params(bottom = False)
for y in range(nse_rcm_pre.shape[0]):
    for x in range(nse_rcm_pre.shape[1]):
        axes[0, 0].text(x + 0.5, y + 0.5, '%.2f' % nse_rcm_pre[y, x],
                 ha="center", va="center", color='k',
                 )
                 
pcm2 = axes[1, 0].pcolormesh(nse_gcm_pre, edgecolors ='k', linewidths = 1.5, norm=norm, cmap='Greys')
axes[1, 0].set_title(u'B) HadGEM2-ES', loc='left', fontweight='bold', fontsize=8)
axes[1, 0].set_xticks(np.arange(nse_gcm_pre.shape[1]) + 0.5)
axes[1, 0].set_yticks(np.arange(nse_gcm_pre.shape[0]) + 0.5)
axes[1, 0].set_xticklabels(xlabels, fontsize=8)
axes[1, 0].set_yticklabels(ylabels, fontsize=8)
clb=fig.colorbar(pcm2, ax=axes[:, 0], extend='max', shrink=0.7)
clb.set_label(u'IOA - Precipitation', rotation=90)
clb.ax.yaxis.set_label_position('right')
for y in range(nse_gcm_pre.shape[0]):
    for x in range(nse_gcm_pre.shape[1]):
        axes[1, 0].text(x + 0.5, y + 0.5, '%.2f' % nse_gcm_pre[y, x],
                 ha="center", va="center", color='k',
                 )
                 
# Second column heatmaps with same colormap
pcm3 = axes[0, 1].pcolormesh(nse_rcm_tas, edgecolors ='k', linewidths = 1.5, norm=norm, cmap='Greys')
axes[0, 1].set_title('C) RegCM4.7', loc='left', fontweight='bold', fontsize=8)
axes[0, 1].set_xticks(np.arange(nse_rcm_tas.shape[1]) + 0.5)
axes[0, 1].set_yticks(np.arange(nse_rcm_tas.shape[0]) + 0.5)
axes[0, 1].set_yticklabels(ylabels, fontsize=8)
axes[0, 1].tick_params(left = False)
axes[0, 1].tick_params(bottom = False)
plt.setp(axes[0, 1].get_xticklabels(), visible=False)
plt.setp(axes[0, 1].get_yticklabels(), visible=False)
for y in range(nse_rcm_tas.shape[0]):
    for x in range(nse_rcm_tas.shape[1]):
        axes[0, 1].text(x + 0.5, y + 0.5, '%.2f' % nse_rcm_tas[y, x],
                 ha="center", va="center", color='k',
                 )
                 
pcm4 = axes[1, 1].pcolormesh(nse_gcm_tas, edgecolors ='k', linewidths = 1.5, norm=norm, cmap='Greys')
axes[1, 1].set_title('D) HadGEM2-ES', loc='left', fontweight='bold', fontsize=8)
axes[1, 1].set_xticks(np.arange(nse_gcm_tas.shape[1]) + 0.5)
axes[1, 1].set_yticks(np.arange(nse_gcm_tas.shape[0]) + 0.5)
axes[1, 1].set_xticklabels(xlabels, fontsize=8)
axes[1, 1].set_yticklabels(ylabels, fontsize=8)
axes[1, 1].tick_params(left = False)
plt.setp(axes[1, 1].get_yticklabels(), visible=False)
clb=fig.colorbar(pcm4, ax=axes[:, 1], extend='max', shrink=0.7)
clb.set_label(u'IOA - Temperature', rotation=90)
clb.ax.yaxis.set_label_position('right')
for y in range(nse_gcm_tas.shape[0]):
    for x in range(nse_gcm_tas.shape[1]):
        axes[1, 1].text(x + 0.5, y + 0.5, '%.2f' % nse_gcm_tas[y, x],
                 ha="center", va="center", color='k',
                 )
                       
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_portrait_diagram_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)	
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()


