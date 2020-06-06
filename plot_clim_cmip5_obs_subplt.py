# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/08/2019"
__description__ = "This script plot climatology graphics from CMIP5 models end OBS basedata"

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from sklearn import metrics
from scipy.stats import norm
from matplotlib.font_manager import FontProperties
from comp_statist_indices import compute_corr, compute_rmse, compute_pbias


def import_cmip5_clim(param, area, model):

	path  = '/home/nice/Documents/ufrn/phd_project/datas/cmip5/hist'
	arq   = '{0}/{1}_{2}_Amon_{3}_historical_r1i1p1_197512-200511.nc'.format(path, param, area,	model)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mdl_data = np.nanmean(np.nanmean(value, axis=1), axis=1)

	mdl_clim = []
	for mon in range(1, 12 + 1):
		mdl = np.nanmean(mdl_data[mon::12], axis=0)
		mdl_clim.append(mdl)
	return mdl_clim


def import_obs_climp(param, area, database):

	path  = '/home/nice/Documents/ufrn/phd_project/datas/obs_data'
	arq   = '{0}/{1}_{2}_{3}_obs_mon_197512-200511.nc'.format(path, param, area, database)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	obs_data = np.nanmean(np.nanmean(value, axis=1), axis=1)

	obs_clim = []
	for mon in range(1, 12 + 1):
		obs = np.nanmean(obs_data[mon::12], axis=0)
		obs_clim.append(obs)
	return obs_clim
	              

# Import cmip5 model end obs database climatology
pre_amz_gcm1 = import_cmip5_clim(u'pr', u'amz', u'BCC-CSM1.1')
tmp_amz_gcm1 = import_cmip5_clim(u'tas', u'amz', u'BCC-CSM1.1')
pre_neb_gcm1 = import_cmip5_clim(u'pr', u'neb', u'BCC-CSM1.1')
tmp_neb_gcm1 = import_cmip5_clim(u'tas', u'neb', u'BCC-CSM1.1')
pre_mato_gcm1 = import_cmip5_clim(u'pr', u'matopiba', u'BCC-CSM1.1')
tmp_mato_gcm1 = import_cmip5_clim(u'tas', u'matopiba', u'BCC-CSM1.1')

pre_amz_gcm2 = import_cmip5_clim(u'pr', u'amz', u'BCC-CSM1.1M')
tmp_amz_gcm2 = import_cmip5_clim(u'tas', u'amz', u'BCC-CSM1.1M')
pre_neb_gcm2 = import_cmip5_clim(u'pr', u'neb', u'BCC-CSM1.1M')
tmp_neb_gcm2 = import_cmip5_clim(u'tas', u'neb', u'BCC-CSM1.1M')
pre_mato_gcm2 = import_cmip5_clim(u'pr', u'matopiba', u'BCC-CSM1.1M')
tmp_mato_gcm2 = import_cmip5_clim(u'tas', u'matopiba', u'BCC-CSM1.1M')

pre_amz_gcm3 = import_cmip5_clim(u'pr', u'amz', u'BNU-ESM')
tmp_amz_gcm3 = import_cmip5_clim(u'tas', u'amz', u'BNU-ESM')
pre_neb_gcm3 = import_cmip5_clim(u'pr', u'neb', u'BNU-ESM')
tmp_neb_gcm3 = import_cmip5_clim(u'tas', u'neb', u'BNU-ESM')
pre_mato_gcm3 = import_cmip5_clim(u'pr', u'matopiba', u'BNU-ESM')
tmp_mato_gcm3 = import_cmip5_clim(u'tas', u'matopiba', u'BNU-ESM')

pre_amz_gcm4 = import_cmip5_clim(u'pr', u'amz', u'CanESM2')
tmp_amz_gcm4 = import_cmip5_clim(u'tas', u'amz', u'CanESM2')
pre_neb_gcm4 = import_cmip5_clim(u'pr', u'neb', u'CanESM2')
tmp_neb_gcm4 = import_cmip5_clim(u'tas', u'neb', u'CanESM2')
pre_mato_gcm4 = import_cmip5_clim(u'pr', u'matopiba', u'CanESM2')
tmp_mato_gcm4 = import_cmip5_clim(u'tas', u'matopiba', u'CanESM2')

pre_amz_gcm5 = import_cmip5_clim(u'pr', u'amz', u'CNRM-CM5')
tmp_amz_gcm5 = import_cmip5_clim(u'tas', u'amz', u'CNRM-CM5')
pre_neb_gcm5 = import_cmip5_clim(u'pr', u'neb', u'CNRM-CM5')
tmp_neb_gcm5 = import_cmip5_clim(u'tas', u'neb', u'CNRM-CM5')
pre_mato_gcm5 = import_cmip5_clim(u'pr', u'matopiba', u'CNRM-CM5')
tmp_mato_gcm5 = import_cmip5_clim(u'tas', u'matopiba', u'CNRM-CM5')

pre_amz_gcm6 = import_cmip5_clim(u'pr', u'amz', u'CSIRO-ACCESS-1')
tmp_amz_gcm6 = import_cmip5_clim(u'tas', u'amz', u'CSIRO-ACCESS-1')
pre_neb_gcm6 = import_cmip5_clim(u'pr', u'neb', u'CSIRO-ACCESS-1')
tmp_neb_gcm6 = import_cmip5_clim(u'tas', u'neb', u'CSIRO-ACCESS-1')
pre_mato_gcm6 = import_cmip5_clim(u'pr', u'matopiba', u'CSIRO-ACCESS-1')
tmp_mato_gcm6 = import_cmip5_clim(u'tas', u'matopiba', u'CSIRO-ACCESS-1')

pre_amz_gcm7 = import_cmip5_clim(u'pr', u'amz', u'CSIRO-ACCESS-3')
tmp_amz_gcm7 = import_cmip5_clim(u'tas', u'amz', u'CSIRO-ACCESS-3')
pre_neb_gcm7 = import_cmip5_clim(u'pr', u'neb', u'CSIRO-ACCESS-3')
tmp_neb_gcm7 = import_cmip5_clim(u'tas', u'neb', u'CSIRO-ACCESS-3')
pre_mato_gcm7 = import_cmip5_clim(u'pr', u'matopiba', u'CSIRO-ACCESS-3')
tmp_mato_gcm7 = import_cmip5_clim(u'tas', u'matopiba', u'CSIRO-ACCESS-3')

pre_amz_gcm8 = import_cmip5_clim(u'pr', u'amz', u'CSIRO-MK36')
tmp_amz_gcm8 = import_cmip5_clim(u'tas', u'amz', u'CSIRO-MK36')
pre_neb_gcm8 = import_cmip5_clim(u'pr', u'neb', u'CSIRO-MK36')
tmp_neb_gcm8 = import_cmip5_clim(u'tas', u'neb', u'CSIRO-MK36')
pre_mato_gcm8 = import_cmip5_clim(u'pr', u'matopiba', u'CSIRO-MK36')
tmp_mato_gcm8 = import_cmip5_clim(u'tas', u'matopiba', u'CSIRO-MK36')

pre_amz_gcm9 = import_cmip5_clim(u'pr', u'amz', u'FIO-ESM')
tmp_amz_gcm9 = import_cmip5_clim(u'tas', u'amz', u'FIO-ESM')
pre_neb_gcm9 = import_cmip5_clim(u'pr', u'neb', u'FIO-ESM')
tmp_neb_gcm9 = import_cmip5_clim(u'tas', u'neb', u'FIO-ESM')
pre_mato_gcm9 = import_cmip5_clim(u'pr', u'matopiba', u'FIO-ESM')
tmp_mato_gcm9 = import_cmip5_clim(u'tas', u'matopiba', u'FIO-ESM')

pre_amz_gcm10 = import_cmip5_clim(u'pr', u'amz', u'GISS-E2-H-CC')
tmp_amz_gcm10 = import_cmip5_clim(u'tas', u'amz', u'GISS-E2-H-CC')
pre_neb_gcm10 = import_cmip5_clim(u'pr', u'neb', u'GISS-E2-H-CC')
tmp_neb_gcm10 = import_cmip5_clim(u'tas', u'neb', u'GISS-E2-H-CC')
pre_mato_gcm10 = import_cmip5_clim(u'pr', u'matopiba', u'GISS-E2-H-CC')
tmp_mato_gcm10 = import_cmip5_clim(u'tas', u'matopiba', u'GISS-E2-H-CC')

pre_amz_gcm11 = import_cmip5_clim(u'pr', u'amz', u'GISS-E2-H')
tmp_amz_gcm11 = import_cmip5_clim(u'tas', u'amz', u'GISS-E2-H')
pre_neb_gcm11 = import_cmip5_clim(u'pr', u'neb', u'GISS-E2-H')
tmp_neb_gcm11 = import_cmip5_clim(u'tas', u'neb', u'GISS-E2-H')
pre_mato_gcm11 = import_cmip5_clim(u'pr', u'matopiba', u'GISS-E2-H')
tmp_mato_gcm11 = import_cmip5_clim(u'tas', u'matopiba', u'GISS-E2-H')

pre_amz_gcm12 = import_cmip5_clim(u'pr', u'amz', u'GISS-E2-R')
tmp_amz_gcm12 = import_cmip5_clim(u'tas', u'amz', u'GISS-E2-R')
pre_neb_gcm12 = import_cmip5_clim(u'pr', u'neb', u'GISS-E2-R')
tmp_neb_gcm12 = import_cmip5_clim(u'tas', u'neb', u'GISS-E2-R')
pre_mato_gcm12 = import_cmip5_clim(u'pr', u'matopiba', u'GISS-E2-R')
tmp_mato_gcm12 = import_cmip5_clim(u'tas', u'matopiba', u'GISS-E2-R')

pre_amz_gcm13 = import_cmip5_clim(u'pr', u'amz', u'HadGEM2-AO')
tmp_amz_gcm13 = import_cmip5_clim(u'tas', u'amz', u'HadGEM2-AO')
pre_neb_gcm13 = import_cmip5_clim(u'pr', u'neb', u'HadGEM2-AO')
tmp_neb_gcm13 = import_cmip5_clim(u'tas', u'neb', u'HadGEM2-AO')
pre_mato_gcm13 = import_cmip5_clim(u'pr', u'matopiba', u'HadGEM2-AO')
tmp_mato_gcm13 = import_cmip5_clim(u'tas', u'matopiba', u'HadGEM2-AO')

pre_amz_gcm14 = import_cmip5_clim(u'pr', u'amz', u'HadGEM2-CC')
tmp_amz_gcm14 = import_cmip5_clim(u'tas', u'amz', u'HadGEM2-CC')
pre_neb_gcm14 = import_cmip5_clim(u'pr', u'neb', u'HadGEM2-CC')
tmp_neb_gcm14 = import_cmip5_clim(u'tas', u'neb', u'HadGEM2-CC')
pre_mato_gcm14 = import_cmip5_clim(u'pr', u'matopiba', u'HadGEM2-CC')
tmp_mato_gcm14 = import_cmip5_clim(u'tas', u'matopiba', u'HadGEM2-CC')

pre_amz_gcm15 = import_cmip5_clim(u'pr', u'amz', u'HadGEM2-ES')
tmp_amz_gcm15 = import_cmip5_clim(u'tas', u'amz', u'HadGEM2-ES')
pre_neb_gcm15 = import_cmip5_clim(u'pr', u'neb', u'HadGEM2-ES')
tmp_neb_gcm15 = import_cmip5_clim(u'tas', u'neb', u'HadGEM2-ES')
pre_mato_gcm15 = import_cmip5_clim(u'pr', u'matopiba', u'HadGEM2-ES')
tmp_mato_gcm15 = import_cmip5_clim(u'tas', u'matopiba', u'HadGEM2-ES')

pre_amz_gcm16 = import_cmip5_clim(u'pr', u'amz', u'INMCM4')
tmp_amz_gcm16 = import_cmip5_clim(u'tas', u'amz', u'INMCM4')
pre_neb_gcm16 = import_cmip5_clim(u'pr', u'neb', u'INMCM4')
tmp_neb_gcm16 = import_cmip5_clim(u'tas', u'neb', u'INMCM4')
pre_mato_gcm16 = import_cmip5_clim(u'pr', u'matopiba', u'INMCM4')
tmp_mato_gcm16 = import_cmip5_clim(u'tas', u'matopiba', u'INMCM4')

pre_amz_gcm17 = import_cmip5_clim(u'pr', u'amz', u'IPSL-CM5A-LR')
tmp_amz_gcm17 = import_cmip5_clim(u'tas', u'amz', u'IPSL-CM5A-LR')
pre_neb_gcm17 = import_cmip5_clim(u'pr', u'neb', u'IPSL-CM5A-LR')
tmp_neb_gcm17 = import_cmip5_clim(u'tas', u'neb', u'IPSL-CM5A-LR')
pre_mato_gcm17 = import_cmip5_clim(u'pr', u'matopiba', u'IPSL-CM5A-LR')
tmp_mato_gcm17 = import_cmip5_clim(u'tas', u'matopiba', u'IPSL-CM5A-LR')

pre_amz_gcm18 = import_cmip5_clim(u'pr', u'amz', u'IPSL-CM5A-MR')
tmp_amz_gcm18 = import_cmip5_clim(u'tas', u'amz', u'IPSL-CM5A-MR')
pre_neb_gcm18 = import_cmip5_clim(u'pr', u'neb', u'IPSL-CM5A-MR')
tmp_neb_gcm18 = import_cmip5_clim(u'tas', u'neb', u'IPSL-CM5A-MR')
pre_mato_gcm18 = import_cmip5_clim(u'pr', u'matopiba', u'IPSL-CM5A-MR')
tmp_mato_gcm18 = import_cmip5_clim(u'tas', u'matopiba', u'IPSL-CM5A-MR')

pre_amz_gcm19 = import_cmip5_clim(u'pr', u'amz', u'IPSL-CM5-LR')
tmp_amz_gcm19 = import_cmip5_clim(u'tas', u'amz', u'IPSL-CM5A-LR')
pre_neb_gcm19 = import_cmip5_clim(u'pr', u'neb', u'IPSL-CM5A-LR')
tmp_neb_gcm19 = import_cmip5_clim(u'tas', u'neb', u'IPSL-CM5A-LR')
pre_mato_gcm19 = import_cmip5_clim(u'pr', u'matopiba', u'IPSL-CM5A-LR')
tmp_mato_gcm19 = import_cmip5_clim(u'tas', u'matopiba', u'IPSL-CM5A-LR')

pre_amz_obs  = import_obs_clim(u'pre', u'amz', u'cru_ts4.02')
tmp_amz_obs  = import_obs_clim(u'tmp', u'amz', u'cru_ts4.02')
pre_neb_obs  = import_obs_clim(u'pre', u'neb', u'cru_ts4.02')
tmp_neb_obs  = import_obs_clim(u'tmp', u'neb', u'cru_ts4.02')
pre_mato_obs  = import_obs_clim(u'pre', u'mato', u'cru_ts4.02')
tmp_mato_obs  = import_obs_clim(u'tmp', u'mato', u'cru_ts4.02')

# Plot model end obs data climatology
fig = plt.figure(figsize=(10, 8))
time = np.arange(0.5, 12 + 0.5)

ax1 = fig.add_subplot(321)
plt_clim1 = ax1.plot(time, mdl1_clim, time, mdl2_clim, time, mdl3_clim,
time,  mdl4_clim, time, mdl5_clim, time, mdl6_clim, time, mdl7_clim,
time, mdl8_clim, time, mdl9_clim, time, mdl10_clim, time, mdl11_clim,
time, mdl12_clim, time, mdl13_clim, time, mdl14_clim, time, mdl15_clim, 
time, mdl16_clim, time, mdl17_clim, time, mdl18_clim, time, mdl19_clim, 
time, mdl20_clim, time, mdl21_clim, time, mdl22_clim, time, mdl23_clim, 
time, mdl24_clim, time, mdl25_clim, time, mdl26_clim, time, mdl27_clim, 
time, mdl28_clim, time, mdl29_clim, time, mdl30_clim, time, mdl31_clim, 
time, mdl32_clim, time, mdl33_clim, time, obs1_clim, time, obs1_clim_p5, 
time, obs1_clim_p95)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 14, 2))
plt.title(u'A)', loc='left', fontweight='bold')
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')


l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35, l36 = plt_clim1
plt.setp(l1,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l2,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l3,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l4,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l5,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l6,  linewidth=1, markeredgewidth=1, color='green')
plt.setp(l7,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l8,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l9,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l10, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l11, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l12, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l13, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l14, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l15, linewidth=1, markeredgewidth=1, color='red')
plt.setp(l16, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l17, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l18, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l19, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l20, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l21, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l22, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l23, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l24, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l25, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l26, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l27, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l28, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l29, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l30, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l31, linewidth=1, markeredgewidth=1, color='gold')
plt.setp(l32, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l33, linewidth=1, markeredgewidth=1, color='blue')
plt.setp(l34, linewidth=1, markeredgewidth=1, color='black')
plt.setp(l35, linewidth=1, markeredgewidth=1, color='slategray')
plt.setp(l36, linewidth=1, markeredgewidth=1, color='slategray')

ax2 = fig.add_subplot(322)
plt_clim2 = ax2.plot(time, mdl1_clim, time, mdl2_clim, time, mdl3_clim,
time,  mdl4_clim, time, mdl5_clim, time, mdl6_clim, time, mdl7_clim,
time, mdl8_clim, time, mdl9_clim, time, mdl10_clim, time, mdl11_clim,
time, mdl12_clim, time, mdl13_clim, time, mdl14_clim, time, mdl15_clim, 
time, mdl16_clim, time, mdl17_clim, time, mdl18_clim, time, mdl19_clim, 
time, mdl20_clim, time, mdl21_clim, time, mdl22_clim, time, mdl23_clim, 
time, mdl24_clim, time, mdl25_clim, time, mdl26_clim, time, mdl27_clim, 
time, mdl28_clim, time, mdl29_clim, time, mdl30_clim, time, mdl31_clim, 
time, mdl32_clim, time, mdl33_clim, time, obs1_clim, time, obs1_clim_p5, 
time, obs1_clim_p95)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(18, 32, 2))
plt.title(u'B)', loc='left', fontweight='bold')
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')

l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35, l36 = plt_clim2
plt.setp(l1,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l2,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l3,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l4,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l5,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l6,  linewidth=1, markeredgewidth=1, color='green')
plt.setp(l7,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l8,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l9,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l10, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l11, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l12, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l13, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l14, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l15, linewidth=1, markeredgewidth=1, color='red')
plt.setp(l16, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l17, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l18, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l19, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l20, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l21, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l22, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l23, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l24, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l25, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l26, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l27, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l28, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l29, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l30, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l31, linewidth=1, markeredgewidth=1, color='gold')
plt.setp(l32, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l33, linewidth=1, markeredgewidth=1, color='blue')
plt.setp(l34, linewidth=1, markeredgewidth=1, color='black')
plt.setp(l35, linewidth=1, markeredgewidth=1, color='slategray')
plt.setp(l36, linewidth=1, markeredgewidth=1, color='slategray')

ax3 = fig.add_subplot(323)
plt_clim3 = ax3.plot(time, mdl1_clim, time, mdl2_clim, time, mdl3_clim,
time,  mdl4_clim, time, mdl5_clim, time, mdl6_clim, time, mdl7_clim,
time, mdl8_clim, time, mdl9_clim, time, mdl10_clim, time, mdl11_clim,
time, mdl12_clim, time, mdl13_clim, time, mdl14_clim, time, mdl15_clim, 
time, mdl16_clim, time, mdl17_clim, time, mdl18_clim, time, mdl19_clim, 
time, mdl20_clim, time, mdl21_clim, time, mdl22_clim, time, mdl23_clim, 
time, mdl24_clim, time, mdl25_clim, time, mdl26_clim, time, mdl27_clim, 
time, mdl28_clim, time, mdl29_clim, time, mdl30_clim, time, mdl31_clim, 
time, mdl32_clim, time, mdl33_clim, time, obs1_clim, time, obs1_clim_p5, 
time, obs1_clim_p95)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 14, 2))
plt.ylabel(u'Rainfall (mm d⁻¹)')
plt.title(u'C)', loc='left', fontweight='bold')
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')

l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35, l36 = plt_clim3
plt.setp(l1,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l2,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l3,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l4,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l5,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l6,  linewidth=1, markeredgewidth=1, color='green')
plt.setp(l7,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l8,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l9,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l10, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l11, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l12, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l13, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l14, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l15, linewidth=1, markeredgewidth=1, color='red')
plt.setp(l16, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l17, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l18, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l19, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l20, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l21, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l22, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l23, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l24, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l25, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l26, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l27, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l28, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l29, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l30, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l31, linewidth=1, markeredgewidth=1, color='gold')
plt.setp(l32, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l33, linewidth=1, markeredgewidth=1, color='blue')
plt.setp(l34, linewidth=1, markeredgewidth=1, color='black')
plt.setp(l35, linewidth=1, markeredgewidth=1, color='slategray')
plt.setp(l36, linewidth=1, markeredgewidth=1, color='slategray')

ax4 = fig.add_subplot(324)
plt_clim4 = ax4.plot(time, mdl1_clim, time, mdl2_clim, time, mdl3_clim,
time,  mdl4_clim, time, mdl5_clim, time, mdl6_clim, time, mdl7_clim,
time, mdl8_clim, time, mdl9_clim, time, mdl10_clim, time, mdl11_clim,
time, mdl12_clim, time, mdl13_clim, time, mdl14_clim, time, mdl15_clim, 
time, mdl16_clim, time, mdl17_clim, time, mdl18_clim, time, mdl19_clim, 
time, mdl20_clim, time, mdl21_clim, time, mdl22_clim, time, mdl23_clim, 
time, mdl24_clim, time, mdl25_clim, time, mdl26_clim, time, mdl27_clim, 
time, mdl28_clim, time, mdl29_clim, time, mdl30_clim, time, mdl31_clim, 
time, mdl32_clim, time, mdl33_clim, time, obs1_clim, time, obs1_clim_p5, 
time, obs1_clim_p95)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(18, 32, 2))
plt.ylabel(u'Temperature (°C)')
plt.title(u'D)', loc='left', fontweight='bold')
ax4.xaxis.grid(True, which='major', linestyle='--')
ax4.yaxis.grid(True, which='major', linestyle='--')

l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35, l36 = plt_clim4
plt.setp(l1,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l2,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l3,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l4,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l5,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l6,  linewidth=1, markeredgewidth=1, color='green')
plt.setp(l7,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l8,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l9,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l10, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l11, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l12, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l13, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l14, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l15, linewidth=1, markeredgewidth=1, color='red')
plt.setp(l16, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l17, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l18, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l19, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l20, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l21, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l22, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l23, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l24, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l25, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l26, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l27, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l28, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l29, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l30, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l31, linewidth=1, markeredgewidth=1, color='gold')
plt.setp(l32, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l33, linewidth=1, markeredgewidth=1, color='blue')
plt.setp(l34, linewidth=1, markeredgewidth=1, color='black')
plt.setp(l35, linewidth=1, markeredgewidth=1, color='slategray')
plt.setp(l36, linewidth=1, markeredgewidth=1, color='slategray')

ax5 = fig.add_subplot(325)
plt_clim5 = ax5.plot(time, mdl1_clim, time, mdl2_clim, time, mdl3_clim,
time,  mdl4_clim, time, mdl5_clim, time, mdl6_clim, time, mdl7_clim,
time, mdl8_clim, time, mdl9_clim, time, mdl10_clim, time, mdl11_clim,
time, mdl12_clim, time, mdl13_clim, time, mdl14_clim, time, mdl15_clim, 
time, mdl16_clim, time, mdl17_clim, time, mdl18_clim, time, mdl19_clim, 
time, mdl20_clim, time, mdl21_clim, time, mdl22_clim, time, mdl23_clim, 
time, mdl24_clim, time, mdl25_clim, time, mdl26_clim, time, mdl27_clim, 
time, mdl28_clim, time, mdl29_clim, time, mdl30_clim, time, mdl31_clim, 
time, mdl32_clim, time, mdl33_clim, time, obs1_clim, time, obs1_clim_p5, 
time, obs1_clim_p95)
plt.xlabel('Months')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 14, 2))
plt.title(u'E)', loc='left', fontweight='bold')
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')

l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35, l36 = plt_clim5
plt.setp(l1,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l2,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l3,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l4,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l5,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l6,  linewidth=1, markeredgewidth=1, color='green')
plt.setp(l7,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l8,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l9,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l10, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l11, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l12, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l13, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l14, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l15, linewidth=1, markeredgewidth=1, color='red')
plt.setp(l16, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l17, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l18, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l19, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l20, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l21, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l22, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l23, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l24, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l25, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l26, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l27, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l28, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l29, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l30, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l31, linewidth=1, markeredgewidth=1, color='gold')
plt.setp(l32, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l33, linewidth=1, markeredgewidth=1, color='blue')
plt.setp(l34, linewidth=1, markeredgewidth=1, color='black')
plt.setp(l35, linewidth=1, markeredgewidth=1, color='slategray')
plt.setp(l36, linewidth=1, markeredgewidth=1, color='slategray')

legend = ('M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 
'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20', 'M21', 'M22', 'M23', 'M24', 'M25', 'M26',
'M27', 'M28', 'M29', 'M30', 'M31', 'M32', 'M33', 'OBS', 'OBS5%', 'OBS95%')
plt.legend(plt_clim5, legend, loc='lower left', bbox_to_anchor=(-0.1, -.9), ncol=9)

ax6 = fig.add_subplot(326)
plt_clim6 = ax6.plot(time, mdl1_clim, time, mdl2_clim, time, mdl3_clim,
time,  mdl4_clim, time, mdl5_clim, time, mdl6_clim, time, mdl7_clim,
time, mdl8_clim, time, mdl9_clim, time, mdl10_clim, time, mdl11_clim,
time, mdl12_clim, time, mdl13_clim, time, mdl14_clim, time, mdl15_clim, 
time, mdl16_clim, time, mdl17_clim, time, mdl18_clim, time, mdl19_clim, 
time, mdl20_clim, time, mdl21_clim, time, mdl22_clim, time, mdl23_clim, 
time, mdl24_clim, time, mdl25_clim, time, mdl26_clim, time, mdl27_clim, 
time, mdl28_clim, time, mdl29_clim, time, mdl30_clim, time, mdl31_clim, 
time, mdl32_clim, time, mdl33_clim, time, obs1_clim, time, obs1_clim_p5, 
time, obs1_clim_p95)
plt.xlabel('Months')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(18, 32, 2))
plt.title(u'F)', loc='left', fontweight='bold')
ax6.xaxis.grid(True, which='major', linestyle='--')
ax6.yaxis.grid(True, which='major', linestyle='--')

l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35, l36 = plt_clim6
plt.setp(l1,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l2,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l3,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l4,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l5,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l6,  linewidth=1, markeredgewidth=1, color='green')
plt.setp(l7,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l8,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l9,  linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l10, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l11, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l12, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l13, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l14, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l15, linewidth=1, markeredgewidth=1, color='red')
plt.setp(l16, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l17, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l18, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l19, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l20, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l21, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l22, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l23, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l24, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l25, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l26, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l27, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l28, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l29, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l30, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l31, linewidth=1, markeredgewidth=1, color='gold')
plt.setp(l32, linewidth=1, markeredgewidth=1, color='gainsboro')
plt.setp(l33, linewidth=1, markeredgewidth=1, color='blue')
plt.setp(l34, linewidth=1, markeredgewidth=1, color='black')
plt.setp(l35, linewidth=1, markeredgewidth=1, color='slategray')
plt.setp(l36, linewidth=1, markeredgewidth=1, color='slategray')

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.20, hspace=0.35)

# Save figure
path_out = '/home/nice'
name_out = 'pyplt_clim_cmip5_cru_1975-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()



