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

	path  = '/home/nice/Documents/dataset/gcm/cmip5'
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


def import_obs_clim(param, area, database):

	path  = '/home/nice/Documents/dataset/obs/gcm'
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

pre_amz_gcm10 = import_cmip5_clim(u'pr', u'amz', u'GISS-E2-H')
tmp_amz_gcm10 = import_cmip5_clim(u'tas', u'amz', u'GISS-E2-H')
pre_neb_gcm10 = import_cmip5_clim(u'pr', u'neb', u'GISS-E2-H')
tmp_neb_gcm10 = import_cmip5_clim(u'tas', u'neb', u'GISS-E2-H')
pre_mato_gcm10 = import_cmip5_clim(u'pr', u'matopiba', u'GISS-E2-H')
tmp_mato_gcm10 = import_cmip5_clim(u'tas', u'matopiba', u'GISS-E2-H')

pre_amz_gcm11 = import_cmip5_clim(u'pr', u'amz', u'GISS-E2-H-CC')
tmp_amz_gcm11 = import_cmip5_clim(u'tas', u'amz', u'GISS-E2-H-CC')
pre_neb_gcm11 = import_cmip5_clim(u'pr', u'neb', u'GISS-E2-H-CC')
tmp_neb_gcm11 = import_cmip5_clim(u'tas', u'neb', u'GISS-E2-H-CC')
pre_mato_gcm11 = import_cmip5_clim(u'pr', u'matopiba', u'GISS-E2-H-CC')
tmp_mato_gcm11 = import_cmip5_clim(u'tas', u'matopiba', u'GISS-E2-H-CC')

pre_amz_gcm12 = import_cmip5_clim(u'pr', u'amz', u'HadGEM2-AO')
tmp_amz_gcm12 = import_cmip5_clim(u'tas', u'amz', u'HadGEM2-AO')
pre_neb_gcm12 = import_cmip5_clim(u'pr', u'neb', u'HadGEM2-AO')
tmp_neb_gcm12 = import_cmip5_clim(u'tas', u'neb', u'HadGEM2-AO')
pre_mato_gcm12 = import_cmip5_clim(u'pr', u'matopiba', u'HadGEM2-AO')
tmp_mato_gcm12 = import_cmip5_clim(u'tas', u'matopiba', u'HadGEM2-AO')

pre_amz_gcm13 = import_cmip5_clim(u'pr', u'amz', u'HadGEM2-CC')
tmp_amz_gcm13 = import_cmip5_clim(u'tas', u'amz', u'HadGEM2-CC')
pre_neb_gcm13 = import_cmip5_clim(u'pr', u'neb', u'HadGEM2-CC')
tmp_neb_gcm13 = import_cmip5_clim(u'tas', u'neb', u'HadGEM2-CC')
pre_mato_gcm13 = import_cmip5_clim(u'pr', u'matopiba', u'HadGEM2-CC')
tmp_mato_gcm13 = import_cmip5_clim(u'tas', u'matopiba', u'HadGEM2-CC')

pre_amz_gcm14 = import_cmip5_clim(u'pr', u'amz', u'HadGEM2-ES')
tmp_amz_gcm14 = import_cmip5_clim(u'tas', u'amz', u'HadGEM2-ES')
pre_neb_gcm14 = import_cmip5_clim(u'pr', u'neb', u'HadGEM2-ES')
tmp_neb_gcm14 = import_cmip5_clim(u'tas', u'neb', u'HadGEM2-ES')
pre_mato_gcm14 = import_cmip5_clim(u'pr', u'matopiba', u'HadGEM2-ES')
tmp_mato_gcm14 = import_cmip5_clim(u'tas', u'matopiba', u'HadGEM2-ES')

pre_amz_gcm15 = import_cmip5_clim(u'pr', u'amz', u'INMCM4')
tmp_amz_gcm15 = import_cmip5_clim(u'tas', u'amz', u'INMCM4')
pre_neb_gcm15 = import_cmip5_clim(u'pr', u'neb', u'INMCM4')
tmp_neb_gcm15 = import_cmip5_clim(u'tas', u'neb', u'INMCM4')
pre_mato_gcm15 = import_cmip5_clim(u'pr', u'matopiba', u'INMCM4')
tmp_mato_gcm15 = import_cmip5_clim(u'tas', u'matopiba', u'INMCM4')

pre_amz_gcm16 = import_cmip5_clim(u'pr', u'amz', u'IPSL-CM5A-LR')
tmp_amz_gcm16 = import_cmip5_clim(u'tas', u'amz', u'IPSL-CM5A-LR')
pre_neb_gcm16 = import_cmip5_clim(u'pr', u'neb', u'IPSL-CM5A-LR')
tmp_neb_gcm16 = import_cmip5_clim(u'tas', u'neb', u'IPSL-CM5A-LR')
pre_mato_gcm16 = import_cmip5_clim(u'pr', u'matopiba', u'IPSL-CM5A-LR')
tmp_mato_gcm16 = import_cmip5_clim(u'tas', u'matopiba', u'IPSL-CM5A-LR')

pre_amz_gcm17 = import_cmip5_clim(u'pr', u'amz', u'IPSL-CM5A-MR')
tmp_amz_gcm17 = import_cmip5_clim(u'tas', u'amz', u'IPSL-CM5A-MR')
pre_neb_gcm17 = import_cmip5_clim(u'pr', u'neb', u'IPSL-CM5A-MR')
tmp_neb_gcm17 = import_cmip5_clim(u'tas', u'neb', u'IPSL-CM5A-MR')
pre_mato_gcm17 = import_cmip5_clim(u'pr', u'matopiba', u'IPSL-CM5A-MR')
tmp_mato_gcm17 = import_cmip5_clim(u'tas', u'matopiba', u'IPSL-CM5A-MR')

pre_amz_gcm18 = import_cmip5_clim(u'pr', u'amz', u'LASG-FGOALS-G2')
tmp_amz_gcm18 = import_cmip5_clim(u'tas', u'amz', u'LASG-FGOALS-G2')
pre_neb_gcm18 = import_cmip5_clim(u'pr', u'neb', u'LASG-FGOALS-G2')
tmp_neb_gcm18 = import_cmip5_clim(u'tas', u'neb', u'LASG-FGOALS-G2')
pre_mato_gcm18 = import_cmip5_clim(u'pr', u'matopiba', u'LASG-FGOALS-G2')
tmp_mato_gcm18 = import_cmip5_clim(u'tas', u'matopiba', u'LASG-FGOALS-G2')

pre_amz_gcm19 = import_cmip5_clim(u'pr', u'amz', u'LASG-FGOALS-S2')
tmp_amz_gcm19 = import_cmip5_clim(u'tas', u'amz', u'LASG-FGOALS-S2')
pre_neb_gcm19 = import_cmip5_clim(u'pr', u'neb', u'LASG-FGOALS-S2')
tmp_neb_gcm19 = import_cmip5_clim(u'tas', u'neb', u'LASG-FGOALS-S2')
pre_mato_gcm19 = import_cmip5_clim(u'pr', u'matopiba', u'LASG-FGOALS-S2')
tmp_mato_gcm19 = import_cmip5_clim(u'tas', u'matopiba', u'LASG-FGOALS-S2')

pre_amz_gcm20 = import_cmip5_clim(u'pr', u'amz', u'MIROC5')
tmp_amz_gcm20 = import_cmip5_clim(u'tas', u'amz', u'MIROC5')
pre_neb_gcm20 = import_cmip5_clim(u'pr', u'neb', u'MIROC5')
tmp_neb_gcm20 = import_cmip5_clim(u'tas', u'neb', u'MIROC5')
pre_mato_gcm20 = import_cmip5_clim(u'pr', u'matopiba', u'MIROC5')
tmp_mato_gcm20 = import_cmip5_clim(u'tas', u'matopiba', u'MIROC5')

pre_amz_gcm21 = import_cmip5_clim(u'pr', u'amz', u'MIROC-ESM')
tmp_amz_gcm21 = import_cmip5_clim(u'tas', u'amz', u'MIROC-ESM')
pre_neb_gcm21 = import_cmip5_clim(u'pr', u'neb', u'MIROC-ESM')
tmp_neb_gcm21 = import_cmip5_clim(u'tas', u'neb', u'MIROC-ESM')
pre_mato_gcm21 = import_cmip5_clim(u'pr', u'matopiba', u'MIROC-ESM')
tmp_mato_gcm21 = import_cmip5_clim(u'tas', u'matopiba', u'MIROC-ESM')

pre_amz_gcm22 = import_cmip5_clim(u'pr', u'amz', u'MIROC-ESM-CHEM')
tmp_amz_gcm22 = import_cmip5_clim(u'tas', u'amz', u'MIROC-ESM-CHEM')
pre_neb_gcm22 = import_cmip5_clim(u'pr', u'neb', u'MIROC-ESM-CHEM')
tmp_neb_gcm22 = import_cmip5_clim(u'tas', u'neb', u'MIROC-ESM-CHEM')
pre_mato_gcm22 = import_cmip5_clim(u'pr', u'matopiba', u'MIROC-ESM-CHEM')
tmp_mato_gcm22 = import_cmip5_clim(u'tas', u'matopiba', u'MIROC-ESM-CHEM')

pre_amz_gcm23 = import_cmip5_clim(u'pr', u'amz', u'MPI-ESM-LR')
tmp_amz_gcm23 = import_cmip5_clim(u'tas', u'amz', u'MPI-ESM-LR')
pre_neb_gcm23 = import_cmip5_clim(u'pr', u'neb', u'MPI-ESM-LR')
tmp_neb_gcm23 = import_cmip5_clim(u'tas', u'neb', u'MPI-ESM-LR')
pre_mato_gcm23 = import_cmip5_clim(u'pr', u'matopiba', u'MPI-ESM-LR')
tmp_mato_gcm23 = import_cmip5_clim(u'tas', u'matopiba', u'MPI-ESM-LR')

pre_amz_gcm24 = import_cmip5_clim(u'pr', u'amz', u'MPI-ESM-MR')
tmp_amz_gcm24 = import_cmip5_clim(u'tas', u'amz', u'MPI-ESM-MR')
pre_neb_gcm24 = import_cmip5_clim(u'pr', u'neb', u'MPI-ESM-MR')
tmp_neb_gcm24 = import_cmip5_clim(u'tas', u'neb', u'MPI-ESM-MR')
pre_mato_gcm24 = import_cmip5_clim(u'pr', u'matopiba', u'MPI-ESM-MR')
tmp_mato_gcm24 = import_cmip5_clim(u'tas', u'matopiba', u'MPI-ESM-MR')

pre_amz_gcm25 = import_cmip5_clim(u'pr', u'amz', u'MRI-CGCM3')
tmp_amz_gcm25 = import_cmip5_clim(u'tas', u'amz', u'MRI-CGCM3')
pre_neb_gcm25 = import_cmip5_clim(u'pr', u'neb', u'MRI-CGCM3')
tmp_neb_gcm25 = import_cmip5_clim(u'tas', u'neb', u'MRI-CGCM3')
pre_mato_gcm25 = import_cmip5_clim(u'pr', u'matopiba', u'MRI-CGCM3')
tmp_mato_gcm25 = import_cmip5_clim(u'tas', u'matopiba', u'MRI-CGCM3')

pre_amz_gcm26 = import_cmip5_clim(u'pr', u'amz', u'NCAR-CCSM4')
tmp_amz_gcm26 = import_cmip5_clim(u'tas', u'amz', u'NCAR-CCSM4')
pre_neb_gcm26 = import_cmip5_clim(u'pr', u'neb', u'NCAR-CCSM4')
tmp_neb_gcm26 = import_cmip5_clim(u'tas', u'neb', u'NCAR-CCSM4')
pre_mato_gcm26 = import_cmip5_clim(u'pr', u'matopiba', u'NCAR-CCSM4')
tmp_mato_gcm26 = import_cmip5_clim(u'tas', u'matopiba', u'NCAR-CCSM4')

pre_amz_gcm27 = import_cmip5_clim(u'pr', u'amz', u'NCAR-CESM1-BGC')
tmp_amz_gcm27 = import_cmip5_clim(u'tas', u'amz', u'NCAR-CESM1-BGC')
pre_neb_gcm27 = import_cmip5_clim(u'pr', u'neb', u'NCAR-CESM1-BGC')
tmp_neb_gcm27 = import_cmip5_clim(u'tas', u'neb', u'NCAR-CESM1-BGC')
pre_mato_gcm27 = import_cmip5_clim(u'pr', u'matopiba', u'NCAR-CESM1-BGC')
tmp_mato_gcm27 = import_cmip5_clim(u'tas', u'matopiba', u'NCAR-CESM1-BGC')

pre_amz_gcm28 = import_cmip5_clim(u'pr', u'amz', u'NCAR-CESM1-CAM5')
tmp_amz_gcm28 = import_cmip5_clim(u'tas', u'amz', u'NCAR-CESM1-CAM5')
pre_neb_gcm28 = import_cmip5_clim(u'pr', u'neb', u'NCAR-CESM1-CAM5')
tmp_neb_gcm28 = import_cmip5_clim(u'tas', u'neb', u'NCAR-CESM1-CAM5')
pre_mato_gcm28 = import_cmip5_clim(u'pr', u'matopiba', u'NCAR-CESM1-CAM5')
tmp_mato_gcm28 = import_cmip5_clim(u'tas', u'matopiba', u'NCAR-CESM1-CAM5')

pre_amz_gcm29 = import_cmip5_clim(u'pr', u'amz', u'NorESM1-M')
tmp_amz_gcm29 = import_cmip5_clim(u'tas', u'amz', u'NorESM1-M')
pre_neb_gcm29 = import_cmip5_clim(u'pr', u'neb', u'NorESM1-M')
tmp_neb_gcm29 = import_cmip5_clim(u'tas', u'neb', u'NorESM1-M')
pre_mato_gcm29 = import_cmip5_clim(u'pr', u'matopiba', u'NorESM1-M')
tmp_mato_gcm29 = import_cmip5_clim(u'tas', u'matopiba', u'NorESM1-M')

pre_amz_gcm30 = import_cmip5_clim(u'pr', u'amz', u'NorESM1-ME')
tmp_amz_gcm30 = import_cmip5_clim(u'tas', u'amz', u'NorESM1-ME')
pre_neb_gcm30 = import_cmip5_clim(u'pr', u'neb', u'NorESM1-ME')
tmp_neb_gcm30 = import_cmip5_clim(u'tas', u'neb', u'NorESM1-ME')
pre_mato_gcm30 = import_cmip5_clim(u'pr', u'matopiba', u'NorESM1-ME')
tmp_mato_gcm30 = import_cmip5_clim(u'tas', u'matopiba', u'NorESM1-ME')

pre_amz_gcm31 = import_cmip5_clim(u'pr', u'amz', u'ensmean_cmip5')
tmp_amz_gcm31 = import_cmip5_clim(u'tas', u'amz', u'ensmean_cmip5')
pre_neb_gcm31 = import_cmip5_clim(u'pr', u'neb', u'ensmean_cmip5')
tmp_neb_gcm31 = import_cmip5_clim(u'tas', u'neb', u'ensmean_cmip5')
pre_mato_gcm31 = import_cmip5_clim(u'pr', u'matopiba', u'ensmean_cmip5')
tmp_mato_gcm31 = import_cmip5_clim(u'tas', u'matopiba', u'ensmean_cmip5')

pre_amz_obs  = import_obs_clim(u'pre', u'amz', u'cru_ts4.02')
tmp_amz_obs  = import_obs_clim(u'tmp', u'amz', u'cru_ts4.02')
pre_neb_obs  = import_obs_clim(u'pre', u'neb', u'cru_ts4.02')
tmp_neb_obs  = import_obs_clim(u'tmp', u'neb', u'cru_ts4.02')
pre_mato_obs  = import_obs_clim(u'pre', u'matopiba', u'cru_ts4.02')
tmp_mato_obs  = import_obs_clim(u'tmp', u'matopiba', u'cru_ts4.02')

# Plot model end obs data climatology
fig = plt.figure(figsize=(10, 8))
time = np.arange(0.5, 12 + 0.5)

ax1 = fig.add_subplot(321)
plt_clim1 = ax1.plot(time, pre_amz_gcm1, time, pre_amz_gcm2, time, pre_amz_gcm3,
time, pre_amz_gcm4, time, pre_amz_gcm5, time, pre_amz_gcm6, time, pre_amz_gcm7,
time, pre_amz_gcm8, time, pre_amz_gcm9, time, pre_amz_gcm10, time, pre_amz_gcm11,
time, pre_amz_gcm12, time, pre_amz_gcm13, time, pre_amz_gcm14, time, pre_amz_gcm15, 
time, pre_amz_gcm16, time, pre_amz_gcm17, time, pre_amz_gcm18, time, pre_amz_gcm19, 
time, pre_amz_gcm20, time, pre_amz_gcm21, time, pre_amz_gcm22, time, pre_amz_gcm23, 
time, pre_amz_gcm24, time, pre_amz_gcm25, time, pre_amz_gcm26, time, pre_amz_gcm27, 
time, pre_amz_gcm28, time, pre_amz_gcm29, time, pre_amz_gcm30, time, pre_amz_gcm31, 
time, pre_amz_obs)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 21, 3))
plt.title(u'A)', loc='left', fontweight='bold')
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32 = plt_clim1
plt.setp(l1)
plt.setp(l2)
plt.setp(l3)
plt.setp(l4)
plt.setp(l5)
plt.setp(l6)
plt.setp(l7)
plt.setp(l8)
plt.setp(l9)
plt.setp(l10)
plt.setp(l11)
plt.setp(l12)
plt.setp(l13)
plt.setp(l14)
plt.setp(l15)
plt.setp(l16)
plt.setp(l17)
plt.setp(l18)
plt.setp(l19)
plt.setp(l20)
plt.setp(l21)
plt.setp(l22)
plt.setp(l23)
plt.setp(l24)
plt.setp(l25)
plt.setp(l26)
plt.setp(l27)
plt.setp(l28)
plt.setp(l29)
plt.setp(l30)
plt.setp(l31, linestyle='dashed', color='black')
plt.setp(l32, color='black')

ax2 = fig.add_subplot(322)
plt_clim2 = ax2.plot(time, tmp_amz_gcm1, time, tmp_amz_gcm2, time, tmp_amz_gcm3,
time, tmp_amz_gcm4, time, tmp_amz_gcm5, time, tmp_amz_gcm6, time, tmp_amz_gcm7,
time, tmp_amz_gcm8, time, tmp_amz_gcm9, time, tmp_amz_gcm10, time, tmp_amz_gcm11,
time, tmp_amz_gcm12, time, tmp_amz_gcm13, time, tmp_amz_gcm14, time, tmp_amz_gcm15, 
time, tmp_amz_gcm16, time, tmp_amz_gcm17, time, tmp_amz_gcm18, time, tmp_amz_gcm19, 
time, tmp_amz_gcm20, time, tmp_amz_gcm21, time, tmp_amz_gcm22, time, tmp_amz_gcm23, 
time, tmp_amz_gcm24, time, tmp_amz_gcm25, time, tmp_amz_gcm26, time, tmp_amz_gcm27, 
time, tmp_amz_gcm28, time, tmp_amz_gcm29, time, tmp_amz_gcm30, time, tmp_amz_gcm31, 
time, tmp_amz_obs)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(20, 38, 3))
plt.title(u'D)', loc='left', fontweight='bold')
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32 = plt_clim2
plt.setp(l1)
plt.setp(l2)
plt.setp(l3)
plt.setp(l4)
plt.setp(l5)
plt.setp(l6)
plt.setp(l7)
plt.setp(l8)
plt.setp(l9)
plt.setp(l10)
plt.setp(l11)
plt.setp(l12)
plt.setp(l13)
plt.setp(l14)
plt.setp(l15)
plt.setp(l16)
plt.setp(l17)
plt.setp(l18)
plt.setp(l19)
plt.setp(l20)
plt.setp(l21)
plt.setp(l22)
plt.setp(l23)
plt.setp(l24)
plt.setp(l25)
plt.setp(l26)
plt.setp(l27)
plt.setp(l28)
plt.setp(l29)
plt.setp(l30)
plt.setp(l31, linestyle='dashed', color='black')
plt.setp(l32, color='black')

ax3 = fig.add_subplot(323)
plt_clim3 = ax3.plot(time, pre_neb_gcm1, time, pre_neb_gcm2, time, pre_neb_gcm3,
time, pre_neb_gcm4, time, pre_neb_gcm5, time, pre_neb_gcm6, time, pre_neb_gcm7,
time, pre_neb_gcm8, time, pre_neb_gcm9, time, pre_neb_gcm10, time, pre_neb_gcm11,
time, pre_neb_gcm12, time, pre_neb_gcm13, time, pre_neb_gcm14, time, pre_neb_gcm15, 
time, pre_neb_gcm16, time, pre_neb_gcm17, time, pre_neb_gcm18, time, pre_neb_gcm19, 
time, pre_neb_gcm20, time, pre_neb_gcm21, time, pre_neb_gcm22, time, pre_neb_gcm23, 
time, pre_neb_gcm24, time, pre_neb_gcm25, time, pre_neb_gcm26, time, pre_neb_gcm27, 
time, pre_neb_gcm28, time, pre_neb_gcm29, time, pre_neb_gcm30, time, pre_neb_gcm31, 
time, pre_neb_obs)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 21, 3))
plt.ylabel(u'Precipitation (mm d⁻¹)')
plt.title(u'B)', loc='left', fontweight='bold')
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32 = plt_clim3
plt.setp(l1)
plt.setp(l2)
plt.setp(l3)
plt.setp(l4)
plt.setp(l5)
plt.setp(l6)
plt.setp(l7)
plt.setp(l8)
plt.setp(l9)
plt.setp(l10)
plt.setp(l11)
plt.setp(l12)
plt.setp(l13)
plt.setp(l14)
plt.setp(l15)
plt.setp(l16)
plt.setp(l17)
plt.setp(l18)
plt.setp(l19)
plt.setp(l20)
plt.setp(l21)
plt.setp(l22)
plt.setp(l23)
plt.setp(l24)
plt.setp(l25)
plt.setp(l26)
plt.setp(l27)
plt.setp(l28)
plt.setp(l29)
plt.setp(l30)
plt.setp(l31, linestyle='dashed', color='black')
plt.setp(l32, color='black')

ax4 = fig.add_subplot(324)
plt_clim4 = ax4.plot(time, tmp_neb_gcm1, time, tmp_neb_gcm2, time, tmp_neb_gcm3,
time, tmp_neb_gcm4, time, tmp_neb_gcm5, time, tmp_neb_gcm6, time, tmp_neb_gcm7,
time, tmp_neb_gcm8, time, tmp_neb_gcm9, time, tmp_neb_gcm10, time, tmp_neb_gcm11,
time, tmp_neb_gcm12, time, tmp_neb_gcm13, time, tmp_neb_gcm14, time, tmp_neb_gcm15, 
time, tmp_neb_gcm16, time, tmp_neb_gcm17, time, tmp_neb_gcm18, time, tmp_neb_gcm19, 
time, tmp_neb_gcm20, time, tmp_neb_gcm21, time, tmp_neb_gcm22, time, tmp_neb_gcm23, 
time, tmp_neb_gcm24, time, tmp_neb_gcm25, time, tmp_neb_gcm26, time, tmp_neb_gcm27, 
time, tmp_neb_gcm28, time, tmp_neb_gcm29, time, tmp_neb_gcm30, time, tmp_neb_gcm31, 
time, tmp_neb_obs)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(20, 38, 3))
plt.ylabel(u'Temperature (°C)')
plt.title(u'E)', loc='left', fontweight='bold')
ax4.xaxis.grid(True, which='major', linestyle='--')
ax4.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32 = plt_clim4
plt.setp(l1)
plt.setp(l2)
plt.setp(l3)
plt.setp(l4)
plt.setp(l5)
plt.setp(l6)
plt.setp(l7)
plt.setp(l8)
plt.setp(l9)
plt.setp(l10)
plt.setp(l11)
plt.setp(l12)
plt.setp(l13)
plt.setp(l14)
plt.setp(l15)
plt.setp(l16)
plt.setp(l17)
plt.setp(l18)
plt.setp(l19)
plt.setp(l20)
plt.setp(l21)
plt.setp(l22)
plt.setp(l23)
plt.setp(l24)
plt.setp(l25)
plt.setp(l26)
plt.setp(l27)
plt.setp(l28)
plt.setp(l29)
plt.setp(l30)
plt.setp(l31, linestyle='dashed', color='black')
plt.setp(l32, color='black')

ax5 = fig.add_subplot(325)
plt_clim5 = ax5.plot(time, pre_mato_gcm1, time, pre_mato_gcm2, time, pre_mato_gcm3,
time, pre_mato_gcm4, time, pre_mato_gcm5, time, pre_mato_gcm6, time, pre_mato_gcm7,
time, pre_mato_gcm8, time, pre_mato_gcm9, time, pre_mato_gcm10, time, pre_mato_gcm11,
time, pre_mato_gcm12, time, pre_mato_gcm13, time, pre_mato_gcm14, time, pre_mato_gcm15, 
time, pre_mato_gcm16, time, pre_mato_gcm17, time, pre_mato_gcm18, time, pre_mato_gcm19, 
time, pre_mato_gcm20, time, pre_mato_gcm21, time, pre_mato_gcm22, time, pre_mato_gcm23, 
time, pre_mato_gcm24, time, pre_mato_gcm25, time, pre_mato_gcm26, time, pre_mato_gcm27, 
time, pre_mato_gcm28, time, pre_mato_gcm29, time, pre_mato_gcm30, time, pre_mato_gcm31, 
time, pre_mato_obs)
plt.xlabel('Months')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(0, 21, 3))
plt.title(u'C)', loc='left', fontweight='bold')
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32 = plt_clim5
plt.setp(l1)
plt.setp(l2)
plt.setp(l3)
plt.setp(l4)
plt.setp(l5)
plt.setp(l6)
plt.setp(l7)
plt.setp(l8)
plt.setp(l9)
plt.setp(l10)
plt.setp(l11)
plt.setp(l12)
plt.setp(l13)
plt.setp(l14)
plt.setp(l15)
plt.setp(l16)
plt.setp(l17)
plt.setp(l18)
plt.setp(l19)
plt.setp(l20)
plt.setp(l21)
plt.setp(l22)
plt.setp(l23)
plt.setp(l24)
plt.setp(l25)
plt.setp(l26)
plt.setp(l27)
plt.setp(l28)
plt.setp(l29)
plt.setp(l30)
plt.setp(l31, linestyle='dashed', color='black')
plt.setp(l32, color='black')

ax6 = fig.add_subplot(326)
plt_clim6 = ax6.plot(time, tmp_mato_gcm1, time, tmp_mato_gcm2, time, tmp_mato_gcm3,
time, tmp_mato_gcm4, time, tmp_mato_gcm5, time, tmp_mato_gcm6, time, tmp_mato_gcm7,
time, tmp_mato_gcm8, time, tmp_mato_gcm9, time, tmp_mato_gcm10, time, tmp_mato_gcm11,
time, tmp_mato_gcm12, time, tmp_mato_gcm13, time, tmp_mato_gcm14, time, tmp_mato_gcm15, 
time, tmp_mato_gcm16, time, tmp_mato_gcm17, time, tmp_mato_gcm18, time, tmp_mato_gcm19, 
time, tmp_mato_gcm20, time, tmp_mato_gcm21, time, tmp_mato_gcm22, time, tmp_mato_gcm23, 
time, tmp_mato_gcm24, time, tmp_mato_gcm25, time, tmp_mato_gcm26, time, tmp_mato_gcm27, 
time, tmp_mato_gcm28, time, tmp_mato_gcm29, time, tmp_mato_gcm30, time, tmp_mato_gcm31, 
time, tmp_mato_obs)
plt.xlabel('Months')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.yticks(np.arange(20, 38, 3))
plt.title(u'F)', loc='left', fontweight='bold')
ax6.xaxis.grid(True, which='major', linestyle='--')
ax6.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32 = plt_clim6
plt.setp(l1)
plt.setp(l2)
plt.setp(l3)
plt.setp(l4)
plt.setp(l5)
plt.setp(l6)
plt.setp(l7)
plt.setp(l8)
plt.setp(l9)
plt.setp(l10)
plt.setp(l11)
plt.setp(l12)
plt.setp(l13)
plt.setp(l14)
plt.setp(l15)
plt.setp(l16)
plt.setp(l17)
plt.setp(l18)
plt.setp(l19)
plt.setp(l20)
plt.setp(l21)
plt.setp(l22)
plt.setp(l23)
plt.setp(l24)
plt.setp(l25)
plt.setp(l26)
plt.setp(l27)
plt.setp(l28)
plt.setp(l29)
plt.setp(l30)
plt.setp(l31, linestyle='dashed', color='black')
plt.setp(l32, color='black')

legend = ['BCC-CSM1.1','BCC-CSM1.1M','BNU-ESM','CanESM2','CNRM-CM5','CSIRO-ACCESS1.0','CSIRO-ACCESS1.3',
'CSIRO-MK36','FIO-ESM','GISS-E2-H','GISS-E2-H-CC','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','INMCM4',
'IPSL-CM5A-LR','IPSL-CM5A-MR','LASG-FGOALS-G2','LASG-FGOALS-S2','MIROC5','MIROC-ESM','MIROC-ESM-CHEM',
'MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NCAR-CCSM4','NCAR-CESM1-BGC','NCAR-CESM1-CAM5','NorESM1-M',
'NorESM1-ME','ensmean_cmip5','CRU']
plt.legend(plt_clim6, legend, loc=(1.019, -0.19))
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.20, hspace=0.35)

# Save figure
path_out = '/home/nice'
name_out = 'pyplt_clim_cmip5_cru_1975-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()



