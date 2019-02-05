# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "01/21/2019"
__description__ = "This script to extract time series from CMIP5 models"


import netCDF4
import numpy as np
# mpl.use('Agg')


def import_cmip5_clim(model):
	
	param = 'tas'
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/cmip5_hist'
	arq   = '{0}/{1}_amz_neb_Amon_{2}_{3}_{4}.nc'.format(path, param,
	model, exp, date)	
	
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


def import_obs_clim(database):
	
	param = 'tmp'
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/obs_data'
	arq   = '{0}/{1}_amz_neb_{2}_obs_mon_{3}.nc'.format(path,
	param, database, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	obs_data = np.nanmean(np.nanmean(value, axis=1), axis=1)

	obs_clim = []
	obs_clim_p5 = []
	obs_clim_p95 = []

	for mon in range(1, 12 + 1):
		
		obs = np.nanmean(obs_data[mon::12], axis=0)
		obs_clim.append(obs)
	
		obs_p5  = norm.ppf(0.05, loc=np.nanmean(obs_data[mon::12], axis=0), scale=np.nanstd(obs_data[mon::12], axis=0))
		obs_clim_p5.append(obs_p5)
		
		obs_p95  = norm.ppf(0.95, loc=np.nanmean(obs_data[mon::12], axis=0), scale=np.nanstd(obs_data[mon::12], axis=0))
		obs_clim_p95.append(obs_p95)
	
	return obs_clim, obs_clim_p5, obs_clim_p95
	
	
# Import model end obs database climatology
model  = u'BCC-CSM1.1'
mdl1_clim = import_cmip5_clim(model)
		
model  = u'BCC-CSM1.1M'
mdl2_clim = import_cmip5_clim(model)

model  = u'BNU-ESM'
mdl3_clim = import_cmip5_clim(model)

model  = u'CanESM2'
mdl4_clim = import_cmip5_clim(model)

model  = u'CNRM-CM5'
mdl5_clim = import_cmip5_clim(model)

model  = u'CSIRO-ACCESS-1'
mdl6_clim = import_cmip5_clim(model)

model  = u'CSIRO-ACCESS-3'
mdl7_clim = import_cmip5_clim(model)

model  = u'CSIRO-MK36'
mdl8_clim = import_cmip5_clim(model)

model  = u'FIO-ESM'
mdl9_clim = import_cmip5_clim(model)

model  = u'GISS-E2-H-CC'
mdl10_clim = import_cmip5_clim(model)

model  = u'GISS-E2-H'
mdl11_clim = import_cmip5_clim(model)

model  = u'GISS-E2-R'
mdl12_clim = import_cmip5_clim(model)

model  = u'HadGEM2-AO'
mdl13_clim = import_cmip5_clim(model)

model  = u'HadGEM2-CC'
mdl14_clim = import_cmip5_clim(model)

model  = u'INMCM4'
mdl15_clim = import_cmip5_clim(model)

model  = u'IPSL-CM5A-LR'
mdl16_clim = import_cmip5_clim(model)

model  = u'IPSL-CM5A-MR'
mdl17_clim = import_cmip5_clim(model)

model  = u'IPSL-CM5B-LR'
mdl18_clim = import_cmip5_clim(model)

model  = u'LASG-FGOALS-G2'
mdl19_clim = import_cmip5_clim(model)

model  = u'LASG-FGOALS-S2'
mdl20_clim = import_cmip5_clim(model)

model  = u'MIROC5'
mdl21_clim = import_cmip5_clim(model)

model  = u'MIROC-ESM-CHEM'
mdl22_clim = import_cmip5_clim(model)

model  = u'MIROC-ESM'
mdl23_clim = import_cmip5_clim(model)

model  = u'MPI-ESM-LR'
mdl24_clim = import_cmip5_clim(model)

model  = u'MPI-ESM-MR'
mdl25_clim = import_cmip5_clim(model)

model  = u'MRI-CGCM3'
mdl26_clim = import_cmip5_clim(model)

model  = u'NCAR-CCSM4'
mdl27_clim = import_cmip5_clim(model)

model  = u'NCAR-CESM1-BGC'
mdl28_clim = import_cmip5_clim(model)

model  = u'NCAR-CESM1-CAM5'
mdl29_clim = import_cmip5_clim(model)

model  = u'NorESM1-ME'
mdl30_clim = import_cmip5_clim(model)

model  = u'NorESM1-M'
mdl31_clim = import_cmip5_clim(model)

model  = u'HadGEM2-ES'
mdl32_clim = import_cmip5_clim(model)

model  = u'ensmean_cmip5'
mdl33_clim = import_cmip5_clim(model)

database  = u'cru_ts4.02'
obs1_clim, obs1_clim_p5 , obs1_clim_p95  = import_obs_clim(database)


