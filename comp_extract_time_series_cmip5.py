# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "01/21/2019"
__description__ = "This script to extract time series from cmip5 models and obs database"

import netCDF4
import numpy as np
# mpl.use('Agg')


def import_cmip5(model):
	
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

	mdl_month = np.nanmean(np.nanmean(value[1:,:,:], axis=1), axis=1)
	mdl_season = np.nanmean(np.nanmean(value[0:360:3,:,:], axis=1), axis=1)
		
	return mdl_month, mdl_season


def import_obs(database):
	
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

	obs_month = np.nanmean(np.nanmean(value[1:,:,:], axis=1), axis=1)
	obs_season = np.nanmean(np.nanmean(value[0:360:3,:,:], axis=1), axis=1)
		
	return obs_month, obs_season
	
	
# Import cmip5 model and obs database
mdl_list = [u'BCC-CSM1.1', u'BCC-CSM1.1M', u'BNU-ESM', u'CanESM2', u'CNRM-CM5',
u'CSIRO-ACCESS-1', u'CSIRO-ACCESS-3', u'CSIRO-MK36', u'FIO-ESM', u'GISS-E2-H-CC',
u'GISS-E2-H', u'GISS-E2-R', u'HadGEM2-AO', u'HadGEM2-CC', u'INMCM4', u'IPSL-CM5A-LR',
u'IPSL-CM5A-MR', u'IPSL-CM5B-LR', u'LASG-FGOALS-G2', u'LASG-FGOALS-S2', u'MIROC5',
u'MIROC-ESM-CHEM', u'MIROC-ESM', u'MPI-ESM-LR', u'MPI-ESM-MR', u'MRI-CGCM3',
u'NCAR-CCSM4', u'NCAR-CESM1-BGC', u'NCAR-CESM1-CAM5', u'NorESM1-ME', u'NorESM1-M',
u'HadGEM2-ES', u'ensmean_cmip5', u'cru_ts4.02']

for mdl in mdl_list:
	print 'model:', mdl
	
	if mdl == 'cru_ts4.02':
		data_month, data_season  = import_obs(mdl)
	
	else:
		data_month, data_season = import_cmip5(mdl)
		
	data_month  = data_month.reshape(359, 1)
	data_season = data_season.reshape(120, 1)
	
	path_out  = '/home/nice/Documentos/ufrn/PhD_project/datas/file_txt'
	
	print 'Save monthly list'
	filename_month = '{0}/temp2m_amz_neb_{1}_monmean_197512-200511.txt'.format(path_out, mdl)	
	np.savetxt(filename_month, (data_month), fmt="%s")
	
	print 'Save seasonal list'
	filename_season = '{0}/temp2m_amz_neb_{1}_seasmean_197512-200511.txt'.format(path_out, mdl)	
	np.savetxt(filename_season, (data_season), fmt="%s")

exit()
		


