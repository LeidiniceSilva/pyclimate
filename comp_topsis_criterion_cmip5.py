# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "06/01/2020"
__description__ = "This script compute criterion TOPSIS from CMIP5 models end OBS basedata"

import os
import netCDF4
import numpy as np
import texttable as tt

from comp_statist_indices import compute_corr
from comp_statist_indices import compute_rmse
from comp_statist_indices import compute_mae
from comp_statist_indices import compute_r2

def import_cmip5(model):
	
	param = 'tas' # pr or tas
	area  = 'amz' # amz or neb
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/phd_project/datas/cmip5_hist'
	arq   = '{0}/{1}_{2}_Amon_{3}_{4}_{5}.nc'.format(path, param, area,	model, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	month_sim = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	return month_sim


def import_obs(database):
	
	param = 'tmp' # pre or tmp
	area  = 'amz' # amz or neb
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/phd_project/datas/obs_data'
	arq   = '{0}/{1}_{2}_{3}_obs_mon_{4}.nc'.format(path, param, area, database, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	month_obs = np.nanmean(np.nanmean(value, axis=1), axis=1) 

	return month_obs
	

models = ['BCC-CSM1.1','BCC-CSM1.1M','BNU-ESM','CanESM2','CNRM-CM5','CSIRO-ACCESS-1','CSIRO-ACCESS-3','CSIRO-MK36',
'FIO-ESM','GISS-E2-H-CC','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','INMCM4','IPSL-CM5A-LR',
'IPSL-CM5A-MR','IPSL-CM5B-LR','LASG-FGOALS-G2','LASG-FGOALS-S2','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR',
'MPI-ESM-MR','MRI-CGCM3','NCAR-CCSM4','NCAR-CESM1-BGC','NCAR-CESM1-CAM5','NorESM1-ME','NorESM1-M','ensmean_cmip5']

metrics_gcm = []

for mdl in models:
	
	# Import cmip5 model end obs database monthly	
	month_gcm = import_cmip5(mdl)
		
	obs  = u'cru_ts4.02'
	month_cru = import_obs(obs)
	
	corr = compute_corr(month_gcm, month_cru)
	r2 = compute_r2(month_gcm, month_cru)
	rmse = compute_rmse(month_gcm, month_cru)
	mae = compute_mae(month_gcm, month_cru)
	
	metrics = [mdl, round((corr),3), round((r2),3), round((rmse),3), round((mae),3)]
	metrics_gcm.append(metrics)

print(metrics_gcm)
exit()

