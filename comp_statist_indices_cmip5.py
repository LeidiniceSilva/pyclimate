# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "02/15/2019"
__description__ = "This script plot compute statistic indices from cmip5 models and obs database"

import os
import netCDF4
import numpy as np
import texttable as tt
from sklearn import metrics

from comp_statist_indices import compute_corr, compute_mae, compute_rmse 
from comp_statist_indices import compute_bias, compute_pbias, compute_apb, compute_effic_coeffic


def import_cmip5(model):
	
	param = 'pr' # pr or tas
	area  = 'neb' # amz or neb
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/cmip5_hist'
	arq   = '{0}/{1}_{2}_Amon_{3}_{4}_{5}.nc'.format(path, param, area,
	model, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mdl_data = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return mdl_data


def import_obs(database):
	
	param = 'pre' # pre or tmp
	area  = 'neb' # amz or neb
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/obs_data'
	arq   = '{0}/{1}_{2}_{3}_obs_mon_{4}.nc'.format(path, param, area, 
	database, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	obs_data = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return obs_data


tab = tt.Texttable(max_width=200)
tab_inform = [[]]
    
mdl_list = ['BCC-CSM1.1','BCC-CSM1.1M','BNU-ESM','CanESM2','CNRM-CM5','CSIRO-ACCESS-1','CSIRO-ACCESS-3','CSIRO-MK36',
'FIO-ESM','GISS-E2-H-CC','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','INMCM4','IPSL-CM5A-LR',
'IPSL-CM5A-MR','IPSL-CM5B-LR','LASG-FGOALS-G2','LASG-FGOALS-S2','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR',
'MPI-ESM-MR','MRI-CGCM3','NCAR-CCSM4','NCAR-CESM1-BGC','NCAR-CESM1-CAM5','NorESM1-ME','NorESM1-M','ensmean_cmip5']
		
for mdl in mdl_list:
	print 'CMIP5 Model:', mdl
	
	# Import cmip5 model and obs database
	mdl_clim = import_cmip5(mdl)
	
	obs  = u'cru_ts4.02'
	obs_clim = import_obs(obs)
	
	# Compute statiscts index from cmip5 model and obs database
	r     = round((np.corrcoef(np.array(mdl_clim), np.array(obs_clim)))[0][1], 3)
	r2    = metrics.r2_score(obs_clim, mdl_clim)
	mae   = metrics.mean_absolute_error(obs_clim, mdl_clim)
	mse   = metrics.mean_squared_error(obs_clim, mdl_clim)
	rmse  = compute_rmse(obs_clim, mdl_clim)
	bias  = compute_bias(mdl_clim, obs_clim)
	pbias = compute_pbias(mdl_clim, obs_clim)
	apb   = compute_apb(mdl_clim, obs_clim)
	effic = compute_effic_coeffic(mdl_clim, obs_clim)
	mae = metrics.median_absolute_error(obs_clim, mdl_clim)
		
	print r, r2, mae, mse, rmse, bias, pbias, apb, effic, Mae
	
	tab_inform.append([mdl, r, r2, mae, mse, rmse, bias, pbias, apb, effic, mae])

tab.add_rows(tab_inform)
tab.set_cols_align(['c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c'])
tab.header([u'CMIP5 Models', u'R', u'R2_score', u'MAE', u'MSE', U'RMSE', U'BIAS', U'PBIAS', U'APB', U'NASH', u'MAE'])
table = str(tab.draw())
	
file_name = 'table_monthly_pre_neb_statist_indices_cmip5_models_1975-2005.asc'
file_save = open(file_name, 'w')
file_save.write(table)
file_save.close()        
exit()
	

		







