# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "01/08/2019"
__description__ = "This script create dictionary monthly with information about the best cmip5 models"

import os
import netCDF4
import numpy as np

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from comp_statist_indices import compute_corr, compute_mae, compute_rmse 
from comp_statist_indices import compute_bias, compute_pbias, compute_apb, compute_effic_coeffic


def import_cmip5_clim(mon, model):
	
	param = 'tas' # pr or tas
	area  = 'amz' # amz or neb
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
	mdl_mon = mdl_data[mon::12]
	
	return mdl_mon


def import_obs_clim(mon, database):
	
	param = 'tmp' # pre or tmp
	area  = 'amz' # amz or neb
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
	obs_mon = obs_data[mon::12]
	
	return obs_mon
	              
	              
monthly_dict = {}
month_dict = {0:'jan', 1:'feb', 2:'mar', 3:'apr', 4:'may', 5:'jun', 6:'jul', 7:'aug', 8:'sep', 9:'oct', 10:'nov', 11:'dec'}

for method in ['r', 'mae', 'rmse', 'bias', 'pbias', 'apb', 'effic']:
	print 'Method:', method

	month_list = []
	for mon in np.arange(0, 11 + 1):
		month = month_dict[mon]
		print 'Month:', month
				  
		# Import cmip5 model end obs database climatology
		model  = u'BCC-CSM1.1'
		mdl1_clim = import_cmip5_clim(mon, model)

		model  = u'BCC-CSM1.1M'
		mdl2_clim = import_cmip5_clim(mon, model)

		model  = u'BNU-ESM'
		mdl3_clim = import_cmip5_clim(mon, model)

		model  = u'CanESM2'
		mdl4_clim = import_cmip5_clim(mon, model)

		model  = u'CNRM-CM5'
		mdl5_clim = import_cmip5_clim(mon, model)

		model  = u'CSIRO-ACCESS-1'
		mdl6_clim = import_cmip5_clim(mon, model)

		model  = u'CSIRO-ACCESS-3'
		mdl7_clim = import_cmip5_clim(mon, model)

		model  = u'CSIRO-MK36'
		mdl8_clim = import_cmip5_clim(mon, model)

		model  = u'FIO-ESM'
		mdl9_clim = import_cmip5_clim(mon, model)

		model  = u'GISS-E2-H-CC'
		mdl10_clim = import_cmip5_clim(mon, model)

		model  = u'GISS-E2-H'
		mdl11_clim = import_cmip5_clim(mon, model)

		model  = u'GISS-E2-R'
		mdl12_clim = import_cmip5_clim(mon, model)

		model  = u'HadGEM2-AO'
		mdl13_clim = import_cmip5_clim(mon, model)

		model  = u'HadGEM2-CC'
		mdl14_clim = import_cmip5_clim(mon, model)
		
		model  = u'HadGEM2-ES'
		mdl15_clim = import_cmip5_clim(mon, model)

		model  = u'INMCM4'
		mdl16_clim = import_cmip5_clim(mon, model)

		model  = u'IPSL-CM5A-LR'
		mdl17_clim = import_cmip5_clim(mon, model)

		model  = u'IPSL-CM5A-MR'
		mdl18_clim = import_cmip5_clim(mon, model)

		model  = u'IPSL-CM5B-LR'
		mdl19_clim = import_cmip5_clim(mon, model)

		model  = u'LASG-FGOALS-G2'
		mdl20_clim = import_cmip5_clim(mon, model)

		model  = u'LASG-FGOALS-S2'
		mdl21_clim = import_cmip5_clim(mon, model)

		model  = u'MIROC5'
		mdl22_clim = import_cmip5_clim(mon, model)

		model  = u'MIROC-ESM-CHEM'
		mdl23_clim = import_cmip5_clim(mon, model)

		model  = u'MIROC-ESM'
		mdl24_clim = import_cmip5_clim(mon, model)

		model  = u'MPI-ESM-LR'
		mdl25_clim = import_cmip5_clim(mon, model)

		model  = u'MPI-ESM-MR'
		mdl26_clim = import_cmip5_clim(mon, model)

		model  = u'MRI-CGCM3'
		mdl27_clim = import_cmip5_clim(mon, model)

		model  = u'NCAR-CCSM4'
		mdl28_clim = import_cmip5_clim(mon, model)

		model  = u'NCAR-CESM1-BGC'
		mdl29_clim = import_cmip5_clim(mon, model)

		model  = u'NCAR-CESM1-CAM5'
		mdl30_clim = import_cmip5_clim(mon, model)

		model  = u'NorESM1-ME'
		mdl31_clim = import_cmip5_clim(mon, model)

		model  = u'NorESM1-M'
		mdl32_clim = import_cmip5_clim(mon, model)

		model  = u'ensmean_cmip5'
		mdl33_clim = import_cmip5_clim(mon, model)

		database  = u'cru_ts4.02'
		obs1_clim = import_obs_clim(mon, database)

		# Compute statiscts index from CMIP5 models
		# Compute Pearson Linear Correlattion Coeficient
		methods_dict = dict(r=[compute_corr(mdl1_clim, obs1_clim),
		compute_corr(mdl2_clim, obs1_clim),
		compute_corr(mdl3_clim, obs1_clim),
		compute_corr(mdl4_clim, obs1_clim),
		compute_corr(mdl5_clim, obs1_clim),
		compute_corr(mdl6_clim, obs1_clim), 
		compute_corr(mdl7_clim, obs1_clim), 
		compute_corr(mdl8_clim, obs1_clim),
		compute_corr(mdl9_clim, obs1_clim), 
		compute_corr(mdl10_clim, obs1_clim), 
		compute_corr(mdl11_clim, obs1_clim),
		compute_corr(mdl12_clim, obs1_clim), 
		compute_corr(mdl13_clim, obs1_clim), 
		compute_corr(mdl14_clim, obs1_clim), 
		compute_corr(mdl15_clim, obs1_clim), 
		compute_corr(mdl16_clim, obs1_clim),
		compute_corr(mdl17_clim, obs1_clim),
		compute_corr(mdl18_clim, obs1_clim), 
		compute_corr(mdl19_clim, obs1_clim),
		compute_corr(mdl20_clim, obs1_clim), 
		compute_corr(mdl21_clim, obs1_clim), 
		compute_corr(mdl22_clim, obs1_clim),
		compute_corr(mdl23_clim, obs1_clim),
		compute_corr(mdl24_clim, obs1_clim),
		compute_corr(mdl25_clim, obs1_clim), 
		compute_corr(mdl26_clim, obs1_clim), 
		compute_corr(mdl27_clim, obs1_clim), 
		compute_corr(mdl28_clim, obs1_clim),
		compute_corr(mdl29_clim, obs1_clim),
		compute_corr(mdl30_clim, obs1_clim),
		compute_corr(mdl31_clim, obs1_clim), 
		compute_corr(mdl32_clim, obs1_clim)],

		mae=[compute_mae(mdl1_clim, obs1_clim),
		compute_mae(mdl2_clim, obs1_clim),
		compute_mae(mdl3_clim, obs1_clim),
		compute_mae(mdl4_clim, obs1_clim),
		compute_mae(mdl5_clim, obs1_clim),
		compute_mae(mdl6_clim, obs1_clim), 
		compute_mae(mdl7_clim, obs1_clim), 
		compute_mae(mdl8_clim, obs1_clim),
		compute_mae(mdl9_clim, obs1_clim), 
		compute_mae(mdl10_clim, obs1_clim), 
		compute_mae(mdl11_clim, obs1_clim),
		compute_mae(mdl12_clim, obs1_clim), 
		compute_mae(mdl13_clim, obs1_clim), 
		compute_mae(mdl14_clim, obs1_clim), 
		compute_mae(mdl15_clim, obs1_clim), 
		compute_mae(mdl16_clim, obs1_clim),
		compute_mae(mdl17_clim, obs1_clim),
		compute_mae(mdl18_clim, obs1_clim), 
		compute_mae(mdl19_clim, obs1_clim),
		compute_mae(mdl20_clim, obs1_clim), 
		compute_mae(mdl21_clim, obs1_clim), 
		compute_mae(mdl22_clim, obs1_clim),
		compute_mae(mdl23_clim, obs1_clim),
		compute_mae(mdl24_clim, obs1_clim),
		compute_mae(mdl25_clim, obs1_clim), 
		compute_mae(mdl26_clim, obs1_clim), 
		compute_mae(mdl27_clim, obs1_clim), 
		compute_mae(mdl28_clim, obs1_clim),
		compute_mae(mdl29_clim, obs1_clim),
		compute_mae(mdl30_clim, obs1_clim),
		compute_mae(mdl31_clim, obs1_clim), 
		compute_mae(mdl32_clim, obs1_clim)],

		rmse=[compute_rmse(mdl1_clim, obs1_clim),
		compute_rmse(mdl2_clim, obs1_clim),
		compute_rmse(mdl3_clim, obs1_clim),
		compute_rmse(mdl4_clim, obs1_clim),
		compute_rmse(mdl5_clim, obs1_clim),
		compute_rmse(mdl6_clim, obs1_clim), 
		compute_rmse(mdl7_clim, obs1_clim), 
		compute_rmse(mdl8_clim, obs1_clim),
		compute_rmse(mdl9_clim, obs1_clim), 
		compute_rmse(mdl10_clim, obs1_clim), 
		compute_rmse(mdl11_clim, obs1_clim),
		compute_rmse(mdl12_clim, obs1_clim), 
		compute_rmse(mdl13_clim, obs1_clim), 
		compute_rmse(mdl14_clim, obs1_clim), 
		compute_rmse(mdl15_clim, obs1_clim), 
		compute_rmse(mdl16_clim, obs1_clim),
		compute_rmse(mdl17_clim, obs1_clim),
		compute_rmse(mdl18_clim, obs1_clim), 
		compute_rmse(mdl19_clim, obs1_clim),
		compute_rmse(mdl20_clim, obs1_clim), 
		compute_rmse(mdl21_clim, obs1_clim), 
		compute_rmse(mdl22_clim, obs1_clim),
		compute_rmse(mdl23_clim, obs1_clim),
		compute_rmse(mdl24_clim, obs1_clim),
		compute_rmse(mdl25_clim, obs1_clim), 
		compute_rmse(mdl26_clim, obs1_clim), 
		compute_rmse(mdl27_clim, obs1_clim), 
		compute_rmse(mdl28_clim, obs1_clim),
		compute_rmse(mdl29_clim, obs1_clim),
		compute_rmse(mdl30_clim, obs1_clim),
		compute_rmse(mdl31_clim, obs1_clim), 
		compute_rmse(mdl32_clim, obs1_clim)],
		
		bias=[compute_bias(mdl1_clim, obs1_clim),
		compute_bias(mdl2_clim, obs1_clim),
		compute_bias(mdl3_clim, obs1_clim),
		compute_bias(mdl4_clim, obs1_clim),
		compute_bias(mdl5_clim, obs1_clim),
		compute_bias(mdl6_clim, obs1_clim), 
		compute_bias(mdl7_clim, obs1_clim), 
		compute_bias(mdl8_clim, obs1_clim),
		compute_bias(mdl9_clim, obs1_clim), 
		compute_bias(mdl10_clim, obs1_clim), 
		compute_bias(mdl11_clim, obs1_clim),
		compute_bias(mdl12_clim, obs1_clim), 
		compute_bias(mdl13_clim, obs1_clim), 
		compute_bias(mdl14_clim, obs1_clim), 
		compute_bias(mdl15_clim, obs1_clim), 
		compute_bias(mdl16_clim, obs1_clim),
		compute_bias(mdl17_clim, obs1_clim),
		compute_bias(mdl18_clim, obs1_clim), 
		compute_bias(mdl19_clim, obs1_clim),
		compute_bias(mdl20_clim, obs1_clim), 
		compute_bias(mdl21_clim, obs1_clim), 
		compute_bias(mdl22_clim, obs1_clim),
		compute_bias(mdl23_clim, obs1_clim),
		compute_bias(mdl24_clim, obs1_clim),
		compute_bias(mdl25_clim, obs1_clim), 
		compute_bias(mdl26_clim, obs1_clim), 
		compute_bias(mdl27_clim, obs1_clim), 
		compute_bias(mdl28_clim, obs1_clim),
		compute_bias(mdl29_clim, obs1_clim),
		compute_bias(mdl30_clim, obs1_clim),
		compute_bias(mdl31_clim, obs1_clim), 
		compute_bias(mdl32_clim, obs1_clim)],
		
		pbias=[compute_pbias(mdl1_clim, obs1_clim),
		compute_pbias(mdl2_clim, obs1_clim),
		compute_pbias(mdl3_clim, obs1_clim),
		compute_pbias(mdl4_clim, obs1_clim),
		compute_pbias(mdl5_clim, obs1_clim),
		compute_pbias(mdl6_clim, obs1_clim), 
		compute_pbias(mdl7_clim, obs1_clim), 
		compute_pbias(mdl8_clim, obs1_clim),
		compute_pbias(mdl9_clim, obs1_clim), 
		compute_pbias(mdl10_clim, obs1_clim), 
		compute_pbias(mdl11_clim, obs1_clim),
		compute_pbias(mdl12_clim, obs1_clim), 
		compute_pbias(mdl13_clim, obs1_clim), 
		compute_pbias(mdl14_clim, obs1_clim), 
		compute_pbias(mdl15_clim, obs1_clim), 
		compute_pbias(mdl16_clim, obs1_clim),
		compute_pbias(mdl17_clim, obs1_clim),
		compute_pbias(mdl18_clim, obs1_clim), 
		compute_pbias(mdl19_clim, obs1_clim),
		compute_pbias(mdl20_clim, obs1_clim), 
		compute_pbias(mdl21_clim, obs1_clim), 
		compute_pbias(mdl22_clim, obs1_clim),
		compute_pbias(mdl23_clim, obs1_clim),
		compute_pbias(mdl24_clim, obs1_clim),
		compute_pbias(mdl25_clim, obs1_clim), 
		compute_pbias(mdl26_clim, obs1_clim), 
		compute_pbias(mdl27_clim, obs1_clim), 
		compute_pbias(mdl28_clim, obs1_clim),
		compute_pbias(mdl29_clim, obs1_clim),
		compute_pbias(mdl30_clim, obs1_clim),
		compute_pbias(mdl31_clim, obs1_clim), 
		compute_pbias(mdl32_clim, obs1_clim)],
		
		apb=[compute_apb(mdl1_clim, obs1_clim),
		compute_apb(mdl2_clim, obs1_clim),
		compute_apb(mdl3_clim, obs1_clim),
		compute_apb(mdl4_clim, obs1_clim),
		compute_apb(mdl5_clim, obs1_clim),
		compute_apb(mdl6_clim, obs1_clim), 
		compute_apb(mdl7_clim, obs1_clim), 
		compute_apb(mdl8_clim, obs1_clim),
		compute_apb(mdl9_clim, obs1_clim), 
		compute_apb(mdl10_clim, obs1_clim), 
		compute_apb(mdl11_clim, obs1_clim),
		compute_apb(mdl12_clim, obs1_clim), 
		compute_apb(mdl13_clim, obs1_clim), 
		compute_apb(mdl14_clim, obs1_clim), 
		compute_apb(mdl15_clim, obs1_clim), 
		compute_apb(mdl16_clim, obs1_clim),
		compute_apb(mdl17_clim, obs1_clim),
		compute_apb(mdl18_clim, obs1_clim), 
		compute_apb(mdl19_clim, obs1_clim),
		compute_apb(mdl20_clim, obs1_clim), 
		compute_apb(mdl21_clim, obs1_clim), 
		compute_apb(mdl22_clim, obs1_clim),
		compute_apb(mdl23_clim, obs1_clim),
		compute_apb(mdl24_clim, obs1_clim),
		compute_apb(mdl25_clim, obs1_clim), 
		compute_apb(mdl26_clim, obs1_clim), 
		compute_apb(mdl27_clim, obs1_clim), 
		compute_apb(mdl28_clim, obs1_clim),
		compute_apb(mdl29_clim, obs1_clim),
		compute_apb(mdl30_clim, obs1_clim),
		compute_apb(mdl31_clim, obs1_clim), 
		compute_apb(mdl32_clim, obs1_clim)],
		
		effic=[compute_effic_coeffic(mdl1_clim, obs1_clim),
		compute_effic_coeffic(mdl2_clim, obs1_clim),
		compute_effic_coeffic(mdl3_clim, obs1_clim),
		compute_effic_coeffic(mdl4_clim, obs1_clim),
		compute_effic_coeffic(mdl5_clim, obs1_clim),
		compute_effic_coeffic(mdl6_clim, obs1_clim), 
		compute_effic_coeffic(mdl7_clim, obs1_clim), 
		compute_effic_coeffic(mdl8_clim, obs1_clim),
		compute_effic_coeffic(mdl9_clim, obs1_clim), 
		compute_effic_coeffic(mdl10_clim, obs1_clim), 
		compute_effic_coeffic(mdl11_clim, obs1_clim),
		compute_effic_coeffic(mdl12_clim, obs1_clim), 
		compute_effic_coeffic(mdl13_clim, obs1_clim), 
		compute_effic_coeffic(mdl14_clim, obs1_clim), 
		compute_effic_coeffic(mdl15_clim, obs1_clim), 
		compute_effic_coeffic(mdl16_clim, obs1_clim),
		compute_effic_coeffic(mdl17_clim, obs1_clim),
		compute_effic_coeffic(mdl18_clim, obs1_clim), 
		compute_effic_coeffic(mdl19_clim, obs1_clim),
		compute_effic_coeffic(mdl20_clim, obs1_clim), 
		compute_effic_coeffic(mdl21_clim, obs1_clim), 
		compute_effic_coeffic(mdl22_clim, obs1_clim),
		compute_effic_coeffic(mdl23_clim, obs1_clim),
		compute_effic_coeffic(mdl24_clim, obs1_clim),
		compute_effic_coeffic(mdl25_clim, obs1_clim), 
		compute_effic_coeffic(mdl26_clim, obs1_clim), 
		compute_effic_coeffic(mdl27_clim, obs1_clim), 
		compute_effic_coeffic(mdl28_clim, obs1_clim),
		compute_effic_coeffic(mdl29_clim, obs1_clim),
		compute_effic_coeffic(mdl30_clim, obs1_clim),
		compute_effic_coeffic(mdl31_clim, obs1_clim), 
		compute_effic_coeffic(mdl32_clim, obs1_clim)])		
		
		mdl_str = ['BCC-CSM1.1','BCC-CSM1.1M','BNU-ESM','CanESM2','CNRM-CM5','CSIRO-ACCESS-1','CSIRO-ACCESS-3','CSIRO-MK36',
		'FIO-ESM','GISS-E2-H-CC','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','INMCM4','IPSL-CM5A-LR',
		'IPSL-CM5A-MR','IPSL-CM5B-LR','LASG-FGOALS-G2','LASG-FGOALS-S2','MIROC','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR',
		'MPI-ESM-MR','MRI-CGCM3','NCAR-CCSM4','NCAR-CESM1-BGC','NCAR-CESM1-CAM5','NorESM1-ME','NorESM1-M','ensmean_cmip5']
			
		statist_list = methods_dict[method]
			
		if method == 'r' or 'effic':
			best_result = mdl_str[np.where(statist_list == np.max(statist_list))[-1][-1]]

		else:
			best_result = mdl_str[np.where(statist_list == np.min(statist_list))[-1][-1]]

		month_list.append(best_result)
		monthly_dict[method] = month_list

	file_name = 'dict_monthly_tmp_amz_best_cmip5_models_1975-2005.py'
	file_save = open(file_name, 'w')
	file_save.write('Dictionary_monthly_tmp_amz_best_method_statist_cmip5_models_1975-2005=' + str(monthly_dict))
	file_save.close()





