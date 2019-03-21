# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "01/08/2019"
__description__ = "This script plot climatology graphics from CMIP5 models end OBS basedata"

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from sklearn import metrics
from comp_statist_indices import compute_pbias
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.font_manager import FontProperties


def import_cmip5_clim(model):
	
	param = 'pr' # pr or tas
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

	mdl_clim = []
	for mon in range(1, 12 + 1):
		mdl = np.nanmean(mdl_data[mon::12], axis=0)
		mdl_clim.append(mdl)
	
	return mdl_clim


def import_obs_clim(database):
	
	param = 'pre' # pre or tmp
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

	obs_clim = []
	for mon in range(1, 12 + 1):
		
		obs = np.nanmean(obs_data[mon::12], axis=0)
		obs_clim.append(obs)
	
	return obs_clim
	              
x = []
y = []
z = []

mdl_list = ['BCC-CSM1.1','BCC-CSM1.1M','BNU-ESM','CanESM2','CNRM-CM5','CSIRO-ACCESS-1','CSIRO-ACCESS-3','CSIRO-MK36',
'FIO-ESM','GISS-E2-H-CC','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','INMCM4','IPSL-CM5A-LR',
'IPSL-CM5A-MR','IPSL-CM5B-LR','LASG-FGOALS-G2','LASG-FGOALS-S2','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR',
'MPI-ESM-MR','MRI-CGCM3','NCAR-CCSM4','NCAR-CESM1-BGC','NCAR-CESM1-CAM5','NorESM1-ME','NorESM1-M','ensmean_cmip5']

for mdl in mdl_list:
	print 'CMIP5 Model:', mdl
	
	# Import cmip5 model end obs database monthly	
	mdl_clim = import_cmip5_clim(mdl)
	
	obs  = u'cru_ts4.02'
	obs_clim = import_obs_clim(obs)
	
	# Compute statiscts index from CMIP5 models
	x.append(metrics.mean_absolute_error(obs_clim, mdl_clim))
	y.append(metrics.mean_squared_error(obs_clim, mdl_clim))
	z.append(compute_pbias(obs_clim, mdl_clim))

fig = pyplot.figure()
ax = Axes3D(fig)

types = ['BCC-CSM1.1','BCC-CSM1.1M','BNU-ESM','CanESM2','CNRM-CM5','CSIRO-ACCESS-1','CSIRO-ACCESS-3','CSIRO-MK36',
'FIO-ESM','GISS-E2-H-CC','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','INMCM4','IPSL-CM5A-LR',
'IPSL-CM5A-MR','IPSL-CM5B-LR','LASG-FGOALS-G2','LASG-FGOALS-S2','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR',
'MPI-ESM-MR','MRI-CGCM3','NCAR-CCSM4','NCAR-CESM1-BGC','NCAR-CESM1-CAM5','NorESM1-ME','NorESM1-M','ensmean_cmip5']

for i,type in enumerate(types):
    
    a = x[i]
    b = y[i]
    c = z[i]
    
    ax.scatter(a, b, c, c='b', s=30, marker='o')
    ax.set_xlabel('MAE')
    ax.set_ylabel('MSE')
    ax.set_zlabel('PBIAS')
    
    plt.title('3D Scatter plot from CMIP5 models')

plt.show()
exit()

