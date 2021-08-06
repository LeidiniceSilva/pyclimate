# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/08/2019"
__description__ = "This script plot scatter plot from cmip5 models and obs database"

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
	
	param = 'tas' # pr or tas
	area  = 'amz' # amz or neb
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/PhD_project/datas/cmip5_hist'
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
	
	param = 'tmp' # pre or tmp
	area  = 'amz' # amz or neb
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/PhD_project/datas/obs_data'
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
	
	# Import cmip5 model and obs database
	mdl_clim = import_cmip5_clim(mdl)
	
	obs  = u'cru_ts4.02'
	obs_clim = import_obs_clim(obs)
	
	# Compute statiscts index from cmip5 models
	x.append(metrics.mean_absolute_error(obs_clim, mdl_clim))
	y.append(metrics.mean_squared_error(obs_clim, mdl_clim))
	z.append(compute_pbias(obs_clim, mdl_clim))

# Plot cmip5 model and obs database
fig = plt.figure()
ax = Axes3D(fig)

a = []
b = []
c = []

colors=plt.matplotlib.cm.Set1(np.linspace(0,1,len(x)))

for i,type in enumerate(mdl_list):
	print(mdl_list[i])
	
	a = (x[i])
	b = (y[i])
	c = (z[i])

	ax.scatter(a, b, c, c=colors[i], s=250,  marker='$%d$' % (i+1))

ax.set_xlabel('EMA', fontsize=12)
ax.set_ylabel('EQM', fontsize=12)
ax.set_zlabel('VIÉS (%)', fontsize=12)

# Choice variable: Rainfall (AMZ and AMZ) or Temperature (AMZ and AMZ) 
out_var    = u'tmp' # pre or tmp
out_area   = u'amz' # amz or neb
area_name  = u'AMZ (Lat:16S 4N, Lon:74W 48W)' # AMZ (Lat:16S 4N, Lon:74W 48W) or NEB (Lat:15S 2N, Lon:46W 34W)

if out_var == 'pre':
	var_name   = u'Precipitação'
else:
	var_name   = u'Temperatura' 
	
plt.title(u'Scatter plot de {0} - {1}  \n CMIP5-hist x CRU-ts4.02 - 1975-2005 (Período de Referência: 1850-2005)'.format(var_name, area_name), fontsize=12, y=0.99)

# Path out to save figure
path_out = '/home/nice'
name_out = 'pyplt_scatter_plot_{0}_{1}_cmip5_cru_1975-2005.png'.format(out_var, out_area)
if not os.path.exists(path_out):
	create_path(path_out)	
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()

