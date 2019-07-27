# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/14/2019"
__description__ = "This script plot boxplot and from CMIP5 models end OBS basedata"

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from scipy.stats import norm
from matplotlib.font_manager import FontProperties


def import_cmip5_clim(model):
	
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


def import_obs_clim(database):
	
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
	
	
# Import cmip5 model end obs database climatology
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

model  = u'HadGEM2-ES'
mdl15_clim = import_cmip5_clim(model)

model  = u'INMCM4'
mdl16_clim = import_cmip5_clim(model)

model  = u'IPSL-CM5A-LR'
mdl17_clim = import_cmip5_clim(model)

model  = u'IPSL-CM5A-MR'
mdl18_clim = import_cmip5_clim(model)

model  = u'IPSL-CM5B-LR'
mdl19_clim = import_cmip5_clim(model)

model  = u'LASG-FGOALS-G2'
mdl20_clim = import_cmip5_clim(model)

model  = u'LASG-FGOALS-S2'
mdl21_clim = import_cmip5_clim(model)

model  = u'MIROC5'
mdl22_clim = import_cmip5_clim(model)

model  = u'MIROC-ESM-CHEM'
mdl23_clim = import_cmip5_clim(model)

model  = u'MIROC-ESM'
mdl24_clim = import_cmip5_clim(model)

model  = u'MPI-ESM-LR'
mdl25_clim = import_cmip5_clim(model)

model  = u'MPI-ESM-MR'
mdl26_clim = import_cmip5_clim(model)

model  = u'MRI-CGCM3'
mdl27_clim = import_cmip5_clim(model)

model  = u'NCAR-CCSM4'
mdl28_clim = import_cmip5_clim(model)

model  = u'NCAR-CESM1-BGC'
mdl29_clim = import_cmip5_clim(model)

model  = u'NCAR-CESM1-CAM5'
mdl30_clim = import_cmip5_clim(model)

model  = u'NorESM1-ME'
mdl31_clim = import_cmip5_clim(model)

model  = u'NorESM1-M'
mdl32_clim = import_cmip5_clim(model)

model  = u'ensmean_cmip5'
mdl33_clim = import_cmip5_clim(model)

database  = u'cru_ts4.02'
obs1_clim = import_obs_clim(database)

# Plot model end obs data boxplot
fig, ax = plt.subplots(figsize=(28, 16))
time = np.arange(1, 35)

data = [mdl1_clim, mdl2_clim, mdl3_clim, mdl4_clim, mdl5_clim, mdl6_clim, mdl7_clim, mdl8_clim, mdl9_clim,
mdl10_clim, mdl11_clim, mdl12_clim, mdl13_clim, mdl14_clim, mdl15_clim, mdl16_clim, mdl17_clim, mdl18_clim, 
mdl19_clim, mdl20_clim, mdl21_clim, mdl22_clim, mdl23_clim, mdl24_clim, mdl25_clim, mdl26_clim, mdl27_clim, 
mdl28_clim, mdl29_clim, mdl30_clim, mdl31_clim, mdl32_clim, mdl33_clim, obs1_clim]

plt_bp = plt.boxplot(data, patch_artist=True, notch=True, bootstrap=10000, vert=1)

# Change outline and fill color
for box in plt_bp['boxes']:
    box.set( color='black', linewidth=2)
    box.set( facecolor = 'gray' )

# Change color and linewidth of the whiskers
for whisker in plt_bp['whiskers']:
    whisker.set(color='black', linewidth=2)

# Change color and linewidth of the caps
for cap in plt_bp['caps']:
    cap.set(color='black', linewidth=2)

# Change color and linewidth of the medians
for median in plt_bp['medians']:
    median.set(color='blue', linewidth=2)

# Change the style of fliers and their fill
for flier in plt_bp['fliers']:
    flier.set(marker='o', color='blue', alpha=1)
   
# choice ariable: Rainfall or Temperature
out_var    = u'pre' # pre or tmp
out_area   = u'neb' # amz or neb
area_name  = u'NEB (Lat:15S 2N, Lon:46W 34W)' # AMZ (Lat:16S 4N, Lon:74W 48W) or NEB (Lat:15S 2N, Lon:46W 34W)

if out_var == 'pre':
	yaxis = np.arange(0, 20, 2)
	var_name   = u'Rainfall'
	label_name = u'Rain (mm/day)' 

else:
	yaxis = np.arange(18, 34, 2)
	var_name   = u'Temperature' 
	label_name = u'TEmperature 2m ($^\circ$C)'
	
fig.suptitle((u'{0} Boxplot - {1} \n CMIP5-hist x CRU-ts4.02 - 1975-2005 (Reference period: 1850-2005)'.format(var_name, area_name)), fontsize=30, y=0.98)

xaxis = ['BCC-CSM1.1','BCC-CSM1.1M','BNU-ESM','CanESM2','CNRM-CM5','CSIRO-ACCESS-1','CSIRO-ACCESS-3','CSIRO-MK36',
'FIO-ESM','GISS-E2-H-CC','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','INMCM4','IPSL-CM5A-LR',
'IPSL-CM5A-MR','IPSL-CM5B-LR','LASG-FGOALS-G2','LASG-FGOALS-S2','MIROC','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR',
'MPI-ESM-MR','MRI-CGCM3','NCAR-CCSM4','NCAR-CESM1-BGC','NCAR-CESM1-CAM5','NorESM1-ME','NorESM1-M','ensmean_cmip5','cru_ts4.02']

plt.xlabel(u'CMIP5-hist and CRU-ts4.02', fontsize=30)
plt.ylabel(u'{}'.format(label_name), fontsize=30) 
plt.xticks(time, xaxis, rotation=90, fontsize=10)
plt.yticks(yaxis, fontsize=30)
plt.tick_params(axis='both', which='major', length=10, width=4, pad=8, labelcolor='black')
ax.yaxis.grid(True, which='major', linestyle='--', linewidth='1.4', zorder=0.6)

path_out = '/home/nice/Documentos/ufrn/PhD_project/results/cmip5'
name_out = 'pyplt_boxplot_{0}_{1}_cmip5_cru_1975-2005.png'.format(out_var, out_area)

if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
plt.show()
exit()
