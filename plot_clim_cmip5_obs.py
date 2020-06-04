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
	              
               
# Import cmip5 model end obs database climatology
model  = u'BCC-CSM1.1'
mdl1_clim = import_cmip5_clim(model)
		
model  = u'BCC-CSM1.1M'
mdl2_clim = import_cmip5_clim(model)

model  = u'CNRM-CM5'
mdl3_clim = import_cmip5_clim(model)

model  = u'CSIRO-ACCESS-3'
mdl4_clim = import_cmip5_clim(model)

model  = u'CSIRO-MK36'
mdl5_clim = import_cmip5_clim(model)

model  = u'FIO-ESM'
mdl6_clim = import_cmip5_clim(model)

model  = u'GISS-E2-H-CC'
mdl7_clim = import_cmip5_clim(model)

model  = u'GISS-E2-H'
mdl8_clim = import_cmip5_clim(model)

model  = u'GISS-E2-R'
mdl9_clim = import_cmip5_clim(model)

model  = u'HadGEM2-AO'
mdl10_clim = import_cmip5_clim(model)

model  = u'HadGEM2-CC'
mdl11_clim = import_cmip5_clim(model)

model  = u'INMCM4'
mdl12_clim = import_cmip5_clim(model)

model  = u'IPSL-CM5A-LR'
mdl13_clim = import_cmip5_clim(model)

model  = u'IPSL-CM5B-LR'
mdl14_clim = import_cmip5_clim(model)

model  = u'LASG-FGOALS-G2'
mdl15_clim = import_cmip5_clim(model)

model  = u'LASG-FGOALS-S2'
mdl16_clim = import_cmip5_clim(model)

model  = u'MIROC-ESM-CHEM'
mdl17_clim = import_cmip5_clim(model)

model  = u'MIROC-ESM'
mdl18_clim = import_cmip5_clim(model)

model  = u'MPI-ESM-LR'
mdl19_clim = import_cmip5_clim(model)

model  = u'MRI-CGCM3'
mdl20_clim = import_cmip5_clim(model)

model  = u'NCAR-CCSM4'
mdl21_clim = import_cmip5_clim(model)

model  = u'NCAR-CESM1-BGC'
mdl22_clim = import_cmip5_clim(model)

model  = u'NCAR-CESM1-CAM5'
mdl23_clim = import_cmip5_clim(model)

model  = u'NorESM1-ME'
mdl24_clim = import_cmip5_clim(model)

model  = u'ensmean_cmip5'
mdl25_clim = import_cmip5_clim(model)

model  = u'BNU-ESM'
mdl26_clim = import_cmip5_clim(model)

model  = u'CanESM2'
mdl27_clim = import_cmip5_clim(model)

model  = u'CSIRO-ACCESS-1'
mdl28_clim = import_cmip5_clim(model)

model  = u'HadGEM2-ES'
mdl29_clim = import_cmip5_clim(model)

model  = u'IPSL-CM5A-MR'
mdl30_clim = import_cmip5_clim(model)

model  = u'MIROC5'
mdl31_clim = import_cmip5_clim(model)

model  = u'MPI-ESM-MR'
mdl32_clim = import_cmip5_clim(model)

model  = u'NorESM1-M'
mdl33_clim = import_cmip5_clim(model)

database  = u'cru_ts4.02'
obs1_clim, obs1_clim_p5 , obs1_clim_p95  = import_obs_clim(database)

# Plot model end obs data climatology
fig, ax = plt.subplots(figsize=(12, 6))
time = np.arange(0.5, 12 + 0.5)

plt_clim = plt.plot(time, mdl1_clim, time, mdl2_clim, time, mdl3_clim,
time,  mdl4_clim, time, mdl5_clim, time, mdl6_clim, time, mdl7_clim,
time, mdl8_clim, time, mdl9_clim, time, mdl10_clim, time, mdl11_clim,
time, mdl12_clim, time, mdl13_clim, time, mdl14_clim, time, mdl15_clim, 
time, mdl16_clim, time, mdl17_clim, time, mdl18_clim, time, mdl19_clim, 
time, mdl20_clim, time, mdl21_clim, time, mdl22_clim, time, mdl23_clim, 
time, mdl24_clim, time, mdl25_clim, time, mdl26_clim, time, mdl27_clim, 
time, mdl28_clim, time, mdl29_clim, time, mdl30_clim, time, mdl31_clim, 
time, mdl32_clim, time, mdl33_clim, time, obs1_clim, time, obs1_clim_p5, 
time, obs1_clim_p95)

l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35, l36 = plt_clim

plt.setp(l1,  linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l2,  linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l3,  linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l4,  linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l5,  linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l6,  linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l7,  linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l8,  linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l9,  linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l10, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l11, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l12, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l13, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l14, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l15, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l16, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l17, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l18, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l19, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l20, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l21, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l22, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l23, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l24, linewidth=2, markeredgewidth=1, color='gainsboro')
plt.setp(l25, linewidth=3, markeredgewidth=1, color='darkgray')
plt.setp(l26, linewidth=3, markeredgewidth=1, color='green')
plt.setp(l27, linewidth=3, markeredgewidth=1, color='orange')
plt.setp(l28, linewidth=3, markeredgewidth=1, color='pink')
plt.setp(l29, linewidth=3, markeredgewidth=1, color='purple')
plt.setp(l30, linewidth=3, markeredgewidth=1, color='brown')
plt.setp(l31, linewidth=3, markeredgewidth=1, color='yellow')
plt.setp(l32, linewidth=3, markeredgewidth=1, color='blue')
plt.setp(l33, linewidth=3, markeredgewidth=1, color='red')
plt.setp(l34, linewidth=3, markeredgewidth=1, color='black')
plt.setp(l35, linewidth=2, markeredgewidth=3, color='slategray')
plt.setp(l36, linewidth=2, markeredgewidth=3, color='slategray')

plt.fill_between(time, obs1_clim_p5, obs1_clim_p95, facecolor='slategray', alpha=0.8, interpolate=True)

# choice variable: Rainfall (AMZ and AMZ) or Temperature (AMZ and AMZ) 
out_var    = u'tmp' # pre or tmp
out_area   = u'amz' # amz or neb
area_name  = u'AMZ (Lat:16S 4N, Lon:74W 48W)' # AMZ (Lat:16S 4N, Lon:74W 48W) or NEB (Lat:15S 2N, Lon:46W 34W)

if out_var == 'pre':
	yaxis = np.arange(0, 14, 2)
	var_name   = u'Precipitação'
	label_name = u'Precipitação (mm/d)' 
	plt.text(9.5, 10.5, u'r2 = {0}'.format(round(r2, 3)), fontsize=12, fontweight='bold')

else:
	yaxis = np.arange(18, 34, 2)
	var_name   = u'Temperatura' 
	label_name = u'Temperatura 2m ($^\circ$C)' 
	plt.text(9.5, 30.5, u'r2 = {0}'.format(round(r2, 3)), fontsize=12, fontweight='bold')

fig.suptitle(u'Ciclo Anual de {0} - {1} \n CMIP5-hist x CRU-ts4.02 - 1975-2005 (Período de Referência: 1850-2005)'.format(var_name, area_name), fontsize=15, y=0.98)

xaxis = ['Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez']

plt.ylabel(u'{0}'.format(label_name), fontsize=12)
plt.xticks(time, xaxis, fontsize=12)
plt.yticks(yaxis, fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=12, length=5, width=2, pad=4, labelcolor='black')

legend = (u'CMIP5-hist', u'ENSMEAN_CMIP5', u'BNU-ESM', u'CanESM2', u'CSIRO-ACCESS-1', u'HadGEM-ES', 
		  u'IPSL-CM5A-MR', u'MIROC5', u'MPI-ESM-MR', u'NorESM1-M', u'CRU', u'CRU_p5%', u'CRU_p95%')    
plt.legend(plt_clim[23:], legend, loc='upper center', bbox_to_anchor=(0.5, -0.07), shadow=True, ncol=7, prop=FontProperties(size=10))
ax.xaxis.grid(True, which='major', linestyle='--', linewidth='1.4', zorder=0.6)
ax.yaxis.grid(True, which='major', linestyle='--', linewidth='1.4', zorder=0.6)
    
path_out = '/home/nice'
name_out = 'pyplt_clim_{0}_{1}_cmip5_cru_1975-2005.png'.format(out_var, out_area)

if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

