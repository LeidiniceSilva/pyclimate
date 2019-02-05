# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "01/08/2019"
__description__ = "This script plot boxplot and time series from CMIP5 models"

import os
import netCDF4
import numpy as np

from pylab import *
from netCDF4 import Dataset
from scipy.stats import norm
from matplotlib.font_manager import FontProperties



def import_cmip5_clim(model):
	
	param = 'pr'
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documentos/cmip5_hist'
	arq   = '{0}/{1}_amz_neb_Amon_{2}_{3}_{4}_seamask.nc'.format(path,
	param, model, exp, date)	
	
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
	
	param = 'pre'
	date  = '197512-200511'

	path  = '/home/nice/Documentos/obs_data'
	arq   = '{0}/{1}_amz_neb_{2}_obs_mon_197512-200511.nc'.format(path,
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

# Plot model end obs data climatology
fig   = plt.figure(figsize=(26, 18))
time = np.arange(1, 12 + 1)

plt_l = plt.plot(time, mdl1_clim, time, mdl2_clim, time, mdl3_clim,
time,  mdl4_clim, time, mdl5_clim, time, mdl6_clim, time, mdl7_clim,
time, mdl8_clim, time, mdl9_clim, time, mdl10_clim, time, mdl11_clim,
time, mdl12_clim, time, mdl13_clim, time, mdl14_clim, time, mdl15_clim, 
time, mdl16_clim, time, mdl17_clim, time, mdl18_clim, time, mdl19_clim, 
time, mdl20_clim, time, mdl21_clim, time, mdl22_clim, time, mdl23_clim, 
time, mdl24_clim, time, mdl25_clim, time, mdl26_clim, time, mdl27_clim, 
time, mdl28_clim, time, mdl29_clim, time, mdl30_clim, time, mdl31_clim, 
time, mdl32_clim, time, mdl33_clim, time, obs1_clim, time, obs1_clim_p5, 
time, obs1_clim_p95)

l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35, l36 = plt_l

plt.setp(l1,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l2,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l3,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l4,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l5,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l6,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l7,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l8,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l9,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l10, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l11, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l12, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l13, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l14, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l15, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l16, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l17, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l18, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l19, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l20, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l21, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l22, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l23, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l24, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l25, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l26, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l27, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l28, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l29, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l30, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l31, linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l32, linewidth=8, markeredgewidth=1, color='gray')
plt.setp(l33, linewidth=8, markeredgewidth=1, color='black')
plt.setp(l34, linewidth=8, markeredgewidth=1, linestyle='--', color='black')
plt.setp(l35, linewidth=8, markeredgewidth=1, linestyle='--', color='lightgray')
plt.setp(l36, linewidth=8, markeredgewidth=1, linestyle='--', color='dimgray')

plt.fill_between(time, obs1_clim_p5, obs1_clim_p95, facecolor='slategray', alpha=0.2, interpolate=True)

fig.suptitle(u'Climatologia de precipitação (1975-2005)', fontsize=30, fontweight='bold', y=0.92)

plt.xlabel(u'Meses', fontsize=30)
plt.ylabel(u'Precipitação (mm/d)', fontsize=30)
xaxis = ['Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez']
plt.xticks(time, xaxis, fontsize=30)
plt.yticks(np.arange(0, 14, 2), fontsize=30)
plt.tick_params(axis='both', which='major', labelsize=30, length=8, width=2, pad=4, labelcolor='black')

legend = (u'CMIP5', u'HadGEM-ES', u'ENSMEAN_CMIP5', u'CRU', u'CRU_p5%', u'CRU_p95%')
font = FontProperties(size=25)		    
plt.legend(plt_l[30:], legend, loc='upper center', bbox_to_anchor=(0.5, -0.07), shadow=True, ncol=6, prop=font)
plt.grid(True, which='major', linestyle='--', linewidth='1.4', zorder=0.6)
			    
path_out = '/home/nice/Documentos/cmip5_hist'
name_out = 'clim_cmip5_obs.png'

if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
plt.show()
exit()

