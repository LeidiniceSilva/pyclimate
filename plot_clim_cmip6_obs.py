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


def import_cmip6_clim(param, model):

	path = '/home/nice/Documentos/data_file/cmip_data/cmip6/historical'
	arq  = '{0}/{1}_Amon_{2}_historical_r1i1p1f1_gn_185001-201412.nc'.format(path, param, model)	

	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:] * 86400
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	mdl = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	mdl_clim = []
	for mon in range(1, 12 + 1):
		mdl_data = np.nanmean(mdl[mon::12], axis=0)
		mdl_clim.append(mdl_data)
	
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
	              
               
# Import cmip6 model end obs database climatology
mdl1_clim = import_cmip6_clim('pr', 'BCC-CSM2-MR')
mdl2_clim = import_cmip6_clim('pr', 'BCC-CSM2-MR')
mdl3_clim = import_cmip6_clim('pr', 'BCC-CSM2-MR')
mdl4_clim = import_cmip6_clim('pr', 'BCC-CSM2-MR')
mdl5_clim = import_cmip6_clim('pr', 'BCC-CSM2-MR')
mdl6_clim = import_cmip6_clim('pr', 'BCC-CSM2-MR')
mdl7_clim = import_cmip6_clim('pr', 'BCC-CSM2-MR')
mdl8_clim = import_cmip6_clim('pr', 'BCC-CSM2-MR')
mdl9_clim = import_cmip6_clim('pr', 'BCC-CSM2-MR')
obs1_clim, obs1_clim_p5 , obs1_clim_p95  = import_obs_clim('cru_ts4.02')

# Compute statiscts index from best CMIP5 models and CRU obs database
r2 = metrics.r2_score(mdl1_clim, obs1_clim)

# Plot model end obs data climatology
fig, ax = plt.subplots()
time = np.arange(0.5, 12 + 0.5)

plt_clim = plt.plot(time, mdl1_clim, time, mdl2_clim, time, mdl3_clim, time, 
mdl4_clim, time, mdl5_clim, time, mdl6_clim, time, mdl7_clim, time, mdl8_clim,
time, mdl9_clim, time, obs1_clim, time, obs1_clim_p5, time, obs1_clim_p95)

l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12 = plt_clim

plt.setp(l1,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l2,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l3,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l4,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l5,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l6,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l7,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l8,  linewidth=3, markeredgewidth=1, color='gainsboro')
plt.setp(l9,  linewidth=3, markeredgewidth=1, color='blue')
plt.setp(l10, linewidth=3, markeredgewidth=1, color='black')
plt.setp(l11, linewidth=3, markeredgewidth=1, color='slategray')
plt.setp(l12, linewidth=3, markeredgewidth=1, color='slategray')

plt.fill_between(time, obs1_clim_p5, obs1_clim_p95, facecolor='slategray', alpha=0.8, interpolate=True)

plt.title(u'CMIP6-hist x CRU-ts4.02 (Reference period: 1850-2005)')
plt.xlabel(u'Meses')
plt.ylabel(u'Precipitação (mm)')
plt.xticks(time, ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
plt.yticks(np.arange(0, 14, 2))
plt.tick_params(axis='both', which='major', labelcolor='black')

legend = (u'ENSMEAN_CMIP5', u'BNU-ESM', u'CanESM2', u'CSIRO-ACCESS-1', u'HadGEM-ES', 
		  u'IPSL-CM5A-LR', u'MIROC5', u'MPI-ESM-LR', u'NorESM1-M', u'CRU', u'CRU_p5%', u'CRU_p95%')    
plt.legend(plt_clim, legend, loc='upper center', bbox_to_anchor=(0.5, -0.07), shadow=True, ncol=6)
ax.xaxis.grid(True, which='major', linestyle='--', linewidth='1.4', zorder=0.6)
ax.yaxis.grid(True, which='major', linestyle='--', linewidth='1.4', zorder=0.6)
    
# Path out to save figure
path_out = '/home/nice/Documentos/ufrn/papers/cmip6/results'
name_out = 'pyplt_clim_pr_cmip6_obs_1979-2014.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

