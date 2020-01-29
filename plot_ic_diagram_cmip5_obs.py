# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "03/25/2019"
__description__ = "This script plot correlation Portrait Diagram from CMIP5 models"

import os
import netCDF4
import numpy as np
import matplotlib
import seaborn as sns 
import matplotlib.pyplot as plt

from comp_statist_indices import compute_corr
from comp_statist_indices import compute_index_agreement


def import_cmip5(model):
	
	param = 'tas' # pr or tas
	area  = 'neb' # amz or neb
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/PhD_project/datas/cmip5_hist'
	arq   = '{0}/{1}_{2}_Amon_{3}_{4}_{5}.nc'.format(path, param, area,	model, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_mdl = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	return annual_mdl


def import_obs(database):
	
	param = 'tmp' # pre or tmp
	area  = 'amz' # amz or neb
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/PhD_project/datas/obs_data'
	arq   = '{0}/{1}_{2}_{3}_obs_mon_{4}.nc'.format(path, param, area, database, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_obs = np.nanmean(np.nanmean(value, axis=1), axis=1) 

	return annual_obs
	

models = ['BCC-CSM1.1','BCC-CSM1.1M','BNU-ESM','CanESM2','CNRM-CM5','CSIRO-ACCESS-1','CSIRO-ACCESS-3','CSIRO-MK36',
'FIO-ESM','GISS-E2-H-CC','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','INMCM4','IPSL-CM5A-LR',
'IPSL-CM5A-MR','IPSL-CM5B-LR','LASG-FGOALS-G2','LASG-FGOALS-S2','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR',
'MPI-ESM-MR','MRI-CGCM3','NCAR-CCSM4','NCAR-CESM1-BGC','NCAR-CESM1-CAM5','NorESM1-ME','NorESM1-M','ensmean_cmip5']


mdl_ic = []

for mdl in models:
	
	# Import cmip5 model end obs database monthly	
	annual_sim = import_cmip5(mdl)
		
	obs  = u'cru_ts4.02'
	annual_cru = import_obs(obs)
	
	icw = compute_index_agreement(annual_sim, annual_cru)
	
	corr = compute_corr(annual_sim, annual_cru)
	
	ic = corr * icw
	
	mdl_ic.append(ic)
	

# Plot cmip5 model end obs database       

fig, ax = plt.subplots()

# Create horizontal bars
plt.barh(np.arange(len(models)), mdl_ic, height=0.7, align='center', alpha=0.8, color='gray', edgecolor='black')
 
# Choice variable: Rainfall (AMZ and AMZ) or Temperature (AMZ and AMZ) 
out_var    = u'tmp' # pre or tmp
out_area   = u'neb' # amz or neb
area_name  = u'NEB (Lat:15S 2N, Lon:46W 34W)' # AMZ (Lat:16S 4N, Lon:74W 48W) or NEB (Lat:15S 2N, Lon:46W 34W)

if out_var == 'pre':
	var_name   = u'Precipitação'
else:
	var_name   = u'Temperatura' 
	
plt.title(u'Índice de Confiança de {0} - {1}  \n CMIP5-hist x CRU-ts4.02 - 1975-2005 (Período de Referência: 1850-2005)'.format(var_name, area_name), fontsize=10, x=0.45, y=0.99)
plt.xticks(np.arange(0, 1.1, 0.1), fontsize=8)
plt.yticks(np.arange(len(models)), models, fontsize=8)
plt.text(0.15, 33, u'Péssimo', fontsize=8)
plt.text(0.41, 33, u'Ruim', fontsize=8)
plt.text(0.52, 33, u'Tolerável', fontsize=8)
plt.text(0.67, 33, u'Bom', fontsize=8)
plt.text(0.86, 33, u'Excelente', fontsize=8)

plt.axvline(0.4, linewidth=1.5, linestyle='dashed', color='gray', zorder=0.5)
plt.axvline(0.51, linewidth=1.5, linestyle='dashed', color='gray', zorder=0.5)
plt.axvline(0.66, linewidth=1.5, linestyle='dashed', color='gray', zorder=0.5)
plt.axvline(0.85, linewidth=1.5, linestyle='dashed', color='gray', zorder=0.5)
plt.axhline(32.5, linewidth=1.5, linestyle='dashed', color='gray', zorder=0.5)

# Save figure
path_out = '/home/nice'
name_out = 'pyplt_ic_diagram_{0}_{1}_cmip5_cru_1975-2005.png'.format(out_var, out_area)

if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


