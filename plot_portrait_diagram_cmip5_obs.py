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

from comp_statist_indices import compute_bias


def import_cmip5(model):
	
	param = 'pr' # pr or tas
	area  = 'amz' # amz or neb
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
	sea_mdl = np.nanmean(np.nanmean(value[0:360:3,:,:], axis=1), axis=1)
	djf_mdl = sea_mdl[0:120:4]
	mam_mdl = sea_mdl[1:120:4]
	jja_mdl = sea_mdl[2:120:4]
	son_mdl = sea_mdl[3:120:4]

	return annual_mdl, djf_mdl, mam_mdl, jja_mdl, son_mdl


def import_obs(database):
	
	param = 'pre' # pre or tmp
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
	sea_obs = np.nanmean(np.nanmean(value[0:360:3,:,:], axis=1), axis=1)
	djf_obs = sea_obs[0:120:4]
	mam_obs = sea_obs[1:120:4]
	jja_obs = sea_obs[2:120:4]
	son_obs = sea_obs[3:120:4]

	return annual_obs, djf_obs, mam_obs, jja_obs, son_obs
	
	
seasons = ['DJF','MAM','JJA','SON', 'Anual']
models = ['BCC-CSM1.1','BCC-CSM1.1M','BNU-ESM','CanESM2','CNRM-CM5','CSIRO-ACCESS-1','CSIRO-ACCESS-3','CSIRO-MK36',
'FIO-ESM','GISS-E2-H-CC','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','INMCM4','IPSL-CM5A-LR',
'IPSL-CM5A-MR','IPSL-CM5B-LR','LASG-FGOALS-G2','LASG-FGOALS-S2','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR',
'MPI-ESM-MR','MRI-CGCM3','NCAR-CCSM4','NCAR-CESM1-BGC','NCAR-CESM1-CAM5','NorESM1-ME','NorESM1-M','ensmean_cmip5']

djf = []
mam = []
jja = []
son = []
annual = []

for mdl in models:
	print('CMIP5 Model:', mdl)
	
	# Import cmip5 model end obs database monthly	
	annual_sim, djf_sim, mam_sim, jja_sim, son_sim= import_cmip5(mdl)
		
	obs  = u'cru_ts4.02'
	annual_cru, djf_cru, mam_cru, jja_cru, son_cru = import_obs(obs)
	
	bias_djf = compute_bias(djf_sim, djf_cru)
	bias_mam = compute_bias(mam_sim, mam_cru)
	bias_jja = compute_bias(jja_sim, jja_cru)
	bias_son = compute_bias(son_sim, son_cru)
	bias_anual = compute_bias(annual_sim, annual_cru)
	
	djf.append(bias_djf)
	mam.append(bias_mam)
	jja.append(bias_jja)
	son.append(bias_son)
	annual.append(bias_anual)
	
harvest = np.array([djf, mam, jja, son, annual])

# Plot cmip5 model end obs database       
fig, ax = plt.subplots(figsize=(12,4))     

# Choice variable: Rainfall (AMZ and AMZ) or Temperature (AMZ and AMZ) 
out_var    = u'pre' # pre or tmp
out_area   = u'amz' # amz or neb
area_name  = u'AMZ (Lat:16S 4N, Lon:74W 48W)' # AMZ (Lat:16S 4N, Lon:74W 48W) or NEB (Lat:15S 2N, Lon:46W 34W)

sns.heatmap(harvest, yticklabels=seasons, xticklabels=models, ax=ax, annot=True, fmt=".1f", annot_kws={'size':8}, vmin=-8, vmax=8, cmap='bwr', linewidths=4)

if out_var == 'pre':
	var_name   = u'Precipitação'
else:
	var_name   = u'Temperatura' 
	
plt.title(u'Viés de {0} - {1}  \n CMIP5-hist x CRU-ts4.02 - 1975-2005 (Período de Referência: 1850-2005)'.format(var_name, area_name), fontsize=12, x=0.45, y=0.80)

# Save figure
path_out = '/home/nice'
name_out = 'pyplt_portrait_diagram_{0}_{1}_cmip5_cru_1975-2005.png'.format(out_var, out_area)

if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


