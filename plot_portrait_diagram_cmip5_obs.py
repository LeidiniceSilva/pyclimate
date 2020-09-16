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
import matplotlib.colors as colors

from comp_statist_indices import compute_bias


def import_cmip5_pr(area, model):
	
	param = 'pr' # pr or tas
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/phd_project/datas/cmip5/hist'
	arq   = '{0}/{1}_{2}_Amon_{3}_{4}_{5}.nc'.format(path, param, area,	model, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_mdl_pr = np.nanmean(np.nanmean(value, axis=1), axis=1) 
	sea_mdl_pr = np.nanmean(np.nanmean(value[0:360:3,:,:], axis=1), axis=1)
	djf_mdl_pr = sea_mdl_pr[0:120:4]
	mam_mdl_pr = sea_mdl_pr[1:120:4]
	jja_mdl_pr = sea_mdl_pr[2:120:4]
	son_mdl_pr = sea_mdl_pr[3:120:4]

	return annual_mdl_pr, djf_mdl_pr, mam_mdl_pr, jja_mdl_pr, son_mdl_pr


def import_cru_pre(area, obs):
	
	param = 'pre' # pre or tmp
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/phd_project/datas/obs_data'
	arq   = '{0}/{1}_{2}_{3}_obs_mon_{4}.nc'.format(path, param, area, obs, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_obs_pre = np.nanmean(np.nanmean(value, axis=1), axis=1) 
	sea_obs_pre = np.nanmean(np.nanmean(value[0:360:3,:,:], axis=1), axis=1)
	djf_obs_pre = sea_obs_pre[0:120:4]
	mam_obs_pre = sea_obs_pre[1:120:4]
	jja_obs_pre = sea_obs_pre[2:120:4]
	son_obs_pre = sea_obs_pre[3:120:4]

	return annual_obs_pre, djf_obs_pre, mam_obs_pre, jja_obs_pre, son_obs_pre


def import_cmip5_tas(area, model):
	
	param = 'tas' # pr or tas
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/phd_project/datas/cmip5/hist'
	arq   = '{0}/{1}_{2}_Amon_{3}_{4}_{5}.nc'.format(path, param, area,	model, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_mdl_tas = np.nanmean(np.nanmean(value, axis=1), axis=1) 
	sea_mdl_tas = np.nanmean(np.nanmean(value[0:360:3,:,:], axis=1), axis=1)
	djf_mdl_tas = sea_mdl_tas[0:120:4]
	mam_mdl_tas = sea_mdl_tas[1:120:4]
	jja_mdl_tas = sea_mdl_tas[2:120:4]
	son_mdl_tas = sea_mdl_tas[3:120:4]

	return annual_mdl_tas, djf_mdl_tas, mam_mdl_tas, jja_mdl_tas, son_mdl_tas


def import_cru_tmp(area, obs):
	
	param = 'tmp' # pre or tmp
	date  = '197512-200511'

	path  = '/home/nice/Documents/ufrn/phd_project/datas/obs_data'
	arq   = '{0}/{1}_{2}_{3}_obs_mon_{4}.nc'.format(path, param, area, obs, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_obs_tmp = np.nanmean(np.nanmean(value, axis=1), axis=1) 
	sea_obs_tmp = np.nanmean(np.nanmean(value[0:360:3,:,:], axis=1), axis=1)
	djf_obs_tmp = sea_obs_tmp[0:120:4]
	mam_obs_tmp = sea_obs_tmp[1:120:4]
	jja_obs_tmp = sea_obs_tmp[2:120:4]
	son_obs_tmp = sea_obs_tmp[3:120:4]

	return annual_obs_tmp, djf_obs_tmp, mam_obs_tmp, jja_obs_tmp, son_obs_tmp
	
	
seasons = ['DJF','MAM','JJA','SON', 'Anual']
areas = ['amz','neb','matopiba']
models = ['BCC-CSM1.1','BCC-CSM1.1M','BNU-ESM','CanESM2','CNRM-CM5','CSIRO-ACCESS-1','CSIRO-ACCESS-3','CSIRO-MK36',
'FIO-ESM','GISS-E2-H-CC','GISS-E2-H','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','INMCM4','IPSL-CM5A-LR',
'IPSL-CM5A-MR','LASG-FGOALS-G2','LASG-FGOALS-S2','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR',
'MPI-ESM-MR','MRI-CGCM3','NCAR-CCSM4','NCAR-CESM1-BGC','NCAR-CESM1-CAM5','NorESM1-ME','NorESM1-M','ensmean_cmip5']

djfap = []
mamap = []
jjaap = []
sonap = []
annualap = []

djfnp = []
mamnp = []
jjanp = []
sonnp = []
annualnp = []

djfmp = []
mammp = []
jjamp = []
sonmp = []
annualmp = []

djfat = []
mamat = []
jjaat = []
sonat = []
annualat = []

djfnt = []
mamnt = []
jjant = []
sonnt = []
annualnt = []

djfmt = []
mammt = []
jjamt = []
sonmt = []
annualmt = []

for model in models:
	print('CMIP5 Model:', model)
	
	# Import cmip5 model end obs database monthly pr
	annual_sim_amz_pr, djf_sim_amz_pr, mam_sim_amz_pr, jja_sim_amz_pr, son_sim_amz_pr = import_cmip5_pr('amz', model)
	annual_sim_neb_pr, djf_sim_neb_pr, mam_sim_neb_pr, jja_sim_neb_pr, son_sim_neb_pr = import_cmip5_pr('neb', model)
	annual_sim_matopiba_pr, djf_sim_matopiba_pr, mam_sim_matopiba_pr, jja_sim_matopiba_pr, son_sim_matopiba_pr = import_cmip5_pr('matopiba', model)

	annual_cru_amz_pre, djf_cru_amz_pre, mam_cru_amz_pre, jja_cru_amz_pre, son_cru_amz_pre = import_cru_pre('amz', u'cru_ts4.02')
	annual_cru_neb_pre, djf_cru_neb_pre, mam_cru_neb_pre, jja_cru_neb_pre, son_cru_neb_pre = import_cru_pre('neb', u'cru_ts4.02')
	annual_cru_matopiba_pre, djf_cru_matopiba_pre, mam_cru_matopiba_pre, jja_cru_matopiba_pre, son_cru_matopiba_pre = import_cru_pre('matopiba', u'cru_ts4.02')

	# Import cmip5 model end obs database monthly tas	
	annual_sim_amz_tas, djf_sim_amz_tas, mam_sim_amz_tas, jja_sim_amz_tas, son_sim_amz_tas = import_cmip5_tas('amz', model)
	annual_sim_neb_tas, djf_sim_neb_tas, mam_sim_neb_tas, jja_sim_neb_tas, son_sim_neb_tas = import_cmip5_tas('neb', model)
	annual_sim_matopiba_tas, djf_sim_matopiba_tas, mam_sim_matopiba_tas, jja_sim_matopiba_tas, son_sim_matopiba_tas = import_cmip5_tas('matopiba', model)

	annual_cru_amz_tmp, djf_cru_amz_tmp, mam_cru_amz_tmp, jja_cru_amz_tmp, son_cru_amz_tmp = import_cru_tmp('amz', u'cru_ts4.02')
	annual_cru_neb_tmp, djf_cru_neb_tmp, mam_cru_neb_tmp, jja_cru_neb_tmp, son_cru_neb_tmp = import_cru_tmp('neb', u'cru_ts4.02')
	annual_cru_matopiba_tmp, djf_cru_matopiba_tmp, mam_cru_matopiba_tmp, jja_cru_matopiba_tmp, son_cru_matopiba_tmp = import_cru_tmp('matopiba', u'cru_ts4.02')

	# Compute bias pr
	bias_djf_amz_pr = compute_bias(djf_sim_amz_pr, djf_cru_amz_pre)
	bias_mam_amz_pr = compute_bias(mam_sim_amz_pr, mam_cru_amz_pre)
	bias_jja_amz_pr = compute_bias(jja_sim_amz_pr, jja_cru_amz_pre)
	bias_son_amz_pr = compute_bias(son_sim_amz_pr, son_cru_amz_pre)
	bias_anual_amz_pr = compute_bias(annual_sim_amz_pr, annual_cru_amz_pre)
	
	bias_djf_neb_pr = compute_bias(djf_sim_neb_pr, djf_cru_neb_pre)
	bias_mam_neb_pr = compute_bias(mam_sim_neb_pr, mam_cru_neb_pre)
	bias_jja_neb_pr = compute_bias(jja_sim_neb_pr, jja_cru_neb_pre)
	bias_son_neb_pr = compute_bias(son_sim_neb_pr, son_cru_neb_pre)
	bias_anual_neb_pr = compute_bias(annual_sim_neb_pr, annual_cru_neb_pre)
	
	bias_djf_matopiba_pr = compute_bias(djf_sim_matopiba_pr, djf_cru_matopiba_pre)
	bias_mam_matopiba_pr = compute_bias(mam_sim_matopiba_pr, mam_cru_matopiba_pre)
	bias_jja_matopiba_pr = compute_bias(jja_sim_matopiba_pr, jja_cru_matopiba_pre)
	bias_son_matopiba_pr = compute_bias(son_sim_matopiba_pr, son_cru_matopiba_pre)
	bias_anual_matopiba_pr = compute_bias(annual_sim_matopiba_pr, annual_cru_matopiba_pre)
	
	djfap.append(bias_djf_amz_pr)
	mamap.append(bias_mam_amz_pr)
	jjaap.append(bias_jja_amz_pr)
	sonap.append(bias_son_amz_pr)
	annualap.append(bias_anual_amz_pr)
	
	djfnp.append(bias_djf_neb_pr)
	mamnp.append(bias_mam_neb_pr)
	jjanp.append(bias_jja_neb_pr)
	sonnp.append(bias_son_neb_pr)
	annualnp.append(bias_anual_neb_pr)
	
	djfmp.append(bias_djf_matopiba_pr)
	mammp.append(bias_mam_matopiba_pr)
	jjamp.append(bias_jja_matopiba_pr)
	sonmp.append(bias_son_matopiba_pr)
	annualmp.append(bias_anual_matopiba_pr)
	
	amz_pr = np.array([djfap, mamap, jjaap, sonap, annualap])
	neb_pr = np.array([djfnp, mamnp, jjanp, sonnp, annualnp])
	matopiba_pr = np.array([djfmp, mammp, jjamp, sonmp, annualmp])
	
	# Compute bias tas
	bias_djf_amz_tas = compute_bias(djf_sim_amz_tas, djf_cru_amz_tmp)
	bias_mam_amz_tas = compute_bias(mam_sim_amz_tas, mam_cru_amz_tmp)
	bias_jja_amz_tas = compute_bias(jja_sim_amz_tas, jja_cru_amz_tmp)
	bias_son_amz_tas = compute_bias(son_sim_amz_tas, son_cru_amz_tmp)
	bias_anual_amz_tas = compute_bias(annual_sim_amz_tas, annual_cru_amz_tmp)
	
	bias_djf_neb_tas = compute_bias(djf_sim_neb_tas, djf_cru_neb_tmp)
	bias_mam_neb_tas = compute_bias(mam_sim_neb_tas, mam_cru_neb_tmp)
	bias_jja_neb_tas = compute_bias(jja_sim_neb_tas, jja_cru_neb_tmp)
	bias_son_neb_tas = compute_bias(son_sim_neb_tas, son_cru_neb_tmp)
	bias_anual_neb_tas = compute_bias(annual_sim_neb_tas, annual_cru_neb_tmp)
	
	bias_djf_matopiba_tas = compute_bias(djf_sim_matopiba_tas, djf_cru_matopiba_tmp)
	bias_mam_matopiba_tas = compute_bias(mam_sim_matopiba_tas, mam_cru_matopiba_tmp)
	bias_jja_matopiba_tas = compute_bias(jja_sim_matopiba_tas, jja_cru_matopiba_tmp)
	bias_son_matopiba_tas = compute_bias(son_sim_matopiba_tas, son_cru_matopiba_tmp)
	bias_anual_matopiba_tas = compute_bias(annual_sim_matopiba_tas, annual_cru_matopiba_tmp)
	
	djfat.append(bias_djf_amz_tas)
	mamat.append(bias_mam_amz_tas)
	jjaat.append(bias_jja_amz_tas)
	sonat.append(bias_son_amz_tas)
	annualat.append(bias_anual_amz_tas)
	
	djfnt.append(bias_djf_neb_tas)
	mamnt.append(bias_mam_neb_tas)
	jjant.append(bias_jja_neb_tas)
	sonnt.append(bias_son_neb_tas)
	annualnt.append(bias_anual_neb_tas)
	
	djfmt.append(bias_djf_matopiba_tas)
	mammt.append(bias_mam_matopiba_tas)
	jjamt.append(bias_jja_matopiba_tas)
	sonmt.append(bias_son_matopiba_tas)
	annualmt.append(bias_anual_matopiba_tas)
	
	amz_tas = np.array([djfat, mamat, jjaat, sonat, annualat])
	neb_tas = np.array([djfnt, mamnt, jjant, sonnt, annualnt])
	matopiba_tas = np.array([djfmt, mammt, jjamt, sonmt, annualmt])


# Plot model end obs data climatology
fig, axes = plt.subplots(nrows=3, ncols=2, constrained_layout=True, figsize=(10, 8))
norm = colors.BoundaryNorm(boundaries=np.linspace(-8, 8, 10), ncolors=256)

xlabels = [u'BCC-CSM1.1',u'',u'BNU-ESM',u'',u'CNRM-CM5',u'',u'CSIRO-ACCESS-3',u'',u'FIO-ESM',u'',
u'GISS-E2-H',u'',u'HadGEM2-CC',u'',u'INMCM4',u'',u'IPSL-CM5A-MR',u'',u'LASG-FGOALS-S2',u'',u'MIROC-ESM-CHEM',
u'',u'MPI-ESM-LR',u'',u'MRI-CGCM3',u'',u'NCAR-CESM1-BGC',u'',u'NorESM1-ME',u'',u'ensmean_cmip5']
ylabels = [u'DJF', u'MAM', u'JJA', u'SON', u'Annual']

# First column heatmaps with same colormap
pcm1 = axes[0, 0].pcolormesh(amz_pr, edgecolors ='k', linewidths = 1, cmap='BrBG')
axes[0, 0].set_title(u'A)', loc='left', fontweight='bold')
axes[0, 0].set_xticks(np.arange(amz_pr.shape[1]) + 0.5)
axes[0, 0].set_yticks(np.arange(amz_pr.shape[0]) + 0.5)
axes[0, 0].set_yticklabels(ylabels)
plt.setp(axes[0, 0].get_xticklabels(),visible=False)

pcm2 = axes[1, 0].pcolormesh(neb_pr, edgecolors ='k', linewidths = 1, norm=norm, cmap='BrBG')
axes[1, 0].set_title(u'C)', loc='left', fontweight='bold')
axes[1, 0].set_xticks(np.arange(neb_pr.shape[1]) + 0.5)
axes[1, 0].set_yticks(np.arange(neb_pr.shape[0]) + 0.5)
axes[1, 0].set_yticklabels(ylabels)
plt.setp(axes[1, 0].get_xticklabels(),visible=False)

pcm3 = axes[2, 0].pcolormesh(matopiba_pr, edgecolors ='k', linewidths = 1, cmap='BrBG')
axes[2, 0].set_title(u'E)', loc='left', fontweight='bold')
axes[2, 0].set_xticks(np.arange(matopiba_pr.shape[1]) + 0.5)
axes[2, 0].set_yticks(np.arange(matopiba_pr.shape[0]) + 0.5)
axes[2, 0].set_xticklabels(xlabels, rotation=90)
axes[2, 0].set_yticklabels(ylabels)
clb=fig.colorbar(pcm2, ax=axes[:, 0], extend='both', shrink=0.7)

# Second column heatmaps with same colormap
pcm4 = axes[0, 1].pcolormesh(amz_tas, edgecolors ='k', linewidths = 1, cmap='coolwarm')
axes[0, 1].set_title('B)', loc='left', fontweight='bold')
axes[0, 1].set_xticks(np.arange(amz_tas.shape[1]) + 0.5)
axes[0, 1].set_yticks(np.arange(amz_tas.shape[0]) + 0.5)
axes[0, 1].set_yticklabels(ylabels)
plt.setp(axes[0, 1].get_xticklabels(),visible=False)

pcm5 = axes[1, 1].pcolormesh(neb_tas, edgecolors ='k', linewidths = 1, norm=norm, cmap='coolwarm')
axes[1, 1].set_title('D)', loc='left', fontweight='bold')
axes[1, 1].set_xticks(np.arange(neb_tas.shape[1]) + 0.5)
axes[1, 1].set_yticks(np.arange(neb_tas.shape[0]) + 0.5)
axes[1, 1].set_yticklabels(ylabels)
plt.setp(axes[1, 1].get_xticklabels(),visible=False)

pcm6 = axes[2, 1].pcolormesh(matopiba_tas, edgecolors ='k', linewidths = 1, cmap='coolwarm')
axes[2, 1].set_title('F)', loc='left', fontweight='bold')
axes[2, 1].set_xticks(np.arange(matopiba_tas.shape[1]) + 0.5)
axes[2, 1].set_yticks(np.arange(matopiba_tas.shape[0]) + 0.5)
axes[2, 1].set_xticklabels(xlabels, rotation=90)
axes[2, 1].set_yticklabels(ylabels)

clb=fig.colorbar(pcm5, ax=axes[:, 1], extend='both', shrink=0.7)
clb.ax.set_title('(°C)')

# Save figure
path_out = '/home/nice'
name_out = 'pyplt_portrait_diagram_cmip5_cru_1975-2005.png'

if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()


