# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
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


def import_sim(exp):
	
	param = 'pr' # pr or tas
	area  = 'amz_neb' # amz or neb
	date  = '2001-2010'

	path  = '/home/nice/Documentos/ufrn/papers/regcm_pbl/datas'
	arq   = '{0}/{1}_{2}_{3}_mon_{4}.nc'.format(path, param, area, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	exp_mdl = np.nanmean(np.nanmean(value, axis=1), axis=1)

	mdl_clim = []
	for mon in range(1, 12 + 1):
		mdl = np.nanmean(exp_mdl[mon::12], axis=0)
		mdl_clim.append(mdl)

	return mdl_clim


def import_obs(obs):
	
	param = 'precip' # precip, pre or tmp
	area  = 'amz_neb' # amz or neb
	date  = '2001-2010'

	path  = '/home/nice/Documentos/ufrn/papers/regcm_pbl/datas'
	arq   = '{0}/{1}_{2}_{3}_mon_{4}.nc'.format(path, param, area, obs, date)	
		
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
exp  = u'regcm_exp1'
exp1_clim = import_sim(exp)
		
exp  = u'regcm_exp2'
exp2_clim = import_sim(exp)

obs  = u'gpcp_v2.2_obs'
obs_clim, obs_clim_p5, obs_clim_p95  = import_obs(obs)

# Plot model end obs data climatology
fig = plt.figure(1, figsize=(34, 24))
time = np.arange(1, 13)
bar_width = .30

# Subplot one
plt.subplot(211)
plt_clim1 = plt.bar(time, exp1_clim, alpha=0.8, color='blue', label='Reg_Exp1', width = 0.25, edgecolor='black')
plt_clim2 = plt.bar(time + .30, exp2_clim, alpha=0.8, color='red', label='Reg_Exp2', width = 0.25, edgecolor='black')
plt_clim3 = plt.bar(time + .60, obs_clim, alpha=0.8, color='black', label='GPCP', width = 0.25, edgecolor='black')

plt.ylabel('Precipitação (mm)',  fontsize=30)
plt.title('Climatologia de Precipitação NAMZ (2001-2010)', fontsize=32,  fontweight='bold')
plt.xticks(time + .30, ('Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez'), fontsize=30)
plt.yticks(np.arange(0, 14, 2), fontsize=30)
plt.legend(loc='upper center', shadow=True, ncol=3, prop=FontProperties(size=30))
plt.tick_params(axis='both', which='major', labelsize=30, length=10, width=4, pad=8, labelcolor='black')

# Subplot two
plt.subplot(212)
plt_clim1 = plt.bar(time, exp2_clim, bar_width, alpha=0.8, color='green', label='Reg_Exp2')
plt_clim2 = plt.bar(time + bar_width, obs_clim, bar_width, alpha=0.8, color='black', label='GPCP')

plt.xlabel('Meses', fontsize=30)
plt.ylabel('Precipitação (mm)', fontsize=30)
plt.xticks(time + .15, ('Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez'), fontsize=30)
plt.yticks(np.arange(0, 14, 2), fontsize=30)
plt.legend(loc='upper center', shadow=True, ncol=2, prop=FontProperties(size=30))
plt.tick_params(axis='both', which='major', labelsize=30, length=10, width=4, pad=8, labelcolor='black')

# Path out to save figure
path_out = '/home/nice/Documentos/ufrn/papers/regcm_pbl/results'
name_out = 'pyplt_clim_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()





# Compute statiscts index from best CMIP5 models and CRU obs database
r2_exp1 = metrics.r2_score(exp1_clim, obs_clim)
r2_exp2 = metrics.r2_score(exp2_clim, obs_clim)

# Plot model end obs data climatology
fig, ax = plt.subplots(figsize=(28, 16))
time = np.arange(1, 13)

plt_clim = plt.plot(time, exp1_clim, time, exp2_clim, time, obs_clim, time, obs_clim_p5, time, obs_clim_p95)

l1, l2, l3, l4, l5 = plt_clim

plt.setp(l1, linewidth=6, markeredgewidth=1, color='blue')
plt.setp(l2, linewidth=6, markeredgewidth=1, color='red')
plt.setp(l3, linewidth=6, markeredgewidth=1, color='black')
plt.setp(l4, linewidth=3, markeredgewidth=3, color='slategray')
plt.setp(l5, linewidth=3, markeredgewidth=3, color='slategray')


plt.fill_between(time, obs_clim_p5, obs_clim_p95, facecolor='slategray', alpha=0.8, interpolate=True)

# choice variable: Rainfall (AMZ and AMZ) or Temperature (AMZ and AMZ) 
out_var    = u'pr' # pre or tmp
out_area   = u'amz_neb' # amz or neb
area_name  = u'AMZ_NEB (Lat:16S 4N, Lon:74W 48W)' # AMZ (Lat:16S 4N, Lon:74W 48W) or NEB (Lat:15S 2N, Lon:46W 34W)

if out_var == 'pr':
	yaxis = np.arange(0, 14, 2)
	var_name   = u'Rainfall'
	label_name = u'Rain (mm/day)' 
	plt.text(9.5, 10.5, u'r2 = {0}'.format(round(r2_exp1, 3)), fontsize=24, fontweight='bold')
	plt.text(9.5, 11.5, u'r2 = {0}'.format(round(r2_exp2, 3)), fontsize=24, fontweight='bold')

else:
	yaxis = np.arange(18, 34, 2)
	var_name   = u'Temperature' 
	label_name = u'Temperature 2m ($^\circ$C)' 
	plt.text(9.5, 10.5, u'r2 = {0}'.format(round(r2_exp1, 3)), fontsize=24, fontweight='bold')
	plt.text(10.5, 10.5, u'r2 = {0}'.format(round(r2_exp2, 3)), fontsize=24, fontweight='bold')
	
fig.suptitle(u'{0} Climatology - {1} (2001-2010)'.format(var_name, area_name), fontsize=30, y=0.98)

xaxis = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

plt.xlabel(u'Months', fontsize=30)
plt.ylabel(u'{0}'.format(label_name), fontsize=30)
plt.xticks(time, xaxis, fontsize=30)
plt.yticks(yaxis, fontsize=30)
plt.tick_params(axis='both', which='major', labelsize=30, length=10, width=4, pad=8, labelcolor='black')

font = FontProperties(size=24)	
legend = (u'Reg_EXP1', u'Reg_EXP2', u'GPCP', u'GPCP_p5%', u'GPCP_p95%')    
plt.legend(plt_clim, legend, loc='upper center', bbox_to_anchor=(0.5, -0.07), shadow=True, ncol=6, prop=font)
ax.xaxis.grid(True, which='major', linestyle='--', linewidth='1.4', zorder=0.6)
ax.yaxis.grid(True, which='major', linestyle='--', linewidth='1.4', zorder=0.6)
    
path_out = '/home/nice/Documentos/ufrn/papers/regcm_exp/results'
name_out = 'pyplt_clim_{0}_{1}_regcm_exp_gpcp_2001-2010.png'.format(out_var, out_area)

if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

