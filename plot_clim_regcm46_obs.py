# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/26/2018"
__description__ = "This script plot climatology graphics from Rec_EXP models end OBS basedata"

import os
import netCDF4
import statistics
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib as mpl

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from sklearn import metrics
from scipy.stats import norm
from matplotlib.font_manager import FontProperties

def import_sim(area, exp):
	
	param = 'pr' # pr or tas
	date  = '2001-2010'

	path  = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
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


def import_obs(area, obs):
	
	param = 'precip' # precip, pre or tmp
	date  = '2001-2010'

	path  = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq   = '{0}/{1}_{2}_{3}_mon_{4}.nc'.format(path, param, area, obs, date)	
		
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
	              
               
# Import regcm exps model end obs database climatology
nam_exp1_clim = import_sim(u'namz', u'regcm_exp1')
sam_exp1_clim = import_sim(u'samz', u'regcm_exp1')
neb_exp1_clim = import_sim(u'neb', u'regcm_exp1')

nam_exp2_clim = import_sim(u'namz', u'regcm_exp2')
sam_exp2_clim = import_sim(u'samz', u'regcm_exp2')
neb_exp2_clim = import_sim(u'neb', u'regcm_exp2')

nam_obs_clim = import_obs(u'namz', u'gpcp_v2.2_obs')
sam_obs_clim = import_obs(u'samz', u'gpcp_v2.2_obs')
neb_obs_clim = import_obs(u'neb', u'gpcp_v2.2_obs')

nam_exp1_median = statistics.median(nam_exp1_clim)
sam_exp1_median = statistics.median(sam_exp1_clim)
neb_exp1_median = statistics.median(neb_exp1_clim)

nam_exp2_median = statistics.median(nam_exp2_clim)
sam_exp2_median = statistics.median(sam_exp2_clim)
neb_exp2_median = statistics.median(neb_exp2_clim)

nam_obs_median = statistics.median(nam_obs_clim)
sam_obs_median = statistics.median(sam_obs_clim)
neb_obs_median = statistics.median(neb_obs_clim)

print(neb_exp1_median)
print(neb_exp2_median)
print(neb_obs_median)

# Plot model end obs data climatology
fig = plt.figure()
time = np.arange(1, 13)
bar_width = .30

# Subplot one
plt.subplot(311)
plt_clim1 = plt.bar(time, nam_exp1_clim, alpha=0.8, color='blue', label='Reg_Exp1', width = 0.25, edgecolor='black')
plt_clim2 = plt.bar(time + .30, nam_exp2_clim, alpha=0.8, color='red', label='Reg_Exp2', width = 0.25, edgecolor='black')
plt_clim3 = plt.bar(time + .60, nam_obs_clim, alpha=0.8, color='black', label='GPCP', width = 0.25, edgecolor='black')
plt.axhline(nam_exp1_median, linewidth=1, linestyle='dashed', color='blue', alpha=0.8)
plt.axhline(nam_exp2_median, linewidth=1, linestyle='dashed', color='red', alpha=0.8)
plt.axhline(nam_obs_median, linewidth=1, linestyle='dashed', color='black', alpha=0.8)
plt.text(0.5, 10, u'A) NAMZ', fontsize=8, fontweight='bold')
plt.xticks(time + .30, ('Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez'))
plt.yticks(np.arange(0, 14, 2))
plt.text(2, 10, 'Reg_Exp1 = 4.06 Reg_Exp2 = 3.78 GPCP = 9.96', fontsize=8)

# Subplot two
plt.subplot(312)
plt_clim1 = plt.bar(time, sam_exp1_clim, alpha=0.8, color='blue', label='Reg_Exp1', width = 0.25, edgecolor='black')
plt_clim2 = plt.bar(time + .30, sam_exp2_clim, alpha=0.8, color='red', label='Reg_Exp2', width = 0.25, edgecolor='black')
plt_clim3 = plt.bar(time + .60, sam_obs_clim, alpha=0.8, color='black', label='GPCP', width = 0.25, edgecolor='black')
plt.axhline(sam_exp1_median, linewidth=1, linestyle='dashed', color='blue', alpha=0.8)
plt.axhline(sam_exp2_median, linewidth=1, linestyle='dashed', color='red', alpha=0.8)
plt.axhline(sam_obs_median, linewidth=1, linestyle='dashed', color='black', alpha=0.8)
plt.ylabel('Precipitação (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.text(0.5, 10, u'B) SAMZ', fontsize=8, fontweight='bold')
plt.xticks(time + .30, ('Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez'))
plt.yticks(np.arange(0, 14, 2))
plt.text(2, 10, 'Reg_Exp1 = 2.38 Reg_Exp2 = 2.74 GPCP = 4.53', fontsize=8)

# Subplot three
plt.subplot(313)
plt_clim1 = plt.bar(time, neb_exp1_clim, alpha=0.8, color='blue', label='Reg_Exp1', width = 0.25, edgecolor='black')
plt_clim2 = plt.bar(time + .30, neb_exp2_clim, alpha=0.8, color='red', label='Reg_Exp2', width = 0.25, edgecolor='black')
plt_clim3 = plt.bar(time + .60, neb_obs_clim, alpha=0.8, color='black', label='GPCP', width = 0.25, edgecolor='black')
plt.axhline(neb_exp1_median, linewidth=1, linestyle='dashed', color='blue', alpha=0.8)
plt.axhline(neb_exp2_median, linewidth=1, linestyle='dashed', color='red', alpha=0.8)
plt.axhline(neb_obs_median, linewidth=1, linestyle='dashed', color='black', alpha=0.8)
plt.xlabel('Meses', fontsize=8, fontweight='bold')
plt.text(0.5, 10, u'C) NEB', fontsize=8, fontweight='bold')
plt.xticks(time + .30, ('Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez'))
plt.yticks(np.arange(0, 14, 2))
plt.text(2, 10, 'Reg_Exp1 = 2.57 Reg_Exp2 = 2.74 GPCP = 2.98', fontsize=8)
plt.legend(loc='upper right', bbox_to_anchor=(1, 0.80), shadow=True, ncol=3, prop=FontProperties(size=8))

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/regcm_pbl/results'
name_out = 'pyplt_clim_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()






