# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot annual climatology from regcm46 and obs database"

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


def import_obs(area, obs):
	
	param = 'precip' # precip, pre or tmp
	date  = '2001-2010'

	path  = '/home/nice/Documents/ufrn/phd_project/papers/paper_rcm_pbl/datas'
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
	

def import_sim(area, exp):
	
	param = 'pr' # pr or tas
	date  = '2001-2010'

	path  = '/home/nice/Documents/ufrn/phd_project/papers/paper_rcm_pbl/datas'
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
	              
               
# Import regcm exps and obs database 
nam_exp1_clim = import_sim(u'namz', u'regcm_exp1')
sam_exp1_clim = import_sim(u'samz', u'regcm_exp1')
neb_exp1_clim = import_sim(u'neb', u'regcm_exp1')

nam_exp2_clim = import_sim(u'namz', u'regcm_exp2')
sam_exp2_clim = import_sim(u'samz', u'regcm_exp2')
neb_exp2_clim = import_sim(u'neb', u'regcm_exp2')

nam_obs_clim = import_obs(u'namz', u'gpcp_v2.3_obs')
sam_obs_clim = import_obs(u'samz', u'gpcp_v2.3_obs')
neb_obs_clim = import_obs(u'neb', u'gpcp_v2.3_obs')

nam_exp1_median = np.nanmean(nam_exp1_clim)
sam_exp1_median = np.nanmean(sam_exp1_clim)
neb_exp1_median = np.nanmean(neb_exp1_clim)

nam_exp2_median = np.nanmean(nam_exp2_clim)
sam_exp2_median = np.nanmean(sam_exp2_clim)
neb_exp2_median = np.nanmean(neb_exp2_clim)

nam_obs_median = np.nanmean(nam_obs_clim)
sam_obs_median = np.nanmean(sam_obs_clim)
neb_obs_median = np.nanmean(neb_obs_clim)

# Plot regcm exps and obs database 
fig = plt.figure()
time = np.arange(1, 13)
bar_width = .30

# Subplot one
ax = fig.add_subplot(3, 1, 1)
plt_clim1 = plt.bar(time, nam_obs_clim, alpha=0.8, color='black', label='GPCP', width = 0.25, edgecolor='black')
plt_clim2 = plt.bar(time + .30, nam_exp1_clim, alpha=0.8, color='blue', label='Reg_Exp1', width = 0.25, edgecolor='black')
plt_clim3 = plt.bar(time + .60, nam_exp2_clim, alpha=0.8, color='red', label='Reg_Exp2', width = 0.25, edgecolor='black')
plt.axhline(nam_exp1_median, linewidth=1, linestyle='dashed', color='blue', alpha=0.8)
plt.axhline(nam_exp2_median, linewidth=1, linestyle='dashed', color='red', alpha=0.8)
plt.axhline(nam_obs_median, linewidth=1, linestyle='dashed', color='black', alpha=0.8)
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time + .30, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Ouc', 'Nov', 'Dec'))
plt.yticks(np.arange(0, 14, 2))
plt.setp(ax.get_xticklabels(), visible=False)
plt.text(7, 7, 'ME = 6.83', fontsize=8, color='black')
plt.text(7, 5, 'ME = 4.46', fontsize=8, color='blue')
plt.text(9, 5, 'ME = 4.22', fontsize=8, color='red')
plt.legend(loc=1, handlelength=0.75, handleheight=0.75, shadow=True, ncol=3, prop=FontProperties(size=8))

# Subplot two
ax = fig.add_subplot(3, 1, 2)
plt_clim1 = plt.bar(time, sam_obs_clim, alpha=0.8, color='black', label='GPCP', width = 0.25, edgecolor='black')
plt_clim2 = plt.bar(time + .30, sam_exp1_clim, alpha=0.8, color='blue', label='Reg_Exp1', width = 0.25, edgecolor='black')
plt_clim3 = plt.bar(time + .60, sam_exp2_clim, alpha=0.8, color='red', label='Reg_Exp2', width = 0.25, edgecolor='black')
plt.axhline(sam_exp1_median, linewidth=1, linestyle='dashed', color='blue', alpha=0.8)
plt.axhline(sam_exp2_median, linewidth=1, linestyle='dashed', color='red', alpha=0.8)
plt.axhline(sam_obs_median, linewidth=1, linestyle='dashed', color='black', alpha=0.8)
plt.ylabel('Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time + .30, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Ouc', 'Nov', 'Dec'))
plt.yticks(np.arange(0, 14, 2))
plt.setp(ax.get_xticklabels(), visible=False)
plt.text(5, 6.5, 'ME = 6.20', fontsize=8, color='black')
plt.text(5, 4, 'ME = 2.83', fontsize=8, color='blue')
plt.text(7, 4, 'ME = 3.15', fontsize=8, color='red')

# Subplot three
ax = fig.add_subplot(3, 1, 3)
plt_clim1 = plt.bar(time, neb_obs_clim, alpha=0.8, color='black', label='GPCP', width = 0.25, edgecolor='black')
plt_clim2 = plt.bar(time + .30, neb_exp1_clim, alpha=0.8, color='blue', label='Reg_Exp1', width = 0.25, edgecolor='black')
plt_clim3 = plt.bar(time + .60, neb_exp2_clim, alpha=0.8, color='red', label='Reg_Exp2', width = 0.25, edgecolor='black')
plt.axhline(neb_exp1_median, linewidth=1, linestyle='dashed', color='blue', alpha=0.8)
plt.axhline(neb_exp2_median, linewidth=1, linestyle='dashed', color='red', alpha=0.8)
plt.axhline(neb_obs_median, linewidth=1, linestyle='dashed', color='black', alpha=0.8)
plt.xlabel('Months', fontsize=8, fontweight='bold')
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time + .30, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Ouc', 'Nov', 'Dec'))
plt.yticks(np.arange(0, 14, 2))
plt.text(5, 4, 'ME = 2.99', fontsize=8, color='black')
plt.text(7, 4, 'ME = 2.66', fontsize=8, color='blue')
plt.text(9, 4, 'ME = 2.60', fontsize=8, color='red')

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/phd_project/papers/paper_rcm_pbl/figs'
name_out = 'pyplt_clim_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()






