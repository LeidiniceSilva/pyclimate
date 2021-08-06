# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/11/2021"
__description__ = "This script plot annual standard desviation from regcm46 and obs database"

import os
import netCDF4
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
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

	sim_annual = np.nanmean(np.nanmean(value[0:120:12,:,:], axis=1), axis=1) 

	return sim_annual


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

	obs_annual = np.nanmean(np.nanmean(value[0:120:12,:,:], axis=1), axis=1) 
	
	return obs_annual
	              
               
# Import regcm exps and obs database 
nam_exp1 = import_sim(u'namz', u'regcm_exp1')
sam_exp1 = import_sim(u'samz', u'regcm_exp1')
neb_exp1 = import_sim(u'neb', u'regcm_exp1')

nam_exp2 = import_sim(u'namz', u'regcm_exp2')
sam_exp2 = import_sim(u'samz', u'regcm_exp2')
neb_exp2 = import_sim(u'neb', u'regcm_exp2')

nam_obs = import_obs(u'namz', u'gpcp_v2.2_obs')
sam_obs = import_obs(u'samz', u'gpcp_v2.2_obs')
neb_obs = import_obs(u'neb', u'gpcp_v2.2_obs')

# Compute mean
nam_exp1_mean = np.mean(nam_exp1)
sam_exp1_mean = np.mean(sam_exp1)
neb_exp1_mean = np.mean(neb_exp1)

nam_exp2_mean = np.mean(nam_exp2)
sam_exp2_mean = np.mean(sam_exp2)
neb_exp2_mean = np.mean(neb_exp2)

nam_obs_mean = np.mean(nam_obs)
sam_obs_mean = np.mean(sam_obs)
neb_obs_mean = np.mean(neb_obs)

# Compute std
nam_exp1_std = np.std(nam_exp1)
sam_exp1_std = np.std(sam_exp1)
neb_exp1_std = np.std(neb_exp1)

nam_exp2_std = np.std(nam_exp2)
sam_exp2_std = np.std(sam_exp2)
neb_exp2_std = np.std(nam_exp2)

nam_obs_std = np.std(nam_obs)
sam_obs_std = np.std(sam_obs)
neb_obs_std = np.std(neb_obs)

# Plot regcm exps and obs database 
fig = plt.figure()
time = np.arange(1, 10 + 1)
lineStyle_exp1={"linestyle":"-", "linewidth":2, "markeredgewidth":2, "elinewidth":2, "capsize":3}
lineStyle_exp2={"linestyle":"-", "linewidth":2, "markeredgewidth":2, "elinewidth":2, "capsize":3}
lineStyle_obs={"linestyle":"-", "linewidth":2, "markeredgewidth":2, "elinewidth":2, "capsize":3}

ax1 = fig.add_subplot(311)
line_exp1 = ax1.errorbar(time, nam_exp1, yerr=nam_exp1_std, **lineStyle_exp1, color='blue', label='Reg_Exp1')
line_exp2 = ax1.errorbar(time, nam_exp2, yerr=nam_exp2_std, **lineStyle_exp2, color='red', label='Reg_Exp2')
line_obs = ax1.errorbar(time, nam_obs, yerr=nam_obs_std, **lineStyle_obs, color='black', label='GPCP')
ax1.set_yticks(np.arange(0, 20, 4))
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'))
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
plt.text(0.75, 14, u'A) NAMZ', fontsize=8, fontweight='bold')
plt.text(2, 14, 'Reg_Exp1 = 4.05(1.53) Reg_Exp2 = 4.01(0.99) GPCP = 7.27(1.66)', fontsize=8)

ax2 = fig.add_subplot(312)
line_exp1 = ax2.errorbar(time, sam_exp1, yerr=sam_exp1_std, **lineStyle_exp1, color='blue', label='Reg_Exp1')
line_exp2 = ax2.errorbar(time, sam_exp2, yerr=sam_exp2_std, **lineStyle_exp2, color='red', label='Reg_Exp2')
line_obs = ax2.errorbar(time, sam_obs, yerr=sam_obs_std, **lineStyle_obs, color='black', label='GPCP')
ax2.set_ylabel('Precipitação (mm d⁻¹)', fontsize=8, fontweight='bold')
ax2.set_yticks(np.arange(0, 20, 4))
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'))
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')
plt.text(0.75, 14, u'B) SAMZ', fontsize=8, fontweight='bold')
plt.text(2, 14, 'Reg_Exp1 = 5.41(0.55) Reg_Exp2 = 4.98(0.55) GPCP = 9.28(1.31)', fontsize=8)

ax3 = fig.add_subplot(313)
line_exp1 = ax3.errorbar(time, neb_exp1, yerr=neb_exp1_std, **lineStyle_exp1, color='blue', label='Reg_Exp1')
line_exp2 = ax3.errorbar(time, neb_exp2, yerr=neb_exp2_std, **lineStyle_exp2, color='red', label='Reg_Exp2')
line_obs = ax3.errorbar(time, neb_obs, yerr=neb_obs_std, **lineStyle_obs, color='black', label='GPCP')
ax3.set_xlabel('Anos', fontsize=8, fontweight='bold')
ax3.set_yticks(np.arange(0, 20, 4))
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'))
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')
plt.text(0.75, 14, u'C) NEB', fontsize=8, fontweight='bold')
plt.legend(loc='upper right', bbox_to_anchor=(1, 0.88), shadow=True, ncol=3, prop=FontProperties(size=8))
plt.text(2, 14, 'Reg_Exp1 = 3.88(1.52) Reg_Exp2 = 4.87(0.99) GPCP = 5.33(2.84)', fontsize=8)

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/regcm_pbl/results'
name_out = 'pyplt_std_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()






