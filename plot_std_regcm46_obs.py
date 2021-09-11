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

	path  = '/home/nice/Documents/ufrn/phd_project/papers/paper_rcm_pbl/datas'
	arq   = '{0}/{1}_{2}_{3}_mon_{4}.nc'.format(path, param, area, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	sim = np.nanmean(np.nanmean(value[0:120:12,:,:], axis=1), axis=1)  

	return sim


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
	obs = np.nanmean(np.nanmean(value[0:120:12,:,:], axis=1), axis=1)  

	return obs
			              
               
# Import regcm exps and obs database 
obs_namz = import_obs(u'namz', u'gpcp_v2.3_obs')
obs_samz = import_obs(u'samz', u'gpcp_v2.3_obs')
obs_neb = import_obs(u'neb', u'gpcp_v2.3_obs')

exp1_namz = import_sim(u'namz', u'regcm_exp1')
exp1_samz = import_sim(u'samz', u'regcm_exp1')
exp1_neb = import_sim(u'neb', u'regcm_exp1')

exp2_namz = import_sim(u'namz', u'regcm_exp2')
exp2_samz = import_sim(u'samz', u'regcm_exp2')
exp2_neb = import_sim(u'neb', u'regcm_exp2')

namz_obs_std = np.std(obs_namz)
samz_obs_std = np.std(obs_samz)
neb_obs_std = np.std(obs_neb)

namz_exp1_std = np.std(exp1_namz)
samz_exp1_std = np.std(exp1_samz)
neb_exp1_std = np.std(exp1_neb)

namz_exp2_std = np.std(exp2_namz)
samz_exp2_std = np.std(exp2_samz)
neb_exp2_std = np.std(exp2_neb)


# generate some data
x = np.arange(0, 10, 0.2)
y = np.sin(x)

# plot it
f, (a0, a1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]})
a0.plot(x, y)
a1.plot(y, x)

f.tight_layout()
f.savefig('grid_figure.pdf')


# Plot regcm exps and obs database 
fig = plt.figure()
time = np.arange(1, 10 + 1)

lineStyle_exp1={"linestyle":"-", "linewidth":2, "markeredgewidth":2, "elinewidth":2, "capsize":3}
lineStyle_exp2={"linestyle":"-", "linewidth":2, "markeredgewidth":2, "elinewidth":2, "capsize":3}
lineStyle_obs={"linestyle":"-", "linewidth":2, "markeredgewidth":2, "elinewidth":2, "capsize":3}

ax1 = fig.add_subplot(311)
line_obs = ax1.errorbar(time, obs_namz, yerr=namz_obs_std, **lineStyle_obs, color='black', label='GPCP')
line_exp1 = ax1.errorbar(time, exp1_namz, yerr=namz_exp1_std, **lineStyle_exp1, color='blue', label='Reg_Exp1')
line_exp2 = ax1.errorbar(time, exp2_namz, yerr=namz_exp2_std, **lineStyle_exp2, color='red', label='Reg_Exp2')
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
ax1.set_yticks(np.arange(0, 20, 4))
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'))
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.legend(loc=1, shadow=True, ncol=3, prop=FontProperties(size=8))
plt.text(1.5, 14, 'SD = 2.05', fontsize=8, color='black')
plt.text(1.5, 11.5, 'SD = 1.69', fontsize=8, color='blue')
plt.text(1.5, 9, 'SD = 1.22', fontsize=8, color='red')

ax2 = fig.add_subplot(312)
line_obs = ax2.errorbar(time, obs_samz, yerr=samz_obs_std, **lineStyle_obs, color='black', label='GPCP')
line_exp1 = ax2.errorbar(time, exp1_samz, yerr=samz_exp1_std, **lineStyle_exp1, color='blue', label='Reg_Exp1')
line_exp2 = ax2.errorbar(time, exp2_samz, yerr=samz_exp2_std, **lineStyle_exp2, color='red', label='Reg_Exp2')
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
ax2.set_ylabel('Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
ax2.set_yticks(np.arange(0, 20, 4))
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'))
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.text(1.5, 14, 'SD = 1.62', fontsize=8, color='black')
plt.text(3.5, 14, 'SD = 0.82', fontsize=8, color='blue')
plt.text(5.5, 14, 'SD = 0.62', fontsize=8, color='red')

ax3 = fig.add_subplot(313)
line_obs = ax3.errorbar(time, obs_neb, yerr=neb_obs_std, **lineStyle_obs, color='black', label='GPCP')
line_exp1 = ax3.errorbar(time, exp1_neb, yerr=neb_exp1_std, **lineStyle_exp1, color='blue', label='Reg_Exp1')
line_exp2 = ax3.errorbar(time, exp2_neb, yerr=neb_exp2_std, **lineStyle_exp2, color='red', label='Reg_Exp2')
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
ax3.set_xlabel('Years', fontsize=8, fontweight='bold')
ax3.set_yticks(np.arange(0, 20, 4))
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'))
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')
plt.text(4.5, 14, 'SD = 3.05', fontsize=8, color='black')
plt.text(6.5, 14, 'SD = 1.56', fontsize=8, color='blue')
plt.text(8.5, 14, 'SD = 1.18', fontsize=8, color='red')

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/regcm_pbl/results'
name_out = 'pyplt_std_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()





