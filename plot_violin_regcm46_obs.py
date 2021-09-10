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

	sim_yr1 = np.nanmean(np.nanmean(value[0:12,:,:], axis=1), axis=1) 
	sim_yr2 = np.nanmean(np.nanmean(value[12:24,:,:], axis=1), axis=1) 
	sim_yr3 = np.nanmean(np.nanmean(value[24:36,:,:], axis=1), axis=1) 
	sim_yr4 = np.nanmean(np.nanmean(value[36:48,:,:], axis=1), axis=1) 
	sim_yr5 = np.nanmean(np.nanmean(value[48:60,:,:], axis=1), axis=1) 
	sim_yr6 = np.nanmean(np.nanmean(value[60:72,:,:], axis=1), axis=1) 
	sim_yr7 = np.nanmean(np.nanmean(value[72:84,:,:], axis=1), axis=1) 
	sim_yr8 = np.nanmean(np.nanmean(value[84:96,:,:], axis=1), axis=1) 
	sim_yr9 = np.nanmean(np.nanmean(value[96:108,:,:], axis=1), axis=1) 
	sim_yr10 = np.nanmean(np.nanmean(value[108:120,:,:], axis=1), axis=1) 

	return sim_yr1, sim_yr2, sim_yr3, sim_yr4, sim_yr5, sim_yr6, sim_yr7, sim_yr8, sim_yr9, sim_yr10 


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

	obs_yr1 = np.nanmean(np.nanmean(value[0:12,:,:], axis=1), axis=1) 
	obs_yr2 = np.nanmean(np.nanmean(value[12:24,:,:], axis=1), axis=1) 
	obs_yr3 = np.nanmean(np.nanmean(value[24:36,:,:], axis=1), axis=1) 
	obs_yr4 = np.nanmean(np.nanmean(value[36:48,:,:], axis=1), axis=1) 
	obs_yr5 = np.nanmean(np.nanmean(value[48:60,:,:], axis=1), axis=1) 
	obs_yr6 = np.nanmean(np.nanmean(value[60:72,:,:], axis=1), axis=1) 
	obs_yr7 = np.nanmean(np.nanmean(value[72:84,:,:], axis=1), axis=1) 
	obs_yr8 = np.nanmean(np.nanmean(value[84:96,:,:], axis=1), axis=1) 
	obs_yr9 = np.nanmean(np.nanmean(value[96:108,:,:], axis=1), axis=1) 
	obs_yr10 = np.nanmean(np.nanmean(value[108:120,:,:], axis=1), axis=1) 

	return obs_yr1, obs_yr2, obs_yr3, obs_yr4, obs_yr5, obs_yr6, obs_yr7, obs_yr8, obs_yr9, obs_yr10 
			              
               
# Import regcm exps and obs database 
obs_yr1_namz, obs_yr2_namz, obs_yr3_namz, obs_yr4_namz, obs_yr5_namz, obs_yr6_namz, obs_yr7_namz, obs_yr8_namz, obs_yr9_namz, obs_yr10_namz = import_obs(u'namz', u'gpcp_v2.3_obs')
obs_yr1_samz, obs_yr2_samz, obs_yr3_samz, obs_yr4_samz, obs_yr5_samz, obs_yr6_samz, obs_yr7_samz, obs_yr8_samz, obs_yr9_samz, obs_yr10_samz = import_obs(u'samz', u'gpcp_v2.3_obs')
obs_yr1_neb, obs_yr2_neb, obs_yr3_neb, obs_yr4_neb, obs_yr5_neb, obs_yr6_neb, obs_yr7_neb, obs_yr8_neb, obs_yr9_neb, obs_yr10_neb = import_obs(u'neb', u'gpcp_v2.3_obs')

exp1_yr1_namz, exp1_yr2_namz, exp1_yr3_namz, exp1_yr4_namz, exp1_yr5_namz, exp1_yr6_namz, exp1_yr7_namz, exp1_yr8_namz, exp1_yr9_namz, exp1_yr10_namz = import_sim(u'namz', u'regcm_exp1')
exp1_yr1_samz, exp1_yr2_samz, exp1_yr3_samz, exp1_yr4_samz, exp1_yr5_samz, exp1_yr6_samz, exp1_yr7_samz, exp1_yr8_samz, exp1_yr9_samz, exp1_yr10_samz = import_sim(u'samz', u'regcm_exp1')
exp1_yr1_neb, exp1_yr2_neb, exp1_yr3_neb, exp1_yr4_neb, exp1_yr5_neb, exp1_yr6_neb, exp1_yr7_neb, exp1_yr8_neb, exp1_yr9_neb, exp1_yr10_neb = import_sim(u'neb', u'regcm_exp1')

exp2_yr1_namz, exp2_yr2_namz, exp2_yr3_namz, exp2_yr4_namz, exp2_yr5_namz, exp2_yr6_namz, exp2_yr7_namz, exp2_yr8_namz, exp2_yr9_namz, exp1_yr10_namz = import_sim(u'namz', u'regcm_exp2')
exp2_yr1_samz, exp2_yr2_samz, exp2_yr3_samz, exp2_yr4_samz, exp2_yr5_samz, exp2_yr6_samz, exp2_yr7_samz, exp2_yr8_samz, exp2_yr9_samz, exp2_yr10_samz = import_sim(u'samz', u'regcm_exp1')
exp2_yr1_neb, exp2_yr2_neb, exp2_yr3_neb, exp2_yr4_neb, exp2_yr5_neb, exp2_yr6_neb, exp2_yr7_neb, exp2_yr8_neb, exp2_yr9_neb, exp2_yr10_neb = import_sim(u'neb', u'regcm_exp1')

obs_namz = [obs_yr1_namz, obs_yr2_namz, obs_yr3_namz, obs_yr4_namz, obs_yr5_namz, obs_yr6_namz, obs_yr7_namz, obs_yr8_namz, obs_yr9_namz, obs_yr10_namz]
exp1_namz = [exp1_yr1_namz, exp1_yr2_namz, exp1_yr3_namz, exp1_yr4_namz, exp1_yr5_namz, exp1_yr6_namz, exp1_yr7_namz, exp1_yr8_namz, exp1_yr9_namz, exp1_yr10_namz]
exp2_namz = [exp2_yr1_namz, exp2_yr2_namz, exp2_yr3_namz, exp2_yr4_namz, exp2_yr5_namz, exp2_yr6_namz, exp2_yr7_namz, exp2_yr8_namz, exp2_yr9_namz, exp1_yr10_namz]

obs_samz = [obs_yr1_samz, obs_yr2_samz, obs_yr3_samz, obs_yr4_samz, obs_yr5_samz, obs_yr6_samz, obs_yr7_samz, obs_yr8_samz, obs_yr9_samz, obs_yr10_samz]
exp1_samz  = [exp1_yr1_samz, exp1_yr2_samz, exp1_yr3_samz, exp1_yr4_samz, exp1_yr5_samz, exp1_yr6_samz, exp1_yr7_samz, exp1_yr8_samz, exp1_yr9_samz, exp1_yr10_samz]
exp2_samz = [exp2_yr1_samz, exp2_yr2_samz, exp2_yr3_samz, exp2_yr4_samz, exp2_yr5_samz, exp2_yr6_samz, exp2_yr7_samz, exp2_yr8_samz, exp2_yr9_samz, exp2_yr10_samz]

obs_neb = [obs_yr1_neb, obs_yr2_neb, obs_yr3_neb, obs_yr4_neb, obs_yr5_neb, obs_yr6_neb, obs_yr7_neb, obs_yr8_neb, obs_yr9_neb, obs_yr10_neb]
exp1_neb  = [exp1_yr1_neb, exp1_yr2_neb, exp1_yr3_neb, exp1_yr4_neb, exp1_yr5_neb, exp1_yr6_neb, exp1_yr7_neb, exp1_yr8_neb, exp1_yr9_neb, exp1_yr10_neb]
exp2_neb = [exp2_yr1_neb, exp2_yr2_neb, exp2_yr3_neb, exp2_yr4_neb, exp2_yr5_neb, exp2_yr6_neb, exp2_yr7_neb, exp2_yr8_neb, exp2_yr9_neb, exp2_yr10_neb]

# Plot regcm exps and obs database 
fig = plt.figure()
time = np.arange(1, 11 + 1)

ax1 = fig.add_subplot(311)
boxplt1 = plt.violinplot(obs_namz, showmeans=False, showmedians=True)
boxplt2 = plt.violinplot(exp1_namz, showmeans=False, showmedians=True)
boxplt3 = plt.violinplot(exp2_namz, showmeans=False, showmedians=True)
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
ax1.set_yticks(np.arange(0, 20, 4))
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'))
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')

ax2 = fig.add_subplot(312)
boxplt1 = plt.violinplot(obs_samz)
boxplt2 = plt.violinplot(exp1_samz)
boxplt3 = plt.violinplot(exp2_samz)
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
ax2.set_ylabel('Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
ax2.set_yticks(np.arange(0, 20, 4))
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'))
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')

ax3 = fig.add_subplot(313)
boxplt1 = plt.violinplot(obs_neb)
boxplt2 = plt.violinplot(exp1_neb)
boxplt3 = plt.violinplot(exp2_neb)
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
ax3.set_xlabel('Years', fontsize=8, fontweight='bold')
ax3.set_yticks(np.arange(0, 20, 4))
plt.xticks(time, ('2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010'))
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/regcm_pbl/results'
name_out = 'pyplt_violin_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
