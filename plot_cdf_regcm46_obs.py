# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script compute CDF functions from Reg and Had models"

import os
import netCDF4
import numpy as np
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.font_manager import FontProperties


def cdf_function(data):
	
	sortedtime = np.sort(data)
	cdf = 1. * np.arange(len(data))/(len(data) - 1)
	
	return sortedtime, cdf


def import_sim(area, exp):
	
	param = 'pr' # pr or tas
	date  = '2001-2010'

	path  = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq   = '{0}/{1}_{2}_{3}_mon_{4}.nc'.format(path, param, area, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	sim = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	return sim


def import_obs(area, obs):
	
	param = 'precip' # precip, pre or tmp
	date  = '2001-2010'

	path  = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq   = '{0}/{1}_{2}_{3}_mon_{4}.nc'.format(path, param, area, obs, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	obs = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1) 
	
	return obs
	
	
# Import regcm exps model end obs database climatology
nam_exp1 = import_sim(u'namz', u'regcm_exp1')
sam_exp1 = import_sim(u'samz', u'regcm_exp1')
neb_exp1 = import_sim(u'neb', u'regcm_exp1')

nam_exp2 = import_sim(u'namz', u'regcm_exp2')
sam_exp2 = import_sim(u'samz', u'regcm_exp2')
neb_exp2 = import_sim(u'neb', u'regcm_exp2')

nam_obs = import_obs(u'namz', u'gpcp_v2.2_obs')
sam_obs = import_obs(u'samz', u'gpcp_v2.2_obs')
neb_obs = import_obs(u'neb', u'gpcp_v2.2_obs')

# Import cdf function
sortedtime_nam_exp1, cdf_nam_exp1 = cdf_function(nam_exp1)
sortedtime_nam_exp2, cdf_nam_exp2 = cdf_function(nam_exp2)
sortedtime_nam_obs, cdf_nam_obs = cdf_function(nam_obs)
sortedtime_sam_exp1, cdf_sam_exp1 = cdf_function(sam_exp1)
sortedtime_sam_exp2, cdf_sam_exp2 = cdf_function(sam_exp2)
sortedtime_sam_obs, cdf_sam_obs = cdf_function(sam_obs)
sortedtime_neb_exp1, cdf_neb_exp1 = cdf_function(neb_exp1)
sortedtime_neb_exp2, cdf_neb_exp2 = cdf_function(neb_exp2)
sortedtime_neb_obs, cdf_neb_obs = cdf_function(neb_obs)

# Plot model end obs data climatology
fig = plt.figure()

ax1 = fig.add_subplot(1, 3, 1)
a=plt.plot(sortedtime_nam_exp1, cdf_nam_exp1, color='blue', label='Reg_Exp1', linestyle='-.', linewidth=1.5)
b=plt.plot(sortedtime_nam_exp2, cdf_nam_exp2,  color='red', label='Reg_Exp2', linestyle='-.', linewidth=1.5)
c=plt.plot(sortedtime_nam_obs, cdf_nam_obs, color='black', label='GPCP', linestyle='-.', linewidth=1.5)
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'A) NAMZ', fontweight='bold')
plt.xticks(np.arange(0, 14, 2))
plt.yticks(np.arange(0, 1.1, 0.1))
plt.ylabel(u'CDF', fontweight='bold')

ax2 = fig.add_subplot(1, 3, 2)
plt.plot(sortedtime_sam_exp1, cdf_sam_exp1, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_sam_exp2, cdf_sam_exp2, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_sam_obs, cdf_sam_obs, label='GPCP', color='black', linestyle='-.', linewidth=1.5)
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'B) SAMZ', fontweight='bold')
plt.xticks(np.arange(0, 14, 2))
plt.yticks(np.arange(0, 1.1, 0.1))
plt.xlabel(u'Precipitação (mm d⁻¹)', fontweight='bold')

ax5 = fig.add_subplot(1, 3, 3)
plt.plot(sortedtime_neb_exp1, cdf_neb_exp1, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_neb_exp2, cdf_neb_exp2, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_neb_obs, cdf_neb_obs,  label='GPCP', color='black', linestyle='-.', linewidth=1.5)
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'C) NEB', fontweight='bold')
plt.xticks(np.arange(0, 14, 2))
plt.yticks(np.arange(0, 1.1, 0.1))
plt.legend(loc='lower right', shadow=True, ncol=1, prop=FontProperties(size=8))

fig.tight_layout()

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/regcm_pbl/results'
name_out = 'pyplt_cdf_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()

