# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script compute cdf function from regcm46 and obs database"

import os
import netCDF4
import numpy as np
import numpy as np
import matplotlib.pyplot as plt

from comp_statist_indices import compute_cdf
from matplotlib.font_manager import FontProperties


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

	season_obs = value[2:120:3,:,:]
	djf_obs = np.nanmean(np.nanmean(season_obs[3:40:4], axis=1), axis=1)
	mam_obs = np.nanmean(np.nanmean(season_obs[0:40:4], axis=1), axis=1)
	jja_obs = np.nanmean(np.nanmean(season_obs[1:40:4], axis=1), axis=1)

	return djf_obs, mam_obs, jja_obs
	

def import_sim(area, exp):
	
	param = 'pr' # pr or tas
	date  = '2001-2010'

	path  = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq   = '{0}/{1}_{2}_{3}_mon_{4}.nc'.format(path, param, area, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_sim = value[2:120:3,:,:]
	djf_sim = np.nanmean(np.nanmean(season_sim[3:40:4], axis=1), axis=1)
	mam_sim = np.nanmean(np.nanmean(season_sim[0:40:4], axis=1), axis=1)
	jja_sim = np.nanmean(np.nanmean(season_sim[1:40:4], axis=1), axis=1)

	return djf_sim, mam_sim, jja_sim
	
	
# Import regcm exps and obs database 
djf_obs_namz, mam_obs_namz, jja_obs_namz = import_obs(u'namz', u'gpcp_v2.3_obs')
djf_exp1_namz, mam_exp1_namz, jja_exp1_namz = import_sim(u'namz', u'regcm_exp1')
djf_exp2_namz, mam_exp2_namz, jja_exp2_namz = import_sim(u'namz', u'regcm_exp2')

djf_obs_samz, mam_obs_samz, jja_obs_samz = import_obs(u'samz', u'gpcp_v2.3_obs')
djf_exp1_samz, mam_exp1_samz, jja_exp1_samz = import_sim(u'samz', u'regcm_exp1')
djf_exp2_samz, mam_exp2_samz, jja_exp2_samz = import_sim(u'samz', u'regcm_exp2')

djf_obs_neb, mam_obs_neb, jja_obs_neb = import_obs(u'neb', u'gpcp_v2.3_obs')
djf_exp1_neb, mam_exp1_neb, jja_exp1_neb = import_sim(u'neb', u'regcm_exp1')
djf_exp2_neb, mam_exp2_neb, jja_exp2_neb = import_sim(u'neb', u'regcm_exp2')

# Import cdf function
sortedtime_djf_obs_namz, cdf_djf_obs_namz = compute_cdf(djf_obs_namz)
sortedtime_djf_obs_samz, cdf_djf_obs_samz = compute_cdf(djf_obs_samz)
sortedtime_djf_obs_neb, cdf_djf_obs_neb = compute_cdf(djf_obs_neb)

sortedtime_mam_obs_namz, cdf_mam_obs_namz = compute_cdf(mam_obs_namz)
sortedtime_mam_obs_samz, cdf_mam_obs_samz = compute_cdf(mam_obs_samz)
sortedtime_mam_obs_neb, cdf_mam_obs_neb = compute_cdf(mam_obs_neb)

sortedtime_jja_obs_namz, cdf_jja_obs_namz = compute_cdf(jja_obs_namz)
sortedtime_jja_obs_samz, cdf_jja_obs_samz = compute_cdf(jja_obs_samz)
sortedtime_jja_obs_neb, cdf_jja_obs_neb = compute_cdf(jja_obs_neb)

sortedtime_djf_exp1_namz, cdf_djf_exp1_namz = compute_cdf(djf_exp1_namz)
sortedtime_djf_exp1_samz, cdf_djf_exp1_samz = compute_cdf(djf_exp1_samz)
sortedtime_djf_exp1_neb, cdf_djf_exp1_neb = compute_cdf(djf_exp1_neb)

sortedtime_mam_exp1_namz, cdf_mam_exp1_namz = compute_cdf(mam_exp1_namz)
sortedtime_mam_exp1_samz, cdf_mam_exp1_samz = compute_cdf(mam_exp1_samz)
sortedtime_mam_exp1_neb, cdf_mam_exp1_neb = compute_cdf(mam_exp1_neb)

sortedtime_jja_exp1_namz, cdf_jja_exp1_namz = compute_cdf(jja_exp1_namz)
sortedtime_jja_exp1_samz, cdf_jja_exp1_samz = compute_cdf(jja_exp1_samz)
sortedtime_jja_exp1_neb, cdf_jja_exp1_neb = compute_cdf(jja_exp1_neb)

sortedtime_djf_exp2_namz, cdf_djf_exp2_namz = compute_cdf(djf_exp2_namz)
sortedtime_djf_exp2_samz, cdf_djf_exp2_samz = compute_cdf(djf_exp2_samz)
sortedtime_djf_exp2_neb, cdf_djf_exp2_neb = compute_cdf(djf_exp2_neb)

sortedtime_mam_exp2_namz, cdf_mam_exp2_namz = compute_cdf(mam_exp2_namz)
sortedtime_mam_exp2_samz, cdf_mam_exp2_samz = compute_cdf(mam_exp2_samz)
sortedtime_mam_exp2_neb, cdf_mam_exp2_neb = compute_cdf(mam_exp2_neb)

sortedtime_jja_exp2_namz, cdf_jja_exp2_namz = compute_cdf(jja_exp2_namz)
sortedtime_jja_exp2_samz, cdf_jja_exp2_samz = compute_cdf(jja_exp2_samz)
sortedtime_jja_exp2_neb, cdf_jja_exp2_neb = compute_cdf(jja_exp2_neb)

# Plot regcm exps and obs database 
fig = plt.figure()

ax1 = fig.add_subplot(3, 3, 1)
plt.plot(sortedtime_djf_obs_namz, cdf_djf_obs_namz, color='black', label='GPCP', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_djf_exp1_namz, cdf_djf_exp1_namz, color='blue', label='Reg_Exp1', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_djf_exp2_namz, cdf_djf_exp2_namz,  color='red', label='Reg_Exp2', linestyle='-.', linewidth=1.5)
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.xticks(np.arange(0, 12, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = fig.add_subplot(3, 3, 2)
plt.plot(sortedtime_mam_obs_namz, cdf_mam_obs_namz, label='GPCP', color='black', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_mam_exp1_namz, cdf_mam_exp1_namz, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_mam_exp2_namz, cdf_mam_exp2_namz, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt.xticks(np.arange(0, 12, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)

ax3 = fig.add_subplot(3, 3, 3)
plt.plot(sortedtime_jja_obs_namz, cdf_jja_obs_namz,  label='GPCP', color='black', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_jja_exp1_namz, cdf_jja_exp1_namz, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_jja_exp2_namz, cdf_jja_exp2_namz, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.xticks(np.arange(0, 12, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)

ax4 = fig.add_subplot(3, 3, 4)
plt.plot(sortedtime_djf_obs_samz, cdf_djf_obs_samz, color='black', label='GPCP', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_djf_exp1_samz, cdf_djf_exp1_samz, color='blue', label='Reg_Exp1', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_djf_exp2_samz, cdf_djf_exp2_samz,  color='red', label='Reg_Exp2', linestyle='-.', linewidth=1.5)
ax4.xaxis.grid(True, which='major', linestyle='--')
ax4.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
plt.xticks(np.arange(0, 12, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.ylabel(u'CDF', fontsize=8, fontweight='bold')
plt.setp(ax4.get_xticklabels(), visible=False)

ax5 = fig.add_subplot(3, 3, 5)
plt.plot(sortedtime_mam_obs_samz, cdf_mam_obs_samz, label='GPCP', color='black', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_mam_exp1_samz, cdf_mam_exp1_samz, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_mam_exp2_samz, cdf_mam_exp2_samz, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
plt.xticks(np.arange(0, 12, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax5.get_xticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)

ax6 = fig.add_subplot(3, 3, 6)
plt.plot(sortedtime_jja_obs_samz, cdf_jja_obs_samz,  label='GPCP', color='black', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_jja_exp1_samz, cdf_jja_exp1_samz, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_jja_exp2_samz, cdf_jja_exp2_samz, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
ax6.xaxis.grid(True, which='major', linestyle='--')
ax6.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
plt.xticks(np.arange(0, 12, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax6.get_xticklabels(), visible=False)
plt.setp(ax6.get_yticklabels(), visible=False)

ax7 = fig.add_subplot(3, 3, 7)
plt.plot(sortedtime_djf_obs_neb, cdf_djf_obs_neb, color='black', label='GPCP', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_djf_exp1_neb, cdf_djf_exp1_neb, color='blue', label='Reg_Exp1', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_djf_exp2_neb, cdf_djf_exp2_neb,  color='red', label='Reg_Exp2', linestyle='-.', linewidth=1.5)
ax7.xaxis.grid(True, which='major', linestyle='--')
ax7.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
plt.xticks(np.arange(0, 12, 2))
plt.yticks(np.arange(0, 1.2, 0.2))

ax8 = fig.add_subplot(3, 3, 8)
plt.plot(sortedtime_mam_obs_neb, cdf_mam_obs_neb, label='GPCP', color='black', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_mam_exp1_neb, cdf_mam_exp1_neb, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_mam_exp2_neb, cdf_mam_exp2_neb, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
ax8.xaxis.grid(True, which='major', linestyle='--')
ax8.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
plt.xticks(np.arange(0, 12, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.xlabel(u'Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.setp(ax8.get_yticklabels(), visible=False)

ax9 = fig.add_subplot(3, 3, 9)
plt.plot(sortedtime_jja_obs_neb, cdf_jja_obs_neb,  label='GPCP', color='black', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_jja_exp1_neb, cdf_jja_exp1_neb, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_jja_exp2_neb, cdf_jja_exp2_neb, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
ax9.xaxis.grid(True, which='major', linestyle='--')
ax9.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
plt.xticks(np.arange(0, 12, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax9.get_yticklabels(), visible=False)
plt.legend(loc='lower right', shadow=True, ncol=1, prop=FontProperties(size=8))

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/regcm_pbl/results'
name_out = 'pyplt_cdf_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()

