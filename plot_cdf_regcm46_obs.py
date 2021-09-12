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

	path  = '/home/nice/Documents/ufrn/phd_project/papers/paper_rcm_pbl/datas'
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

	path  = '/home/nice/Documents/ufrn/phd_project/papers/paper_rcm_pbl/datas'
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

# Compute bias
bias_djf_exp1_namz = djf_exp1_namz - djf_obs_namz 
bias_djf_exp1_samz = djf_exp1_samz - djf_obs_samz 
bias_djf_exp1_neb = djf_exp1_neb - djf_obs_neb 
bias_djf_exp2_namz = djf_exp2_namz - djf_obs_namz 
bias_djf_exp2_samz = djf_exp2_samz - djf_obs_samz 
bias_djf_exp2_neb = djf_exp2_neb - djf_obs_neb 

bias_mam_exp1_namz = mam_exp1_namz - mam_obs_namz 
bias_mam_exp1_samz = mam_exp1_samz - mam_obs_samz 
bias_mam_exp1_neb = mam_exp1_neb - mam_obs_neb 
bias_mam_exp2_namz = mam_exp2_namz - mam_obs_namz 
bias_mam_exp2_samz = mam_exp2_samz - mam_obs_samz 
bias_mam_exp2_neb = mam_exp2_neb - mam_obs_neb

bias_jja_exp1_namz = jja_exp1_namz - jja_obs_namz 
bias_jja_exp1_samz = jja_exp1_samz - jja_obs_samz 
bias_jja_exp1_neb = jja_exp1_neb - jja_obs_neb 
bias_jja_exp2_namz = jja_exp2_namz - jja_obs_namz 
bias_jja_exp2_samz = jja_exp2_samz - jja_obs_samz 
bias_jja_exp2_neb = jja_exp2_neb - jja_obs_neb

# Import cdf function
sortedtime_djf_exp1_namz, cdf_djf_exp1_namz = compute_cdf(bias_djf_exp1_namz)
sortedtime_djf_exp1_samz, cdf_djf_exp1_samz = compute_cdf(bias_djf_exp1_samz)
sortedtime_djf_exp1_neb, cdf_djf_exp1_neb = compute_cdf(bias_djf_exp1_neb)

sortedtime_mam_exp1_namz, cdf_mam_exp1_namz = compute_cdf(bias_mam_exp1_namz)
sortedtime_mam_exp1_samz, cdf_mam_exp1_samz = compute_cdf(bias_mam_exp1_samz)
sortedtime_mam_exp1_neb, cdf_mam_exp1_neb = compute_cdf(bias_mam_exp1_neb)

sortedtime_jja_exp1_namz, cdf_jja_exp1_namz = compute_cdf(bias_jja_exp1_namz)
sortedtime_jja_exp1_samz, cdf_jja_exp1_samz = compute_cdf(bias_jja_exp1_samz)
sortedtime_jja_exp1_neb, cdf_jja_exp1_neb = compute_cdf(bias_jja_exp1_neb)

sortedtime_djf_exp2_namz, cdf_djf_exp2_namz = compute_cdf(bias_djf_exp2_namz)
sortedtime_djf_exp2_samz, cdf_djf_exp2_samz = compute_cdf(bias_djf_exp2_samz)
sortedtime_djf_exp2_neb, cdf_djf_exp2_neb = compute_cdf(bias_djf_exp2_neb)

sortedtime_mam_exp2_namz, cdf_mam_exp2_namz = compute_cdf(bias_mam_exp2_namz)
sortedtime_mam_exp2_samz, cdf_mam_exp2_samz = compute_cdf(bias_mam_exp2_samz)
sortedtime_mam_exp2_neb, cdf_mam_exp2_neb = compute_cdf(bias_mam_exp2_neb)

sortedtime_jja_exp2_namz, cdf_jja_exp2_namz = compute_cdf(jja_exp2_namz)
sortedtime_jja_exp2_samz, cdf_jja_exp2_samz = compute_cdf(jja_exp2_samz)
sortedtime_jja_exp2_neb, cdf_jja_exp2_neb = compute_cdf(jja_exp2_neb)

# Plot regcm exps and obs database 
fig = plt.figure()

ax1 = fig.add_subplot(3, 3, 1)
plt.axvspan(-8, 0, color='red', alpha=0.2)
plt.axvspan(0, 4, color='blue', alpha=0.2)
plt.plot(sortedtime_djf_exp1_namz, cdf_djf_exp1_namz, color='blue', label='Reg_Exp1', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_djf_exp2_namz, cdf_djf_exp2_namz,  color='red', label='Reg_Exp2', linestyle='-.', linewidth=1.5)
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-8, 4)
plt.ylim(0, 1)
plt.xticks(np.arange(-8, 6, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax1.get_xticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')
plt.text(-7, 0.2, '20%', color='black', fontsize=8)
plt.text(-7, 0.8, '80%', color='black', fontsize=8)
plt.text(-3, 0.8, 'Dry', color='black', fontsize=8)
plt.text(1, 0.8, 'Wet', color='black', fontsize=8)
plt.legend(loc=6, handlelength=0.50, handleheight=0.50, shadow=True, ncol=1, prop=FontProperties(size=6))

ax2 = fig.add_subplot(3, 3, 2)
plt.axvspan(-8, 0, color='red', alpha=0.2)
plt.axvspan(0, 4, color='blue', alpha=0.2)
plt.plot(sortedtime_mam_exp1_namz, cdf_mam_exp1_namz, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_mam_exp2_namz, cdf_mam_exp2_namz, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-8, 4)
plt.ylim(0, 1)
plt.xticks(np.arange(-8, 6, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')

ax3 = fig.add_subplot(3, 3, 3)
plt.axvspan(-8, 0, color='red', alpha=0.2)
plt.axvspan(0, 4, color='blue', alpha=0.2)
plt.plot(sortedtime_jja_exp1_namz, cdf_jja_exp1_namz, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_jja_exp2_namz, cdf_jja_exp2_namz, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-8, 4)
plt.ylim(0, 1)
plt.xticks(np.arange(-8, 6, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')

ax4 = fig.add_subplot(3, 3, 4)
plt.axvspan(-8, 0, color='red', alpha=0.2)
plt.axvspan(0, 4, color='blue', alpha=0.2)
plt.plot(sortedtime_djf_exp1_samz, cdf_djf_exp1_samz, color='blue', label='Reg_Exp1', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_djf_exp2_samz, cdf_djf_exp2_samz,  color='red', label='Reg_Exp2', linestyle='-.', linewidth=1.5)
ax4.xaxis.grid(True, which='major', linestyle='--')
ax4.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-8, 4)
plt.ylim(0, 1)
plt.xticks(np.arange(-8, 6, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.ylabel(u'CDF', fontsize=8, fontweight='bold')
plt.setp(ax4.get_xticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')

ax5 = fig.add_subplot(3, 3, 5)
plt.axvspan(-8, 0, color='red', alpha=0.2)
plt.axvspan(0, 4, color='blue', alpha=0.2)
plt.plot(sortedtime_mam_exp1_samz, cdf_mam_exp1_samz, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_mam_exp2_samz, cdf_mam_exp2_samz, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-8, 4)
plt.ylim(0, 1)
plt.xticks(np.arange(-8, 6, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax5.get_xticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')

ax6 = fig.add_subplot(3, 3, 6)
plt.axvspan(-8, 0, color='red', alpha=0.2)
plt.axvspan(0, 4, color='blue', alpha=0.2)
plt.plot(sortedtime_jja_exp1_samz, cdf_jja_exp1_samz, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_jja_exp2_samz, cdf_jja_exp2_samz, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
ax6.xaxis.grid(True, which='major', linestyle='--')
ax6.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-8, 4)
plt.ylim(0, 1)
plt.xticks(np.arange(-8, 6, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax6.get_xticklabels(), visible=False)
plt.setp(ax6.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')

ax7 = fig.add_subplot(3, 3, 7)
plt.axvspan(-8, 0, color='red', alpha=0.2)
plt.axvspan(0, 4, color='blue', alpha=0.2)
plt.plot(sortedtime_djf_exp1_neb, cdf_djf_exp1_neb, color='blue', label='Reg_Exp1', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_djf_exp2_neb, cdf_djf_exp2_neb,  color='red', label='Reg_Exp2', linestyle='-.', linewidth=1.5)
ax7.xaxis.grid(True, which='major', linestyle='--')
ax7.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-8, 4)
plt.ylim(0, 1)
plt.xticks(np.arange(-8, 6, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')

ax8 = fig.add_subplot(3, 3, 8)
plt.axvspan(-8, 0, color='red', alpha=0.2)
plt.axvspan(0, 4, color='blue', alpha=0.2)
plt.plot(sortedtime_mam_exp1_neb, cdf_mam_exp1_neb, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_mam_exp2_neb, cdf_mam_exp2_neb, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
ax8.xaxis.grid(True, which='major', linestyle='--')
ax8.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-8, 4)
plt.ylim(0, 1)
plt.xticks(np.arange(-8, 6, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.xlabel(u'Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.setp(ax8.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')

ax9 = fig.add_subplot(3, 3, 9)
plt.axvspan(-8, 0, color='red', alpha=0.2)
plt.axvspan(0, 4, color='blue', alpha=0.2)
plt.plot(sortedtime_jja_exp1_neb, cdf_jja_exp1_neb, label='Reg_Exp1', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_jja_exp2_neb, cdf_jja_exp2_neb, label='Reg_Exp2', color='red', linestyle='-.', linewidth=1.5)
ax9.xaxis.grid(True, which='major', linestyle='--')
ax9.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-8, 4)
plt.ylim(0, 1)
plt.xticks(np.arange(-8, 6, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax9.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/phd_project/papers/paper_rcm_pbl/figs'
name_out = 'pyplt_cdf_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
