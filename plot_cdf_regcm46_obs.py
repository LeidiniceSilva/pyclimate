# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script compute cdf function from regcm46 and obs database"

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

from comp_statist_indices import compute_cdf
from matplotlib.font_manager import FontProperties


def import_obs(area, obs, time):
	
	param = 'pr' 
	date  = '2001-2005'

	path  = '/home/nice/Documentos/dataset/obs/reg_pbl'
	arq   = '{0}/{1}_{2}_{3}_{4}_{5}.nc'.format(path, param, area, obs, time, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables['cmorph'][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean_obs = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return lat, lon, mean_obs
	
	
def import_sim(area, exp, time):
	
	param = 'pr' 
	date  = '2001-2005'

	path  = '/home/nice/Documentos/dataset/rcm/reg_pbl'
	arq   = '{0}/{1}_{2}_{3}_{4}_{5}.nc'.format(path, param, area, exp, time, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean_sim = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return lat, lon, mean_sim
	
		
# Import regcm exps and obs database 
djf_obs_namz = import_obs(u'namz', u'cmorph_v1_obs', 'djf')
djf_exp1_namz = import_sim(u'namz', u'regcm_exp1', 'djf')
djf_exp2_namz = import_sim(u'namz', u'regcm_exp2', 'djf')

jja_obs_namz = import_obs(u'namz', u'cmorph_v1_obs', 'jja')
jja_exp1_namz = import_sim(u'namz', u'regcm_exp1', 'jja')
jja_exp2_namz = import_sim(u'namz', u'regcm_exp2', 'jja')

ann_obs_namz = import_obs(u'namz', u'cmorph_v1_obs', 'ann')
ann_exp1_namz = import_sim(u'namz', u'regcm_exp1', 'ann')
ann_exp2_namz = import_sim(u'namz', u'regcm_exp2', 'ann')

djf_obs_samz = import_obs(u'samz', u'cmorph_v1_obs', 'djf')
djf_exp1_samz = import_sim(u'samz', u'regcm_exp1', 'djf')
djf_exp2_samz = import_sim(u'samz', u'regcm_exp2', 'djf')

jja_obs_samz = import_obs(u'samz', u'cmorph_v1_obs', 'jja')
jja_exp1_samz = import_sim(u'samz', u'regcm_exp1', 'jja')
jja_exp2_samz = import_sim(u'samz', u'regcm_exp2', 'jja')

ann_obs_samz = import_obs(u'samz', u'cmorph_v1_obs', 'ann')
ann_exp1_samz = import_sim(u'samz', u'regcm_exp1', 'ann')
ann_exp2_samz = import_sim(u'samz', u'regcm_exp2', 'ann')

djf_obs_neb = import_obs(u'neb', u'cmorph_v1_obs', 'djf')
djf_exp1_neb = import_sim(u'neb', u'regcm_exp1', 'djf')
djf_exp2_neb = import_sim(u'neb', u'regcm_exp2', 'djf')

jja_obs_neb = import_obs(u'neb', u'cmorph_v1_obs', 'jja')
jja_exp1_neb = import_sim(u'neb', u'regcm_exp1', 'jja')
jja_exp2_neb = import_sim(u'neb', u'regcm_exp2', 'jja')

ann_obs_neb = import_obs(u'neb', u'cmorph_v1_obs', 'ann')
ann_exp1_neb = import_sim(u'neb', u'regcm_exp1', 'ann')
ann_exp2_neb = import_sim(u'neb', u'regcm_exp2', 'ann')

# Compute bias
bias_djf_exp1_namz = djf_exp1_namz[2] - djf_obs_namz[2] 
bias_djf_exp1_samz = djf_exp1_samz[2] - djf_obs_samz[2] 
bias_djf_exp1_neb = djf_exp1_neb[2] - djf_obs_neb[2] 
bias_djf_exp2_namz = djf_exp2_namz[2] - djf_obs_namz[2] 
bias_djf_exp2_samz = djf_exp2_samz[2] - djf_obs_samz[2] 
bias_djf_exp2_neb = djf_exp2_neb[2] - djf_obs_neb[2] 

bias_jja_exp1_namz = jja_exp1_namz[2] - jja_obs_namz[2] 
bias_jja_exp1_samz = jja_exp1_samz[2] - jja_obs_samz[2] 
bias_jja_exp1_neb = jja_exp1_neb[2] - jja_obs_neb[2] 
bias_jja_exp2_namz = jja_exp2_namz[2] - jja_obs_namz[2] 
bias_jja_exp2_samz = jja_exp2_samz[2] - jja_obs_samz[2] 
bias_jja_exp2_neb = jja_exp2_neb[2] - jja_obs_neb[2]

bias_ann_exp1_namz = ann_exp1_namz[2] - ann_obs_namz[2] 
bias_ann_exp1_samz = ann_exp1_samz[2] - ann_obs_samz[2] 
bias_ann_exp1_neb = ann_exp1_neb[2] - ann_obs_neb[2] 
bias_ann_exp2_namz = ann_exp2_namz[2] - ann_obs_namz[2] 
bias_ann_exp2_samz = ann_exp2_samz[2] - ann_obs_samz[2] 
bias_ann_exp2_neb = ann_exp2_neb[2] - ann_obs_neb[2]

# Import cdf function
sortedtime_djf_exp1_namz, cdf_djf_exp1_namz = compute_cdf(bias_djf_exp1_namz)
sortedtime_djf_exp1_samz, cdf_djf_exp1_samz = compute_cdf(bias_djf_exp1_samz)
sortedtime_djf_exp1_neb, cdf_djf_exp1_neb = compute_cdf(bias_djf_exp1_neb)

sortedtime_jja_exp1_namz, cdf_jja_exp1_namz = compute_cdf(bias_jja_exp1_namz)
sortedtime_jja_exp1_samz, cdf_jja_exp1_samz = compute_cdf(bias_jja_exp1_samz)
sortedtime_jja_exp1_neb, cdf_jja_exp1_neb = compute_cdf(bias_jja_exp1_neb)

sortedtime_ann_exp1_namz, cdf_ann_exp1_namz = compute_cdf(bias_ann_exp1_namz)
sortedtime_ann_exp1_samz, cdf_ann_exp1_samz = compute_cdf(bias_ann_exp1_samz)
sortedtime_ann_exp1_neb, cdf_ann_exp1_neb = compute_cdf(bias_ann_exp1_neb)

sortedtime_djf_exp2_namz, cdf_djf_exp2_namz = compute_cdf(bias_djf_exp2_namz)
sortedtime_djf_exp2_samz, cdf_djf_exp2_samz = compute_cdf(bias_djf_exp2_samz)
sortedtime_djf_exp2_neb, cdf_djf_exp2_neb = compute_cdf(bias_djf_exp2_neb)

sortedtime_jja_exp2_namz, cdf_jja_exp2_namz = compute_cdf(bias_jja_exp2_namz)
sortedtime_jja_exp2_samz, cdf_jja_exp2_samz = compute_cdf(bias_jja_exp2_samz)
sortedtime_jja_exp2_neb, cdf_jja_exp2_neb = compute_cdf(bias_jja_exp2_neb)

sortedtime_ann_exp2_namz, cdf_ann_exp2_namz = compute_cdf(bias_ann_exp2_namz)
sortedtime_ann_exp2_samz, cdf_ann_exp2_samz = compute_cdf(bias_ann_exp2_samz)
sortedtime_ann_exp2_neb, cdf_ann_exp2_neb = compute_cdf(bias_ann_exp2_neb)

# Plot regcm exps and obs database 
fig = plt.figure()

ax1 = fig.add_subplot(3, 3, 1)
plt.axvspan(-6, 0, color='red', alpha=0.2)
plt.axvspan(0, 6, color='blue', alpha=0.2)
plt.plot(sortedtime_djf_exp1_namz, cdf_djf_exp1_namz, color='blue', label='Reg_Holtslag', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_djf_exp2_namz, cdf_djf_exp2_namz,  color='red', label='Reg_UW-PBL', linestyle='-.', linewidth=1.5)
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'A) NAMZ DJF', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-6, 3)
plt.ylim(0, 1)
plt.xticks(np.arange(-6, 4, 1), fontsize=8)
plt.yticks(np.arange(0, 1.2, 0.2), fontsize=8)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')
plt.text(-5.8, 0.2, '20%', color='black', fontsize=8)
plt.text(-5.8, 0.8, '80%', color='black', fontsize=8)
plt.text(-3, 0.8, 'Dry', color='black', fontsize=8)
plt.text(1, 0.8, 'Wet', color='black', fontsize=8)
plt.legend(loc=5, handlelength=0.50, handleheight=0.50, shadow=True, ncol=1, prop=FontProperties(size=6))
plt.setp(ax1.get_xticklabels(), visible=False)	

ax2 = fig.add_subplot(3, 3, 2)
plt.axvspan(-6, 0, color='red', alpha=0.2)
plt.axvspan(0, 6, color='blue', alpha=0.2)
plt.plot(sortedtime_jja_exp1_namz, cdf_jja_exp1_namz, label='Reg_Holtslag', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_jja_exp2_namz, cdf_jja_exp2_namz, label='Reg_UW-PBL', color='red', linestyle='-.', linewidth=1.5)
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'B) NAMZ JJA', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-6, 3)
plt.ylim(0, 1)
plt.xticks(np.arange(-6, 4, 1), fontsize=8)
plt.yticks(np.arange(0, 1.2, 0.2), fontsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')
plt.setp(ax2.get_xticklabels(), visible=False)	

ax3 = fig.add_subplot(3, 3, 3)
plt.axvspan(-6, 0, color='red', alpha=0.2)
plt.axvspan(0, 6, color='blue', alpha=0.2)
plt.plot(sortedtime_ann_exp1_namz, cdf_ann_exp1_namz, label='Reg_Holtslag', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_ann_exp2_namz, cdf_ann_exp2_namz, label='Reg_UW-PBL', color='red', linestyle='-.', linewidth=1.5)
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'C) NAMZ ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-6, 3)
plt.ylim(0, 1)
plt.xticks(np.arange(-6, 4, 1), fontsize=8)
plt.yticks(np.arange(0, 1.2, 0.2), fontsize=8)
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')

ax4 = fig.add_subplot(3, 3, 4)
plt.axvspan(-6, 0, color='red', alpha=0.2)
plt.axvspan(0, 6, color='blue', alpha=0.2)
plt.plot(sortedtime_djf_exp1_samz, cdf_djf_exp1_samz, color='blue', label='Reg_Holtslag', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_djf_exp2_samz, cdf_djf_exp2_samz,  color='red', label='Reg_UW-PBL', linestyle='-.', linewidth=1.5)
ax4.xaxis.grid(True, which='major', linestyle='--')
ax4.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'D) SAMZ DJF', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-6, 3)
plt.ylim(0, 1)
plt.xticks(np.arange(-6, 4, 1), fontsize=8)
plt.yticks(np.arange(0, 1.2, 0.2), fontsize=8)
plt.ylabel(u'CDF', fontsize=8, fontweight='bold')
plt.setp(ax4.get_xticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')
plt.setp(ax4.get_xticklabels(), visible=False)	

ax5 = fig.add_subplot(3, 3, 5)
plt.axvspan(-6, 0, color='red', alpha=0.2)
plt.axvspan(0, 6, color='blue', alpha=0.2)
plt.plot(sortedtime_jja_exp1_samz, cdf_jja_exp1_samz, label='Reg_Holtslag', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_jja_exp2_samz, cdf_jja_exp2_samz, label='Reg_UW-PBL', color='red', linestyle='-.', linewidth=1.5)
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'E) SAMZ JJA', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-6, 3)
plt.ylim(0, 1)
plt.xticks(np.arange(-6, 4, 1), fontsize=8)
plt.yticks(np.arange(0, 1.2, 0.2), fontsize=8)
plt.setp(ax5.get_xticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')
plt.setp(ax5.get_xticklabels(), visible=False)	

ax6 = fig.add_subplot(3, 3, 6)
plt.axvspan(-6, 0, color='red', alpha=0.2)
plt.axvspan(0, 6, color='blue', alpha=0.2)
plt.plot(sortedtime_ann_exp1_samz, cdf_ann_exp1_samz, label='Reg_Holtslag', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_ann_exp2_samz, cdf_ann_exp2_samz, label='Reg_UW-PBL', color='red', linestyle='-.', linewidth=1.5)
ax6.xaxis.grid(True, which='major', linestyle='--')
ax6.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'F) SAMZ ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-6, 3)
plt.ylim(0, 1)
plt.xticks(np.arange(-6, 4, 1), fontsize=8)
plt.yticks(np.arange(0, 1.2, 0.2), fontsize=8)
plt.setp(ax6.get_xticklabels(), visible=False)
plt.setp(ax6.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')

ax7 = fig.add_subplot(3, 3, 7)
plt.axvspan(-6, 0, color='red', alpha=0.2)
plt.axvspan(0, 6, color='blue', alpha=0.2)
plt.plot(sortedtime_djf_exp1_neb, cdf_djf_exp1_neb, color='blue', label='Reg_Holtslag', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_djf_exp2_neb, cdf_djf_exp2_neb,  color='red', label='Reg_UW-PBL', linestyle='-.', linewidth=1.5)
ax7.xaxis.grid(True, which='major', linestyle='--')
ax7.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'G) NEB DJF', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-6, 3)
plt.ylim(0, 1)
plt.xticks(np.arange(-6, 4, 1), fontsize=8)
plt.yticks(np.arange(0, 1.2, 0.2), fontsize=8)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')

ax8 = fig.add_subplot(3, 3, 8)
plt.axvspan(-6, 0, color='red', alpha=0.2)
plt.axvspan(0, 6, color='blue', alpha=0.2)
plt.plot(sortedtime_jja_exp1_neb, cdf_jja_exp1_neb, label='Reg_Holtslag', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_jja_exp2_neb, cdf_jja_exp2_neb, label='Reg_UW-PBL', color='red', linestyle='-.', linewidth=1.5)
ax8.xaxis.grid(True, which='major', linestyle='--')
ax8.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'H) NEB JJA', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-6, 3)
plt.ylim(0, 1)
plt.xticks(np.arange(-6, 4, 1), fontsize=8)
plt.yticks(np.arange(0, 1.2, 0.2), fontsize=8)
plt.xlabel(u'Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.setp(ax8.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')

ax9 = fig.add_subplot(3, 3, 9)
plt.axvspan(-6, 0, color='red', alpha=0.2)
plt.axvspan(0, 6, color='blue', alpha=0.2)
plt.plot(sortedtime_ann_exp1_neb, cdf_ann_exp1_neb, label='Reg_Holtslag', color='blue', linestyle='-.', linewidth=1.5)
plt.plot(sortedtime_ann_exp2_neb, cdf_ann_exp2_neb, label='Reg_UW-PBL', color='red', linestyle='-.', linewidth=1.5)
ax9.xaxis.grid(True, which='major', linestyle='--')
ax9.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'I) NEB ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlim(-6, 3)
plt.ylim(0, 1)
plt.xticks(np.arange(-6, 4, 1), fontsize=8)
plt.yticks(np.arange(0, 1.2, 0.2), fontsize=8)
plt.setp(ax9.get_yticklabels(), visible=False)
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0.8, linewidth=1., linestyle='dashed', color='black')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_cdf_pr_regcm_pbl_obs_2001-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
