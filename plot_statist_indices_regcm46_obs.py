# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot taylor diagram from regcm46 and obs database"

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from comp_statist_indices import compute_corr, compute_rmse
from matplotlib.font_manager import FontProperties


def import_obs(area, obs, time):
	
	param = 'precip' 
	date  = '2001-2010'

	path  = '/home/nice/Documents/dataset/obs/reg_pbl'
	arq   = '{0}/{1}_{2}_{3}_{4}_{5}.nc'.format(path, param, area, obs, time, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean_obs = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return lat, lon, mean_obs
	
	
def import_sim(area, exp, time):
	
	param = 'pr' 
	date  = '2001-2010'

	path  = '/home/nice/Documents/dataset/rcm/reg_pbl'
	arq   = '{0}/{1}_{2}_{3}_{4}_{5}.nc'.format(path, param, area, exp, time, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean_sim = np.nanmean(np.nanmean(value, axis=1), axis=1)

	return lat, lon, mean_sim
	
		
# Import regcm exps and obs database 
djf_obs_namz = import_obs(u'namz', u'gpcp_v2.3_obs', 'djf')
djf_exp1_namz = import_sim(u'namz', u'regcm_exp1', 'djf')
djf_exp2_namz = import_sim(u'namz', u'regcm_exp2', 'djf')

jja_obs_namz = import_obs(u'namz', u'gpcp_v2.3_obs', 'jja')
jja_exp1_namz = import_sim(u'namz', u'regcm_exp1', 'jja')
jja_exp2_namz = import_sim(u'namz', u'regcm_exp2', 'jja')

ann_obs_namz = import_obs(u'namz', u'gpcp_v2.3_obs', 'ann')
ann_exp1_namz = import_sim(u'namz', u'regcm_exp1', 'ann')
ann_exp2_namz = import_sim(u'namz', u'regcm_exp2', 'ann')

djf_obs_samz = import_obs(u'samz', u'gpcp_v2.3_obs', 'djf')
djf_exp1_samz = import_sim(u'samz', u'regcm_exp1', 'djf')
djf_exp2_samz = import_sim(u'samz', u'regcm_exp2', 'djf')

jja_obs_samz = import_obs(u'samz', u'gpcp_v2.3_obs', 'jja')
jja_exp1_samz = import_sim(u'samz', u'regcm_exp1', 'jja')
jja_exp2_samz = import_sim(u'samz', u'regcm_exp2', 'jja')

ann_obs_samz = import_obs(u'samz', u'gpcp_v2.3_obs', 'ann')
ann_exp1_samz = import_sim(u'samz', u'regcm_exp1', 'ann')
ann_exp2_samz = import_sim(u'samz', u'regcm_exp2', 'ann')

djf_obs_neb = import_obs(u'neb', u'gpcp_v2.3_obs', 'djf')
djf_exp1_neb = import_sim(u'neb', u'regcm_exp1', 'djf')
djf_exp2_neb = import_sim(u'neb', u'regcm_exp2', 'djf')

jja_obs_neb = import_obs(u'neb', u'gpcp_v2.3_obs', 'jja')
jja_exp1_neb = import_sim(u'neb', u'regcm_exp1', 'jja')
jja_exp2_neb = import_sim(u'neb', u'regcm_exp2', 'jja')

ann_obs_neb = import_obs(u'neb', u'gpcp_v2.3_obs', 'ann')
ann_exp1_neb = import_sim(u'neb', u'regcm_exp1', 'ann')
ann_exp2_neb = import_sim(u'neb', u'regcm_exp2', 'ann')

# Compute correlation
corr_djf_exp1_namz = compute_corr(djf_obs_namz[2], djf_exp1_namz[2]) 
corr_djf_exp1_samz = compute_corr(djf_obs_samz[2], djf_exp1_samz[2]) 
corr_djf_exp1_neb  = compute_corr(djf_obs_neb[2], djf_exp1_neb[2])
corr_djf_exp2_namz = compute_corr(djf_obs_namz[2], djf_exp2_namz[2]) 
corr_djf_exp2_samz = compute_corr(djf_obs_samz[2], djf_exp2_samz[2]) 
corr_djf_exp2_neb  = compute_corr(djf_obs_neb[2], djf_exp2_neb[2]) 

corr_jja_exp1_namz = compute_corr(jja_obs_namz[2], jja_exp1_namz[2]) 
corr_jja_exp1_samz = compute_corr(jja_obs_samz[2], jja_exp1_samz[2]) 
corr_jja_exp1_neb  = compute_corr(jja_obs_neb[2], jja_exp1_neb[2])
corr_jja_exp2_namz = compute_corr(jja_obs_namz[2], jja_exp2_namz[2]) 
corr_jja_exp2_samz = compute_corr(jja_obs_samz[2], jja_exp2_samz[2]) 
corr_jja_exp2_neb  = compute_corr(jja_obs_neb[2], jja_exp2_neb[2]) 

corr_ann_exp1_namz = compute_corr(ann_obs_namz[2], ann_exp1_namz[2]) 
corr_ann_exp1_samz = compute_corr(ann_obs_samz[2], ann_exp1_samz[2]) 
corr_ann_exp1_neb  = compute_corr(ann_obs_neb[2], ann_exp1_neb[2])
corr_ann_exp2_namz = compute_corr(ann_obs_namz[2], ann_exp2_namz[2]) 
corr_ann_exp2_samz = compute_corr(ann_obs_samz[2], ann_exp2_samz[2]) 
corr_ann_exp2_neb  = compute_corr(ann_obs_neb[2], ann_exp2_neb[2]) 

# Compute correlation
rmse_djf_exp1_namz = compute_rmse(djf_exp1_namz[2], djf_obs_namz[2]) 
rmse_djf_exp1_samz = compute_rmse(djf_exp1_samz[2], djf_obs_samz[2]) 
rmse_djf_exp1_neb  = compute_rmse(djf_exp1_neb[2], djf_obs_neb[2])
rmse_djf_exp2_namz = compute_rmse(djf_exp2_namz[2], djf_obs_namz[2]) 
rmse_djf_exp2_samz = compute_rmse(djf_exp2_samz[2], djf_obs_samz[2]) 
rmse_djf_exp2_neb  = compute_rmse(djf_exp2_neb[2], djf_obs_neb[2]) 

rmse_jja_exp1_namz = compute_rmse(jja_exp1_namz[2], jja_obs_namz[2]) 
rmse_jja_exp1_samz = compute_rmse(jja_exp1_samz[2], jja_obs_samz[2]) 
rmse_jja_exp1_neb  = compute_rmse(jja_exp1_neb[2], jja_obs_neb[2]) 
rmse_jja_exp2_namz = compute_rmse(jja_exp2_namz[2], jja_obs_namz[2]) 
rmse_jja_exp2_samz = compute_rmse(jja_exp2_samz[2], jja_obs_samz[2]) 
rmse_jja_exp2_neb  = compute_rmse(jja_exp2_neb[2], jja_obs_neb[2])

rmse_ann_exp1_namz = compute_rmse(ann_exp1_namz[2], ann_obs_namz[2]) 
rmse_ann_exp1_samz = compute_rmse(ann_exp1_samz[2], ann_obs_samz[2]) 
rmse_ann_exp1_neb  = compute_rmse(ann_exp1_neb[2], ann_obs_neb[2]) 
rmse_ann_exp2_namz = compute_rmse(ann_exp2_namz[2], ann_obs_namz[2]) 
rmse_ann_exp2_samz = compute_rmse(ann_exp2_samz[2], ann_obs_samz[2]) 
rmse_ann_exp2_neb  = compute_rmse(ann_exp2_neb[2], ann_obs_neb[2])

# Plot models and obs database 
fig = plt.figure()

ax = fig.add_subplot(2, 2, 1)
plt.plot(corr_djf_exp1_namz, rmse_djf_exp1_namz,  marker=r"$\mathcircled{1}$", markersize=10, color='red')
plt.plot(corr_djf_exp2_namz, rmse_djf_exp2_namz,  marker=r"$\mathcircled{1}$", markersize=10, color='blue')
plt.plot(corr_jja_exp1_namz, rmse_jja_exp1_namz,  marker=r"$\mathcircled{2}$", markersize=10, color='red')
plt.plot(corr_jja_exp2_namz, rmse_jja_exp2_namz,  marker=r"$\mathcircled{2}$", markersize=10, color='blue')
plt.plot(corr_ann_exp1_namz, rmse_ann_exp1_namz,  marker=r"$\mathcircled{3}$", markersize=10, color='red')
plt.plot(corr_ann_exp2_namz, rmse_ann_exp2_namz,  marker=r"$\mathcircled{3}$", markersize=10, color='blue')
plt.title(u'A) NAMZ', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'RMSE', fontsize=8, fontweight='bold')
plt.xlim(-1, 1)
plt.ylim(0, 6)
plt.yticks(np.arange(0, 8, 1), fontsize=8)
plt.xticks(np.arange(-1, 1.2, 0.2), fontsize=8)
plt.axvline(0, linewidth=1., linestyle='--', color='black')
plt.axhline(2, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(2, 2, 2)
plt.plot(corr_djf_exp1_samz, rmse_djf_exp1_samz,  marker=r"$\mathcircled{1}$", markersize=10, color='red')
plt.plot(corr_djf_exp2_samz, rmse_djf_exp2_samz,  marker=r"$\mathcircled{1}$", markersize=10, color='blue')
plt.plot(corr_jja_exp1_samz, rmse_jja_exp1_samz,  marker=r"$\mathcircled{2}$", markersize=10, color='red')
plt.plot(corr_jja_exp2_samz, rmse_jja_exp2_samz,  marker=r"$\mathcircled{2}$", markersize=10, color='blue')
plt.plot(corr_ann_exp1_samz, rmse_ann_exp1_samz,  marker=r"$\mathcircled{3}$", markersize=10, color='red')
plt.plot(corr_ann_exp2_samz, rmse_ann_exp2_samz,  marker=r"$\mathcircled{3}$", markersize=10, color='blue')
plt.title(u'B) SAMZ', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'PCC', fontsize=8, fontweight='bold')
plt.xlim(-1, 1)
plt.ylim(0, 6)
plt.yticks(np.arange(0, 8, 1), fontsize=8)
plt.xticks(np.arange(-1, 1.2, 0.2), fontsize=8)
plt.axvline(0, linewidth=1., linestyle='--', color='black')
plt.axhline(2, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(2, 2, 3)
plt.plot(corr_djf_exp1_neb, rmse_djf_exp1_neb,  marker=r"$\mathcircled{1}$", markersize=10, color='red')
plt.plot(corr_djf_exp2_neb, rmse_djf_exp2_neb,  marker=r"$\mathcircled{1}$", markersize=10, color='blue')
plt.plot(corr_jja_exp1_neb, rmse_jja_exp1_neb,  marker=r"$\mathcircled{2}$", markersize=10, color='red')
plt.plot(corr_jja_exp2_neb, rmse_jja_exp2_neb,  marker=r"$\mathcircled{2}$", markersize=10, color='blue')
plt.plot(corr_ann_exp1_neb, rmse_ann_exp1_neb,  marker=r"$\mathcircled{3}$", markersize=10, color='red')
plt.plot(corr_ann_exp2_neb, rmse_ann_exp2_neb,  marker=r"$\mathcircled{3}$", markersize=10, color='blue')
plt.title(u'C) NEB', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'PCC', fontsize=8, fontweight='bold')
plt.ylabel(u'RMSE', fontsize=8, fontweight='bold')
plt.xlim(-1, 1)
plt.ylim(0, 6)
plt.yticks(np.arange(0, 8, 1), fontsize=8)
plt.xticks(np.arange(-1, 1.2, 0.2), fontsize=8)
plt.axvline(0, linewidth=1., linestyle='--', color='black')
plt.axhline(2, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

Reg_H_DJF  = mlines.Line2D([], [], color='red',  marker=r"$\mathcircled{1}$", linestyle='None', markersize=10, label='DJF Reg_H')
Reg_UW_DJF = mlines.Line2D([], [], color='blue', marker=r"$\mathcircled{1}$", linestyle='None', markersize=10, label='DJF Reg_UW')
Reg_H_MAM  = mlines.Line2D([], [], color='red',  marker=r"$\mathcircled{2}$", linestyle='None', markersize=10, label='MAM Reg_H')
Reg_UW_MAM = mlines.Line2D([], [], color='blue', marker=r"$\mathcircled{2}$", linestyle='None', markersize=10, label='MAM Reg_UW')
Reg_H_JJA  = mlines.Line2D([], [], color='red',  marker=r"$\mathcircled{3}$", linestyle='None', markersize=10, label='JJA Reg_H')
Reg_UW_JJA = mlines.Line2D([], [], color='blue', marker=r"$\mathcircled{3}$", linestyle='None', markersize=10, label='JJA Reg_UW')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, shadow=True, handles=[Reg_H_DJF, Reg_UW_DJF, Reg_H_MAM, Reg_UW_MAM, Reg_H_JJA, Reg_UW_JJA])
         	
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_statist_indices_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
