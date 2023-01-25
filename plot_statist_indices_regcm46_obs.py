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

import scipy.stats as st
from comp_statist_indices import compute_corr, compute_rmse
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

# Compute correlation
corr_djf_exp1_namz = st.pearsonr(djf_obs_namz[2], djf_exp1_namz[2]) 
corr_djf_exp1_samz = st.pearsonr(djf_obs_samz[2], djf_exp1_samz[2]) 
corr_djf_exp1_neb  = st.pearsonr(djf_obs_neb[2], djf_exp1_neb[2])
corr_djf_exp2_namz = st.pearsonr(djf_obs_namz[2], djf_exp2_namz[2]) 
corr_djf_exp2_samz = st.pearsonr(djf_obs_samz[2], djf_exp2_samz[2]) 
corr_djf_exp2_neb  = st.pearsonr(djf_obs_neb[2], djf_exp2_neb[2]) 

corr_jja_exp1_namz = st.pearsonr(jja_obs_namz[2], jja_exp1_namz[2]) 
corr_jja_exp1_samz = st.pearsonr(jja_obs_samz[2], jja_exp1_samz[2]) 
corr_jja_exp1_neb  = st.pearsonr(jja_obs_neb[2], jja_exp1_neb[2])
corr_jja_exp2_namz = st.pearsonr(jja_obs_namz[2], jja_exp2_namz[2]) 
corr_jja_exp2_samz = st.pearsonr(jja_obs_samz[2], jja_exp2_samz[2]) 
corr_jja_exp2_neb  = st.pearsonr(jja_obs_neb[2], jja_exp2_neb[2]) 

corr_ann_exp1_namz = st.pearsonr(ann_obs_namz[2], ann_exp1_namz[2]) 
corr_ann_exp1_samz = st.pearsonr(ann_obs_samz[2], ann_exp1_samz[2]) 
corr_ann_exp1_neb  = st.pearsonr(ann_obs_neb[2], ann_exp1_neb[2])
corr_ann_exp2_namz = st.pearsonr(ann_obs_namz[2], ann_exp2_namz[2]) 
corr_ann_exp2_samz = st.pearsonr(ann_obs_samz[2], ann_exp2_samz[2]) 
corr_ann_exp2_neb  = st.pearsonr(ann_obs_neb[2], ann_exp2_neb[2]) 

print(st.pearsonr(djf_obs_namz[2], djf_exp1_namz[2])) 
print(st.pearsonr(djf_obs_samz[2], djf_exp1_samz[2])) 
print(st.pearsonr(djf_obs_neb[2], djf_exp1_neb[2]))
print(st.pearsonr(djf_obs_namz[2], djf_exp2_namz[2])) 
print(st.pearsonr(djf_obs_samz[2], djf_exp2_samz[2])) 
print(st.pearsonr(djf_obs_neb[2], djf_exp2_neb[2])) 

print(st.pearsonr(jja_obs_namz[2], jja_exp1_namz[2])) 
print(st.pearsonr(jja_obs_samz[2], jja_exp1_samz[2])) 
print(st.pearsonr(jja_obs_neb[2], jja_exp1_neb[2]))
print(st.pearsonr(jja_obs_namz[2], jja_exp2_namz[2])) 
print(st.pearsonr(jja_obs_samz[2], jja_exp2_samz[2])) 
print(st.pearsonr(jja_obs_neb[2], jja_exp2_neb[2])) 

print(st.pearsonr(ann_obs_namz[2], ann_exp1_namz[2])) 
print(st.pearsonr(ann_obs_samz[2], ann_exp1_samz[2]))
print(st.pearsonr(ann_obs_neb[2], ann_exp1_neb[2]))
print(st.pearsonr(ann_obs_namz[2], ann_exp2_namz[2])) 
print(st.pearsonr(ann_obs_samz[2], ann_exp2_samz[2]))
print(st.pearsonr(ann_obs_neb[2], ann_exp2_neb[2]))

# Compute rmse
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
fig = plt.figure(figsize=(8,3))

ax = fig.add_subplot(1, 3, 1)
plt.axvspan(0, 3, ymin = 0.50, ymax = 1, color='gray', alpha=0.4)
plt.axvspan(3, 6, ymin = 0.0, ymax = 0.50, color='gray', alpha=0.4)
plt.plot(rmse_djf_exp1_namz, corr_djf_exp1_namz[0],  marker='o', markersize=10, markeredgecolor='black', color='red')
plt.plot(rmse_djf_exp2_namz, corr_djf_exp2_namz[0],  marker='o', markersize=10, markeredgecolor='black', color='blue')
plt.plot(rmse_jja_exp1_namz, corr_jja_exp1_namz[0],  marker='s', markersize=10, markeredgecolor='black', color='red')
plt.plot(rmse_jja_exp2_namz, corr_jja_exp2_namz[0],  marker='s', markersize=10, markeredgecolor='black', color='blue')
plt.plot(rmse_ann_exp1_namz, corr_ann_exp1_namz[0],  marker='^', markersize=10, markeredgecolor='black', color='red')
plt.plot(rmse_ann_exp2_namz, corr_ann_exp2_namz[0],  marker='^', markersize=10, markeredgecolor='black', color='blue')
plt.title(u'A) NAMZ', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'RMSE (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.ylabel(u'PCC', fontsize=8, fontweight='bold')
plt.xlim(0, 6)
plt.ylim(-1, 1)
plt.xticks(np.arange(0, 7, 1), fontsize=8)
plt.yticks(np.arange(-1, 1.2, 0.2), fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.axvline(3, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(1, 3, 2)
plt.axvspan(0, 3, ymin = 0.50, ymax = 1, color='gray', alpha=0.4)
plt.axvspan(3, 6, ymin = 0.0, ymax = 0.50, color='gray', alpha=0.4)
plt.plot(rmse_djf_exp1_samz, corr_djf_exp1_samz[0],  marker='o', markersize=10, markeredgecolor='black', color='red')
plt.plot(rmse_djf_exp2_samz, corr_djf_exp2_samz[0],  marker='o', markersize=10, markeredgecolor='black', color='blue')
plt.plot(rmse_jja_exp1_samz, corr_jja_exp1_samz[0],  marker='s', markersize=10, markeredgecolor='black', color='red')
plt.plot(rmse_jja_exp2_samz, corr_jja_exp2_samz[0],  marker='s', markersize=10, markeredgecolor='black', color='blue')
plt.plot(rmse_ann_exp1_samz, corr_ann_exp1_samz[0],  marker='^', markersize=10, markeredgecolor='black', color='red')
plt.plot(rmse_ann_exp2_samz, corr_ann_exp2_samz[0],  marker='^', markersize=10, markeredgecolor='black', color='blue')
plt.title(u'B) SAMZ', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'RMSE (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.xlim(0, 6)
plt.ylim(-1, 1)
plt.xticks(np.arange(0, 7, 1), fontsize=8)
plt.yticks(np.arange(-1, 1.2, 0.2), fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.axvline(3, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(1, 3, 3)
plt.axvspan(0, 3, ymin = 0.50, ymax = 1, color='gray', alpha=0.4)
plt.axvspan(3, 6, ymin = 0.0, ymax = 0.50, color='gray', alpha=0.4)
plt.plot(rmse_djf_exp1_neb, corr_djf_exp1_neb[0],  marker='o', markersize=10, markeredgecolor='black', color='red')
plt.plot(rmse_djf_exp2_neb, corr_djf_exp2_neb[0],  marker='o', markersize=10, markeredgecolor='black', color='blue')
plt.plot(rmse_jja_exp1_neb, corr_jja_exp1_neb[0],  marker='s', markersize=10, markeredgecolor='black', color='red')
plt.plot(rmse_jja_exp2_neb, corr_jja_exp2_neb[0],  marker='s', markersize=10, markeredgecolor='black', color='blue')
plt.plot(rmse_ann_exp1_neb, corr_ann_exp1_neb[0],  marker='^', markersize=10, markeredgecolor='black', color='red')
plt.plot(rmse_ann_exp2_neb, corr_ann_exp2_neb[0],  marker='^', markersize=10, markeredgecolor='black', color='blue')
plt.title(u'C) NEB', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'RMSE (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.xlim(0, 6)
plt.ylim(-1, 1)
plt.xticks(np.arange(0, 7, 1), fontsize=8)
plt.yticks(np.arange(-1, 1.2, 0.2), fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.axvline(3, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

Reg_H_DJF  = mlines.Line2D([], [], color='red',  marker='o', linestyle='None', markersize=10, markeredgecolor='black', label='DJF Reg_Holtslag')
Reg_UW_DJF = mlines.Line2D([], [], color='blue', marker='o', linestyle='None', markersize=10, markeredgecolor='black', label='DJF Reg_UW-PBL')
Reg_H_JJA  = mlines.Line2D([], [], color='red',  marker='s', linestyle='None', markersize=10, markeredgecolor='black', label='JJA Reg_Holtslag')
Reg_UW_JJA = mlines.Line2D([], [], color='blue', marker='s', linestyle='None', markersize=10, markeredgecolor='black', label='JJA Reg_UW-PBL')
Reg_H_ANN  = mlines.Line2D([], [], color='red',  marker='^', linestyle='None', markersize=10, markeredgecolor='black', label='ANN Reg_Holtslag')
Reg_UW_ANN = mlines.Line2D([], [], color='blue', marker='^', linestyle='None', markersize=10, markeredgecolor='black', label='ANN Reg_UW-PBL')
plt.legend(loc=8, fontsize=7, shadow=True, handles=[Reg_H_DJF, Reg_UW_DJF, Reg_H_JJA, Reg_UW_JJA, Reg_H_ANN, Reg_UW_ANN])
         	
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_statist_indices_pr_regcm_pbl_obs_2001-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
