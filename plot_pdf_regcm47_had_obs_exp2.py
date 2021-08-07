# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/05/2021"
__description__ = "This script plot pdf function from regcm47 and hadgem models and obs database"

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
from scipy.stats import norm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from comp_statist_indices import compute_pdf


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/rcm_exp1'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	ann_obs = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1) 
	season_obs = np.nanmean(np.nanmean(value[2:239:3,:,:], axis=1), axis=1)
	djf_obs = season_obs[3:80:4]
	jja_obs = season_obs[1:80:4]

	return ann_obs, djf_obs, jja_obs
	
	
def import_rcm_exp1(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp1/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	ann_rcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1),axis=1) 
	season_rcm = np.nanmean(np.nanmean(value[2:239:3,:,:], axis=1), axis=1)
	djf_rcm = season_rcm[3:80:4]
	jja_rcm = season_rcm[1:80:4]

	return ann_rcm, djf_rcm, jja_rcm


def import_rcm_exp2(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp2/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	ann_rcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1),axis=1) 
	season_rcm = np.nanmean(np.nanmean(value[2:239:3,:,:], axis=1), axis=1)
	djf_rcm = season_rcm[3:80:4]
	jja_rcm = season_rcm[1:80:4]

	return ann_rcm, djf_rcm, jja_rcm
	
	
def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/rcm_exp1/hist'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	ann_gcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1) 
	season_gcm = np.nanmean(np.nanmean(value[2:239:3,:,:], axis=1), axis=1)
	djf_gcm = season_gcm[3:80:4]
	jja_gcm = season_gcm[1:80:4]

	return ann_gcm, djf_gcm, jja_gcm


def compute_bias(sim, obs):
	
	x = len(sim)
	
	bias = []
	for mon in range(0, x):
		diff = sim[mon] - obs[mon]
		bias.append(diff)
	
	return bias
	
	              
# Import models and obs database 
# Precipitation
pre_ann_cru_samz, pre_djf_cru_samz, pre_jja_cru_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
pre_ann_rcm_exp1_samz, pre_djf_rcm_exp1_samz, pre_jja_rcm_exp1_samz = import_rcm_exp1('pr', 'samz', 'hist', '1986-2005')
pre_ann_rcm_exp2_samz, pre_djf_rcm_exp2_samz, pre_jja_rcm_exp2_samz = import_rcm_exp2('pr', 'samz', 'hist', '1986-2005')
pre_ann_gcm_samz, pre_djf_gcm_samz, pre_jja_gcm_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')

pre_ann_cru_eneb, pre_djf_cru_eneb, pre_jja_cru_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
pre_ann_rcm_exp1_eneb, pre_djf_rcm_exp1_eneb, pre_jja_rcm_exp1_eneb = import_rcm_exp1('pr', 'eneb', 'hist', '1986-2005')
pre_ann_rcm_exp2_eneb, pre_djf_rcm_exp2_eneb, pre_jja_rcm_exp2_eneb = import_rcm_exp2('pr', 'eneb', 'hist', '1986-2005')
pre_ann_gcm_eneb, pre_djf_gcm_eneb, pre_jja_gcm_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')

pre_ann_cru_matopiba, pre_djf_cru_matopiba, pre_jja_cru_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')
pre_ann_rcm_exp1_matopiba, pre_djf_rcm_exp1_matopiba, pre_jja_rcm_exp1_matopiba = import_rcm_exp1('pr', 'matopiba', 'hist', '1986-2005')
pre_ann_rcm_exp2_matopiba, pre_djf_rcm_exp2_matopiba, pre_jja_rcm_exp2_matopiba = import_rcm_exp2('pr', 'matopiba', 'hist', '1986-2005')
pre_ann_gcm_matopiba, pre_djf_gcm_matopiba, pre_jja_gcm_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

# Temperature
tas_ann_cru_samz, tas_djf_cru_samz, tas_jja_cru_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
tas_ann_rcm_exp1_samz, tas_djf_rcm_exp1_samz, tas_jja_rcm_exp1_samz = import_rcm_exp1('tas', 'samz', 'hist', '1986-2005')
tas_ann_rcm_exp2_samz, tas_djf_rcm_exp2_samz, tas_jja_rcm_exp2_samz = import_rcm_exp2('tas', 'samz', 'hist', '1986-2005')
tas_ann_gcm_samz, tas_djf_gcm_samz, tas_jja_gcm_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')

tas_ann_cru_eneb, tas_djf_cru_eneb, tas_jja_cru_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
tas_ann_rcm_exp1_eneb, tas_djf_rcm_exp1_eneb, tas_jja_rcm_exp1_eneb = import_rcm_exp1('tas', 'eneb', 'hist', '1986-2005')
tas_ann_rcm_exp2_eneb, tas_djf_rcm_exp2_eneb, tas_jja_rcm_exp2_eneb = import_rcm_exp2('tas', 'eneb', 'hist', '1986-2005')
tas_ann_gcm_eneb, tas_djf_gcm_eneb, tas_jja_gcm_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')

tas_ann_cru_matopiba, tas_djf_cru_matopiba, tas_jja_cru_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')
tas_ann_rcm_exp1_matopiba, tas_djf_rcm_exp1_matopiba, tas_jja_rcm_exp1_matopiba = import_rcm_exp1('tas', 'matopiba', 'hist', '1986-2005')
tas_ann_rcm_exp2_matopiba, tas_djf_rcm_exp2_matopiba, tas_jja_rcm_exp2_matopiba = import_rcm_exp2('tas', 'matopiba', 'hist', '1986-2005')
tas_ann_gcm_matopiba, tas_djf_gcm_matopiba, tas_jja_gcm_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

# Compute bias
# Precipitation
bias_pre_djf_rcm_exp1_samz = compute_bias(pre_djf_rcm_exp1_samz, pre_djf_cru_samz)
bias_pre_djf_rcm_exp2_samz = compute_bias(pre_djf_rcm_exp2_samz, pre_djf_cru_samz)
bias_pre_djf_gcm_samz = compute_bias(pre_djf_gcm_samz, pre_djf_cru_samz)

bias_pre_jja_rcm_exp1_samz = compute_bias(pre_jja_rcm_exp1_samz, pre_jja_cru_samz)
bias_pre_jja_rcm_exp2_samz = compute_bias(pre_jja_rcm_exp2_samz, pre_jja_cru_samz)
bias_pre_jja_gcm_samz = compute_bias(pre_jja_gcm_samz, pre_jja_cru_samz)

bias_pre_ann_rcm_exp1_samz = compute_bias(pre_ann_rcm_exp1_samz, pre_ann_cru_samz)
bias_pre_ann_rcm_exp2_samz = compute_bias(pre_ann_rcm_exp2_samz, pre_ann_cru_samz)
bias_pre_ann_gcm_samz = compute_bias(pre_ann_gcm_samz, pre_ann_cru_samz)

bias_pre_djf_rcm_exp1_eneb = compute_bias(pre_djf_rcm_exp1_eneb, pre_djf_cru_eneb)
bias_pre_djf_rcm_exp2_eneb = compute_bias(pre_djf_rcm_exp2_eneb, pre_djf_cru_eneb)
bias_pre_djf_gcm_eneb = compute_bias(pre_djf_gcm_eneb, pre_djf_cru_eneb)

bias_pre_jja_rcm_exp1_eneb = compute_bias(pre_jja_rcm_exp1_eneb, pre_jja_cru_samz)
bias_pre_jja_rcm_exp2_eneb = compute_bias(pre_jja_rcm_exp2_eneb, pre_jja_cru_samz)
bias_pre_jja_gcm_eneb = compute_bias(pre_jja_gcm_eneb, pre_jja_cru_samz)

bias_pre_ann_rcm_exp1_eneb = compute_bias(pre_ann_rcm_exp1_eneb, pre_ann_cru_eneb)
bias_pre_ann_rcm_exp2_eneb = compute_bias(pre_ann_rcm_exp2_eneb, pre_ann_cru_eneb)
bias_pre_ann_gcm_eneb = compute_bias(pre_ann_gcm_eneb, pre_ann_cru_eneb)

bias_pre_djf_rcm_exp1_matopiba = compute_bias(pre_djf_rcm_exp1_matopiba, pre_djf_cru_matopiba)
bias_pre_djf_rcm_exp2_matopiba = compute_bias(pre_djf_rcm_exp2_matopiba, pre_djf_cru_matopiba)
bias_pre_djf_gcm_matopiba = compute_bias(pre_djf_gcm_matopiba, pre_djf_cru_eneb)

bias_pre_jja_rcm_exp1_matopiba = compute_bias(pre_jja_rcm_exp1_matopiba, pre_jja_cru_matopiba)
bias_pre_jja_rcm_exp2_matopiba = compute_bias(pre_jja_rcm_exp2_matopiba, pre_jja_cru_matopiba)
bias_pre_jja_gcm_matopiba = compute_bias(pre_jja_gcm_matopiba, pre_jja_cru_matopiba)

bias_pre_ann_rcm_exp1_matopiba = compute_bias(pre_ann_rcm_exp1_matopiba, pre_ann_cru_matopiba)
bias_pre_ann_rcm_exp2_matopiba = compute_bias(pre_ann_rcm_exp2_matopiba, pre_ann_cru_matopiba)
bias_pre_ann_gcm_matopiba = compute_bias(pre_ann_gcm_matopiba, pre_ann_cru_matopiba)

# Temperature
#~ bias_tas_reg_exp1_cru_samz = compute_bias(np.nanmean(mon_tas_reg_exp1_samz, axis=1), mon_tas_cru_samz)
#~ bias_tas_reg_exp1_cru_eneb = compute_bias(np.nanmean(mon_tas_reg_exp1_eneb, axis=1), mon_tas_cru_eneb)
#~ bias_tas_reg_exp1_cru_matopiba = compute_bias(np.nanmean(mon_tas_reg_exp1_matopiba, axis=1), mon_tas_cru_matopiba)

#~ bias_tas_reg_exp2_cru_samz = compute_bias(np.nanmean(mon_tas_reg_exp2_samz, axis=1), mon_tas_cru_samz)
#~ bias_tas_reg_exp2_cru_eneb = compute_bias(np.nanmean(mon_tas_reg_exp2_eneb, axis=1), mon_tas_cru_eneb)
#~ bias_tas_reg_exp2_cru_matopiba = compute_bias(np.nanmean(mon_tas_reg_exp2_matopiba, axis=1), mon_tas_cru_matopiba)

#~ bias_tas_had_cru_samz = compute_bias(mon_tas_had_samz, mon_tas_cru_samz)
#~ bias_tas_had_cru_eneb = compute_bias(mon_tas_had_eneb, mon_tas_cru_eneb)
#~ bias_tas_had_cru_matopiba = compute_bias(mon_tas_had_matopiba, mon_tas_cru_matopiba)

# Calculate PDF function
# Precipitation
xpdf_pre_djf_rcm_exp1_samz, pdf_pre_djf_rcm_exp1_samz = compute_pdf(bias_pre_djf_rcm_exp1_samz)
xpdf_pre_djf_rcm_exp2_samz, pdf_pre_djf_rcm_exp2_samz = compute_pdf(bias_pre_djf_rcm_exp2_samz)
xpdf_pre_djf_gcm_samz, pdf_pre_djf_gcm_samz = compute_pdf(bias_pre_djf_gcm_samz)

xpdf_pre_djf_rcm_exp1_eneb, pdf_pre_djf_rcm_exp1_eneb = compute_pdf(bias_pre_djf_rcm_exp1_eneb)
xpdf_pre_djf_rcm_exp2_eneb, pdf_pre_djf_rcm_exp2_eneb = compute_pdf(bias_pre_djf_rcm_exp2_eneb)
xpdf_pre_djf_gcm_eneb, pdf_pre_djf_gcm_eneb = compute_pdf(bias_pre_djf_gcm_eneb)

xpdf_pre_djf_rcm_exp1_matopiba, pdf_pre_djf_rcm_exp1_matopiba = compute_pdf(bias_pre_djf_rcm_exp1_matopiba)
xpdf_pre_djf_rcm_exp2_matopiba, pdf_pre_djf_rcm_exp2_matopiba = compute_pdf(bias_pre_djf_rcm_exp2_matopiba)
xpdf_pre_djf_gcm_matopiba, pdf_pre_djf_gcm_matopiba = compute_pdf(bias_pre_djf_gcm_matopiba)

# Temperature
#~ xpdf_tas_cru_samz, pdf_tas_cru_samz = compute_pdf(mon_tas_cru_samz)
#~ xpdf_tas_cru_eneb, pdf_tas_cru_eneb = compute_pdf(mon_tas_cru_eneb)
#~ xpdf_tas_cru_matopiba, pdf_tas_cru_matopiba = compute_pdf(mon_tas_cru_matopiba)

#~ xpdf_tas_reg_samz, pdf_tas_reg_samz = compute_pdf(np.nanmean(mon_tas_reg_samz, axis=0))
#~ xpdf_tas_reg_eneb, pdf_tas_reg_eneb = compute_pdf(np.nanmean(mon_tas_reg_eneb, axis=0))
#~ xpdf_tas_reg_matopiba, pdf_tas_reg_matopiba = compute_pdf(np.nanmean(mon_tas_reg_matopiba, axis=0))

#~ xpdf_tas_had_samz, pdf_tas_had_samz = compute_pdf(mon_tas_had_samz)
#~ xpdf_tas_had_eneb, pdf_tas_had_eneb = compute_pdf(mon_tas_had_eneb)
#~ xpdf_tas_had_matopiba, pdf_tas_had_matopiba = compute_pdf(mon_tas_had_matopiba)

# Plot models and obs database 
fig = plt.figure()

ax1 = fig.add_subplot(3, 2, 1)
pdf_line1 = ax1.plot(xpdf_pre_djf_rcm_exp1_samz, pdf_pre_djf_rcm_exp1_samz, color='blue', label='RegCM4.7_EXP1', linestyle='-',linewidth=1.5) 
pdf_line2 = ax1.plot(xpdf_pre_djf_rcm_exp2_samz, pdf_pre_djf_rcm_exp2_samz, color='red', label='RegCM4.7_EXP2', linestyle='-', linewidth=1.5) 
pdf_line3 = ax1.plot(xpdf_pre_djf_gcm_samz, pdf_pre_djf_gcm_samz, color='gray', label='HadGEM2-ES', linestyle='-', linewidth=1.5) 
plt.axvline(0, color='black', linestyle='-', linewidth=1.5)
plt.axvline(-2, color='black', linestyle='-', linewidth=1.5)
plt.axvline(2, color='black', linestyle='-', linewidth=1.5)
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
ax1.set_xlim(-5, 5)
ax1.set_ylim(0, 0.6)
ax1.set_xticks(np.arange(-5, 6, 1))
ax1.set_yticks(np.arange(0, 0.8, 0.2))
ax1.xaxis.set_tick_params(labelsize=8)
ax1.yaxis.set_tick_params(labelsize=8)
plt.setp(ax1.get_xticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')
plt.legend(fontsize=6, loc=1, shadow=True, ncol=1)

ax2 = fig.add_subplot(3, 2, 2)
pdf_line1 = ax2.plot(xpdf_pre_djf_rcm_exp1_samz, pdf_pre_djf_rcm_exp1_samz, color='blue', label='RegCM4.7_EXP1', linestyle='-',linewidth=1.5) 
pdf_line2 = ax2.plot(xpdf_pre_djf_rcm_exp2_samz, pdf_pre_djf_rcm_exp2_samz, color='red', label='RegCM4.7_EXP2', linestyle='-', linewidth=1.5) 
pdf_line3 = ax2.plot(xpdf_pre_djf_gcm_samz, pdf_pre_djf_gcm_samz, color='gray', label='HadGEM2-ES', linestyle='-', linewidth=1.5) 
plt.axvline(0, color='black', linestyle='-', linewidth=1.5)
plt.axvline(-2, color='black', linestyle='-', linewidth=1.5)
plt.axvline(2, color='black', linestyle='-', linewidth=1.5)
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
ax2.set_xlim(-5, 5)
ax2.set_ylim(0, 0.6)
ax2.set_xticks(np.arange(-5, 6, 1))
ax2.set_yticks(np.arange(0, 0.8, 0.2))
ax2.xaxis.set_tick_params(labelsize=8)
ax2.yaxis.set_tick_params(labelsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)	
plt.setp(ax2.get_yticklabels(), visible=False)		     	     
plt.grid(True, which='major', linestyle='--')
	 
ax3 = fig.add_subplot(3, 2, 3)
pdf_line1 = ax3.plot(xpdf_pre_djf_rcm_exp1_eneb, pdf_pre_djf_rcm_exp1_eneb, color='blue', label='RegCM4.7_EXP1', linestyle='-',linewidth=1.5) 
pdf_line2 = ax3.plot(xpdf_pre_djf_rcm_exp2_eneb, pdf_pre_djf_rcm_exp2_eneb, color='red', label='RegCM4.7_EXP2', linestyle='-', linewidth=1.5) 
pdf_line3 = ax3.plot(xpdf_pre_djf_gcm_eneb, pdf_pre_djf_gcm_eneb, color='gray', label='HadGEM2-ES', linestyle='-', linewidth=1.5) 
plt.axvline(0, color='black', linestyle='-', linewidth=1.5)
plt.axvline(-2, color='black', linestyle='-', linewidth=1.5)
plt.axvline(2, color='black', linestyle='-', linewidth=1.5)
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
ax3.set_xlim(-5, 5)
ax3.set_ylim(0, 0.6)
ax3.set_xticks(np.arange(-5, 6, 1))
ax3.set_yticks(np.arange(0, 0.8, 0.2))
ax3.xaxis.set_tick_params(labelsize=8)
ax3.yaxis.set_tick_params(labelsize=8)
plt.setp(ax3.get_xticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')
	 
ax4 = fig.add_subplot(3, 2, 4)
pdf_line1 = ax4.plot(xpdf_pre_djf_rcm_exp1_samz, pdf_pre_djf_rcm_exp1_samz, color='blue', label='RegCM4.7_EXP1', linestyle='-',linewidth=1.5) 
pdf_line2 = ax4.plot(xpdf_pre_djf_rcm_exp2_samz, pdf_pre_djf_rcm_exp2_samz, color='red', label='RegCM4.7_EXP2', linestyle='-', linewidth=1.5) 
pdf_line3 = ax4.plot(xpdf_pre_djf_gcm_samz, pdf_pre_djf_gcm_samz, color='gray', label='HadGEM2-ES', linestyle='-', linewidth=1.5) 
plt.axvline(0, color='black', linestyle='-', linewidth=1.5)
plt.axvline(-2, color='black', linestyle='-', linewidth=1.5)
plt.axvline(2, color='black', linestyle='-', linewidth=1.5)
plt.title(u'E)', loc='left', fontweight='bold', fontsize=8)
ax4.set_xlim(-5, 5)
ax4.set_ylim(0, 0.6)
ax4.set_xticks(np.arange(-5, 6, 1))
ax4.set_yticks(np.arange(0, 0.8, 0.2))
ax4.xaxis.set_tick_params(labelsize=8)
ax4.yaxis.set_tick_params(labelsize=8)
plt.setp(ax4.get_xticklabels(), visible=False)		
plt.setp(ax4.get_yticklabels(), visible=False)		          
plt.grid(True, which='major', linestyle='--')
		 
ax5 = fig.add_subplot(3, 2, 5)
pdf_line1 = ax5.plot(xpdf_pre_djf_rcm_exp1_matopiba, pdf_pre_djf_rcm_exp1_matopiba, color='blue', label='RegCM4.7_EXP1', linestyle='-',linewidth=1.5) 
pdf_line2 = ax5.plot(xpdf_pre_djf_rcm_exp2_matopiba, pdf_pre_djf_rcm_exp2_matopiba, color='red', label='RegCM4.7_EXP2', linestyle='-', linewidth=1.5) 
pdf_line3 = ax5.plot(xpdf_pre_djf_gcm_matopiba, pdf_pre_djf_gcm_matopiba, color='gray', label='HadGEM2-ES', linestyle='-', linewidth=1.5) 
plt.axvline(0, color='black', linestyle='-', linewidth=1.5)
plt.axvline(-2, color='black', linestyle='-', linewidth=1.5)
plt.axvline(2, color='black', linestyle='-', linewidth=1.5)
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
ax5.set_xlim(-5, 5)
ax5.set_ylim(0, 0.6)
ax5.set_xticks(np.arange(-5, 6, 1))
ax5.set_yticks(np.arange(0, 0.8, 0.2))
ax5.xaxis.set_tick_params(labelsize=8)
ax5.yaxis.set_tick_params(labelsize=8)
plt.grid(True, which='major', linestyle='--')

ax6 = fig.add_subplot(3, 2, 6)
pdf_line1 = ax6.plot(xpdf_pre_djf_rcm_exp1_samz, pdf_pre_djf_rcm_exp1_samz, color='blue', label='RegCM4.7_EXP1', linestyle='-',linewidth=1.5) 
pdf_line2 = ax6.plot(xpdf_pre_djf_rcm_exp2_samz, pdf_pre_djf_rcm_exp2_samz, color='red', label='RegCM4.7_EXP2', linestyle='-', linewidth=1.5) 
pdf_line3 = ax6.plot(xpdf_pre_djf_gcm_samz, pdf_pre_djf_gcm_samz, color='gray', label='HadGEM2-ES', linestyle='-', linewidth=1.5) 
plt.axvline(0, color='black', linestyle='-', linewidth=1.5)
plt.axvline(-2, color='black', linestyle='-', linewidth=1.5)
plt.axvline(2, color='black', linestyle='-', linewidth=1.5)
plt.title(u'F)', loc='left', fontweight='bold', fontsize=8)
ax6.set_xlim(-5, 5)
ax6.set_ylim(0, 0.6)
ax6.set_xticks(np.arange(-5, 6, 1))
ax6.set_yticks(np.arange(0, 0.8, 0.2))
ax6.xaxis.set_tick_params(labelsize=8)
ax6.yaxis.set_tick_params(labelsize=8)
plt.setp(ax6.get_yticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_pdf_reg_exp2.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()






