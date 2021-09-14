# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
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
	
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
		
	return value
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp1/hist'
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	return value


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/rcm_exp1/hist'
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	return value

	              
# Import models and obs database 
# Precipitation
mon_pre_cru_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
mon_pre_cru_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
mon_pre_cru_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')

mon_pre_reg_samz = import_rcm('pr', 'samz', 'hist', '1986-2005')
mon_pre_reg_eneb = import_rcm('pr', 'eneb', 'hist', '1986-2005')
mon_pre_reg_matopiba = import_rcm('pr', 'matopiba', 'hist', '1986-2005')

mon_pre_had_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')
mon_pre_had_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')
mon_pre_had_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

# Temperature
mon_tas_cru_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
mon_tas_cru_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
mon_tas_cru_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')

mon_tas_reg_samz = import_rcm('tas', 'samz', 'hist', '1986-2005')
mon_tas_reg_eneb = import_rcm('tas', 'eneb', 'hist', '1986-2005')
mon_tas_reg_matopiba = import_rcm('tas', 'matopiba', 'hist', '1986-2005')

mon_tas_had_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')
mon_tas_had_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')
mon_tas_had_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

# Calculate PDF function
# Precipitation
xpdf_pre_cru_samz, pdf_pre_cru_samz = compute_pdf(mon_pre_cru_samz)
xpdf_pre_cru_eneb, pdf_pre_cru_eneb = compute_pdf(mon_pre_cru_eneb)
xpdf_pre_cru_matopiba, pdf_pre_cru_matopiba = compute_pdf(mon_pre_cru_matopiba)

xpdf_pre_reg_samz, pdf_pre_reg_samz = compute_pdf(mon_pre_reg_samz)
xpdf_pre_reg_eneb, pdf_pre_reg_eneb = compute_pdf(mon_pre_reg_eneb)
xpdf_pre_reg_matopiba, pdf_pre_reg_matopiba = compute_pdf(mon_pre_reg_matopiba)

xpdf_pre_had_samz, pdf_pre_had_samz = compute_pdf(mon_pre_had_samz)
xpdf_pre_had_eneb, pdf_pre_had_eneb = compute_pdf(mon_pre_had_eneb)
xpdf_pre_had_matopiba, pdf_pre_had_matopiba = compute_pdf(mon_pre_had_matopiba)

# Temperature
xpdf_tas_cru_samz, pdf_tas_cru_samz = compute_pdf(mon_tas_cru_samz)
xpdf_tas_cru_eneb, pdf_tas_cru_eneb = compute_pdf(mon_tas_cru_eneb)
xpdf_tas_cru_matopiba, pdf_tas_cru_matopiba = compute_pdf(mon_tas_cru_matopiba)

xpdf_tas_reg_samz, pdf_tas_reg_samz = compute_pdf(np.nanmean(mon_tas_reg_samz, axis=0))
xpdf_tas_reg_eneb, pdf_tas_reg_eneb = compute_pdf(np.nanmean(mon_tas_reg_eneb, axis=0))
xpdf_tas_reg_matopiba, pdf_tas_reg_matopiba = compute_pdf(np.nanmean(mon_tas_reg_matopiba, axis=0))

xpdf_tas_had_samz, pdf_tas_had_samz = compute_pdf(mon_tas_had_samz)
xpdf_tas_had_eneb, pdf_tas_had_eneb = compute_pdf(mon_tas_had_eneb)
xpdf_tas_had_matopiba, pdf_tas_had_matopiba = compute_pdf(mon_tas_had_matopiba)

# Plot models and obs database 
fig = plt.figure()

ax1 = fig.add_subplot(3, 2, 1)
pdf_line1 = ax1.plot(xpdf_pre_cru_samz, pdf_pre_cru_samz, color='black', label='CRU', linestyle='-', linewidth=1.5) 
pdf_line2 = ax1.plot(xpdf_pre_reg_samz, pdf_pre_reg_samz, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5)
pdf_line3 = ax1.plot(xpdf_pre_had_samz, pdf_pre_had_samz, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_cru_samz), color='black', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_reg_samz), color='gray', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_had_samz), color='dimgray', linestyle='--', linewidth=1.5)
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
ax1.set_xlim(0, 14)
ax1.set_ylim(0, 0.4)
ax1.set_xticks(np.arange(0, 16, 2))
ax1.set_yticks(np.arange(0, 0.6, 0.2))
ax1.xaxis.set_tick_params(labelsize=8)
ax1.yaxis.set_tick_params(labelsize=8)
plt.setp(ax1.get_xticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')
plt.legend(fontsize=6, loc=1, shadow=True, ncol=1)

ax2 = fig.add_subplot(3, 2, 2)
pdf_line1 = ax2.plot(xpdf_tas_cru_samz, pdf_tas_cru_samz, color='black', label='CRU', linestyle='-', linewidth=1.5)
pdf_line2 = ax2.plot(xpdf_tas_reg_samz, pdf_tas_reg_samz, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5)
pdf_line3 = ax2.plot(xpdf_tas_had_samz, pdf_tas_had_samz, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_cru_samz), color='black', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_reg_samz), color='gray', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_had_samz), color='dimgray', linestyle='--', linewidth=1.5)
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
ax2.set_xlim(20, 30)
ax2.set_ylim(0, 0.8)
ax2.set_xticks(np.arange(20, 32, 2))
ax2.set_yticks(np.arange(0, 1, 0.2))
ax2.xaxis.set_tick_params(labelsize=8)
ax2.yaxis.set_tick_params(labelsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)	     
plt.grid(True, which='major', linestyle='--')
	 
ax3 = fig.add_subplot(3, 2, 3)
pdf_line1 = ax3.plot(xpdf_pre_cru_eneb, pdf_pre_cru_eneb, color='black', label='CRU', linestyle='-', linewidth=1.5)
pdf_line2 = ax3.plot(xpdf_pre_reg_eneb, pdf_pre_reg_eneb, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5)
pdf_line3 = ax3.plot(xpdf_pre_had_eneb, pdf_pre_had_eneb, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_cru_eneb), color='black', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_reg_eneb), color='gray', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_had_eneb), color='dimgray', linestyle='--', linewidth=1.5)
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
ax3.set_ylabel(u'PDF', fontsize=8)
ax3.set_xlim(0, 14)
ax3.set_ylim(0, 0.4)
ax3.set_xticks(np.arange(0, 16, 2))
ax3.set_yticks(np.arange(0, 0.6, 0.2))
ax3.xaxis.set_tick_params(labelsize=8)
ax3.yaxis.set_tick_params(labelsize=8)
plt.setp(ax3.get_xticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')
	 
ax4 = fig.add_subplot(3, 2, 4)
pdf_line1 = ax4.plot(xpdf_tas_cru_eneb, pdf_tas_cru_eneb, color='black', label='CRU', linestyle='-', linewidth=1.5)
pdf_line2 = ax4.plot(xpdf_tas_reg_eneb, pdf_tas_reg_eneb, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5)
pdf_line3 = ax4.plot(xpdf_tas_had_eneb, pdf_tas_had_eneb, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_cru_eneb), color='black', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_reg_eneb), color='gray', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_had_eneb), color='dimgray', linestyle='--', linewidth=1.5)
plt.title(u'E)', loc='left', fontweight='bold', fontsize=8)
ax4.set_ylabel(u'PDF', fontsize=8)
ax4.set_xlim(20, 30)
ax4.set_ylim(0, 0.8)
ax4.set_xticks(np.arange(20, 32, 2))
ax4.set_yticks(np.arange(0, 1, 0.2))
ax4.xaxis.set_tick_params(labelsize=8)
ax4.yaxis.set_tick_params(labelsize=8)
plt.setp(ax4.get_xticklabels(), visible=False)	     
plt.grid(True, which='major', linestyle='--')
		 
ax5 = fig.add_subplot(3, 2, 5)
pdf_line1 = ax5.plot(xpdf_pre_cru_matopiba, pdf_pre_cru_matopiba, color='black', label='CRU', linestyle='-', linewidth=1.5)
pdf_line2 = ax5.plot(xpdf_pre_reg_matopiba, pdf_pre_reg_matopiba, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5)
pdf_line3 = ax5.plot(xpdf_pre_had_matopiba, pdf_pre_had_matopiba, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_cru_matopiba), color='black', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_reg_matopiba), color='gray', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_had_matopiba), color='dimgray', linestyle='--', linewidth=1.5)
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'Precipitation (mm d⁻¹)', fontsize=8)
ax5.set_xlim(0, 14)
ax5.set_ylim(0, 0.4)
ax5.set_xticks(np.arange(0, 16, 2))
ax5.set_yticks(np.arange(0, 0.6, 0.2))
ax5.xaxis.set_tick_params(labelsize=8)
ax5.yaxis.set_tick_params(labelsize=8)
plt.grid(True, which='major', linestyle='--')

ax6 = fig.add_subplot(3, 2, 6)
pdf_line1 = ax6.plot(xpdf_tas_cru_matopiba, pdf_tas_cru_matopiba, color='black', label='CRU', linestyle='-', linewidth=1.5)
pdf_line2 = ax6.plot(xpdf_tas_reg_matopiba, pdf_tas_reg_matopiba, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5)
pdf_line3 = ax6.plot(xpdf_tas_had_matopiba, pdf_tas_had_matopiba, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_cru_matopiba), color='black', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_reg_matopiba), color='gray', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_had_matopiba), color='dimgray', linestyle='--', linewidth=1.5)
plt.title(u'F)', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'Temperature (°C)', fontsize=8)
ax6.set_xlim(20, 30)
ax6.set_ylim(0, 0.8)
ax6.set_xticks(np.arange(20, 32, 2))
ax6.set_yticks(np.arange(0, 1, 0.2))
ax6.xaxis.set_tick_params(labelsize=8)
ax6.yaxis.set_tick_params(labelsize=8)
plt.grid(True, which='major', linestyle='--')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_pdf_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()


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
	mam_obs = season_obs[0:80:4]
	jja_obs = season_obs[1:80:4]
	son_obs = season_obs[2:80:4]

	return ann_obs, djf_obs, mam_obs, jja_obs, son_obs 
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp1/hist'
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	ann_rcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1),axis=1) 
	season_rcm = np.nanmean(np.nanmean(value[2:239:3,:,:], axis=1), axis=1)
	djf_rcm = season_rcm[3:80:4]
	mam_rcm = season_rcm[0:80:4]
	jja_rcm = season_rcm[1:80:4]
	son_rcm = season_rcm[2:80:4]

	return ann_rcm, djf_rcm, mam_rcm, jja_rcm, son_rcm


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
	mam_gcm = season_gcm[0:80:4]
	jja_gcm = season_gcm[1:80:4]
	son_gcm = season_gcm[2:80:4]

	return ann_gcm, djf_gcm, mam_gcm, jja_gcm, son_gcm

	              
# Import models and obs database 
# Precipitation
ann_obs_pre_samz, djf_obs_pre_samz, mam_obs_pre_samz, jja_obs_pre_samz, son_obs_pre_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
ann_rcm_pre_samz, djf_rcm_pre_samz, mam_rcm_pre_samz, jja_rcm_pre_samz, son_rcm_pre_samz = import_rcm('pr', 'samz', 'hist', '1986-2005')
ann_gcm_pre_samz, djf_gcm_pre_samz, mam_gcm_pre_samz, jja_gcm_pre_samz, son_gcm_pre_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')

ann_obs_pre_eneb, djf_obs_pre_eneb, mam_obs_pre_eneb, jja_obs_pre_eneb, son_obs_pre_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
ann_rcm_pre_eneb, djf_rcm_pre_eneb, mam_rcm_pre_eneb, jja_rcm_pre_eneb, son_rcm_pre_eneb = import_rcm('pr', 'eneb', 'hist', '1986-2005')
ann_gcm_pre_eneb, djf_gcm_pre_eneb, mam_gcm_pre_eneb, jja_gcm_pre_eneb, son_gcm_pre_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')

ann_obs_pre_matopiba, djf_obs_pre_matopiba, mam_obs_pre_matopiba, jja_obs_pre_matopiba, son_obs_pre_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')
ann_rcm_pre_matopiba, djf_rcm_pre_matopiba, mam_rcm_pre_matopiba, jja_rcm_pre_matopiba, son_rcm_pre_matopiba = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
ann_gcm_pre_matopiba, djf_gcm_pre_matopiba, mam_gcm_pre_matopiba, jja_gcm_pre_matopiba, son_gcm_pre_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

# Temperature
ann_obs_tas_samz, djf_obs_tas_samz, mam_obs_tas_samz, jja_obs_tas_samz, son_obs_tas_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
ann_rcm_tas_samz, djf_rcm_tas_samz, mam_rcm_tas_samz, jja_rcm_tas_samz, son_rcm_tas_samz = import_rcm('tas', 'samz', 'hist', '1986-2005')
ann_gcm_tas_samz, djf_gcm_tas_samz, mam_gcm_tas_samz, jja_gcm_tas_samz, son_gcm_tas_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')

ann_obs_tas_eneb, djf_obs_tas_eneb, mam_obs_tas_eneb, jja_obs_tas_eneb, son_obs_tas_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
ann_rcm_tas_eneb, djf_rcm_tas_eneb, mam_rcm_tas_eneb, jja_rcm_tas_eneb, son_rcm_tas_eneb = import_rcm('tas', 'eneb', 'hist', '1986-2005')
ann_gcm_tas_eneb, djf_gcm_tas_eneb, mam_gcm_tas_eneb, jja_gcm_tas_eneb, son_gcm_tas_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')

ann_obs_tas_matopiba, djf_obs_tas_matopiba, mam_obs_tas_matopiba, jja_obs_tas_matopiba, son_obs_tas_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')
ann_rcm_tas_matopiba, djf_rcm_tas_matopiba, mam_rcm_tas_matopiba, jja_rcm_tas_matopiba, son_rcm_tas_matopiba = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
ann_gcm_tas_matopiba, djf_gcm_tas_matopiba, mam_gcm_tas_matopiba, jja_gcm_tas_matopiba, son_gcm_tas_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

# Compute bias
# Precipitation
bias_djf_rcm_pre_samz = djf_rcm_pre_samz - djf_obs_pre_samz
bias_mam_rcm_pre_samz = mam_rcm_pre_samz - mam_obs_pre_samz
bias_jja_rcm_pre_samz = jja_rcm_pre_samz - jja_obs_pre_samz
bias_son_rcm_pre_samz = son_rcm_pre_samz - son_obs_pre_samz
bias_ann_rcm_pre_samz = ann_rcm_pre_samz - ann_obs_pre_samz

bias_djf_gcm_pre_samz = djf_gcm_pre_samz - djf_obs_pre_samz
bias_mam_gcm_pre_samz = mam_gcm_pre_samz - mam_obs_pre_samz
bias_jja_gcm_pre_samz = jja_gcm_pre_samz - jja_obs_pre_samz
bias_son_gcm_pre_samz = son_gcm_pre_samz - son_obs_pre_samz
bias_ann_gcm_pre_samz = ann_gcm_pre_samz - ann_obs_pre_samz

bias_djf_rcm_pre_eneb = djf_rcm_pre_eneb - djf_obs_pre_eneb
bias_mam_rcm_pre_eneb = mam_rcm_pre_eneb - mam_obs_pre_eneb
bias_jja_rcm_pre_eneb = jja_rcm_pre_eneb - jja_obs_pre_eneb
bias_son_rcm_pre_eneb = son_rcm_pre_eneb - son_obs_pre_eneb
bias_ann_rcm_pre_eneb = ann_rcm_pre_eneb - ann_obs_pre_eneb

bias_djf_gcm_pre_eneb = djf_gcm_pre_eneb - djf_obs_pre_eneb
bias_mam_gcm_pre_eneb = mam_gcm_pre_eneb - mam_obs_pre_eneb
bias_jja_gcm_pre_eneb = jja_gcm_pre_eneb - jja_obs_pre_eneb
bias_son_gcm_pre_eneb = son_gcm_pre_eneb - son_obs_pre_eneb
bias_ann_gcm_pre_eneb = ann_gcm_pre_eneb - ann_obs_pre_eneb

bias_djf_rcm_pre_matopiba = djf_rcm_pre_matopiba - djf_obs_pre_matopiba
bias_mam_rcm_pre_matopiba = mam_rcm_pre_matopiba - mam_obs_pre_matopiba
bias_jja_rcm_pre_matopiba = jja_rcm_pre_matopiba - jja_obs_pre_matopiba
bias_son_rcm_pre_matopiba = son_rcm_pre_matopiba - son_obs_pre_matopiba
bias_ann_rcm_pre_matopiba = ann_rcm_pre_matopiba - ann_obs_pre_matopiba

bias_djf_gcm_pre_matopiba = djf_gcm_pre_matopiba - djf_obs_pre_matopiba
bias_mam_gcm_pre_matopiba = mam_gcm_pre_matopiba - mam_obs_pre_matopiba
bias_jja_gcm_pre_matopiba = jja_gcm_pre_matopiba - jja_obs_pre_matopiba
bias_son_gcm_pre_matopiba = son_gcm_pre_matopiba - son_obs_pre_matopiba
bias_ann_gcm_pre_matopiba = ann_gcm_pre_matopiba - ann_obs_pre_matopiba

# Calculate PDF function
xpdf_ann_rcm_pre_samz, pdf_ann_rcm_pre_samz = compute_pdf(bias_ann_rcm_pre_samz)
xpdf_djf_rcm_pre_samz, pdf_djf_rcm_pre_samz = compute_pdf(bias_djf_rcm_pre_samz)
xpdf_mam_rcm_pre_samz, pdf_mam_rcm_pre_samz = compute_pdf(bias_mam_rcm_pre_samz)
xpdf_jja_rcm_pre_samz, pdf_jja_rcm_pre_samz = compute_pdf(bias_jja_rcm_pre_samz)
xpdf_son_rcm_pre_samz, pdf_son_rcm_pre_samz = compute_pdf(bias_son_rcm_pre_samz)

xpdf_ann_gcm_pre_samz, pdf_ann_gcm_pre_samz = compute_pdf(bias_ann_gcm_pre_samz)
xpdf_djf_gcm_pre_samz, pdf_djf_gcm_pre_samz = compute_pdf(bias_djf_gcm_pre_samz)
xpdf_mam_gcm_pre_samz, pdf_mam_gcm_pre_samz = compute_pdf(bias_mam_gcm_pre_samz)
xpdf_jja_gcm_pre_samz, pdf_jja_gcm_pre_samz = compute_pdf(bias_jja_gcm_pre_samz)
xpdf_son_gcm_pre_samz, pdf_son_gcm_pre_samz = compute_pdf(bias_son_gcm_pre_samz)

xpdf_ann_rcm_pre_eneb, pdf_ann_rcm_pre_eneb = compute_pdf(bias_ann_rcm_pre_eneb)
xpdf_djf_rcm_pre_eneb, pdf_djf_rcm_pre_eneb = compute_pdf(bias_djf_rcm_pre_eneb)
xpdf_mam_rcm_pre_eneb, pdf_mam_rcm_pre_eneb = compute_pdf(bias_mam_rcm_pre_eneb)
xpdf_jja_rcm_pre_eneb, pdf_jja_rcm_pre_eneb = compute_pdf(bias_jja_rcm_pre_eneb)
xpdf_son_rcm_pre_eneb, pdf_son_rcm_pre_eneb = compute_pdf(bias_son_rcm_pre_eneb)

xpdf_ann_gcm_pre_eneb, pdf_ann_gcm_pre_eneb = compute_pdf(bias_ann_gcm_pre_eneb)
xpdf_djf_gcm_pre_eneb, pdf_djf_gcm_pre_eneb = compute_pdf(bias_djf_gcm_pre_eneb)
xpdf_mam_gcm_pre_eneb, pdf_mam_gcm_pre_eneb = compute_pdf(bias_mam_gcm_pre_eneb)
xpdf_jja_gcm_pre_eneb, pdf_jja_gcm_pre_eneb = compute_pdf(bias_jja_gcm_pre_eneb)
xpdf_son_gcm_pre_eneb, pdf_son_gcm_pre_eneb = compute_pdf(bias_son_gcm_pre_eneb)

xpdf_ann_rcm_pre_matopiba, pdf_ann_rcm_pre_matopiba = compute_pdf(bias_ann_rcm_pre_matopiba)
xpdf_djf_rcm_pre_matopiba, pdf_djf_rcm_pre_matopiba = compute_pdf(bias_djf_rcm_pre_matopiba)
xpdf_mam_rcm_pre_matopiba, pdf_mam_rcm_pre_matopiba = compute_pdf(bias_mam_rcm_pre_matopiba)
xpdf_jja_rcm_pre_matopiba, pdf_jja_rcm_pre_matopiba = compute_pdf(bias_jja_rcm_pre_matopiba)
xpdf_son_rcm_pre_matopiba, pdf_son_rcm_pre_matopiba = compute_pdf(bias_son_rcm_pre_matopiba)

xpdf_ann_gcm_pre_matopiba, pdf_ann_gcm_pre_matopiba = compute_pdf(bias_ann_gcm_pre_matopiba)
xpdf_djf_gcm_pre_matopiba, pdf_djf_gcm_pre_matopiba = compute_pdf(bias_djf_gcm_pre_matopiba)
xpdf_mam_gcm_pre_matopiba, pdf_mam_gcm_pre_matopiba = compute_pdf(bias_mam_gcm_pre_matopiba)
xpdf_jja_gcm_pre_matopiba, pdf_jja_gcm_pre_matopiba = compute_pdf(bias_jja_gcm_pre_matopiba)
xpdf_son_gcm_pre_matopiba, pdf_son_gcm_pre_matopiba = compute_pdf(bias_son_gcm_pre_matopiba)

# Plot models and obs database 
fig = plt.figure()

ax = fig.add_subplot(3, 5, 1)
ax.plot(xpdf_djf_rcm_pre_samz, pdf_djf_rcm_pre_samz, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_djf_gcm_pre_samz, pdf_djf_gcm_pre_samz, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_xticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 2)
ax.plot(xpdf_mam_rcm_pre_samz, pdf_mam_rcm_pre_samz, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_mam_gcm_pre_samz, pdf_mam_gcm_pre_samz, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_xticklabels(), visible=False)		     
plt.setp(ax.get_yticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 3)
ax.plot(xpdf_jja_rcm_pre_samz, pdf_jja_rcm_pre_samz, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_jja_gcm_pre_samz, pdf_jja_gcm_pre_samz, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_xticklabels(), visible=False)		     
plt.setp(ax.get_yticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 4)
ax.plot(xpdf_son_rcm_pre_samz, pdf_son_rcm_pre_samz, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_son_gcm_pre_samz, pdf_son_gcm_pre_samz, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_xticklabels(), visible=False)		     
plt.setp(ax.get_yticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 5)
ax.plot(xpdf_ann_rcm_pre_samz, pdf_ann_rcm_pre_samz, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_ann_gcm_pre_samz, pdf_ann_gcm_pre_samz, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'E)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_xticklabels(), visible=False)		     
plt.setp(ax.get_yticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 6)
ax.plot(xpdf_djf_rcm_pre_samz, pdf_djf_rcm_pre_eneb, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_djf_gcm_pre_samz, pdf_djf_gcm_pre_eneb, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'F)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_xticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 7)
ax.plot(xpdf_mam_rcm_pre_samz, pdf_mam_rcm_pre_eneb, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_mam_gcm_pre_samz, pdf_mam_gcm_pre_eneb, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'G)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_xticklabels(), visible=False)		     
plt.setp(ax.get_yticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 8)
ax.plot(xpdf_jja_rcm_pre_samz, pdf_jja_rcm_pre_eneb, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_jja_gcm_pre_samz, pdf_jja_gcm_pre_eneb, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'H)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_xticklabels(), visible=False)		     
plt.setp(ax.get_yticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 9)
ax.plot(xpdf_son_rcm_pre_samz, pdf_son_rcm_pre_eneb, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_son_gcm_pre_samz, pdf_son_gcm_pre_eneb, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'I)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_xticklabels(), visible=False)		     
plt.setp(ax.get_yticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 10)
ax.plot(xpdf_ann_rcm_pre_samz, pdf_ann_rcm_pre_eneb, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_ann_gcm_pre_samz, pdf_ann_gcm_pre_eneb, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'J)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_xticklabels(), visible=False)		     
plt.setp(ax.get_yticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 11)
ax.plot(xpdf_djf_rcm_pre_samz, pdf_djf_rcm_pre_matopiba, color='gray', label='CRU', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_djf_gcm_pre_samz, pdf_djf_gcm_pre_matopiba, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'K)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 12)
ax.plot(xpdf_mam_rcm_pre_samz, pdf_mam_rcm_pre_matopiba, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_mam_gcm_pre_samz, pdf_mam_gcm_pre_matopiba, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'L)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_yticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 13)
ax.plot(xpdf_jja_rcm_pre_samz, pdf_jja_rcm_pre_matopiba, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_jja_gcm_pre_samz, pdf_jja_gcm_pre_matopiba, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'M)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_yticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 14)
ax.plot(xpdf_son_rcm_pre_samz, pdf_son_rcm_pre_matopiba, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_son_gcm_pre_samz, pdf_son_gcm_pre_matopiba, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'N)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_yticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 5, 15)
ax.plot(xpdf_ann_rcm_pre_samz, pdf_ann_rcm_pre_matopiba, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5) 
ax.plot(xpdf_ann_gcm_pre_samz, pdf_ann_gcm_pre_matopiba, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5) 
plt.title(u'O)', loc='left', fontweight='bold', fontsize=8)
plt.xlim(-8, 8)
plt.ylim(0, 0.6)
plt.xticks(np.arange(-8, 10, 4))
plt.yticks(np.arange(0, 0.8, 0.2))
plt.axvline(0, color='black', linestyle='--', linewidth=1.)
plt.setp(ax.get_yticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_pdf_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()


