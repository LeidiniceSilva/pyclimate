# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/05/2021"
__description__ = "This script plot portrait diagram from Reg and Had models and obs database"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from matplotlib.font_manager import FontProperties
from comp_statist_indices import compute_diso


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

	              
# Import regcm exps model end obs database climatology
# Precipitation
ann_cru_pre_samz, djf_cru_pre_samz, jja_cru_pre_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
ann_rcm_exp1_pre_samz, djf_rcm_exp1_pre_samz, jja_rcm_exp1_pre_samz = import_rcm_exp1('pr', 'samz', 'hist', '1986-2005')
ann_rcm_exp2_pre_samz, djf_rcm_exp2_pre_samz, jja_rcm_exp2_pre_samz = import_rcm_exp2('pr', 'samz', 'hist', '1986-2005')
ann_gcm_pre_samz, djf_gcm_pre_samz, jja_gcm_pre_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')

ann_cru_pre_eneb, djf_cru_pre_eneb, jja_cru_pre_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
ann_rcm_exp1_pre_eneb, djf_rcm_exp1_pre_eneb, jja_rcm_exp1_pre_eneb = import_rcm_exp1('pr', 'eneb', 'hist', '1986-2005')
ann_rcm_exp2_pre_eneb, djf_rcm_exp2_pre_eneb, jja_rcm_exp2_pre_eneb = import_rcm_exp2('pr', 'eneb', 'hist', '1986-2005')
ann_gcm_pre_eneb, djf_gcm_pre_eneb, jja_gcm_pre_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')

ann_cru_pre_matopiba, djf_cru_pre_matopiba, jja_cru_pre_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')
ann_rcm_exp1_pre_matopiba, djf_rcm_exp1_pre_matopiba, jja_rcm_exp1_pre_matopiba = import_rcm_exp1('pr', 'matopiba', 'hist', '1986-2005')
ann_rcm_exp2_pre_matopiba, djf_rcm_exp2_pre_matopiba, jja_rcm_exp2_pre_matopiba = import_rcm_exp2('pr', 'matopiba', 'hist', '1986-2005')
ann_gcm_pre_matopiba, djf_gcm_pre_matopiba, jja_gcm_pre_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

# Temperature
ann_cru_tas_samz, djf_cru_tas_samz, jja_cru_tas_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
ann_rcm_exp1_tas_samz, djf_rcm_exp1_tas_samz, jja_rcm_exp1_tas_samz = import_rcm_exp1('tas', 'samz', 'hist', '1986-2005')
ann_rcm_exp2_tas_samz, djf_rcm_exp2_tas_samz, jja_rcm_exp2_tas_samz = import_rcm_exp2('tas', 'samz', 'hist', '1986-2005')
ann_gcm_tas_samz, djf_gcm_tas_samz, jja_gcm_tas_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')

ann_cru_tas_eneb, djf_cru_tas_eneb, jja_cru_tas_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
ann_rcm_exp1_tas_eneb, djf_rcm_exp1_tas_eneb, jja_rcm_exp1_tas_eneb = import_rcm_exp1('tas', 'eneb', 'hist', '1986-2005')
ann_rcm_exp2_tas_eneb, djf_rcm_exp2_tas_eneb, jja_rcm_exp2_tas_eneb = import_rcm_exp2('tas', 'eneb', 'hist', '1986-2005')
ann_gcm_tas_eneb, djf_gcm_tas_eneb, jja_gcm_tas_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')

ann_cru_tas_matopiba, djf_cru_tas_matopiba, jja_cru_tas_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')
ann_rcm_exp1_tas_matopiba, djf_rcm_exp1_tas_matopiba, jja_rcm_exp1_tas_matopiba = import_rcm_exp1('tas', 'matopiba', 'hist', '1986-2005')
ann_rcm_exp2_tas_matopiba, djf_rcm_exp2_tas_matopiba, jja_rcm_exp2_tas_matopiba = import_rcm_exp2('tas', 'matopiba', 'hist', '1986-2005')
ann_gcm_tas_matopiba, djf_gcm_tas_matopiba, jja_gcm_tas_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

# Compute Nash–Sutcliffe Efficient Coefficient
# Precipitation
# RegCM
diso_djf_rcm_exp1_pre_samz = compute_diso(djf_rcm_exp1_pre_samz, djf_cru_pre_samz)
diso_jja_rcm_exp1_pre_samz = compute_diso(jja_rcm_exp1_pre_samz, jja_cru_pre_samz)
diso_ann_rcm_exp1_pre_samz = compute_diso(ann_rcm_exp1_pre_samz, ann_cru_pre_samz)
diso_djf_rcm_exp1_pre_eneb = compute_diso(djf_rcm_exp1_pre_eneb, djf_cru_pre_eneb)
diso_jja_rcm_exp1_pre_eneb = compute_diso(jja_rcm_exp1_pre_eneb, jja_cru_pre_eneb)
diso_ann_rcm_exp1_pre_eneb = compute_diso(ann_rcm_exp1_pre_eneb, ann_cru_pre_eneb)
diso_djf_rcm_exp1_pre_matopiba = compute_diso(djf_rcm_exp1_pre_matopiba, djf_cru_pre_matopiba)
diso_jja_rcm_exp1_pre_matopiba = compute_diso(jja_rcm_exp1_pre_matopiba, jja_cru_pre_matopiba)
diso_ann_rcm_exp1_pre_matopiba = compute_diso(ann_rcm_exp1_pre_matopiba, ann_cru_pre_matopiba)

diso_djf_rcm_exp2_pre_samz = compute_diso(djf_rcm_exp2_pre_samz, djf_cru_pre_samz)
diso_jja_rcm_exp2_pre_samz = compute_diso(jja_rcm_exp2_pre_samz, jja_cru_pre_samz)
diso_ann_rcm_exp2_pre_samz = compute_diso(ann_rcm_exp2_pre_samz, ann_cru_pre_samz)
diso_djf_rcm_exp2_pre_eneb = compute_diso(djf_rcm_exp2_pre_eneb, djf_cru_pre_eneb)
diso_jja_rcm_exp2_pre_eneb = compute_diso(jja_rcm_exp2_pre_eneb, jja_cru_pre_eneb)
diso_ann_rcm_exp2_pre_eneb = compute_diso(ann_rcm_exp2_pre_eneb, ann_cru_pre_eneb)
diso_djf_rcm_exp2_pre_matopiba = compute_diso(djf_rcm_exp2_pre_matopiba, djf_cru_pre_matopiba)
diso_jja_rcm_exp2_pre_matopiba = compute_diso(jja_rcm_exp2_pre_matopiba, jja_cru_pre_matopiba)
diso_ann_rcm_exp2_pre_matopiba = compute_diso(ann_rcm_exp2_pre_matopiba, ann_cru_pre_matopiba)

# HadGEM2-ES
diso_djf_gcm_pre_samz = compute_diso(djf_gcm_pre_samz, djf_cru_pre_samz)
diso_jja_gcm_pre_samz = compute_diso(jja_gcm_pre_samz, jja_cru_pre_samz)
diso_ann_gcm_pre_samz = compute_diso(ann_gcm_pre_samz, ann_cru_pre_samz)
diso_djf_gcm_pre_eneb = compute_diso(djf_gcm_pre_eneb, djf_cru_pre_eneb)
diso_jja_gcm_pre_eneb = compute_diso(jja_gcm_pre_eneb, jja_cru_pre_eneb)
diso_ann_gcm_pre_eneb = compute_diso(ann_gcm_pre_eneb, ann_cru_pre_eneb)
diso_djf_gcm_pre_matopiba = compute_diso(djf_gcm_pre_matopiba, djf_cru_pre_matopiba)
diso_jja_gcm_pre_matopiba = compute_diso(jja_gcm_pre_matopiba, jja_cru_pre_matopiba)
diso_ann_gcm_pre_matopiba = compute_diso(ann_gcm_pre_matopiba, ann_cru_pre_matopiba)

# Temperature
# RegCM
diso_djf_rcm_exp1_tas_samz = compute_diso(np.nanmean(djf_rcm_exp1_tas_samz, axis=1), djf_cru_tas_samz)
diso_jja_rcm_exp1_tas_samz = compute_diso(np.nanmean(jja_rcm_exp1_tas_samz, axis=1), jja_cru_tas_samz)
diso_ann_rcm_exp1_tas_samz = compute_diso(np.nanmean(ann_rcm_exp1_tas_samz, axis=1), ann_cru_tas_samz)
diso_djf_rcm_exp1_tas_eneb = compute_diso(np.nanmean(djf_rcm_exp1_tas_eneb, axis=1), djf_cru_tas_eneb)
diso_jja_rcm_exp1_tas_eneb = compute_diso(np.nanmean(jja_rcm_exp1_tas_eneb, axis=1), jja_cru_tas_eneb)
diso_ann_rcm_exp1_tas_eneb = compute_diso(np.nanmean(ann_rcm_exp1_tas_eneb, axis=1), ann_cru_tas_eneb)
diso_djf_rcm_exp1_tas_matopiba = compute_diso(np.nanmean(djf_rcm_exp1_tas_matopiba, axis=1), djf_cru_tas_matopiba)
diso_jja_rcm_exp1_tas_matopiba = compute_diso(np.nanmean(jja_rcm_exp1_tas_matopiba, axis=1), jja_cru_tas_matopiba)
diso_ann_rcm_exp1_tas_matopiba = compute_diso(np.nanmean(ann_rcm_exp1_tas_matopiba, axis=1), ann_cru_tas_matopiba)

diso_djf_rcm_exp2_tas_samz = compute_diso(np.nanmean(djf_rcm_exp2_tas_samz, axis=1), djf_cru_tas_samz)
diso_jja_rcm_exp2_tas_samz = compute_diso(np.nanmean(jja_rcm_exp2_tas_samz, axis=1), jja_cru_tas_samz)
diso_ann_rcm_exp2_tas_samz = compute_diso(np.nanmean(ann_rcm_exp2_tas_samz, axis=1), ann_cru_tas_samz)
diso_djf_rcm_exp2_tas_eneb = compute_diso(np.nanmean(djf_rcm_exp2_tas_eneb, axis=1), djf_cru_tas_eneb)
diso_jja_rcm_exp2_tas_eneb = compute_diso(np.nanmean(jja_rcm_exp2_tas_eneb, axis=1), jja_cru_tas_eneb)
diso_ann_rcm_exp2_tas_eneb = compute_diso(np.nanmean(ann_rcm_exp2_tas_eneb, axis=1), ann_cru_tas_eneb)
diso_djf_rcm_exp2_tas_matopiba = compute_diso(np.nanmean(djf_rcm_exp2_tas_matopiba, axis=1), djf_cru_tas_matopiba)
diso_jja_rcm_exp2_tas_matopiba = compute_diso(np.nanmean(jja_rcm_exp2_tas_matopiba, axis=1), jja_cru_tas_matopiba)
diso_ann_rcm_exp2_tas_matopiba = compute_diso(np.nanmean(ann_rcm_exp2_tas_matopiba, axis=1), ann_cru_tas_matopiba)

# HadGEM2-ES
diso_djf_gcm_tas_samz = compute_diso(djf_gcm_tas_samz, djf_cru_tas_samz)
diso_jja_gcm_tas_samz = compute_diso(jja_gcm_tas_samz, jja_cru_tas_samz)
diso_ann_gcm_tas_samz = compute_diso(ann_gcm_tas_samz, ann_cru_tas_samz)
diso_djf_gcm_tas_eneb = compute_diso(djf_gcm_tas_eneb, djf_cru_tas_eneb)
diso_jja_gcm_tas_eneb = compute_diso(jja_gcm_tas_eneb, jja_cru_tas_eneb)
diso_ann_gcm_tas_eneb = compute_diso(ann_gcm_tas_eneb, ann_cru_tas_eneb)
diso_djf_gcm_tas_matopiba = compute_diso(djf_gcm_tas_matopiba, djf_cru_tas_matopiba)
diso_jja_gcm_tas_matopiba = compute_diso(jja_gcm_tas_matopiba, jja_cru_tas_matopiba)
diso_ann_gcm_tas_matopiba = compute_diso(ann_gcm_tas_matopiba, ann_cru_tas_matopiba)

# Precipitation
diso_samz_pre = np.array([[diso_djf_rcm_exp1_pre_samz, diso_jja_rcm_exp1_pre_samz, diso_ann_rcm_exp1_pre_samz],
[diso_djf_rcm_exp2_pre_samz, diso_jja_rcm_exp2_pre_samz, diso_ann_rcm_exp2_pre_samz], 
[diso_djf_gcm_pre_samz, diso_jja_gcm_pre_samz, diso_ann_gcm_pre_samz]])

diso_eneb_pre = np.array([[diso_djf_rcm_exp1_pre_eneb, diso_jja_rcm_exp1_pre_eneb, diso_ann_rcm_exp1_pre_eneb],
[diso_djf_rcm_exp2_pre_eneb, diso_jja_rcm_exp2_pre_eneb, diso_ann_rcm_exp2_pre_eneb], 
[diso_djf_gcm_pre_eneb, diso_jja_gcm_pre_eneb, diso_ann_gcm_pre_eneb]])

diso_matopiba_pre = np.array([[diso_djf_rcm_exp1_pre_matopiba, diso_jja_rcm_exp1_pre_matopiba, diso_ann_rcm_exp1_pre_matopiba],
[diso_djf_rcm_exp2_pre_matopiba, diso_jja_rcm_exp2_pre_matopiba, diso_ann_rcm_exp2_pre_matopiba], 
[diso_djf_gcm_pre_matopiba, diso_jja_gcm_pre_matopiba, diso_ann_gcm_pre_matopiba]])

# Temperature
diso_samz_tas = np.array([[diso_djf_rcm_exp1_tas_samz, diso_jja_rcm_exp1_tas_samz, diso_ann_rcm_exp1_tas_samz],
[diso_djf_rcm_exp2_tas_samz, diso_jja_rcm_exp2_tas_samz, diso_ann_rcm_exp2_tas_samz], 
[diso_djf_gcm_tas_samz, diso_jja_gcm_tas_samz, diso_ann_gcm_tas_samz]])

diso_eneb_tas = np.array([[diso_djf_rcm_exp1_tas_eneb, diso_jja_rcm_exp1_tas_eneb, diso_ann_rcm_exp1_tas_eneb],
[diso_djf_rcm_exp2_tas_eneb, diso_jja_rcm_exp2_tas_eneb, diso_ann_rcm_exp2_tas_eneb], 
[diso_djf_gcm_tas_eneb, diso_jja_gcm_tas_eneb, diso_ann_gcm_tas_eneb]])

diso_matopiba_tas = np.array([[diso_djf_rcm_exp1_tas_matopiba, diso_jja_rcm_exp1_tas_matopiba, diso_ann_rcm_exp1_tas_matopiba],
[diso_djf_rcm_exp2_tas_matopiba, diso_jja_rcm_exp2_tas_matopiba, diso_ann_rcm_exp2_tas_matopiba], 
[diso_djf_gcm_tas_matopiba, diso_jja_gcm_tas_matopiba, diso_ann_gcm_tas_matopiba]])

# Plot model end obs data climatology
fig = plt.figure()

norm = colors.BoundaryNorm(boundaries=np.arange(0, 2, 0.2), ncolors=256)
color_map = plt.cm.get_cmap('RdYlGn')
reversed_color_map = color_map.reversed()

xlabels = [u'DJF', u'JJA', u'ANN']
ylabels = [u'HadGEM2-ES', u'RegCM4.7_EXP2', u'RegCM4.7_EXP1']

ax1 = fig.add_subplot(3, 2, 1)
pcm1 = ax1.pcolormesh(diso_samz_pre[::-1], edgecolors ='k', linewidths = 1.5, norm=norm, cmap=reversed_color_map)
ax1.set_title(u'A)', loc='left', fontweight='bold', fontsize=8)
ax1.set_xticks(np.arange(diso_samz_pre.shape[1]) + 0.5)
ax1.set_yticks(np.arange(diso_samz_pre.shape[0]) + 0.5)
ax1.set_xticklabels(xlabels, fontsize=8)
ax1.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax1.get_xticklabels(), visible=False)
cbar = fig.colorbar(pcm1, extend='max')
cbar.ax.tick_params(labelsize=8)   
for y in range(diso_samz_pre.shape[0]):
    for x in range(diso_samz_pre.shape[1]):
        plt.text(x + 0.5, y + 0.5, '%.2f' % diso_samz_pre[::-1][y, x],
                 ha="center", va="center", color='k',
                 )
                 
ax2 = fig.add_subplot(3, 2, 2)                 
pcm2 = ax2.pcolormesh(diso_samz_tas[::-1], edgecolors ='k', linewidths = 1.5, norm=norm, cmap=reversed_color_map)
ax2.set_title('D)', loc='left', fontweight='bold', fontsize=8)
ax2.set_xticks(np.arange(diso_samz_tas.shape[1]) + 0.5)
ax2.set_yticks(np.arange(diso_samz_tas.shape[0]) + 0.5)
ax2.set_xticklabels(xlabels, fontsize=8)
ax2.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
cbar = fig.colorbar(pcm2, extend='max')
cbar.ax.tick_params(labelsize=8)
for y in range(diso_samz_tas.shape[0]):
    for x in range(diso_samz_tas.shape[1]):
        plt.text(x + 0.5, y + 0.5, '%.2f' % diso_samz_tas[::-1][y, x],
                 ha="center", va="center", color='k',
                 )   

ax3 = fig.add_subplot(3, 2, 3)                 
pcm3 = ax3.pcolormesh(diso_eneb_pre[::-1], edgecolors ='k', linewidths = 1.5, norm=norm, cmap=reversed_color_map)
ax3.set_title(u'B)', loc='left', fontweight='bold', fontsize=8)
ax3.set_xticks(np.arange(diso_eneb_pre.shape[1]) + 0.5)
ax3.set_yticks(np.arange(diso_eneb_pre.shape[0]) + 0.5)
ax3.set_xticklabels(xlabels, fontsize=8)
ax3.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax3.get_xticklabels(), visible=False)
cbar = fig.colorbar(pcm3, extend='max')
cbar.set_label('Precipitation (mm d⁻¹)', rotation=90, fontsize=8)
cbar.ax.yaxis.set_label_position('right')
cbar.ax.tick_params(labelsize=8)  
for y in range(diso_eneb_pre.shape[0]):
    for x in range(diso_eneb_pre.shape[1]):
        plt.text(x + 0.5, y + 0.5, '%.2f' % diso_eneb_pre[::-1][y, x],
                 ha="center", va="center", color='k',
                 ) 
                 
ax4 = fig.add_subplot(3, 2, 4)
pcm4 = ax4.pcolormesh(diso_eneb_tas[::-1], edgecolors ='k', linewidths = 1.5, norm=norm, cmap=reversed_color_map)
ax4.set_title('E)', loc='left', fontweight='bold', fontsize=8)
ax4.set_xticks(np.arange(diso_eneb_tas.shape[1]) + 0.5)
ax4.set_yticks(np.arange(diso_eneb_tas.shape[0]) + 0.5)
ax4.set_xticklabels(xlabels, fontsize=8)
ax4.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax4.get_xticklabels(), visible=False)
plt.setp(ax4.get_yticklabels(), visible=False)
cbar = fig.colorbar(pcm4, extend='max')
cbar.set_label('Temperature (°C)', rotation=90, fontsize=8)
cbar.ax.yaxis.set_label_position('right')
cbar.ax.tick_params(labelsize=8)  
for y in range(diso_eneb_tas.shape[0]):
    for x in range(diso_eneb_tas.shape[1]):
        plt.text(x + 0.5, y + 0.5, '%.2f' % diso_eneb_tas[::-1][y, x],
                 ha="center", va="center", color='k',
                 ) 
                 
ax5 = fig.add_subplot(3, 2, 5)                 
pcm5 = ax5.pcolormesh(diso_matopiba_pre[::-1], edgecolors ='k', linewidths = 1.5, norm=norm, cmap=reversed_color_map)
ax5.set_title(u'C)', loc='left', fontweight='bold', fontsize=8)
ax5.set_xticks(np.arange(diso_matopiba_pre.shape[1]) + 0.5)
ax5.set_yticks(np.arange(diso_matopiba_pre.shape[0]) + 0.5)
ax5.set_xticklabels(xlabels, fontsize=8)
ax5.set_yticklabels(ylabels, fontsize=8)
cbar = fig.colorbar(pcm5, extend='max')
cbar.ax.tick_params(labelsize=8) 
for y in range(diso_matopiba_pre.shape[0]):
    for x in range(diso_matopiba_pre.shape[1]):
        plt.text(x + 0.5, y + 0.5, '%.2f' % diso_matopiba_pre[::-1][y, x],
                 ha="center", va="center", color='k',
                 ) 
                 
ax6 = fig.add_subplot(3, 2, 6)
pcm6 = ax6.pcolormesh(diso_matopiba_tas[::-1], edgecolors ='k', linewidths = 1.5, norm=norm, cmap=reversed_color_map)
ax6.set_title('F)', loc='left', fontweight='bold', fontsize=8)
ax6.set_xticks(np.arange(diso_matopiba_tas.shape[1]) + 0.5)
ax6.set_yticks(np.arange(diso_matopiba_tas.shape[0]) + 0.5)
ax6.set_xticklabels(xlabels, fontsize=8)
ax6.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax6.get_yticklabels(), visible=False)
cbar = fig.colorbar(pcm6, extend='max')
cbar.ax.tick_params(labelsize=8)  
for y in range(diso_matopiba_tas.shape[0]):
    for x in range(diso_matopiba_tas.shape[1]):
        plt.text(x + 0.5, y + 0.5, '%.2f' % diso_matopiba_tas[::-1][y, x],
                 ha="center", va="center", color='k',
                 ) 
     
# Save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_portrait_diagram_reg_exp2.png'

if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()


