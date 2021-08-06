# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot boxplot from regcm47 and hadgem models to rcp"

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
from sklearn import metrics
from scipy.stats import norm
from matplotlib.font_manager import FontProperties


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/rcm'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return obs
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return gcm

	              
# Import models 
pre_cru_samz_hist = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
pre_reg_samz_hist = import_rcm('pr', 'samz', 'hist', '1986-2005')
pre_had_samz_hist = import_gcm('pr', 'samz', 'hist', '1986-2005')
pre_reg_samz_rcp26 = import_rcm('pr', 'samz', 'rcp26', '2080-2099')
pre_reg_samz_rcp85 = import_rcm('pr', 'samz', 'rcp85', '2080-2099')
pre_had_samz_rcp26 = import_gcm('pr', 'samz', 'rcp26', '2080-2099')
pre_had_samz_rcp85 = import_gcm('pr', 'samz', 'rcp85', '2080-2099')
pre_samz = [pre_cru_samz_hist, pre_reg_samz_hist, pre_had_samz_hist, pre_reg_samz_rcp26, pre_had_samz_rcp26, pre_reg_samz_rcp85, pre_had_samz_rcp85]

pre_cru_eneb_hist = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
pre_reg_eneb_hist = import_rcm('pr', 'eneb', 'hist', '1986-2005')
pre_had_eneb_hist = import_gcm('pr', 'eneb', 'hist', '1986-2005')
pre_reg_eneb_rcp26 = import_rcm('pr', 'eneb', 'rcp26', '2080-2099')
pre_reg_eneb_rcp85 = import_rcm('pr', 'eneb', 'rcp85', '2080-2099')
pre_had_eneb_rcp26 = import_gcm('pr', 'eneb', 'rcp26', '2080-2099')
pre_had_eneb_rcp85 = import_gcm('pr', 'eneb', 'rcp85', '2080-2099')
pre_eneb = [pre_cru_eneb_hist, pre_reg_eneb_hist, pre_had_eneb_hist, pre_reg_eneb_rcp26, pre_had_eneb_rcp26, pre_reg_eneb_rcp85, pre_had_eneb_rcp85]

pre_cru_matopiba_hist = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')
pre_reg_matopiba_hist = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
pre_had_matopiba_hist = import_gcm('pr', 'matopiba', 'hist', '1986-2005')
pre_reg_matopiba_rcp26 = import_rcm('pr', 'matopiba', 'rcp26', '2080-2099')
pre_reg_matopiba_rcp85 = import_rcm('pr', 'matopiba', 'rcp85', '2080-2099')
pre_had_matopiba_rcp26 = import_gcm('pr', 'matopiba', 'rcp26', '2080-2099')
pre_had_matopiba_rcp85 = import_gcm('pr', 'matopiba', 'rcp85', '2080-2099')
pre_matopiba = [pre_cru_matopiba_hist, pre_reg_matopiba_hist, pre_had_matopiba_hist, pre_reg_matopiba_rcp26, pre_had_matopiba_rcp26, pre_reg_matopiba_rcp85, pre_had_matopiba_rcp85]

tmp_cru_samz_hist = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
tas_reg_samz_hist = import_rcm('tas', 'samz', 'hist', '1986-2005')
tas_had_samz_hist = import_gcm('tas', 'samz', 'hist', '1986-2005')
tas_reg_samz_rcp26 = import_rcm('tas', 'samz', 'rcp26', '2080-2099')
tas_reg_samz_rcp85 = import_rcm('tas', 'samz', 'rcp85', '2080-2099')
tas_had_samz_rcp26 = import_gcm('tas', 'samz', 'rcp26', '2080-2099')
tas_had_samz_rcp85 = import_gcm('tas', 'samz', 'rcp85', '2080-2099')
tas_samz = [tmp_cru_samz_hist, np.nanmean(tas_reg_samz_hist, axis=1), tas_had_samz_hist, np.nanmean(tas_reg_samz_rcp26, axis=1), tas_had_samz_rcp26, np.nanmean(tas_reg_samz_rcp85, axis=1), tas_had_samz_rcp85]

tmp_cru_eneb_hist = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
tas_reg_eneb_hist = import_rcm('tas', 'eneb', 'hist', '1986-2005')
tas_had_eneb_hist = import_gcm('tas', 'eneb', 'hist', '1986-2005')
tas_reg_eneb_rcp26 = import_rcm('tas', 'eneb', 'rcp26', '2080-2099')
tas_reg_eneb_rcp85 = import_rcm('tas', 'eneb', 'rcp85', '2080-2099')
tas_had_eneb_rcp26 = import_gcm('tas', 'eneb', 'rcp26', '2080-2099')
tas_had_eneb_rcp85 = import_gcm('tas', 'eneb', 'rcp85', '2080-2099')
tas_eneb = [tmp_cru_eneb_hist, np.nanmean(tas_reg_eneb_hist, axis=1), tas_had_eneb_hist, np.nanmean(tas_reg_eneb_rcp26, axis=1), tas_had_eneb_rcp26, np.nanmean(tas_reg_eneb_rcp85, axis=1), tas_had_eneb_rcp85]

tmp_cru_matopiba_hist = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')
tas_reg_matopiba_hist = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
tas_had_matopiba_hist = import_gcm('tas', 'matopiba', 'hist', '1986-2005')
tas_reg_matopiba_rcp26 = import_rcm('tas', 'matopiba', 'rcp26', '2080-2099')
tas_reg_matopiba_rcp85 = import_rcm('tas', 'matopiba', 'rcp85', '2080-2099')
tas_had_matopiba_rcp26 = import_gcm('tas', 'matopiba', 'rcp26', '2080-2099')
tas_had_matopiba_rcp85 = import_gcm('tas', 'matopiba', 'rcp85', '2080-2099')
tas_matopiba = [tmp_cru_matopiba_hist, np.nanmean(tas_reg_matopiba_hist, axis=1), tas_had_matopiba_hist, np.nanmean(tas_reg_matopiba_rcp26, axis=1), tas_had_matopiba_rcp26, np.nanmean(tas_reg_matopiba_rcp85, axis=1), tas_had_matopiba_rcp85]

# Plot models 
fig = plt.figure()
time = np.arange(1, 8)
colors = ['black', 'gray', 'dimgray', 'blue', 'blue', 'red', 'red']

ax1 = fig.add_subplot(3, 2, 1)
box = plt.boxplot(pre_samz, notch=True, patch_artist=True)
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
for median in box['medians']:
    median.set(color='black', linewidth=1)
for flier in box['fliers']:
    flier.set(marker='+', color='black', alpha=1)
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8) 
plt.xticks(time, (u'CRU', u'RegCM4.7', u'HadGEM2-ES', u'RegCM4.7', u'HadGEM2-ES', u'RegCM4.7', u'HadGEM2-ES'), fontsize=6)
plt.ylim(0, 20)
plt.yticks(np.arange(0, 24, 4), fontsize=8)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.axhspan(16, 20, alpha=0.5, color='silver')
plt.axhline(16, lw=1., linestyle='-', color='black')
plt.axvline(3.5, lw=1., linestyle='-', color='black')
plt.axvline(5.5, lw=1., linestyle='-', color='black')
plt.text(0.54, 17., '         Hist            RCP2.6     RCP8.5', fontweight='bold', fontsize=8)
plt.grid(True, which='major', linestyle='--')
                  
ax2 = fig.add_subplot(3, 2, 2)
box = plt.boxplot(tas_samz, notch=True, patch_artist=True)
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
for median in box['medians']:
    median.set(color='black', linewidth=1)
for flier in box['fliers']:
    flier.set(marker='+', color='black', alpha=1)
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, (u'CRU', u'RegCM4.7', u'HadGEM2-ES', u'RegCM4.7', u'HadGEM2-ES', u'RegCM4.7', u'HadGEM2-ES'), fontsize=6)
plt.ylim(22, 42)
plt.yticks(np.arange(22, 46, 4), fontsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.axhspan(38, 42, alpha=0.5, color='silver')
plt.axhline(38, lw=1., linestyle='-', color='black')
plt.axvline(3.5, lw=1., linestyle='-', color='black', alpha=2)
plt.axvline(5.5, lw=1., linestyle='-', color='black', alpha=2)
plt.text(0.54, 39., '         Hist            RCP2.6     RCP8.5', fontweight='bold', fontsize=8)
plt.grid(True, which='major', linestyle='--')
                      
ax3 = fig.add_subplot(3, 2, 3)
box = plt.boxplot(pre_eneb, notch=True, patch_artist=True)
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
for median in box['medians']:
    median.set(color='black', linewidth=1)
for flier in box['fliers']:
    flier.set(marker='+', color='black', alpha=1)
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Precipitation (mm d⁻¹)', fontsize=8)
plt.xticks(time, (u'CRU', u'RegCM4.7', u'HadGEM2-ES', u'RegCM4.7', u'HadGEM2-ES', u'RegCM4.7', u'HadGEM2-ES'), fontsize=6)
plt.ylim(0, 20)
plt.yticks(np.arange(0, 24, 4), fontsize=8)
plt.setp(ax3.get_xticklabels(), visible=False)
plt.axhspan(16, 20, alpha=0.5, color='silver')
plt.axhline(16, lw=1., linestyle='-', color='black')
plt.axvline(3.5, lw=1., linestyle='-', color='black', alpha=2)
plt.axvline(5.5, lw=1., linestyle='-', color='black', alpha=2)
plt.text(0.54, 17., '         Hist            RCP2.6     RCP8.5', fontweight='bold', fontsize=8)
plt.grid(True, which='major', linestyle='--')
                               
ax4 = fig.add_subplot(3, 2, 4)
box = plt.boxplot(tas_eneb, notch=True, patch_artist=True)
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
for median in box['medians']:
    median.set(color='black', linewidth=1)
for flier in box['fliers']:
    flier.set(marker='+', color='black', alpha=1)
plt.title(u'E)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Temperature (°C)', fontsize=8)
plt.xticks(time, (u'CRU', u'RegCM4.7', u'HadGEM2-ES', u'RegCM4.7', u'HadGEM2-ES', u'RegCM4.7', u'HadGEM2-ES'), fontsize=8)
plt.ylim(22, 42)
plt.yticks(np.arange(22, 44, 4), fontsize=8)
plt.setp(ax4.get_xticklabels(), visible=False)
plt.axhspan(38, 42, alpha=0.5, color='silver')
plt.axhline(38, lw=1., linestyle='-', color='black')
plt.axvline(3.5, lw=1., linestyle='-', color='black', alpha=2)
plt.axvline(5.5, lw=1., linestyle='-', color='black', alpha=2)
plt.text(0.54, 39., '         Hist            RCP2.6     RCP8.5', fontweight='bold', fontsize=8)
plt.grid(True, which='major', linestyle='--')
              
ax5 = fig.add_subplot(3, 2, 5)
box = plt.boxplot(pre_matopiba, notch=True, patch_artist=True)
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
for median in box['medians']:
    median.set(color='black', linewidth=1)
for flier in box['fliers']:
    flier.set(marker='+', color='black', alpha=1)
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.xlabel('Dataset', fontsize=8)
plt.xticks(time, (u'CRU', u'RegCM4.7', u'HadGEM2-ES', u'RegCM4.7', u'HadGEM2-ES', u'RegCM4.7', u'HadGEM2-ES'), rotation=90, fontsize=8)
plt.ylim(0, 20)
plt.yticks(np.arange(0, 24, 4), fontsize=8)
plt.axhspan(16, 20, alpha=0.5, color='silver')
plt.axhline(16, lw=1., linestyle='-', color='black')
plt.axvline(3.5, lw=1., linestyle='-', color='black', alpha=2)
plt.axvline(5.5, lw=1., linestyle='-', color='black', alpha=2)
plt.text(0.54, 17.,  '         Hist            RCP2.6     RCP8.5', fontweight='bold', fontsize=8)
plt.grid(True, which='major', linestyle='--')
                              
ax6 = fig.add_subplot(3, 2, 6)
box = plt.boxplot(tas_matopiba, notch=True, patch_artist=True)
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
for median in box['medians']:
    median.set(color='black', linewidth=1)
for flier in box['fliers']:
    flier.set(marker='+', color='black', alpha=1)
plt.title(u'F)', loc='left', fontweight='bold', fontsize=8)
plt.xlabel('Dataset', fontsize=8)
plt.xticks(time, (u'CRU', u'RegCM4.7', u'HadGEM2-ES', u'RegCM4.7', u'HadGEM2-ES', u'RegCM4.7', u'HadGEM2-ES'), rotation=90, fontsize=8)
plt.ylim(22, 42)
plt.yticks(np.arange(22, 44, 4), fontsize=8)
plt.axhspan(38, 42, alpha=0.5, color='silver')
plt.axhline(38, lw=1., linestyle='-', color='black')
plt.axvline(3.5, lw=1., linestyle='-', color='black', alpha=2)
plt.axvline(5.5, lw=1., linestyle='-', color='black', alpha=2)
plt.text(0.54, 39.,  '         Hist            RCP2.6     RCP8.5', fontweight='bold', fontsize=8)
plt.grid(True, which='major', linestyle='--')
             
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_boxplot_reg_had_rcp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()


