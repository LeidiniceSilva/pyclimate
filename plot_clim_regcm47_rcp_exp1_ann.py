# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot annual climatology from regcm47 and hadgem models to rcp"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from scipy.stats import t
from scipy.stats import norm
from comp_statist_indices import compute_relative_change
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
	

def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp1/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	rcm = []
		
	for mon in range(0, 11 + 1):
		rcm.append(np.nanmean(value[mon::12], axis=0))

	return rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp1/{0}'.format(exp)
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]	
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	gcm = []
	
	for mon in range(0, 11 + 1):
		gcm.append(np.nanmean(value[mon::12], axis=0))
		
	return gcm
	

def comp_diff_rcp_hist_pre(rcp, hist):
	
	rc_rcp_hist = []
	for mon in range(0, 12):
		p1 = (rcp[mon] - hist[mon]) / hist[mon]
		p2 = p1 * 100
		rc_rcp_hist.append(p2)
	
	return rc_rcp_hist

    
def comp_diff_rcp_hist_tas(rcp, hist):
	
	diff_rcp_hist = []
	for mon in range(0, 12):
		diff = rcp[mon] - hist[mon]
		diff_rcp_hist.append(diff)
	
	return diff_rcp_hist
	
	      
# Import models 
# Precipitation
pre_reg_samz_hist = import_rcm('pr', 'samz', 'hist', '1986-2005')
pre_had_samz_hist = import_gcm('pr', 'samz', 'hist', '1986-2005')
pre_reg_eneb_hist = import_rcm('pr', 'eneb', 'hist', '1986-2005')
pre_had_eneb_hist = import_gcm('pr', 'eneb', 'hist', '1986-2005')
pre_reg_matopiba_hist = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
pre_had_matopiba_hist = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

pre_reg_samz_rcp26 = import_rcm('pr', 'samz', 'rcp26', '2080-2099')
pre_had_samz_rcp26 = import_gcm('pr', 'samz', 'rcp26', '2080-2099')
pre_reg_eneb_rcp26 = import_rcm('pr', 'eneb', 'rcp26', '2080-2099')
pre_had_eneb_rcp26 = import_gcm('pr', 'eneb', 'rcp26', '2080-2099')
pre_reg_matopiba_rcp26 = import_rcm('pr', 'matopiba', 'rcp26', '2080-2099')
pre_had_matopiba_rcp26 = import_gcm('pr', 'matopiba', 'rcp26', '2080-2099')

pre_reg_samz_rcp85 = import_rcm('pr', 'samz', 'rcp85', '2080-2099')
pre_had_samz_rcp85 = import_gcm('pr', 'samz', 'rcp85', '2080-2099')
pre_reg_eneb_rcp85 = import_rcm('pr', 'eneb', 'rcp85', '2080-2099')
pre_had_eneb_rcp85 = import_gcm('pr', 'eneb', 'rcp85', '2080-2099')
pre_reg_matopiba_rcp85 = import_rcm('pr', 'matopiba', 'rcp85', '2080-2099')
pre_had_matopiba_rcp85 = import_gcm('pr', 'matopiba', 'rcp85', '2080-2099')

# Temperature
tas_reg_samz_hist = import_rcm('tas', 'samz', 'hist', '1986-2005')
tas_had_samz_hist = import_gcm('tas', 'samz', 'hist', '1986-2005')
tas_reg_eneb_hist = import_rcm('tas', 'eneb', 'hist', '1986-2005')
tas_had_eneb_hist = import_gcm('tas', 'eneb', 'hist', '1986-2005')
tas_reg_matopiba_hist = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
tas_had_matopiba_hist = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

tas_reg_samz_rcp26 = import_rcm('tas', 'samz', 'rcp26', '2080-2099')
tas_had_samz_rcp26 = import_gcm('tas', 'samz', 'rcp26', '2080-2099')
tas_reg_eneb_rcp26 = import_rcm('tas', 'eneb', 'rcp26', '2080-2099')
tas_had_eneb_rcp26 = import_gcm('tas', 'eneb', 'rcp26', '2080-2099')
tas_reg_matopiba_rcp26 = import_rcm('tas', 'matopiba', 'rcp26', '2080-2099')
tas_had_matopiba_rcp26 = import_gcm('tas', 'matopiba', 'rcp26', '2080-2099')

tas_reg_samz_rcp85 = import_rcm('tas', 'samz', 'rcp85', '2080-2099')
tas_had_samz_rcp85 = import_gcm('tas', 'samz', 'rcp85', '2080-2099')
tas_reg_eneb_rcp85 = import_rcm('tas', 'eneb', 'rcp85', '2080-2099')
tas_had_eneb_rcp85 = import_gcm('tas', 'eneb', 'rcp85', '2080-2099')
tas_reg_matopiba_rcp85 = import_rcm('tas', 'matopiba', 'rcp85', '2080-2099')
tas_had_matopiba_rcp85 = import_gcm('tas', 'matopiba', 'rcp85', '2080-2099')

# Compute difference rcp - hist
diff_pre_reg_samz_rcp26_hist = comp_diff_rcp_hist_pre(pre_reg_samz_rcp26, pre_reg_samz_hist)
diff_pre_reg_samz_rcp85_hist = comp_diff_rcp_hist_pre(pre_reg_samz_rcp85, pre_reg_samz_hist)
diff_pre_had_samz_rcp26_hist = comp_diff_rcp_hist_pre(pre_had_samz_rcp26, pre_had_samz_hist)
diff_pre_had_samz_rcp85_hist = comp_diff_rcp_hist_pre(pre_had_samz_rcp85, pre_had_samz_hist)

diff_pre_reg_eneb_rcp26_hist = comp_diff_rcp_hist_pre(pre_reg_eneb_rcp26, pre_reg_eneb_hist)
diff_pre_reg_eneb_rcp85_hist = comp_diff_rcp_hist_pre(pre_reg_eneb_rcp85, pre_reg_eneb_hist)
diff_pre_had_eneb_rcp26_hist = comp_diff_rcp_hist_pre(pre_had_eneb_rcp26, pre_had_eneb_hist)
diff_pre_had_eneb_rcp85_hist = comp_diff_rcp_hist_pre(pre_had_eneb_rcp85, pre_had_eneb_hist)

diff_pre_reg_matopiba_rcp26_hist = comp_diff_rcp_hist_pre(pre_reg_matopiba_rcp26, pre_reg_matopiba_hist)
diff_pre_reg_matopiba_rcp85_hist = comp_diff_rcp_hist_pre(pre_reg_matopiba_rcp85, pre_reg_matopiba_hist)
diff_pre_had_matopiba_rcp26_hist = comp_diff_rcp_hist_pre(pre_had_matopiba_rcp26, pre_had_matopiba_hist)
diff_pre_had_matopiba_rcp85_hist = comp_diff_rcp_hist_pre(pre_had_matopiba_rcp85, pre_had_matopiba_hist)

# Temperature
diff_tas_reg_samz_rcp26_hist = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_samz_rcp26, axis=1), np.nanmean(tas_reg_samz_hist, axis=1))
diff_tas_reg_samz_rcp85_hist = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_samz_rcp85, axis=1), np.nanmean(tas_reg_samz_hist, axis=1))
diff_tas_had_samz_rcp26_hist = comp_diff_rcp_hist_tas(tas_had_samz_rcp26, tas_had_samz_hist)
diff_tas_had_samz_rcp85_hist = comp_diff_rcp_hist_tas(tas_had_samz_rcp85, tas_had_samz_hist)

diff_tas_reg_eneb_rcp26_hist = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_eneb_rcp26, axis=1), np.nanmean(tas_reg_eneb_hist, axis=1))
diff_tas_reg_eneb_rcp85_hist = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_eneb_rcp85, axis=1), np.nanmean(tas_reg_eneb_hist, axis=1))
diff_tas_had_eneb_rcp26_hist = comp_diff_rcp_hist_tas(tas_had_eneb_rcp26, tas_had_eneb_hist)
diff_tas_had_eneb_rcp85_hist = comp_diff_rcp_hist_tas(tas_had_eneb_rcp85, tas_had_eneb_hist)

diff_tas_reg_matopiba_rcp26_hist = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_matopiba_rcp26, axis=1), np.nanmean(tas_reg_matopiba_hist, axis=1))
diff_tas_reg_matopiba_rcp85_hist = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_matopiba_rcp85, axis=1), np.nanmean(tas_reg_matopiba_hist, axis=1))
diff_tas_had_matopiba_rcp26_hist = comp_diff_rcp_hist_tas(tas_had_matopiba_rcp26, tas_had_matopiba_hist)
diff_tas_had_matopiba_rcp85_hist = comp_diff_rcp_hist_tas(tas_had_matopiba_rcp85, tas_had_matopiba_hist)

# Plot models
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)

ax = fig.add_subplot(3, 2, 1)
ax.plot(time, diff_pre_reg_samz_rcp26_hist, label='RegCM4.7 RCP2.6',   markerfacecolor='lightgray', markeredgecolor='black', color='white', marker='o')
ax.plot(time, diff_pre_reg_samz_rcp85_hist, label='RegCM4.7 RCP8.5',   markerfacecolor='dimgray',   markeredgecolor='black', color='white', marker='o')
ax.plot(time, diff_pre_had_samz_rcp26_hist, label='HadGEM2-ES RCP2.6', markerfacecolor='lightgray', markeredgecolor='black', color='white', marker='s')
ax.plot(time, diff_pre_had_samz_rcp85_hist, label='HadGEM2-ES RCP8.5', markerfacecolor='dimgray',   markeredgecolor='black', color='white', marker='s')
plt.title(u'A) SAMZ', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(-80, 80)
plt.yticks(np.arange(-80, 100, 20), fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)
plt.legend(loc=9, handlelength=0.50, handleheight=0.50, shadow=True, ncol=2, fontsize=5.5)

ax = fig.add_subplot(3, 2, 2)
ax.plot(time, diff_tas_reg_samz_rcp26_hist, markerfacecolor='lightgray', markeredgecolor='black', color='white', marker='o')
ax.plot(time, diff_tas_reg_samz_rcp85_hist, markerfacecolor='dimgray',   markeredgecolor='black', color='white', marker='o')
ax.plot(time, diff_tas_had_samz_rcp26_hist, markerfacecolor='lightgray', markeredgecolor='black', color='white', marker='s')
ax.plot(time, diff_tas_had_samz_rcp85_hist, markerfacecolor='dimgray',   markeredgecolor='black', color='white', marker='s')
plt.title(u'D) SAMZ', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(-8, 8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)

ax = fig.add_subplot(3, 2, 3)
ax.plot(time, diff_pre_reg_eneb_rcp26_hist, label='RegCM4.7 RCP2.6',   markerfacecolor='lightgray', markeredgecolor='black', color='white', marker='o')
ax.plot(time, diff_pre_reg_eneb_rcp85_hist, label='RegCM4.7 RCP8.5',   markerfacecolor='dimgray',   markeredgecolor='black', color='white', marker='o')
ax.plot(time, diff_pre_had_eneb_rcp26_hist, label='HadGEM2-ES RCP2.6', markerfacecolor='lightgray', markeredgecolor='black', color='white', marker='s')
ax.plot(time, diff_pre_had_eneb_rcp85_hist, label='HadGEM2-ES RCP8.5', markerfacecolor='dimgray',   markeredgecolor='black', color='white', marker='s')
plt.title(u'B) ENEB', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Precipitation change(%)', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(-80, 80)
plt.yticks(np.arange(-80, 100, 20), fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)

ax = fig.add_subplot(3, 2, 4)
ax.plot(time, diff_tas_reg_eneb_rcp26_hist, markerfacecolor='lightgray', markeredgecolor='black', color='white', marker='o')
ax.plot(time, diff_tas_reg_eneb_rcp85_hist, markerfacecolor='dimgray',   markeredgecolor='black', color='white', marker='o')
ax.plot(time, diff_tas_had_eneb_rcp26_hist, markerfacecolor='lightgray', markeredgecolor='black', color='white', marker='s')
ax.plot(time, diff_tas_had_eneb_rcp85_hist, markerfacecolor='dimgray',   markeredgecolor='black', color='white', marker='s')
plt.title(u'E) ENEB', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Temperature change (°C)', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(-8, 8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.setp(ax.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)

ax = fig.add_subplot(3, 2, 5)
ax.plot(time, diff_pre_reg_matopiba_rcp26_hist, label='RegCM4.7 RCP2.6',   markerfacecolor='lightgray', markeredgecolor='black', color='white', marker='o')
ax.plot(time, diff_pre_reg_matopiba_rcp85_hist, label='RegCM4.7 RCP8.5',   markerfacecolor='dimgray',   markeredgecolor='black', color='white', marker='o')
ax.plot(time, diff_pre_had_matopiba_rcp26_hist, label='HadGEM2-ES RCP2.6', markerfacecolor='lightgray', markeredgecolor='black', color='white', marker='s')
ax.plot(time, diff_pre_had_matopiba_rcp85_hist, label='HadGEM2-ES RCP8.5', markerfacecolor='dimgray',   markeredgecolor='black', color='white', marker='s')
plt.title(u'C) MATOPIBA', loc='left', fontweight='bold', fontsize=8)
plt.xlabel('Months', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(-80, 80)
plt.yticks(np.arange(-80, 100, 20), fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

ax = fig.add_subplot(3, 2, 6)
ax.plot(time, diff_tas_reg_matopiba_rcp26_hist, markerfacecolor='lightgray', markeredgecolor='black', color='white', marker='o')
ax.plot(time, diff_tas_reg_matopiba_rcp85_hist, markerfacecolor='dimgray',   markeredgecolor='black', color='white', marker='o')
ax.plot(time, diff_tas_had_matopiba_rcp26_hist, markerfacecolor='lightgray', markeredgecolor='black', color='white', marker='s')
ax.plot(time, diff_tas_had_matopiba_rcp85_hist, markerfacecolor='dimgray',   markeredgecolor='black', color='white', marker='s')
plt.title(u'F) MATOPIBA', loc='left', fontweight='bold', fontsize=8)
plt.xlabel('Months', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.ylim(-8, 8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.axhline(0, linewidth=1., linestyle='--', color='black')
plt.grid(True, which='major', linestyle='--')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_ann_diff_reg_had_rcp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()






