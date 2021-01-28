# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot annual cycle graphics from Reg and Had models end obs database"

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


def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	rcm = []
	for mon in range(1, 12 + 1):
		rcm.append(np.nanmean(value[mon::12], axis=0))
	
	return rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]	
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	gcm = []
	for mon in range(1, 12 + 1):
		gcm.append(np.nanmean(value[mon::12], axis=0))
	
	return gcm
	
 
def comp_diff_rcp_hist(rcp, hist):
	
	diff_rcp_hist = []
	for mon in range(0, 12):
		diff = rcp[mon] - hist[mon]
		diff_rcp_hist.append(diff)
	
	return diff_rcp_hist
	
	             
# Import regcm exps model end obs database climatology
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
	
diff_pre_reg_samz_rcp26_hist = comp_diff_rcp_hist(pre_reg_samz_rcp26, pre_reg_samz_hist)
diff_pre_reg_samz_rcp85_hist = comp_diff_rcp_hist(pre_reg_samz_rcp85, pre_reg_samz_hist)
diff_pre_had_samz_rcp26_hist = comp_diff_rcp_hist(pre_had_samz_rcp26, pre_had_samz_hist)
diff_pre_had_samz_rcp85_hist = comp_diff_rcp_hist(pre_had_samz_rcp85, pre_had_samz_hist)

diff_pre_reg_eneb_rcp26_hist = comp_diff_rcp_hist(pre_reg_eneb_rcp26, pre_reg_eneb_hist)
diff_pre_reg_eneb_rcp85_hist = comp_diff_rcp_hist(pre_reg_eneb_rcp85, pre_reg_eneb_hist)
diff_pre_had_eneb_rcp26_hist = comp_diff_rcp_hist(pre_had_eneb_rcp26, pre_had_eneb_hist)
diff_pre_had_eneb_rcp85_hist = comp_diff_rcp_hist(pre_had_eneb_rcp85, pre_had_eneb_hist)

diff_pre_reg_matopiba_rcp26_hist = comp_diff_rcp_hist(pre_reg_matopiba_rcp26, pre_reg_matopiba_hist)
diff_pre_reg_matopiba_rcp85_hist = comp_diff_rcp_hist(pre_reg_matopiba_rcp85, pre_reg_matopiba_hist)
diff_pre_had_matopiba_rcp26_hist = comp_diff_rcp_hist(pre_had_matopiba_rcp26, pre_had_matopiba_hist)
diff_pre_had_matopiba_rcp85_hist = comp_diff_rcp_hist(pre_had_matopiba_rcp85, pre_had_matopiba_hist)

# Temperature
diff_tas_reg_samz_rcp26_hist = comp_diff_rcp_hist(np.nanmean(tas_reg_samz_rcp26, axis=1), np.nanmean(tas_reg_samz_hist, axis=1))
diff_tas_reg_samz_rcp85_hist = comp_diff_rcp_hist(np.nanmean(tas_reg_samz_rcp85, axis=1), np.nanmean(tas_reg_samz_hist, axis=1))
diff_tas_had_samz_rcp26_hist = comp_diff_rcp_hist(tas_had_samz_rcp26, tas_had_samz_hist)
diff_tas_had_samz_rcp85_hist = comp_diff_rcp_hist(tas_had_samz_rcp85, tas_had_samz_hist)

diff_tas_reg_eneb_rcp26_hist = comp_diff_rcp_hist(np.nanmean(tas_reg_eneb_rcp26, axis=1), np.nanmean(tas_reg_eneb_hist, axis=1))
diff_tas_reg_eneb_rcp85_hist = comp_diff_rcp_hist(np.nanmean(tas_reg_eneb_rcp85, axis=1), np.nanmean(tas_reg_eneb_hist, axis=1))
diff_tas_had_eneb_rcp26_hist = comp_diff_rcp_hist(tas_had_eneb_rcp26, tas_had_eneb_hist)
diff_tas_had_eneb_rcp85_hist = comp_diff_rcp_hist(tas_had_eneb_rcp85, tas_had_eneb_hist)

diff_tas_reg_matopiba_rcp26_hist = comp_diff_rcp_hist(np.nanmean(tas_reg_matopiba_rcp26, axis=1), np.nanmean(tas_reg_matopiba_hist, axis=1))
diff_tas_reg_matopiba_rcp85_hist = comp_diff_rcp_hist(np.nanmean(tas_reg_matopiba_rcp85, axis=1), np.nanmean(tas_reg_matopiba_hist, axis=1))
diff_tas_had_matopiba_rcp26_hist = comp_diff_rcp_hist(tas_had_matopiba_rcp26, tas_had_matopiba_hist)
diff_tas_had_matopiba_rcp85_hist = comp_diff_rcp_hist(tas_had_matopiba_rcp85, tas_had_matopiba_hist)

# Plot model end obs data climatology
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)

ax1 = fig.add_subplot(3, 2, 1)
plt.axhline(0, linewidth=1., linestyle='-', color='black')
annual_cycle1 = ax1.plot(time, diff_pre_reg_samz_rcp26_hist, time, diff_pre_reg_samz_rcp85_hist, 
time, diff_pre_had_samz_rcp26_hist, time, diff_pre_had_samz_rcp85_hist)
plt.title(u'A) SAMZ', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
ax1.set_ylim(-2, 2)
plt.yticks(np.arange(-2, 3, 1))
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4 = annual_cycle1
plt.setp(l1, linewidth=1, color='blue', marker='o', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l2, linewidth=1, color='red', marker='o', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l3, linewidth=1, color='blue', marker='s', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l4, linewidth=1, color='red', marker='s', markeredgecolor='black', markersize=6, linestyle='dashed')
legend = ['Reg RCP26-Hist', 'Reg RCP85-Hist', 'Had RCP26-Hist', 'Had RCP85-Hist']
plt.legend(annual_cycle1, legend, loc=9, shadow=True, ncol=2, fontsize=6)

ax2 = fig.add_subplot(3, 2, 2)
plt.axhline(0, linewidth=1., linestyle='-', color='black')
annual_cycle2 = ax2.plot(time, diff_tas_reg_samz_rcp26_hist, time, diff_tas_reg_samz_rcp85_hist, 
time, diff_tas_had_samz_rcp26_hist, time, diff_tas_had_samz_rcp85_hist)
plt.title(u'D) SAMZ', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
ax2.set_ylim(-1, 9)
plt.yticks(np.arange(-1, 10, 2))
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4 = annual_cycle2
plt.setp(l1, linewidth=1, color='blue', marker='o', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l2, linewidth=1, color='red', marker='o', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l3, linewidth=1, color='blue', marker='s', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l4, linewidth=1, color='red', marker='s', markeredgecolor='black', markersize=6, linestyle='dashed')

ax3 = fig.add_subplot(3, 2, 3)
plt.axhline(0, linewidth=1., linestyle='-', color='black')
annual_cycle3 = ax3.plot(time, diff_pre_reg_eneb_rcp26_hist, time, diff_pre_reg_eneb_rcp85_hist, 
time, diff_pre_had_eneb_rcp26_hist, time, diff_pre_had_eneb_rcp85_hist)
plt.title(u'C) ENEB', fontweight='bold')
plt.ylabel(u'Precipitation change (mm d⁻¹)', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
ax3.set_ylim(-2, 2)
plt.yticks(np.arange(-2, 3, 1))
plt.setp(ax3.get_xticklabels(), visible=False)
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4 = annual_cycle3
plt.setp(l1, linewidth=1, color='blue', marker='o', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l2, linewidth=1, color='red', marker='o', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l3, linewidth=1, color='blue', marker='s', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l4, linewidth=1, color='red', marker='s', markeredgecolor='black', markersize=6, linestyle='dashed')

ax4 = fig.add_subplot(3, 2, 4)
plt.axhline(0, linewidth=1., linestyle='-', color='black')
annual_cycle4 = ax4.plot(time, diff_tas_reg_eneb_rcp26_hist, time, diff_tas_reg_eneb_rcp85_hist, 
time, diff_tas_had_eneb_rcp26_hist, time, diff_tas_had_eneb_rcp85_hist)
plt.title(u'E) ENEB', fontweight='bold')
plt.ylabel(u'Temperature change (°C)', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
ax4.set_ylim(-1, 9)
plt.yticks(np.arange(-1, 10, 2))
plt.setp(ax4.get_xticklabels(), visible=False)
ax4.xaxis.grid(True, which='major', linestyle='--')
ax4.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4 = annual_cycle4
plt.setp(l1, linewidth=1, color='blue', marker='o', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l2, linewidth=1, color='red', marker='o', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l3, linewidth=1, color='blue', marker='s', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l4, linewidth=1, color='red', marker='s', markeredgecolor='black', markersize=6, linestyle='dashed')

ax5 = fig.add_subplot(3, 2, 5)
plt.axhline(0, linewidth=1., linestyle='-', color='black')
annual_cycle5 = ax5.plot(time, diff_pre_reg_matopiba_rcp26_hist, time, diff_pre_reg_matopiba_rcp85_hist, 
time, diff_pre_had_matopiba_rcp26_hist, time, diff_pre_had_matopiba_rcp85_hist)
plt.title(u'C) MATOPIBA', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.xlabel('Months', fontweight='bold')
ax5.set_ylim(-2, 2)
plt.yticks(np.arange(-2, 3, 1))
plt.setp(ax2.get_xticklabels(), visible=False)
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4 = annual_cycle5
plt.setp(l1, linewidth=1, color='blue', marker='o', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l2, linewidth=1, color='red', marker='o', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l3, linewidth=1, color='blue', marker='s', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l4, linewidth=1, color='red', marker='s', markeredgecolor='black', markersize=6, linestyle='dashed')

ax6 = fig.add_subplot(3, 2, 6)
plt.axhline(0, linewidth=1., linestyle='-', color='black')
annual_cycle6 = ax6.plot(time, diff_tas_reg_matopiba_rcp26_hist, time, diff_tas_reg_matopiba_rcp85_hist, 
time, diff_tas_had_matopiba_rcp26_hist, time, diff_tas_had_matopiba_rcp85_hist)
plt.title(u'F) MATOPIBA', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.xlabel('Months', fontweight='bold')
ax6.set_ylim(-1, 9)
plt.yticks(np.arange(-1, 10, 2))
ax6.xaxis.grid(True, which='major', linestyle='--')
ax6.yaxis.grid(True, which='major', linestyle='--')
l1, l2, l3, l4 = annual_cycle6
plt.setp(l1, linewidth=1, color='blue', marker='o', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l2, linewidth=1, color='red', marker='o', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l3, linewidth=1, color='blue', marker='s', markeredgecolor='black', markersize=6, linestyle='dashed')
plt.setp(l4, linewidth=1, color='red', marker='s', markeredgecolor='black', markersize=6, linestyle='dashed')

fig.tight_layout()
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.35, hspace=0.35)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_annual_cycle_reg_had_rcp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()






