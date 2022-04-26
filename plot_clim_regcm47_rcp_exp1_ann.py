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
	ci_p5 = []
	ci_p95 = []
		
	for mon in range(1, 12 + 1):
		rcm.append(np.nanmean(value[mon::12], axis=0))
		
		rcm_p5  = norm.ppf(0.05, loc=np.nanmean(value[mon::12], axis=0), scale=np.nanstd(value[mon::12], axis=0))
		ci_p5.append(rcm_p5)
		
		rcm_p95  = norm.ppf(0.95, loc=np.nanmean(value[mon::12], axis=0), scale=np.nanstd(value[mon::12], axis=0))
		ci_p95.append(rcm_p95)
	
	return rcm, ci_p5, ci_p95


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp1/{0}'.format(exp)
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]	
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	gcm = []
	ci_p5 = []
	ci_p95 = []
	
	for mon in range(1, 12 + 1):
		gcm.append(np.nanmean(value[mon::12], axis=0))
		
		gcm_p5  = norm.ppf(0.05, loc=np.nanmean(value[mon::12], axis=0), scale=np.nanstd(value[mon::12], axis=0))
		ci_p5.append(gcm_p5)
		
		gcm_p95  = norm.ppf(0.95, loc=np.nanmean(value[mon::12], axis=0), scale=np.nanstd(value[mon::12], axis=0))
		ci_p95.append(gcm_p95)
	
	return gcm, ci_p5, ci_p95
	

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
pre_reg_samz_hist, pre_reg_samz_hist_p5, pre_reg_samz_hist_p95 = import_rcm('pr', 'samz', 'hist', '1986-2005')
pre_had_samz_hist, pre_had_samz_hist_p5, pre_had_samz_hist_p95 = import_gcm('pr', 'samz', 'hist', '1986-2005')
pre_reg_eneb_hist, pre_reg_eneb_hist_p5, pre_reg_eneb_hist_p95 = import_rcm('pr', 'eneb', 'hist', '1986-2005')
pre_had_eneb_hist, pre_had_eneb_hist_p5, pre_had_eneb_hist_p95 = import_gcm('pr', 'eneb', 'hist', '1986-2005')
pre_reg_matopiba_hist, pre_reg_matopiba_hist_p5, pre_reg_matopiba_hist_p95 = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
pre_had_matopiba_hist, pre_had_matopiba_hist_p5, pre_had_matopiba_hist_p95 = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

pre_reg_samz_rcp26, pre_reg_samz_rcp26_p5, pre_reg_samz_rcp26_p95 = import_rcm('pr', 'samz', 'rcp26', '2080-2099')
pre_had_samz_rcp26, pre_had_samz_rcp26_p5, pre_had_samz_rcp26_p95 = import_gcm('pr', 'samz', 'rcp26', '2080-2099')
pre_reg_eneb_rcp26, pre_reg_eneb_rcp26_p5, pre_reg_eneb_rcp26_p95 = import_rcm('pr', 'eneb', 'rcp26', '2080-2099')
pre_had_eneb_rcp26, pre_had_eneb_rcp26_p5, pre_had_eneb_rcp26_p95 = import_gcm('pr', 'eneb', 'rcp26', '2080-2099')
pre_reg_matopiba_rcp26, pre_reg_matopiba_rcp26_p5, pre_reg_matopiba_rcp26_p95 = import_rcm('pr', 'matopiba', 'rcp26', '2080-2099')
pre_had_matopiba_rcp26, pre_had_matopiba_rcp26_p5, pre_had_matopiba_rcp26_p95 = import_gcm('pr', 'matopiba', 'rcp26', '2080-2099')

pre_reg_samz_rcp85, pre_reg_samz_rcp85_p5, pre_reg_samz_rcp85_p95 = import_rcm('pr', 'samz', 'rcp85', '2080-2099')
pre_had_samz_rcp85, pre_had_samz_rcp85_p5, pre_had_samz_rcp85_p95 = import_gcm('pr', 'samz', 'rcp85', '2080-2099')
pre_reg_eneb_rcp85, pre_reg_eneb_rcp85_p5, pre_reg_eneb_rcp85_p95 = import_rcm('pr', 'eneb', 'rcp85', '2080-2099')
pre_had_eneb_rcp85, pre_had_eneb_rcp85_p5, pre_had_eneb_rcp85_p95 = import_gcm('pr', 'eneb', 'rcp85', '2080-2099')
pre_reg_matopiba_rcp85,pre_reg_matopiba_rcp85_p5, pre_reg_matopiba_rcp85_p95 = import_rcm('pr', 'matopiba', 'rcp85', '2080-2099')
pre_had_matopiba_rcp85, pre_had_matopiba_rcp85_p5, pre_had_matopiba_rcp85_p95 = import_gcm('pr', 'matopiba', 'rcp85', '2080-2099')

# Temperature
tas_reg_samz_hist, tas_reg_samz_hist_p5, tas_reg_samz_hist_p95 = import_rcm('tas', 'samz', 'hist', '1986-2005')
tas_had_samz_hist, tas_had_samz_hist_p5, tas_had_samz_hist_p95 = import_gcm('tas', 'samz', 'hist', '1986-2005')
tas_reg_eneb_hist, tas_reg_eneb_hist_p5, tas_reg_eneb_hist_p95 = import_rcm('tas', 'eneb', 'hist', '1986-2005')
tas_had_eneb_hist, tas_had_eneb_hist_p5, tas_had_eneb_hist_p95 = import_gcm('tas', 'eneb', 'hist', '1986-2005')
tas_reg_matopiba_hist, tas_reg_matopiba_hist_p5, tas_reg_matopiba_hist_p95 = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
tas_had_matopiba_hist, tas_had_matopiba_hist_p5, tas_had_matopiba_hist_p95 = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

tas_reg_samz_rcp26, tas_reg_samz_rcp26_p5, tas_reg_samz_rcp26_p95 = import_rcm('tas', 'samz', 'rcp26', '2080-2099')
tas_had_samz_rcp26, tas_had_samz_rcp26_p5, tas_had_samz_rcp26_p95 = import_gcm('tas', 'samz', 'rcp26', '2080-2099')
tas_reg_eneb_rcp26, tas_reg_eneb_rcp26_p5, tas_reg_eneb_rcp26_p95 = import_rcm('tas', 'eneb', 'rcp26', '2080-2099')
tas_had_eneb_rcp26, tas_had_eneb_rcp26_p5, tas_had_eneb_rcp26_p95 = import_gcm('tas', 'eneb', 'rcp26', '2080-2099')
tas_reg_matopiba_rcp26, tas_reg_matopiba_rcp26_p5, tas_reg_matopiba_rcp26_p95 = import_rcm('tas', 'matopiba', 'rcp26', '2080-2099')
tas_had_matopiba_rcp26, tas_had_matopiba_rcp26_p5, tas_had_matopiba_rcp26_p95 = import_gcm('tas', 'matopiba', 'rcp26', '2080-2099')

tas_reg_samz_rcp85, tas_reg_samz_rcp85_p5, tas_reg_samz_rcp85_p95 = import_rcm('tas', 'samz', 'rcp85', '2080-2099')
tas_had_samz_rcp85, tas_had_samz_rcp85_p5, tas_had_samz_rcp85_p95 = import_gcm('tas', 'samz', 'rcp85', '2080-2099')
tas_reg_eneb_rcp85, tas_reg_eneb_rcp85_p5, tas_reg_eneb_rcp85_p95 = import_rcm('tas', 'eneb', 'rcp85', '2080-2099')
tas_had_eneb_rcp85, tas_had_eneb_rcp85_p5, tas_had_eneb_rcp85_p95 = import_gcm('tas', 'eneb', 'rcp85', '2080-2099')
tas_reg_matopiba_rcp85, tas_reg_matopiba_rcp85_p5, tas_reg_matopiba_rcp85_p95 = import_rcm('tas', 'matopiba', 'rcp85', '2080-2099')
tas_had_matopiba_rcp85, tas_had_matopiba_rcp85_p5, tas_had_matopiba_rcp85_p95 = import_gcm('tas', 'matopiba', 'rcp85', '2080-2099')

# Compute difference rcp - hist
diff_pre_reg_samz_rcp26_hist = comp_diff_rcp_hist_pre(pre_reg_samz_rcp26, pre_reg_samz_hist)
diff_pre_reg_samz_rcp85_hist = comp_diff_rcp_hist_pre(pre_reg_samz_rcp85, pre_reg_samz_hist)
diff_pre_had_samz_rcp26_hist = comp_diff_rcp_hist_pre(pre_had_samz_rcp26, pre_had_samz_hist)
diff_pre_had_samz_rcp85_hist = comp_diff_rcp_hist_pre(pre_had_samz_rcp85, pre_had_samz_hist)
diff_pre_reg_samz_rcp26_hist_p5 = comp_diff_rcp_hist_pre(pre_reg_samz_rcp26_p5, pre_reg_samz_hist_p5)
diff_pre_reg_samz_rcp85_hist_p5 = comp_diff_rcp_hist_pre(pre_reg_samz_rcp85_p5, pre_reg_samz_hist_p5)
diff_pre_had_samz_rcp26_hist_p5 = comp_diff_rcp_hist_pre(pre_had_samz_rcp26_p5, pre_had_samz_hist_p5)
diff_pre_had_samz_rcp85_hist_p5 = comp_diff_rcp_hist_pre(pre_had_samz_rcp85_p5, pre_had_samz_hist_p5)
diff_pre_reg_samz_rcp26_hist_p95 = comp_diff_rcp_hist_pre(pre_reg_samz_rcp26_p95, pre_reg_samz_hist_p95)
diff_pre_reg_samz_rcp85_hist_p95 = comp_diff_rcp_hist_pre(pre_reg_samz_rcp85_p95, pre_reg_samz_hist_p95)
diff_pre_had_samz_rcp26_hist_p95 = comp_diff_rcp_hist_pre(pre_had_samz_rcp26_p95, pre_had_samz_hist_p95)
diff_pre_had_samz_rcp85_hist_p95 = comp_diff_rcp_hist_pre(pre_had_samz_rcp85_p95, pre_had_samz_hist_p95)

diff_pre_reg_eneb_rcp26_hist = comp_diff_rcp_hist_pre(pre_reg_eneb_rcp26, pre_reg_eneb_hist)
diff_pre_reg_eneb_rcp85_hist = comp_diff_rcp_hist_pre(pre_reg_eneb_rcp85, pre_reg_eneb_hist)
diff_pre_had_eneb_rcp26_hist = comp_diff_rcp_hist_pre(pre_had_eneb_rcp26, pre_had_eneb_hist)
diff_pre_had_eneb_rcp85_hist = comp_diff_rcp_hist_pre(pre_had_eneb_rcp85, pre_had_eneb_hist)
diff_pre_reg_eneb_rcp26_hist_p5 = comp_diff_rcp_hist_pre(pre_reg_eneb_rcp26_p5, pre_reg_eneb_hist_p5)
diff_pre_reg_eneb_rcp85_hist_p5 = comp_diff_rcp_hist_pre(pre_reg_eneb_rcp85_p5, pre_reg_eneb_hist_p5)
diff_pre_had_eneb_rcp26_hist_p5 = comp_diff_rcp_hist_pre(pre_had_eneb_rcp26_p5, pre_had_eneb_hist_p5)
diff_pre_had_eneb_rcp85_hist_p5 = comp_diff_rcp_hist_pre(pre_had_eneb_rcp85_p5, pre_had_eneb_hist_p5)
diff_pre_reg_eneb_rcp26_hist_p95 = comp_diff_rcp_hist_pre(pre_reg_eneb_rcp26_p95, pre_reg_eneb_hist_p95)
diff_pre_reg_eneb_rcp85_hist_p95 = comp_diff_rcp_hist_pre(pre_reg_eneb_rcp85_p95, pre_reg_eneb_hist_p95)
diff_pre_had_eneb_rcp26_hist_p95 = comp_diff_rcp_hist_pre(pre_had_eneb_rcp26_p95, pre_had_eneb_hist_p95)
diff_pre_had_eneb_rcp85_hist_p95 = comp_diff_rcp_hist_pre(pre_had_eneb_rcp85_p95, pre_had_eneb_hist_p95)

diff_pre_reg_matopiba_rcp26_hist = comp_diff_rcp_hist_pre(pre_reg_matopiba_rcp26, pre_reg_matopiba_hist)
diff_pre_reg_matopiba_rcp85_hist = comp_diff_rcp_hist_pre(pre_reg_matopiba_rcp85, pre_reg_matopiba_hist)
diff_pre_had_matopiba_rcp26_hist = comp_diff_rcp_hist_pre(pre_had_matopiba_rcp26, pre_had_matopiba_hist)
diff_pre_had_matopiba_rcp85_hist = comp_diff_rcp_hist_pre(pre_had_matopiba_rcp85, pre_had_matopiba_hist)
diff_pre_reg_matopiba_rcp26_hist_p5 = comp_diff_rcp_hist_pre(pre_reg_matopiba_rcp26_p5, pre_reg_matopiba_hist_p5)
diff_pre_reg_matopiba_rcp85_hist_p5 = comp_diff_rcp_hist_pre(pre_reg_matopiba_rcp85_p5, pre_reg_matopiba_hist_p5)
diff_pre_had_matopiba_rcp26_hist_p5 = comp_diff_rcp_hist_pre(pre_had_matopiba_rcp26_p5, pre_had_matopiba_hist_p5)
diff_pre_had_matopiba_rcp85_hist_p5 = comp_diff_rcp_hist_pre(pre_had_matopiba_rcp85_p5, pre_had_matopiba_hist_p5)
diff_pre_reg_matopiba_rcp26_hist_p95 = comp_diff_rcp_hist_pre(pre_reg_matopiba_rcp26_p95, pre_reg_matopiba_hist_p95)
diff_pre_reg_matopiba_rcp85_hist_p95 = comp_diff_rcp_hist_pre(pre_reg_matopiba_rcp85_p95, pre_reg_matopiba_hist_p95)
diff_pre_had_matopiba_rcp26_hist_p95 = comp_diff_rcp_hist_pre(pre_had_matopiba_rcp26_p95, pre_had_matopiba_hist_p95)
diff_pre_had_matopiba_rcp85_hist_p95 = comp_diff_rcp_hist_pre(pre_had_matopiba_rcp85_p95, pre_had_matopiba_hist_p95)

# Temperature
diff_tas_reg_samz_rcp26_hist = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_samz_rcp26, axis=1), np.nanmean(tas_reg_samz_hist, axis=1))
diff_tas_reg_samz_rcp85_hist = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_samz_rcp85, axis=1), np.nanmean(tas_reg_samz_hist, axis=1))
diff_tas_had_samz_rcp26_hist = comp_diff_rcp_hist_tas(tas_had_samz_rcp26, tas_had_samz_hist)
diff_tas_had_samz_rcp85_hist = comp_diff_rcp_hist_tas(tas_had_samz_rcp85, tas_had_samz_hist)
diff_tas_reg_samz_rcp26_hist_p5 = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_samz_rcp26_p5, axis=1), np.nanmean(tas_reg_samz_hist_p5, axis=1))
diff_tas_reg_samz_rcp85_hist_p5 = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_samz_rcp85_p5, axis=1), np.nanmean(tas_reg_samz_hist_p5, axis=1))
diff_tas_had_samz_rcp26_hist_p5 = comp_diff_rcp_hist_tas(tas_had_samz_rcp26_p5, tas_had_samz_hist_p5)
diff_tas_had_samz_rcp85_hist_p5 = comp_diff_rcp_hist_tas(tas_had_samz_rcp85_p5, tas_had_samz_hist_p5)
diff_tas_reg_samz_rcp26_hist_p95 = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_samz_rcp26_p95, axis=1), np.nanmean(tas_reg_samz_hist_p95, axis=1))
diff_tas_reg_samz_rcp85_hist_p95 = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_samz_rcp85_p95, axis=1), np.nanmean(tas_reg_samz_hist_p95, axis=1))
diff_tas_had_samz_rcp26_hist_p95 = comp_diff_rcp_hist_tas(tas_had_samz_rcp26_p95, tas_had_samz_hist_p95)
diff_tas_had_samz_rcp85_hist_p95 = comp_diff_rcp_hist_tas(tas_had_samz_rcp85_p95, tas_had_samz_hist_p95)

diff_tas_reg_eneb_rcp26_hist = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_eneb_rcp26, axis=1), np.nanmean(tas_reg_eneb_hist, axis=1))
diff_tas_reg_eneb_rcp85_hist = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_eneb_rcp85, axis=1), np.nanmean(tas_reg_eneb_hist, axis=1))
diff_tas_had_eneb_rcp26_hist = comp_diff_rcp_hist_tas(tas_had_eneb_rcp26, tas_had_eneb_hist)
diff_tas_had_eneb_rcp85_hist = comp_diff_rcp_hist_tas(tas_had_eneb_rcp85, tas_had_eneb_hist)
diff_tas_reg_eneb_rcp26_hist_p5 = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_eneb_rcp26_p5, axis=1), np.nanmean(tas_reg_eneb_hist_p5, axis=1))
diff_tas_reg_eneb_rcp85_hist_p5 = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_eneb_rcp85_p5, axis=1), np.nanmean(tas_reg_eneb_hist_p5, axis=1))
diff_tas_had_eneb_rcp26_hist_p5 = comp_diff_rcp_hist_tas(tas_had_eneb_rcp26_p5, tas_had_eneb_hist_p5)
diff_tas_had_eneb_rcp85_hist_p5 = comp_diff_rcp_hist_tas(tas_had_eneb_rcp85_p5, tas_had_eneb_hist_p5)
diff_tas_reg_eneb_rcp26_hist_p95 = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_eneb_rcp26_p95, axis=1), np.nanmean(tas_reg_eneb_hist_p95, axis=1))
diff_tas_reg_eneb_rcp85_hist_p95 = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_eneb_rcp85_p95, axis=1), np.nanmean(tas_reg_eneb_hist_p95, axis=1))
diff_tas_had_eneb_rcp26_hist_p95 = comp_diff_rcp_hist_tas(tas_had_eneb_rcp26_p95, tas_had_eneb_hist_p95)
diff_tas_had_eneb_rcp85_hist_p95 = comp_diff_rcp_hist_tas(tas_had_eneb_rcp85_p95, tas_had_eneb_hist_p95)

diff_tas_reg_matopiba_rcp26_hist = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_matopiba_rcp26, axis=1), np.nanmean(tas_reg_matopiba_hist, axis=1))
diff_tas_reg_matopiba_rcp85_hist = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_matopiba_rcp85, axis=1), np.nanmean(tas_reg_matopiba_hist, axis=1))
diff_tas_had_matopiba_rcp26_hist = comp_diff_rcp_hist_tas(tas_had_matopiba_rcp26, tas_had_matopiba_hist)
diff_tas_had_matopiba_rcp85_hist = comp_diff_rcp_hist_tas(tas_had_matopiba_rcp85, tas_had_matopiba_hist)
diff_tas_reg_matopiba_rcp26_hist_p5 = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_matopiba_rcp26_p5, axis=1), np.nanmean(tas_reg_matopiba_hist_p5, axis=1))
diff_tas_reg_matopiba_rcp85_hist_p5 = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_matopiba_rcp85_p5, axis=1), np.nanmean(tas_reg_matopiba_hist_p5, axis=1))
diff_tas_had_matopiba_rcp26_hist_p5 = comp_diff_rcp_hist_tas(tas_had_matopiba_rcp26_p5, tas_had_matopiba_hist_p5)
diff_tas_had_matopiba_rcp85_hist_p5 = comp_diff_rcp_hist_tas(tas_had_matopiba_rcp85_p5, tas_had_matopiba_hist_p5)
diff_tas_reg_matopiba_rcp26_hist_p95 = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_matopiba_rcp26_p95, axis=1), np.nanmean(tas_reg_matopiba_hist_p95, axis=1))
diff_tas_reg_matopiba_rcp85_hist_p95 = comp_diff_rcp_hist_tas(np.nanmean(tas_reg_matopiba_rcp85_p95, axis=1), np.nanmean(tas_reg_matopiba_hist_p95, axis=1))
diff_tas_had_matopiba_rcp26_hist_p95 = comp_diff_rcp_hist_tas(tas_had_matopiba_rcp26_p95, tas_had_matopiba_hist_p95)
diff_tas_had_matopiba_rcp85_hist_p95 = comp_diff_rcp_hist_tas(tas_had_matopiba_rcp85_p95, tas_had_matopiba_hist_p95)

# Plot models
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)

ax1 = fig.add_subplot(3, 2, 1)
annual_cycle1 = ax1.plot(time, diff_pre_reg_samz_rcp26_hist, time, diff_pre_reg_samz_rcp85_hist, 
time, diff_pre_had_samz_rcp26_hist, time, diff_pre_had_samz_rcp85_hist)
plt.title(u'A) SAMZ', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
ax1.set_ylim(-80, 80)
plt.yticks(np.arange(-80, 100, 20), fontsize=8)
plt.setp(ax1.get_xticklabels(), visible=False)
l1, l2, l3, l4 = annual_cycle1
plt.setp(l1, linewidth=1.5, color='blue', linestyle='-')
plt.setp(l2, linewidth=1.5, color='red', linestyle='-')
plt.setp(l3, linewidth=1.5, color='blue', linestyle='--')
plt.setp(l4, linewidth=1.5, color='red', linestyle='--')
plt.legend(annual_cycle1, ['RegCM4.7 RCP2.6', 'RegCM4.7 RCP8.5', 'HadGEM2-ES RCP2.6', 'HadGEM2-ES RCP8.5'], loc=9, shadow=True, ncol=2, fontsize=5.5)
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axvline(3.5, linewidth=1., linestyle='-', color='black')
plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
plt.fill_between(time, diff_pre_reg_samz_rcp26_hist_p5, diff_pre_reg_samz_rcp26_hist_p95, facecolor='blue', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_pre_had_samz_rcp26_hist_p5, diff_pre_had_samz_rcp26_hist_p95, facecolor='blue', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_pre_reg_samz_rcp85_hist_p5, diff_pre_reg_samz_rcp85_hist_p95, facecolor='red', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_pre_had_samz_rcp85_hist_p5, diff_pre_had_samz_rcp85_hist_p95, facecolor='red', alpha=0.3, interpolate=True)	
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)

ax2 = fig.add_subplot(3, 2, 2)
annual_cycle2 = ax2.plot(time, diff_tas_reg_samz_rcp26_hist, time, diff_tas_reg_samz_rcp85_hist, 
time, diff_tas_had_samz_rcp26_hist, time, diff_tas_had_samz_rcp85_hist)
plt.title(u'D) SAMZ', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
ax2.set_ylim(-8, 8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)
l1, l2, l3, l4 = annual_cycle2
plt.setp(l1, linewidth=1.5, color='blue', linestyle='-')
plt.setp(l2, linewidth=1.5, color='red', linestyle='-')
plt.setp(l3, linewidth=1.5, color='blue', linestyle='--')
plt.setp(l4, linewidth=1.5, color='red', linestyle='--')
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axvline(3.5, linewidth=1., linestyle='-', color='black')
plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
plt.fill_between(time, diff_tas_reg_samz_rcp26_hist_p5, diff_tas_reg_samz_rcp26_hist_p95, facecolor='blue', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_tas_had_samz_rcp26_hist_p5, diff_tas_had_samz_rcp26_hist_p95, facecolor='blue', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_tas_reg_samz_rcp85_hist_p5, diff_tas_reg_samz_rcp85_hist_p95, facecolor='red', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_tas_had_samz_rcp85_hist_p5, diff_tas_had_samz_rcp85_hist_p95, facecolor='red', alpha=0.3, interpolate=True)	
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)

ax3 = fig.add_subplot(3, 2, 3)
annual_cycle3 = ax3.plot(time, diff_pre_reg_eneb_rcp26_hist, time, diff_pre_reg_eneb_rcp85_hist, 
time, diff_pre_had_eneb_rcp26_hist, time, diff_pre_had_eneb_rcp85_hist)
plt.title(u'B) ENEB', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Mudança de precipitação (%)', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
ax3.set_ylim(-80, 80)
plt.yticks(np.arange(-80, 100, 20), fontsize=8)
plt.setp(ax3.get_xticklabels(), visible=False)
l1, l2, l3, l4 = annual_cycle3
plt.setp(l1, linewidth=1.5, color='blue', linestyle='-')
plt.setp(l2, linewidth=1.5, color='red', linestyle='-')
plt.setp(l3, linewidth=1.5, color='blue', linestyle='--')
plt.setp(l4, linewidth=1.5, color='red', linestyle='--')
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axvline(3.5, linewidth=1., linestyle='-', color='black')
plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
plt.fill_between(time, diff_pre_reg_eneb_rcp26_hist_p5, diff_pre_reg_eneb_rcp26_hist_p95, facecolor='blue', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_pre_had_eneb_rcp26_hist_p5, diff_pre_had_eneb_rcp26_hist_p95, facecolor='blue', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_pre_reg_eneb_rcp85_hist_p5, diff_pre_reg_eneb_rcp85_hist_p95, facecolor='red', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_pre_had_eneb_rcp85_hist_p5, diff_pre_had_eneb_rcp85_hist_p95, facecolor='red', alpha=0.3, interpolate=True)	
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)

ax4 = fig.add_subplot(3, 2, 4)
annual_cycle4 = ax4.plot(time, diff_tas_reg_eneb_rcp26_hist, time, diff_tas_reg_eneb_rcp85_hist, 
time, diff_tas_had_eneb_rcp26_hist, time, diff_tas_had_eneb_rcp85_hist)
plt.title(u'E) ENEB', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Mudança de temperatura (°C)', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
ax4.set_ylim(-8, 8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
plt.setp(ax4.get_xticklabels(), visible=False)
l1, l2, l3, l4 = annual_cycle4
plt.setp(l1, linewidth=1.5, color='blue', linestyle='-')
plt.setp(l2, linewidth=1.5, color='red', linestyle='-')
plt.setp(l3, linewidth=1.5, color='blue', linestyle='--')
plt.setp(l4, linewidth=1.5, color='red', linestyle='--')
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axvline(3.5, linewidth=1., linestyle='-', color='black')
plt.axvline(8.5, linewidth=1., linestyle='-', color='black')	
plt.fill_between(time, diff_tas_reg_eneb_rcp26_hist_p5, diff_tas_reg_eneb_rcp26_hist_p95, facecolor='blue', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_tas_had_eneb_rcp26_hist_p5, diff_tas_had_eneb_rcp26_hist_p95, facecolor='blue', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_tas_reg_eneb_rcp85_hist_p5, diff_tas_reg_eneb_rcp85_hist_p95, facecolor='red', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_tas_had_eneb_rcp85_hist_p5, diff_tas_had_eneb_rcp85_hist_p95, facecolor='red', alpha=0.3, interpolate=True)	
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)

ax5 = fig.add_subplot(3, 2, 5)
annual_cycle5 = ax5.plot(time, diff_pre_reg_matopiba_rcp26_hist, time, diff_pre_reg_matopiba_rcp85_hist, 
time, diff_pre_had_matopiba_rcp26_hist, time, diff_pre_had_matopiba_rcp85_hist)
plt.title(u'C) MATOPIBA', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.xlabel('Meses', fontsize=8)
ax5.set_ylim(-80, 80)
plt.yticks(np.arange(-80, 100, 20), fontsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)
l1, l2, l3, l4 = annual_cycle5
plt.setp(l1, linewidth=1.5, color='blue', linestyle='-')
plt.setp(l2, linewidth=1.5, color='red', linestyle='-')
plt.setp(l3, linewidth=1.5, color='blue', linestyle='--')
plt.setp(l4, linewidth=1.5, color='red', linestyle='--')
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axvline(3.5, linewidth=1., linestyle='-', color='black')
plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
plt.fill_between(time, diff_pre_reg_matopiba_rcp26_hist_p5, diff_pre_reg_matopiba_rcp26_hist_p95, facecolor='blue', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_pre_had_matopiba_rcp26_hist_p5, diff_pre_had_matopiba_rcp26_hist_p95, facecolor='blue', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_pre_reg_matopiba_rcp85_hist_p5, diff_pre_reg_matopiba_rcp85_hist_p95, facecolor='red', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_pre_had_matopiba_rcp85_hist_p5, diff_pre_had_matopiba_rcp85_hist_p95, facecolor='red', alpha=0.3, interpolate=True)	
plt.grid(True, which='major', linestyle='--')

ax6 = fig.add_subplot(3, 2, 6)
annual_cycle6 = ax6.plot(time, diff_tas_reg_matopiba_rcp26_hist, time, diff_tas_reg_matopiba_rcp85_hist, 
time, diff_tas_had_matopiba_rcp26_hist, time, diff_tas_had_matopiba_rcp85_hist)
plt.title(u'F) MATOPIBA', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.xlabel('Meses', fontsize=8)
ax6.set_ylim(-8, 8)
plt.yticks(np.arange(-8, 10, 2), fontsize=8)
l1, l2, l3, l4 = annual_cycle6
plt.setp(l1, linewidth=1.5, color='blue', linestyle='-')
plt.setp(l2, linewidth=1.5, color='red', linestyle='-')
plt.setp(l3, linewidth=1.5, color='blue', linestyle='--')
plt.setp(l4, linewidth=1.5, color='red', linestyle='--')
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axvline(3.5, linewidth=1., linestyle='-', color='black')
plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
plt.fill_between(time, diff_tas_reg_matopiba_rcp26_hist_p5, diff_tas_reg_matopiba_rcp26_hist_p95, facecolor='blue', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_tas_had_matopiba_rcp26_hist_p5, diff_tas_had_matopiba_rcp26_hist_p95, facecolor='blue', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_tas_reg_matopiba_rcp85_hist_p5, diff_tas_reg_matopiba_rcp85_hist_p95, facecolor='red', alpha=0.3, interpolate=True)	
plt.fill_between(time, diff_tas_had_matopiba_rcp85_hist_p5, diff_tas_had_matopiba_rcp85_hist_p95, facecolor='red', alpha=0.3, interpolate=True)	
plt.grid(True, which='major', linestyle='--')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_annual_cycle_reg_had_rcp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()






