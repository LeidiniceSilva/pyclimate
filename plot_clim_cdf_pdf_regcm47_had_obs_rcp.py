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
from scipy.stats import norm
from scipy.stats import t
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def cdf_function(data):
	
	x = np.linspace(np.min(data), np.max(data))
	y = np.nanmean(x)
	z = np.nanstd(x)
	cdf = norm.cdf(x,y,z)
	
	return x, cdf


def pdf_function(data):

	x = np.linspace(np.min(data), np.max(data))
	y = np.nanmean(x)
	z = np.nanstd(x)
	pdf = norm.pdf(x,y,z)
	
	return x, pdf
	

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

# Calculate CDF function
xcdf_pre_reg_samz_rcp26_hist, cdf_pre_reg_samz_rcp26_hist = cdf_function(diff_pre_reg_samz_rcp26_hist)
xcdf_pre_reg_samz_rcp85_hist, cdf_pre_reg_samz_rcp85_hist = cdf_function(diff_pre_reg_samz_rcp85_hist)
xcdf_pre_reg_eneb_rcp26_hist, cdf_pre_reg_eneb_rcp26_hist = cdf_function(diff_pre_reg_eneb_rcp26_hist)
xcdf_pre_reg_eneb_rcp85_hist, cdf_pre_reg_eneb_rcp85_hist = cdf_function(diff_pre_reg_eneb_rcp85_hist)
xcdf_pre_reg_matopiba_rcp26_hist, cdf_pre_reg_matopiba_rcp26_hist = cdf_function(diff_pre_reg_matopiba_rcp26_hist)
xcdf_pre_reg_matopiba_rcp85_hist, cdf_pre_reg_matopiba_rcp85_hist = cdf_function(diff_pre_reg_matopiba_rcp85_hist)
xcdf_pre_had_samz_rcp26_hist, cdf_pre_had_samz_rcp26_hist = cdf_function(diff_pre_had_samz_rcp26_hist)
xcdf_pre_had_samz_rcp85_hist, cdf_pre_had_samz_rcp85_hist = cdf_function(diff_pre_had_samz_rcp85_hist)
xcdf_pre_had_eneb_rcp26_hist, cdf_pre_had_eneb_rcp26_hist = cdf_function(diff_pre_had_eneb_rcp26_hist)
xcdf_pre_had_eneb_rcp85_hist, cdf_pre_had_eneb_rcp85_hist = cdf_function(diff_pre_had_eneb_rcp85_hist)
xcdf_pre_had_matopiba_rcp26_hist, cdf_pre_had_matopiba_rcp26_hist = cdf_function(diff_pre_had_matopiba_rcp26_hist)
xcdf_pre_had_matopiba_rcp85_hist, cdf_pre_had_matopiba_rcp85_hist = cdf_function(diff_pre_had_matopiba_rcp85_hist)

xcdf_tas_reg_samz_rcp26_hist, cdf_tas_reg_samz_rcp26_hist = cdf_function(diff_tas_reg_samz_rcp26_hist)
xcdf_tas_reg_samz_rcp85_hist, cdf_tas_reg_samz_rcp85_hist = cdf_function(diff_tas_reg_samz_rcp85_hist)
xcdf_tas_reg_eneb_rcp26_hist, cdf_tas_reg_eneb_rcp26_hist = cdf_function(diff_tas_reg_eneb_rcp26_hist)
xcdf_tas_reg_eneb_rcp85_hist, cdf_tas_reg_eneb_rcp85_hist = cdf_function(diff_tas_reg_eneb_rcp85_hist)
xcdf_tas_reg_matopiba_rcp26_hist, cdf_tas_reg_matopiba_rcp26_hist = cdf_function(diff_tas_reg_matopiba_rcp26_hist)
xcdf_tas_reg_matopiba_rcp85_hist, cdf_tas_reg_matopiba_rcp85_hist = cdf_function(diff_tas_reg_matopiba_rcp85_hist)
xcdf_tas_had_samz_rcp26_hist, cdf_tas_had_samz_rcp26_hist = cdf_function(diff_tas_had_samz_rcp26_hist)
xcdf_tas_had_samz_rcp85_hist, cdf_tas_had_samz_rcp85_hist = cdf_function(diff_tas_had_samz_rcp85_hist)
xcdf_tas_had_eneb_rcp26_hist, cdf_tas_had_eneb_rcp26_hist = cdf_function(diff_tas_had_eneb_rcp26_hist)
xcdf_tas_had_eneb_rcp85_hist, cdf_tas_had_eneb_rcp85_hist = cdf_function(diff_tas_had_eneb_rcp85_hist)
xcdf_tas_had_matopiba_rcp26_hist, cdf_tas_had_matopiba_rcp26_hist = cdf_function(diff_tas_had_matopiba_rcp26_hist)
xcdf_tas_had_matopiba_rcp85_hist, cdf_tas_had_matopiba_rcp85_hist = cdf_function(diff_tas_had_matopiba_rcp85_hist)

# Calculate PDF function
xpdf_pre_reg_samz_rcp26_hist, pdf_pre_reg_samz_rcp26_hist = pdf_function(diff_pre_reg_samz_rcp26_hist)
xpdf_pre_reg_samz_rcp85_hist, pdf_pre_reg_samz_rcp85_hist = pdf_function(diff_pre_reg_samz_rcp85_hist)
xpdf_pre_reg_eneb_rcp26_hist, pdf_pre_reg_eneb_rcp26_hist = pdf_function(diff_pre_reg_eneb_rcp26_hist)
xpdf_pre_reg_eneb_rcp85_hist, pdf_pre_reg_eneb_rcp85_hist = pdf_function(diff_pre_reg_eneb_rcp85_hist)
xpdf_pre_reg_matopiba_rcp26_hist, pdf_pre_reg_matopiba_rcp26_hist = pdf_function(diff_pre_reg_matopiba_rcp26_hist)
xpdf_pre_reg_matopiba_rcp85_hist, pdf_pre_reg_matopiba_rcp85_hist = pdf_function(diff_pre_reg_matopiba_rcp85_hist)
xpdf_pre_had_samz_rcp26_hist, pdf_pre_had_samz_rcp26_hist = pdf_function(diff_pre_had_samz_rcp26_hist)
xpdf_pre_had_samz_rcp85_hist, pdf_pre_had_samz_rcp85_hist = pdf_function(diff_pre_had_samz_rcp85_hist)
xpdf_pre_had_eneb_rcp26_hist, pdf_pre_had_eneb_rcp26_hist = pdf_function(diff_pre_had_eneb_rcp26_hist)
xpdf_pre_had_eneb_rcp85_hist, pdf_pre_had_eneb_rcp85_hist = pdf_function(diff_pre_had_eneb_rcp85_hist)
xpdf_pre_had_matopiba_rcp26_hist, pdf_pre_had_matopiba_rcp26_hist = pdf_function(diff_pre_had_matopiba_rcp26_hist)
xpdf_pre_had_matopiba_rcp85_hist, pdf_pre_had_matopiba_rcp85_hist = pdf_function(diff_pre_had_matopiba_rcp85_hist)

xpdf_tas_reg_samz_rcp26_hist, pdf_tas_reg_samz_rcp26_hist = pdf_function(diff_tas_reg_samz_rcp26_hist)
xpdf_tas_reg_samz_rcp85_hist, pdf_tas_reg_samz_rcp85_hist = pdf_function(diff_tas_reg_samz_rcp85_hist)
xpdf_tas_reg_eneb_rcp26_hist, pdf_tas_reg_eneb_rcp26_hist = pdf_function(diff_tas_reg_eneb_rcp26_hist,)
xpdf_tas_reg_eneb_rcp85_hist, pdf_tas_reg_eneb_rcp85_hist = pdf_function(diff_tas_reg_eneb_rcp85_hist)
xpdf_tas_reg_matopiba_rcp26_hist, pdf_pre_reg_matopiba_rcp26_hist = pdf_function(diff_tas_reg_matopiba_rcp26_hist)
xpdf_tas_reg_matopiba_rcp85_hist, pdf_pre_reg_matopiba_rcp85_hist = pdf_function(diff_tas_reg_matopiba_rcp85_hist)
xpdf_tas_had_samz_rcp26_hist, pdf_tas_had_samz_rcp26_hist = pdf_function(diff_tas_had_samz_rcp26_hist)
xpdf_tas_had_samz_rcp85_hist, pdf_tas_had_samz_rcp85_hist = pdf_function(diff_tas_had_samz_rcp85_hist)
xpdf_tas_had_eneb_rcp26_hist, pdf_tas_had_eneb_rcp26_hist = pdf_function(diff_tas_had_eneb_rcp26_hist)
xpdf_tas_had_eneb_rcp85_hist, pdf_tas_had_eneb_rcp85_hist = pdf_function(diff_tas_had_eneb_rcp85_hist)
xpdf_tas_had_matopiba_rcp26_hist, pdf_tas_had_matopiba_rcp26_hist = pdf_function(diff_tas_had_matopiba_rcp26_hist)
xpdf_tas_had_matopiba_rcp85_hist, pdf_tas_had_matopiba_rcp85_hist = pdf_function(diff_tas_had_matopiba_rcp85_hist)

# Plot model end obs data climatology
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)

ax1 = fig.add_subplot(3, 2, 1)
annual_cycle1 = ax1.plot(time, diff_pre_reg_samz_rcp26_hist, time, diff_pre_reg_samz_rcp85_hist, 
time, diff_pre_had_samz_rcp26_hist, time, diff_pre_had_samz_rcp85_hist)
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
ax1.set_ylim(-4, 4)
plt.yticks(np.arange(-4, 6, 2), fontsize=8)
plt.setp(ax1.get_xticklabels(), visible=False)
l1, l2, l3, l4 = annual_cycle1
plt.setp(l1, linewidth=1.5, color='blue', linestyle='-')
plt.setp(l2, linewidth=1.5, color='red', linestyle='-')
plt.setp(l3, linewidth=1.5, color='blue', linestyle='--')
plt.setp(l4, linewidth=1.5, color='red', linestyle='--')
plt.grid(True, which='major', linestyle='--')
plt.legend(annual_cycle1, ['Reg RCP2.6-Hist', 'Reg RCP8.5-Hist', 'Had RCP2.6-Hist', 'Had RCP8.5-Hist'], loc=9, shadow=True, ncol=2, fontsize=6)
plt.fill_between(time, -1, 1, facecolor='slategray', alpha=0.3, interpolate=True)	
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axvline(3.5, linewidth=1., linestyle='-', color='black')
plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
		 
#~ ax1_inset = inset_axes(ax1, width='30%', height='30%', loc=9)
#~ line1, = ax1_inset.plot(xcdf_pre_reg_samz_rcp26_hist, cdf_pre_reg_samz_rcp26_hist, color='blue', label='Reg RCP2.6', linestyle='-', linewidth=1.5)
#~ line1, = ax1_inset.plot(xcdf_pre_reg_samz_rcp85_hist, cdf_pre_reg_samz_rcp85_hist, color='red', label='Reg RCP8.5', linestyle='-', linewidth=1.5)
#~ line1, = ax1_inset.plot(xcdf_pre_had_samz_rcp26_hist, cdf_pre_had_samz_rcp26_hist, color='blue', label='Had RCP2.6', linestyle='--', linewidth=1.5)
#~ line1, = ax1_inset.plot(xcdf_pre_had_samz_rcp85_hist, cdf_pre_had_samz_rcp85_hist, color='red', label='Had RCP8.5', linestyle='--', linewidth=1.5)
#~ ax1_inset.set_ylabel(u'CDF', fontsize=6)
#~ ax1_inset.set_xlabel(u'(mm d⁻¹)', fontsize=6)
#~ ax1_inset.set_xlim(-2, 2)
#~ ax1_inset.set_ylim(0, 1)
#~ ax1_inset.set_xticks(np.arange(-2, 3, 1))
#~ ax1_inset.set_yticks(np.arange(0, 1.5, 0.5))
#~ ax1_inset.xaxis.set_tick_params(labelsize=6)
#~ ax1_inset.yaxis.set_tick_params(labelsize=6)
#~ plt.axvline(0, linewidth=1., linestyle='-', color='black')
#~ plt.grid(True, which='major', linestyle='--')

ax2 = fig.add_subplot(3, 2, 2)
annual_cycle2 = ax2.plot(time, diff_tas_reg_samz_rcp26_hist, time, diff_tas_reg_samz_rcp85_hist, 
time, diff_tas_had_samz_rcp26_hist, time, diff_tas_had_samz_rcp85_hist)
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
ax2.set_ylim(0, 9)
plt.yticks(np.arange(0, 10, 3), fontsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)
l1, l2, l3, l4 = annual_cycle2
plt.setp(l1, linewidth=1.5, color='blue', linestyle='-')
plt.setp(l2, linewidth=1.5, color='red', linestyle='-')
plt.setp(l3, linewidth=1.5, color='blue', linestyle='--')
plt.setp(l4, linewidth=1.5, color='red', linestyle='--')
plt.grid(True, which='major', linestyle='--')
plt.fill_between(time, 0.5, 2.5, facecolor='slategray', alpha=0.3, interpolate=True)
plt.fill_between(time, 5, 9, facecolor='slategray', alpha=0.3, interpolate=True)
plt.axvline(3.5, linewidth=1., linestyle='-', color='black')
plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
		 		 
#~ ax2_inset = inset_axes(ax2, width='30%', height='30%', loc=10)
#~ line2, = ax2_inset.plot(xcdf_tas_reg_samz_rcp26_hist, cdf_tas_reg_samz_rcp26_hist, color='blue', label='Reg RCP2.6', linestyle='-', linewidth=1.5)
#~ line2, = ax2_inset.plot(xcdf_tas_reg_samz_rcp85_hist, cdf_tas_reg_samz_rcp85_hist, color='red', label='Reg RCP8.5', linestyle='-', linewidth=1.5)
#~ line2, = ax2_inset.plot(xcdf_tas_had_samz_rcp26_hist, cdf_tas_had_samz_rcp26_hist, color='blue', label='Had RCP2.6', linestyle='--', linewidth=1.5)
#~ line2, = ax2_inset.plot(xcdf_tas_had_samz_rcp85_hist, cdf_tas_had_samz_rcp85_hist, color='red', label='Had RCP8.5', linestyle='--', linewidth=1.5)
#~ ax2_inset.set_ylabel(u'CDF', fontsize=6)
#~ ax2_inset.set_xlabel(u'(°C)', fontsize=6)
#~ ax2_inset.set_xlim(0, 9)
#~ ax2_inset.set_ylim(0, 1)
#~ ax2_inset.set_xticks(np.arange(0, 10, 3))
#~ ax2_inset.set_yticks(np.arange(0, 1.5, 0.5))
#~ ax2_inset.xaxis.set_tick_params(labelsize=6)
#~ ax2_inset.yaxis.set_tick_params(labelsize=6)
#~ plt.grid(True, which='major', linestyle='--')

ax3 = fig.add_subplot(3, 2, 3)
annual_cycle3 = ax3.plot(time, diff_pre_reg_eneb_rcp26_hist, time, diff_pre_reg_eneb_rcp85_hist, 
time, diff_pre_had_eneb_rcp26_hist, time, diff_pre_had_eneb_rcp85_hist)
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Precipitation change (mm d⁻¹)', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
ax3.set_ylim(-4, 4)
plt.yticks(np.arange(-4, 6, 2), fontsize=8)
plt.setp(ax3.get_xticklabels(), visible=False)
l1, l2, l3, l4 = annual_cycle3
plt.setp(l1, linewidth=1.5, color='blue', linestyle='-')
plt.setp(l2, linewidth=1.5, color='red', linestyle='-')
plt.setp(l3, linewidth=1.5, color='blue', linestyle='--')
plt.setp(l4, linewidth=1.5, color='red', linestyle='--')
plt.grid(True, which='major', linestyle='--')
plt.fill_between(time, -1, 1, facecolor='slategray', alpha=0.3, interpolate=True)
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axvline(3.5, linewidth=1., linestyle='-', color='black')
plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
		 
#~ ax3_inset = inset_axes(ax3, width='30%', height='30%', loc=9)
#~ line3, = ax3_inset.plot(xcdf_pre_reg_eneb_rcp26_hist, cdf_pre_reg_eneb_rcp26_hist, color='blue', label='Reg RCP2.6', linestyle='-', linewidth=1.5)
#~ line3, = ax3_inset.plot(xcdf_pre_reg_eneb_rcp85_hist, cdf_pre_reg_eneb_rcp85_hist, color='red', label='Reg RCP8.5', linestyle='-', linewidth=1.5)
#~ line3, = ax3_inset.plot(xcdf_pre_had_eneb_rcp26_hist, cdf_pre_had_eneb_rcp26_hist, color='blue', label='Had RCP2.6', linestyle='--', linewidth=1.5)
#~ line3, = ax3_inset.plot(xcdf_pre_had_eneb_rcp85_hist, cdf_pre_had_eneb_rcp85_hist, color='red', label='Had RCP8.5', linestyle='--', linewidth=1.5)
#~ ax3_inset.set_xlim(-2, 2)
#~ ax3_inset.set_ylim(0, 1)
#~ ax3_inset.set_xticks(np.arange(-2, 3, 1))
#~ ax3_inset.set_yticks(np.arange(0, 1.5, 0.5))
#~ ax3_inset.xaxis.set_tick_params(labelsize=6)
#~ ax3_inset.yaxis.set_tick_params(labelsize=6)
#~ plt.axvline(0, linewidth=1., linestyle='-', color='black')
#~ plt.grid(True, which='major', linestyle='--')

ax4 = fig.add_subplot(3, 2, 4)
annual_cycle4 = ax4.plot(time, diff_tas_reg_eneb_rcp26_hist, time, diff_tas_reg_eneb_rcp85_hist, 
time, diff_tas_had_eneb_rcp26_hist, time, diff_tas_had_eneb_rcp85_hist)
plt.title(u'E)', loc='left', fontweight='bold', fontsize=8)
plt.ylabel(u'Temperature change (°C)', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
ax4.set_ylim(0, 9)
plt.yticks(np.arange(0, 10, 3), fontsize=8)
plt.setp(ax4.get_xticklabels(), visible=False)
l1, l2, l3, l4 = annual_cycle4
plt.setp(l1, linewidth=1.5, color='blue', linestyle='-')
plt.setp(l2, linewidth=1.5, color='red', linestyle='-')
plt.setp(l3, linewidth=1.5, color='blue', linestyle='--')
plt.setp(l4, linewidth=1.5, color='red', linestyle='--')
plt.grid(True, which='major', linestyle='--')
plt.fill_between(time, 0.5, 2.5, facecolor='slategray', alpha=0.3, interpolate=True)
plt.fill_between(time, 3.5, 5.5, facecolor='slategray', alpha=0.3, interpolate=True)
plt.axvline(3.5, linewidth=1., linestyle='-', color='black')
plt.axvline(8.5, linewidth=1., linestyle='-', color='black')	
	
#~ ax4_inset = inset_axes(ax4, width='30%', height='30%', loc=9)
#~ line4, = ax4_inset.plot(xcdf_tas_reg_eneb_rcp26_hist, cdf_tas_reg_eneb_rcp26_hist, color='blue', label='Reg RCP2.6', linestyle='-', linewidth=1.5)
#~ line4, = ax4_inset.plot(xcdf_tas_reg_eneb_rcp85_hist, cdf_tas_reg_eneb_rcp85_hist, color='red', label='Reg RCP8.5', linestyle='-', linewidth=1.5)
#~ line4, = ax4_inset.plot(xcdf_tas_had_eneb_rcp26_hist, cdf_tas_had_eneb_rcp26_hist, color='blue', label='Had RCP2.6', linestyle='--', linewidth=1.5)
#~ line4, = ax4_inset.plot(xcdf_tas_had_eneb_rcp85_hist, cdf_tas_had_eneb_rcp85_hist, color='red', label='Had RCP8.5', linestyle='--', linewidth=1.5)
#~ ax4_inset.set_xlim(0, 9)
#~ ax4_inset.set_ylim(0, 1)
#~ ax4_inset.set_xticks(np.arange(0, 10, 3))
#~ ax4_inset.set_yticks(np.arange(0, 1.5, 0.5))
#~ ax4_inset.xaxis.set_tick_params(labelsize=6)
#~ ax4_inset.yaxis.set_tick_params(labelsize=6)
#~ plt.grid(True, which='major', linestyle='--')

ax5 = fig.add_subplot(3, 2, 5)
annual_cycle5 = ax5.plot(time, diff_pre_reg_matopiba_rcp26_hist, time, diff_pre_reg_matopiba_rcp85_hist, 
time, diff_pre_had_matopiba_rcp26_hist, time, diff_pre_had_matopiba_rcp85_hist)
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.xlabel('Months', fontsize=8)
ax5.set_ylim(-4, 4)
plt.yticks(np.arange(-4, 6, 2), fontsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)
l1, l2, l3, l4 = annual_cycle5
plt.setp(l1, linewidth=1.5, color='blue', linestyle='-')
plt.setp(l2, linewidth=1.5, color='red', linestyle='-')
plt.setp(l3, linewidth=1.5, color='blue', linestyle='--')
plt.setp(l4, linewidth=1.5, color='red', linestyle='--')
plt.grid(True, which='major', linestyle='--')
plt.fill_between(time, -1, 1, facecolor='slategray', alpha=0.3, interpolate=True)
plt.axhline(0, linewidth=1., linestyle='-', color='black')
plt.axvline(3.5, linewidth=1., linestyle='-', color='black')
plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
		 
#~ ax5_inset = inset_axes(ax5, width='30%', height='30%', loc=9)
#~ line5, = ax5_inset.plot(xcdf_pre_reg_matopiba_rcp26_hist, cdf_pre_reg_matopiba_rcp26_hist, color='blue', label='Reg RCP2.6', linestyle='-', linewidth=1.5)
#~ line5, = ax5_inset.plot(xcdf_pre_reg_matopiba_rcp85_hist, cdf_pre_reg_matopiba_rcp85_hist, color='red', label='Reg RCP8.5', linestyle='-', linewidth=1.5)
#~ line5, = ax5_inset.plot(xcdf_pre_had_matopiba_rcp26_hist, cdf_pre_had_matopiba_rcp26_hist, color='blue', label='Had RCP2.6', linestyle='--', linewidth=1.5)
#~ line5, = ax5_inset.plot(xcdf_pre_had_matopiba_rcp85_hist, cdf_pre_had_matopiba_rcp85_hist, color='red', label='Had RCP8.5', linestyle='--', linewidth=1.5)
#~ ax5_inset.set_xlim(-2, 2)
#~ ax5_inset.set_ylim(0, 1)
#~ ax5_inset.set_xticks(np.arange(-2, 3, 1))
#~ ax5_inset.set_yticks(np.arange(0, 1.5, 0.5))
#~ ax5_inset.xaxis.set_tick_params(labelsize=6)
#~ ax5_inset.yaxis.set_tick_params(labelsize=6)
#~ plt.axvline(0, linewidth=1., linestyle='-', color='black')
#~ plt.grid(True, which='major', linestyle='--')

ax6 = fig.add_subplot(3, 2, 6)
annual_cycle6 = ax6.plot(time, diff_tas_reg_matopiba_rcp26_hist, time, diff_tas_reg_matopiba_rcp85_hist, 
time, diff_tas_had_matopiba_rcp26_hist, time, diff_tas_had_matopiba_rcp85_hist)
plt.title(u'F)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.xlabel('Months', fontsize=8)
ax6.set_ylim(0, 9)
plt.yticks(np.arange(0, 10, 3), fontsize=8)
l1, l2, l3, l4 = annual_cycle6
plt.setp(l1, linewidth=1.5, color='blue', linestyle='-')
plt.setp(l2, linewidth=1.5, color='red', linestyle='-')
plt.setp(l3, linewidth=1.5, color='blue', linestyle='--')
plt.setp(l4, linewidth=1.5, color='red', linestyle='--')
plt.grid(True, which='major', linestyle='--')
plt.fill_between(time, 0.5, 2.5, facecolor='slategray', alpha=0.3, interpolate=True)
plt.fill_between(time, 4.5, 7, facecolor='slategray', alpha=0.3, interpolate=True)
plt.axvline(3.5, linewidth=1., linestyle='-', color='black')
plt.axvline(8.5, linewidth=1., linestyle='-', color='black')
		
#~ ax6_inset = inset_axes(ax6, width='30%', height='30%', loc=10)
#~ line6, = ax6_inset.plot(xcdf_tas_reg_matopiba_rcp26_hist, cdf_tas_reg_matopiba_rcp26_hist, color='blue', label='Reg RCP2.6', linestyle='-', linewidth=1.5)
#~ line6, = ax6_inset.plot(xcdf_tas_reg_matopiba_rcp85_hist, cdf_tas_reg_matopiba_rcp85_hist, color='red', label='Reg RCP8.5', linestyle='-', linewidth=1.5)
#~ line6, = ax6_inset.plot(xcdf_tas_had_matopiba_rcp26_hist, cdf_tas_had_matopiba_rcp26_hist, color='blue', label='Had RCP2.6', linestyle='--', linewidth=1.5)
#~ line6, = ax6_inset.plot(xcdf_tas_had_matopiba_rcp85_hist, cdf_tas_had_matopiba_rcp85_hist, color='red', label='Had RCP8.5', linestyle='--', linewidth=1.5)
#~ ax6_inset.set_xlim(0, 9)
#~ ax6_inset.set_ylim(0, 1)
#~ ax6_inset.set_xticks(np.arange(0, 10, 3))
#~ ax6_inset.set_yticks(np.arange(0, 1.5, 0.5))
#~ ax6_inset.xaxis.set_tick_params(labelsize=6)
#~ ax6_inset.yaxis.set_tick_params(labelsize=6)
#~ plt.grid(True, which='major', linestyle='--')

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.35, hspace=0.35)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_annual_cycle_cdf_reg_had_rcp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()






