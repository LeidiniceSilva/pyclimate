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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from comp_statist_indices import compute_pdf
	

def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return value


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]	
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return value


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
	
	
# Import regcm exps model end obs database climatology
# Precipitation
pre_reg_samz_hist = import_rcm('pr', 'samz', 'hist', '1986-2005')
pre_reg_eneb_hist = import_rcm('pr', 'eneb', 'hist', '1986-2005')
pre_reg_matopiba_hist = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
pre_had_samz_hist = import_gcm('pr', 'samz', 'hist', '1986-2005')
pre_had_eneb_hist= import_gcm('pr', 'eneb', 'hist', '1986-2005')
pre_had_matopiba_hist = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

pre_reg_samz_rcp26 = import_rcm('pr', 'samz', 'rcp26', '2080-2099')
pre_reg_eneb_rcp26 = import_rcm('pr', 'eneb', 'rcp26', '2080-2099')
pre_reg_matopiba_rcp26 = import_rcm('pr', 'matopiba', 'rcp26', '2080-2099')
pre_had_samz_rcp26 = import_gcm('pr', 'samz', 'rcp26', '2080-2099')
pre_had_eneb_rcp26 = import_gcm('pr', 'eneb', 'rcp26', '2080-2099')
pre_had_matopiba_rcp26 = import_gcm('pr', 'matopiba', 'rcp26', '2080-2099')

pre_reg_samz_rcp85 = import_rcm('pr', 'samz', 'rcp85', '2080-2099')
pre_reg_eneb_rcp85 = import_rcm('pr', 'eneb', 'rcp85', '2080-2099')
pre_reg_matopiba_rcp85 = import_rcm('pr', 'matopiba', 'rcp85', '2080-2099')
pre_had_samz_rcp85 = import_gcm('pr', 'samz', 'rcp85', '2080-2099')
pre_had_eneb_rcp85 = import_gcm('pr', 'eneb', 'rcp85', '2080-2099')
pre_had_matopiba_rcp85 = import_gcm('pr', 'matopiba', 'rcp85', '2080-2099')

# Temperature
tas_reg_samz_hist = import_rcm('tas', 'samz', 'hist', '1986-2005')
tas_reg_eneb_hist = import_rcm('tas', 'eneb', 'hist', '1986-2005')
tas_reg_matopiba_hist = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
tas_had_samz_hist = import_gcm('tas', 'samz', 'hist', '1986-2005')
tas_had_eneb_hist = import_gcm('tas', 'eneb', 'hist', '1986-2005')
tas_had_matopiba_hist = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

tas_reg_samz_rcp26 = import_rcm('tas', 'samz', 'rcp26', '2080-2099')
tas_reg_eneb_rcp26 = import_rcm('tas', 'eneb', 'rcp26', '2080-2099')
tas_reg_matopiba_rcp26 = import_rcm('tas', 'matopiba', 'rcp26', '2080-2099')
tas_had_samz_rcp26 = import_gcm('tas', 'samz', 'rcp26', '2080-2099')
tas_had_eneb_rcp26 = import_gcm('tas', 'eneb', 'rcp26', '2080-2099')
tas_had_matopiba_rcp26 = import_gcm('tas', 'matopiba', 'rcp26', '2080-2099')

tas_reg_samz_rcp85 = import_rcm('tas', 'samz', 'rcp85', '2080-2099')
tas_reg_eneb_rcp85 = import_rcm('tas', 'eneb', 'rcp85', '2080-2099')
tas_reg_matopiba_rcp85 = import_rcm('tas', 'matopiba', 'rcp85', '2080-2099')
tas_had_samz_rcp85 = import_gcm('tas', 'samz', 'rcp85', '2080-2099')
tas_had_eneb_rcp85 = import_gcm('tas', 'eneb', 'rcp85', '2080-2099')
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

# Calculate PDF function
# Precipitation
xpdf_pre_reg_samz_rcp26_hist, pdf_pre_reg_samz_rcp26_hist = compute_pdf(diff_pre_reg_samz_rcp26_hist)
xpdf_pre_reg_eneb_rcp26_hist, pdf_pre_reg_eneb_rcp26_hist = compute_pdf(diff_pre_reg_eneb_rcp26_hist)
xpdf_pre_reg_matopiba_rcp26_hist, pdf_pre_reg_matopiba_rcp26_hist = compute_pdf(diff_pre_reg_matopiba_rcp26_hist)

xpdf_pre_reg_samz_rcp85_hist, pdf_pre_reg_samz_rcp85_hist = compute_pdf(diff_pre_reg_samz_rcp85_hist)
xpdf_pre_reg_eneb_rcp85_hist, pdf_pre_reg_eneb_rcp85_hist = compute_pdf(diff_pre_reg_eneb_rcp85_hist)
xpdf_pre_reg_matopiba_rcp85_hist, pdf_pre_reg_matopiba_rcp85_hist = compute_pdf(diff_pre_reg_matopiba_rcp85_hist)

xpdf_pre_had_samz_rcp26_hist, pdf_pre_had_samz_rcp26_hist = compute_pdf(diff_pre_had_samz_rcp26_hist)
xpdf_pre_had_eneb_rcp26_hist, pdf_pre_had_eneb_rcp26_hist = compute_pdf(diff_pre_had_eneb_rcp26_hist)
xpdf_pre_had_matopiba_rcp26_hist, pdf_pre_had_matopiba_rcp26_hist = compute_pdf(diff_pre_had_matopiba_rcp26_hist)

xpdf_pre_had_samz_rcp85_hist, pdf_pre_had_samz_rcp85_hist = compute_pdf(diff_pre_had_samz_rcp85_hist)
xpdf_pre_had_eneb_rcp85_hist, pdf_pre_had_eneb_rcp85_hist = compute_pdf(diff_pre_had_eneb_rcp85_hist)
xpdf_pre_had_matopiba_rcp85_hist, pdf_pre_had_matopiba_rcp85_hist = compute_pdf(diff_pre_had_matopiba_rcp85_hist)

# Temperature
xpdf_tas_reg_samz_rcp26_hist, pdf_tas_reg_samz_rcp26_hist = compute_pdf(diff_tas_reg_samz_rcp26_hist)
xpdf_tas_reg_eneb_rcp26_hist, pdf_tas_reg_eneb_rcp26_hist = compute_pdf(diff_tas_reg_eneb_rcp26_hist)
xpdf_tas_reg_matopiba_rcp26_hist, pdf_tas_reg_matopiba_rcp26_hist = compute_pdf(diff_tas_reg_matopiba_rcp26_hist)

xpdf_tas_reg_samz_rcp85_hist, pdf_tas_reg_samz_rcp85_hist = compute_pdf(diff_tas_reg_samz_rcp85_hist)
xpdf_tas_reg_eneb_rcp85_hist, pdf_tas_reg_eneb_rcp85_hist = compute_pdf(diff_tas_reg_eneb_rcp85_hist)
xpdf_tas_reg_matopiba_rcp85_hist, pdf_tas_reg_matopiba_rcp85_hist = compute_pdf(diff_tas_reg_matopiba_rcp85_hist)

xpdf_tas_had_samz_rcp26_hist, pdf_tas_had_samz_rcp26_hist = compute_pdf(diff_tas_had_samz_rcp26_hist)
xpdf_tas_had_eneb_rcp26_hist, pdf_tas_had_eneb_rcp26_hist = compute_pdf(diff_tas_had_eneb_rcp26_hist)
xpdf_tas_had_matopiba_rcp26_hist, pdf_tas_had_matopiba_rcp26_hist = compute_pdf(diff_tas_had_matopiba_rcp26_hist)

xpdf_tas_had_samz_rcp85_hist, pdf_tas_had_samz_rcp85_hist = compute_pdf(diff_tas_had_samz_rcp85_hist)
xpdf_tas_had_eneb_rcp85_hist, pdf_tas_had_eneb_rcp85_hist = compute_pdf(diff_tas_had_eneb_rcp85_hist)
xpdf_tas_had_matopiba_rcp85_hist, pdf_tas_had_matopiba_rcp85_hist = compute_pdf(diff_tas_had_matopiba_rcp85_hist)

# Plot model end obs data climatology
fig = plt.figure()

ax1 = fig.add_subplot(3, 2, 1)
pdf_line1 = ax1.plot(xpdf_pre_reg_samz_rcp26_hist, pdf_pre_reg_samz_rcp26_hist, color='blue', label='RegCM4.7 RCP2.6', linestyle='-', linewidth=1.5)
pdf_line3 = ax1.plot(xpdf_pre_reg_samz_rcp85_hist, pdf_pre_reg_samz_rcp85_hist, color='red', label='RegCM4.7 RCP8.5', linestyle='-', linewidth=1.5)
pdf_line2 = ax1.plot(xpdf_pre_had_samz_rcp26_hist, pdf_pre_reg_samz_rcp26_hist, color='blue', label='HadGEM2-ES RCP2.6', linestyle='--', linewidth=1.5)
pdf_line4 = ax1.plot(xpdf_pre_had_samz_rcp85_hist, pdf_pre_reg_samz_rcp85_hist, color='red', label='HadGEM2-ES RCP8.5', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_reg_samz_rcp26_hist), color='blue', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_had_samz_rcp26_hist), color='blue', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_reg_samz_rcp85_hist), color='red', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_had_samz_rcp85_hist), color='red', linestyle='--', linewidth=1.5)
plt.axvline(0, color='black', linestyle='-', linewidth=1.5)
plt.axvspan(-100, -10, alpha=0.3, color='red')
plt.axvspan(10, 300, alpha=0.3, color='blue')
plt.title(u'A)', loc='left', fontweight='bold', fontsize=8)
ax1.set_xlim(-100, 300)
ax1.set_ylim(0, 0.03)
ax1.set_xticks(np.arange(-100, 350, 50))
ax1.set_yticks(np.arange(0, 0.04, 0.01))
ax1.xaxis.set_tick_params(labelsize=8)
ax1.yaxis.set_tick_params(labelsize=8)
plt.setp(ax1.get_xticklabels(), visible=False)	     
plt.grid(True, which='major', linestyle='--')
plt.legend(fontsize=6, loc=1, shadow=True, ncol=1)

ax2 = fig.add_subplot(3, 2, 2)
pdf_line1 = ax2.plot(xpdf_tas_reg_samz_rcp26_hist, pdf_tas_reg_samz_rcp26_hist, color='blue', label='RegCM4.7 RCP2.6', linestyle='-', linewidth=1.5)
pdf_line2 = ax2.plot(xpdf_tas_had_samz_rcp26_hist, pdf_tas_reg_samz_rcp26_hist, color='blue', label='HadGEM2-ES RCP2.6', linestyle='--', linewidth=1.5)
pdf_line3 = ax2.plot(xpdf_tas_reg_samz_rcp85_hist, pdf_tas_reg_samz_rcp85_hist, color='red', label='RegCM4.7 RCP8.5', linestyle='-', linewidth=1.5)
pdf_line4 = ax2.plot(xpdf_tas_had_samz_rcp85_hist, pdf_tas_reg_samz_rcp85_hist, color='red', label='HadGEM2-ES RCP8.5', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_reg_samz_rcp26_hist), color='blue', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_had_samz_rcp26_hist), color='blue', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_reg_samz_rcp85_hist), color='red', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_had_samz_rcp85_hist), color='red', linestyle='--', linewidth=1.5)
plt.axvline(0, color='black', linestyle='-', linewidth=1.5)
plt.axvspan(-2, -1, alpha=0.3, color='blue')
plt.axvspan(1, 10, alpha=0.3, color='red')
plt.title(u'D)', loc='left', fontweight='bold', fontsize=8)
ax2.set_xlim(-2, 10)
ax2.set_ylim(0, 0.8)
ax2.set_xticks(np.arange(-2, 12, 2))
ax2.set_yticks(np.arange(0, 1, 0.2))
ax2.xaxis.set_tick_params(labelsize=8)
ax2.yaxis.set_tick_params(labelsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)	     
plt.grid(True, which='major', linestyle='--')

ax3 = fig.add_subplot(3, 2, 3)
pdf_line1 = ax3.plot(xpdf_pre_reg_eneb_rcp26_hist, pdf_pre_reg_eneb_rcp26_hist, color='blue', label='RegCM4.7 RCP2.6', linestyle='-', linewidth=1.5)
pdf_line2 = ax3.plot(xpdf_pre_had_eneb_rcp26_hist, pdf_pre_reg_eneb_rcp26_hist, color='blue', label='HadGEM2-ES RCP2.6', linestyle='--', linewidth=1.5)
pdf_line3 = ax3.plot(xpdf_pre_reg_eneb_rcp85_hist, pdf_pre_reg_eneb_rcp85_hist, color='red', label='RegCM4.7 RCP8.5', linestyle='-', linewidth=1.5)
pdf_line4 = ax3.plot(xpdf_pre_had_eneb_rcp85_hist, pdf_pre_reg_eneb_rcp85_hist, color='red', label='HadGEM2-ES RCP8.5', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_reg_eneb_rcp26_hist), color='blue', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_had_eneb_rcp26_hist), color='blue', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_reg_eneb_rcp85_hist), color='red', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_had_eneb_rcp85_hist), color='red', linestyle='--', linewidth=1.5)
plt.axvline(0, color='black', linestyle='-', linewidth=1.5)
plt.axvspan(-100, -10, alpha=0.3, color='red')
plt.axvspan(10, 300, alpha=0.3, color='blue')
plt.title(u'B)', loc='left', fontweight='bold', fontsize=8)
ax3.set_ylabel(u'PDF', fontsize=8)
ax3.set_xlim(-100, 300)
ax3.set_ylim(0, 0.01)
ax3.set_xticks(np.arange(-100, 350, 50))
ax3.set_yticks(np.arange(0, 0.02, 0.005))
ax3.xaxis.set_tick_params(labelsize=8)
ax3.yaxis.set_tick_params(labelsize=8)
plt.setp(ax3.get_xticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')
	 
ax4 = fig.add_subplot(3, 2, 4)
pdf_line1 = ax4.plot(xpdf_tas_reg_eneb_rcp26_hist, pdf_tas_reg_eneb_rcp26_hist, color='blue', label='RegCM4.7 RCP2.6', linestyle='-', linewidth=1.5)
pdf_line2 = ax4.plot(xpdf_tas_had_eneb_rcp26_hist, pdf_tas_reg_eneb_rcp26_hist, color='blue', label='HadGEM2-ES RCP2.6', linestyle='--', linewidth=1.5)
pdf_line3 = ax4.plot(xpdf_tas_reg_eneb_rcp85_hist, pdf_tas_reg_eneb_rcp85_hist, color='red', label='RegCM4.7 RCP8.5', linestyle='-', linewidth=1.5)
pdf_line4 = ax4.plot(xpdf_tas_had_eneb_rcp85_hist, pdf_tas_reg_eneb_rcp85_hist, color='red', label='HadGEM2-ES RCP8.5', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_reg_eneb_rcp26_hist), color='blue', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_had_eneb_rcp26_hist), color='blue', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_reg_eneb_rcp85_hist), color='red', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_had_eneb_rcp85_hist), color='red', linestyle='--', linewidth=1.5)
plt.axvline(0, color='black', linestyle='-', linewidth=1.5)
plt.axvspan(-2, -1, alpha=0.3, color='blue')
plt.axvspan(1, 10, alpha=0.3, color='red')
plt.title(u'E)', loc='left', fontweight='bold', fontsize=8)
ax4.set_ylabel(u'PDF', fontsize=8)
ax4.set_xlim(-2, 10)
ax4.set_ylim(0, 0.8)
ax4.set_xticks(np.arange(-2, 12, 2))
ax4.set_yticks(np.arange(0, 1, 0.2))
ax4.xaxis.set_tick_params(labelsize=8)
ax4.yaxis.set_tick_params(labelsize=8)
plt.setp(ax4.get_xticklabels(), visible=False)	     
plt.grid(True, which='major', linestyle='--')
		 
ax5 = fig.add_subplot(3, 2, 5)
pdf_line1 = ax5.plot(xpdf_pre_reg_matopiba_rcp26_hist, pdf_pre_reg_matopiba_rcp26_hist, color='blue', label='RegCM4.7 RCP2.6', linestyle='-', linewidth=1.5)
pdf_line2 = ax5.plot(xpdf_pre_had_matopiba_rcp26_hist, pdf_pre_reg_matopiba_rcp26_hist, color='blue', label='HadGEM2-ES RCP2.6', linestyle='--', linewidth=1.5)
pdf_line3 = ax5.plot(xpdf_pre_reg_matopiba_rcp85_hist, pdf_pre_reg_matopiba_rcp85_hist, color='red', label='RegCM4.7 RCP8.5', linestyle='-', linewidth=1.5)
pdf_line4 = ax5.plot(xpdf_pre_had_matopiba_rcp85_hist, pdf_pre_reg_matopiba_rcp85_hist, color='red', label='HadGEM2-ES RCP8.5', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_reg_matopiba_rcp26_hist), color='blue', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_had_matopiba_rcp26_hist), color='blue', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_reg_matopiba_rcp85_hist), color='red', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_had_matopiba_rcp85_hist), color='red', linestyle='--', linewidth=1.5)
plt.axvline(0, color='black', linestyle='-', linewidth=1.5)
plt.axvspan(-100, -10, alpha=0.3, color='red')
plt.axvspan(10, 300, alpha=0.3, color='blue')
plt.title(u'C)', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'Precipitation change (%)', fontsize=8)
ax5.set_xlim(-100, 300)
ax5.set_ylim(0, 0.01)
ax5.set_xticks(np.arange(-100, 350, 50))
ax5.set_yticks(np.arange(0, 0.02, 0.005))
ax5.xaxis.set_tick_params(labelsize=8)
ax5.yaxis.set_tick_params(labelsize=8)
plt.grid(True, which='major', linestyle='--')

ax6 = fig.add_subplot(3, 2, 6)
pdf_line1 = ax6.plot(xpdf_tas_reg_matopiba_rcp26_hist, pdf_tas_reg_matopiba_rcp26_hist, color='blue', label='RegCM4.7 RCP2.6', linestyle='-', linewidth=1.5)
pdf_line2 = ax6.plot(xpdf_tas_had_matopiba_rcp26_hist, pdf_tas_reg_matopiba_rcp26_hist, color='blue', label='HadGEM2-ES RCP2.6', linestyle='--', linewidth=1.5)
pdf_line3 = ax6.plot(xpdf_tas_reg_matopiba_rcp85_hist, pdf_tas_reg_matopiba_rcp85_hist, color='red', label='RegCM4.7 RCP8.5', linestyle='-', linewidth=1.5)
pdf_line4 = ax6.plot(xpdf_tas_had_matopiba_rcp85_hist, pdf_tas_reg_matopiba_rcp85_hist, color='red', label='HadGEM2-ES RCP8.5', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_reg_matopiba_rcp26_hist), color='blue', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_had_matopiba_rcp26_hist), color='blue', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_reg_matopiba_rcp85_hist), color='red', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_had_matopiba_rcp85_hist), color='red', linestyle='--', linewidth=1.5)
plt.axvline(0, color='black', linestyle='-', linewidth=1.5)
plt.axvspan(-2, -1, alpha=0.3, color='blue')
plt.axvspan(1, 10, alpha=0.3, color='red')
plt.title(u'F)', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'Temperature change (°C)', fontsize=8)
ax6.set_xlim(-2, 10)
ax6.set_ylim(0, 0.8)
ax6.set_xticks(np.arange(-2, 12, 2))
ax6.set_yticks(np.arange(0, 1, 0.2))
ax6.xaxis.set_tick_params(labelsize=8)
ax6.yaxis.set_tick_params(labelsize=8)
plt.grid(True, which='major', linestyle='--')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_pdf_reg_had_rcp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')

plt.show()
exit()






