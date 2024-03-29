# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot pdf function from regcm47 and hadgem models and obs database"

import os
import netCDF4
import statistics
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# mpl.use('agg')

from pylab import *
from netCDF4 import Dataset
from scipy.stats import norm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from comp_statist_indices import compute_pdf


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/reg_exp1'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
		
	return value
	
	
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp1/hist'
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	return value


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp1/hist'
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

xpdf_tas_reg_samz, pdf_tas_reg_samz = compute_pdf(np.nanmean(mon_tas_reg_samz, axis=1))
xpdf_tas_reg_eneb, pdf_tas_reg_eneb = compute_pdf(np.nanmean(mon_tas_reg_eneb, axis=1))
xpdf_tas_reg_matopiba, pdf_tas_reg_matopiba = compute_pdf(np.nanmean(mon_tas_reg_matopiba, axis=1))

xpdf_tas_had_samz, pdf_tas_had_samz = compute_pdf(mon_tas_had_samz)
xpdf_tas_had_eneb, pdf_tas_had_eneb = compute_pdf(mon_tas_had_eneb)
xpdf_tas_had_matopiba, pdf_tas_had_matopiba = compute_pdf(mon_tas_had_matopiba)

# Plot models and obs database 
fig = plt.figure()

xpdf_pre_cru_samz_mean=np.nanmean(xpdf_pre_cru_samz)
xpdf_pre_reg_samz_mean=np.nanmean(xpdf_pre_reg_samz)
xpdf_pre_had_samz_mean=np.nanmean(xpdf_pre_had_samz)

xpdf_pre_cru_eneb_mean=np.nanmean(xpdf_pre_cru_eneb)
xpdf_pre_reg_eneb_mean=np.nanmean(xpdf_pre_reg_eneb)
xpdf_pre_had_eneb_mean=np.nanmean(xpdf_pre_had_eneb)

xpdf_pre_cru_matopiba_mean=np.nanmean(xpdf_pre_cru_matopiba)
xpdf_pre_reg_matopiba_mean=np.nanmean(xpdf_pre_reg_matopiba)
xpdf_pre_had_matopiba_mean=np.nanmean(xpdf_pre_had_matopiba)

xpdf_tas_cru_samz_mean=np.nanmean(xpdf_tas_cru_samz)
xpdf_tas_reg_samz_mean=np.nanmean(xpdf_tas_reg_samz)
xpdf_tas_had_samz_mean=np.nanmean(xpdf_tas_had_samz)

xpdf_tas_cru_eneb_mean=np.nanmean(xpdf_tas_cru_eneb)
xpdf_tas_reg_eneb_mean=np.nanmean(xpdf_tas_reg_eneb)
xpdf_tas_had_eneb_mean=np.nanmean(xpdf_tas_had_eneb)

xpdf_tas_cru_matopiba_mean=np.nanmean(xpdf_tas_cru_matopiba)
xpdf_tas_reg_matopiba_mean=np.nanmean(xpdf_tas_reg_matopiba)
xpdf_tas_had_matopiba_mean=np.nanmean(xpdf_tas_had_matopiba)

ax1 = fig.add_subplot(3, 2, 1)
pdf_line1 = ax1.plot(xpdf_pre_cru_samz, pdf_pre_cru_samz, color='black', label='CRU', linestyle='-', linewidth=1.5) 
pdf_line2 = ax1.plot(xpdf_pre_reg_samz, pdf_pre_reg_samz, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5)
pdf_line3 = ax1.plot(xpdf_pre_had_samz, pdf_pre_had_samz, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_cru_samz), color='black', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_reg_samz), color='gray', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_had_samz), color='dimgray', linestyle='--', linewidth=1.5)
plt.title(u'A) SAMZ', loc='left', fontweight='bold', fontsize=8)
ax1.set_xlim(0, 14)
ax1.set_ylim(0, 0.4)
ax1.set_xticks(np.arange(0, 16, 2))
ax1.set_yticks(np.arange(0, 0.5, 0.1))
ax1.xaxis.set_tick_params(labelsize=8)
ax1.yaxis.set_tick_params(labelsize=8)
plt.setp(ax1.get_xticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)
plt.legend(fontsize=6, loc=1, shadow=True, ncol=1)
#~ plt.text(0.5, 0.30, 'CRU = %.2f'%(xpdf_pre_cru_samz_mean), fontweight='bold', fontsize=6, color='black')
#~ plt.text(0.5, 0.25, 'RegCM4.7 = %.2f'%(xpdf_pre_reg_samz_mean), fontweight='bold', fontsize=6, color='gray')
#~ plt.text(0.5, 0.20, 'HadGEM2-ES = %.2f'%(xpdf_pre_had_samz_mean), fontweight='bold', fontsize=6, color='dimgray')

ax2 = fig.add_subplot(3, 2, 2)
pdf_line1 = ax2.plot(xpdf_tas_cru_samz, pdf_tas_cru_samz, color='black', label='CRU', linestyle='-', linewidth=1.5)
pdf_line2 = ax2.plot(xpdf_tas_reg_samz, pdf_tas_reg_samz, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5)
pdf_line3 = ax2.plot(xpdf_tas_had_samz, pdf_tas_had_samz, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_cru_samz), color='black', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_reg_samz), color='gray', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_had_samz), color='dimgray', linestyle='--', linewidth=1.5)
plt.title(u'D) SAMZ', loc='left', fontweight='bold', fontsize=8)
ax2.set_xlim(20, 34)
ax2.set_ylim(0, 0.4)
ax2.set_xticks(np.arange(20, 36, 2))
ax2.set_yticks(np.arange(0, 0.5, 0.1))
ax2.xaxis.set_tick_params(labelsize=8)
ax2.yaxis.set_tick_params(labelsize=8)
plt.setp(ax2.get_xticklabels(), visible=False)	     
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)
#~ plt.text(20.5, 0.7, 'CRU=26.46', fontweight='bold', fontsize=6, color='black')
#~ plt.text(20.5, 0.5, 'RegCM4.7=27.48', fontweight='bold', fontsize=6, color='gray')
#~ plt.text(20.5, 0.3, 'HadGEM2-ES=24.94', fontweight='bold', fontsize=6, color='dimgray')
	 
ax3 = fig.add_subplot(3, 2, 3)
pdf_line1 = ax3.plot(xpdf_pre_cru_eneb, pdf_pre_cru_eneb, color='black', label='CRU', linestyle='-', linewidth=1.5)
pdf_line2 = ax3.plot(xpdf_pre_reg_eneb, pdf_pre_reg_eneb, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5)
pdf_line3 = ax3.plot(xpdf_pre_had_eneb, pdf_pre_had_eneb, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_cru_eneb), color='black', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_reg_eneb), color='gray', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_had_eneb), color='dimgray', linestyle='--', linewidth=1.5)
plt.title(u'B) ENEB', loc='left', fontweight='bold', fontsize=8)
ax3.set_ylabel(u'PDF', fontsize=8)
ax3.set_xlim(0, 14)
ax3.set_ylim(0, 0.4)
ax3.set_xticks(np.arange(0, 16, 2))
ax3.set_yticks(np.arange(0, 0.5, 0.1))
ax3.xaxis.set_tick_params(labelsize=8)
ax3.yaxis.set_tick_params(labelsize=8)
plt.setp(ax3.get_xticklabels(), visible=False)		     
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)
#~ plt.text(8, 0.35, 'CRU=4.27', fontweight='bold', fontsize=6, color='black')
#~ plt.text(8, 0.25, 'RegCM4.7=6.93', fontweight='bold', fontsize=6, color='gray')
#~ plt.text(8, 0.15, 'HadGEM2-ES=6.32', fontweight='bold', fontsize=6, color='dimgray')
 
ax4 = fig.add_subplot(3, 2, 4)
pdf_line1 = ax4.plot(xpdf_tas_cru_eneb, pdf_tas_cru_eneb, color='black', label='CRU', linestyle='-', linewidth=1.5)
pdf_line2 = ax4.plot(xpdf_tas_reg_eneb, pdf_tas_reg_eneb, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5)
pdf_line3 = ax4.plot(xpdf_tas_had_eneb, pdf_tas_had_eneb, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_cru_eneb), color='black', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_reg_eneb), color='gray', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_had_eneb), color='dimgray', linestyle='--', linewidth=1.5)
plt.title(u'E) ENEB', loc='left', fontweight='bold', fontsize=8)
ax4.set_ylabel(u'PDF', fontsize=8)
ax4.set_xlim(20, 34)
ax4.set_ylim(0, 0.4)
ax4.set_xticks(np.arange(20, 36, 2))
ax4.set_yticks(np.arange(0, 0.5, 0.1))
ax4.xaxis.set_tick_params(labelsize=8)
ax4.yaxis.set_tick_params(labelsize=8)
plt.setp(ax4.get_xticklabels(), visible=False)	     
plt.grid(True, which='major', linestyle='--')
plt.tick_params(bottom = False)
#~ plt.text(20.5, 0.7, 'CRU=25.46', fontweight='bold', fontsize=6, color='black')
#~ plt.text(20.5, 0.5, 'RegCM4.7=26.69', fontweight='bold', fontsize=6, color='gray')
#~ plt.text(20.5, 0.3, 'HadGEM2-ES=24.34', fontweight='bold', fontsize=6, color='dimgray')
	 	 
ax5 = fig.add_subplot(3, 2, 5)
pdf_line1 = ax5.plot(xpdf_pre_cru_matopiba, pdf_pre_cru_matopiba, color='black', label='CRU', linestyle='-', linewidth=1.5)
pdf_line2 = ax5.plot(xpdf_pre_reg_matopiba, pdf_pre_reg_matopiba, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5)
pdf_line3 = ax5.plot(xpdf_pre_had_matopiba, pdf_pre_had_matopiba, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_cru_matopiba), color='black', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_reg_matopiba), color='gray', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_pre_had_matopiba), color='dimgray', linestyle='--', linewidth=1.5)
plt.title(u'C) MATOPIBA', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'Precipitation (mm d⁻¹)', fontsize=8)
ax5.set_xlim(0, 14)
ax5.set_ylim(0, 0.4)
ax5.set_xticks(np.arange(0, 16, 2))
ax5.set_yticks(np.arange(0, 0.5, 0.1))
ax5.xaxis.set_tick_params(labelsize=8)
ax5.yaxis.set_tick_params(labelsize=8)
plt.grid(True, which='major', linestyle='--')
#~ plt.text(8, 0.3, 'CRU=7.12', fontweight='bold', fontsize=6, color='black')
#~ plt.text(8, 0.2, 'RegCM4.7=5.39', fontweight='bold', fontsize=6, color='gray')
#~ plt.text(8, 0.1, 'HadGEM2-ES=6.97', fontweight='bold', fontsize=6, color='dimgray')
 
ax6 = fig.add_subplot(3, 2, 6)
pdf_line1 = ax6.plot(xpdf_tas_cru_matopiba, pdf_tas_cru_matopiba, color='black', label='CRU', linestyle='-', linewidth=1.5)
pdf_line2 = ax6.plot(xpdf_tas_reg_matopiba, pdf_tas_reg_matopiba, color='gray', label='RegCM4.7', linestyle='-', linewidth=1.5)
pdf_line3 = ax6.plot(xpdf_tas_had_matopiba, pdf_tas_had_matopiba, color='dimgray', label='HadGEM2-ES', linestyle='--', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_cru_matopiba), color='black', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_reg_matopiba), color='gray', linestyle='-', linewidth=1.5)
plt.axvline(np.nanmean(xpdf_tas_had_matopiba), color='dimgray', linestyle='--', linewidth=1.5)
plt.title(u'F) MATOPIBA', loc='left', fontweight='bold', fontsize=8)
plt.xlabel(u'Temperature (°C)', fontsize=8)
ax6.set_xlim(20, 34)
ax6.set_ylim(0, 0.4)
ax6.set_xticks(np.arange(20, 36, 2))
ax6.set_yticks(np.arange(0, 0.5, 0.1))
ax6.xaxis.set_tick_params(labelsize=8)
ax6.yaxis.set_tick_params(labelsize=8)
plt.grid(True, which='major', linestyle='--')
#~ plt.text(20.5, 0.7, 'CRU=26.90', fontweight='bold', fontsize=6, color='black')
#~ plt.text(20.5, 0.5, 'RegCM4.7=28.67', fontweight='bold', fontsize=6, color='gray')
#~ plt.text(20.5, 0.3, 'HadGEM2-ES=24.87', fontweight='bold', fontsize=6, color='dimgray')
	 	
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_pdf_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()

