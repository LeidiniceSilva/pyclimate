# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script compute CDF functions from Reg and Had models"

import os
import scipy.stats
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import norm


def cdf_function(data):
	
	x = np.sort(data)
	cdf = 1. * np.arange(len(data))/(len(data) - 1)
	
	return x, cdf


def pdf_function(data):

	x = np.linspace(np.min(data), np.max(data))
	y = np.nanmean(x)
	z = np.nanstd(x)
	pdf = norm.pdf(x,y,z)
	
	return x, pdf


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
	
	path = '/home/nice/Documents/dataset/gcm/hist'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return gcm
	
	
# Import regcm exps model end obs database climatology
pre_cru_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
pre_reg_samz = import_rcm('pr', 'samz', 'hist', '1986-2005')
pre_had_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')

pre_cru_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
pre_reg_eneb = import_rcm('pr', 'eneb', 'hist', '1986-2005')
pre_had_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')

pre_cru_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')
pre_reg_matopiba = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
pre_had_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

tas_cru_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
tas_reg_samz = import_rcm('tas', 'samz', 'hist', '1986-2005')
tas_had_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')

tas_cru_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
tas_reg_eneb = import_rcm('tas', 'eneb', 'hist', '1986-2005')
tas_had_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')

tas_cru_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')
tas_reg_matopiba = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
tas_had_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

# Calculate CFD function
xcdf_pre_cru_samz, cdf_pre_cru_samz = cdf_function(pre_cru_samz)
xcdf_pre_reg_samz, cdf_pre_reg_samz = cdf_function(pre_reg_samz)
xcdf_pre_had_samz, cdf_pre_had_samz = cdf_function(pre_had_samz)
xcdf_pre_cru_eneb, cdf_pre_cru_eneb = cdf_function(pre_cru_eneb)
xcdf_pre_reg_eneb, cdf_pre_reg_eneb = cdf_function(pre_reg_eneb)
xcdf_pre_had_eneb, cdf_pre_had_eneb = cdf_function(pre_had_eneb)
xcdf_pre_cru_matopiba, cdf_pre_cru_matopiba = cdf_function(pre_cru_matopiba)
xcdf_pre_reg_matopiba, cdf_pre_reg_matopiba = cdf_function(pre_reg_matopiba)
xcdf_pre_had_matopiba, cdf_pre_had_matopiba = cdf_function(pre_had_matopiba)

xcdf_tas_cru_samz, cdf_tas_cru_samz = cdf_function(tas_cru_samz)
xcdf_tas_reg_samz, cdf_tas_reg_samz = cdf_function(np.nanmean(tas_reg_samz, axis=0))
xcdf_tas_had_samz, cdf_tas_had_samz = cdf_function(tas_had_samz)
xcdf_tas_cru_eneb, cdf_tas_cru_eneb = cdf_function(tas_cru_eneb)
xcdf_tas_reg_eneb, cdf_tas_reg_eneb = cdf_function(np.nanmean(tas_reg_eneb, axis=0))
xcdf_tas_had_eneb, cdf_tas_had_eneb = cdf_function(tas_had_eneb)
xcdf_tas_cru_matopiba, cdf_tas_cru_matopiba = cdf_function(tas_cru_matopiba)
xcdf_tas_reg_matopiba, cdf_tas_reg_matopiba = cdf_function(np.nanmean(tas_reg_matopiba, axis=0))
xcdf_tas_had_matopiba, cdf_tas_had_matopiba = cdf_function(tas_had_matopiba)

# Calculate PDF function
xpdf_pre_cru_samz, pdf_pre_cru_samz = pdf_function(pre_cru_samz)
xpdf_pre_reg_samz, pdf_pre_reg_samz = pdf_function(pre_reg_samz)
xpdf_pre_had_samz, pdf_pre_had_samz = pdf_function(pre_had_samz)
xpdf_pre_cru_eneb, pdf_pre_cru_eneb = pdf_function(pre_cru_eneb)
xpdf_pre_reg_eneb, pdf_pre_reg_eneb = pdf_function(pre_reg_eneb)
xpdf_pre_had_eneb, pdf_pre_had_eneb = pdf_function(pre_had_eneb)
xpdf_pre_cru_matopiba, pdf_pre_cru_matopiba = pdf_function(pre_cru_matopiba)
xpdf_pre_reg_matopiba, pdf_pre_reg_matopiba = pdf_function(pre_reg_matopiba)
xpdf_pre_had_matopiba, pdf_pre_had_matopiba = pdf_function(pre_had_matopiba)

xpdf_tas_cru_samz, pdf_tas_cru_samz = pdf_function(tas_cru_samz)
xpdf_tas_reg_samz, pdf_tas_reg_samz = pdf_function(np.nanmean(tas_reg_samz, axis=0))
xpdf_tas_had_samz, pdf_tas_had_samz = pdf_function(tas_had_samz)
xpdf_tas_cru_eneb, pdf_tas_cru_eneb = pdf_function(tas_cru_eneb)
xpdf_tas_reg_eneb, pdf_tas_reg_eneb = pdf_function(np.nanmean(tas_reg_eneb, axis=0))
xpdf_tas_had_eneb, pdf_tas_had_eneb = pdf_function(tas_had_eneb)
xpdf_tas_cru_matopiba, pdf_tas_cru_matopiba = pdf_function(tas_cru_matopiba)
xpdf_tas_reg_matopiba, pdf_tas_reg_matopiba = pdf_function(np.nanmean(tas_reg_matopiba, axis=0))
xpdf_tas_had_matopiba, pdf_tas_had_matopiba = pdf_function(tas_had_matopiba)

# Plot model end obs data climatology
fig = plt.figure()

ax1 = fig.add_subplot(3, 4, 1)
plt.plot(xcdf_pre_cru_samz, cdf_pre_cru_samz, color='black', linestyle='--', linewidth=1.5)
plt.plot(xcdf_pre_reg_samz, cdf_pre_reg_samz, color='green', linestyle='--', linewidth=1.5)
plt.plot(xcdf_pre_had_samz, cdf_pre_had_samz, color='magenta', linestyle='--', linewidth=1.5)
plt.title(u'A)', fontweight='bold')
plt.xticks(np.arange(0, 18, 3))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = fig.add_subplot(3, 4, 2)
plt.plot(xpdf_pre_cru_samz, pdf_pre_cru_samz, color='black', linestyle='--', linewidth=1.5)
plt.plot(xpdf_pre_reg_samz, pdf_pre_reg_samz, color='green', linestyle='--', linewidth=1.5)
plt.plot(xpdf_pre_had_samz, pdf_pre_had_samz, color='magenta', linestyle='--', linewidth=1.5)
plt.title(u'D)', fontweight='bold')
plt.xticks(np.arange(0, 18, 3))
plt.setp(ax2.get_xticklabels(), visible=False)

ax3 = fig.add_subplot(3, 4, 3)
plt.plot(xcdf_tas_cru_samz, cdf_tas_cru_samz, color='black', linestyle='--', linewidth=1.5)
plt.plot(xcdf_tas_reg_samz, cdf_tas_reg_samz, color='green', linestyle='--', linewidth=1.5)
plt.plot(xcdf_tas_had_samz, cdf_tas_had_samz, color='magenta', linestyle='--', linewidth=1.5)
plt.title(u'G)', fontweight='bold')
plt.xticks(np.arange(20, 32, 3))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax3.get_xticklabels(), visible=False)

ax4 = fig.add_subplot(3, 4, 4)
plt.plot(xpdf_tas_cru_samz, pdf_tas_cru_samz, color='black', linestyle='--', linewidth=1.5)
plt.plot(xpdf_tas_reg_samz, pdf_tas_reg_samz, color='green', linestyle='--', linewidth=1.5)
plt.plot(xpdf_tas_had_samz, pdf_tas_had_samz, color='magenta', linestyle='--', linewidth=1.5)
plt.title(u'J)', fontweight='bold')
plt.xticks(np.arange(20, 32, 3))
plt.setp(ax4.get_xticklabels(), visible=False)

ax5 = fig.add_subplot(3, 4, 5)
plt.plot(xcdf_pre_cru_eneb, cdf_pre_cru_eneb, color='black', linestyle='--', linewidth=1.5)
plt.plot(xcdf_pre_reg_eneb, cdf_pre_reg_eneb, color='green', linestyle='--', linewidth=1.5)
plt.plot(xcdf_pre_had_eneb, cdf_pre_had_eneb, color='magenta', linestyle='--', linewidth=1.5)
plt.title(u'B)', fontweight='bold')
plt.xticks(np.arange(0, 18, 3))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax5.get_xticklabels(), visible=False)
plt.ylabel(u'CDF', fontweight='bold')

ax6 = fig.add_subplot(3, 4, 6)
plt.plot(xpdf_pre_cru_eneb, pdf_pre_cru_eneb, color='black', linestyle='--', linewidth=1.5)
plt.plot(xpdf_pre_reg_eneb, pdf_pre_reg_eneb, color='green', linestyle='--', linewidth=1.5)
plt.plot(xpdf_pre_had_eneb, pdf_pre_had_eneb, color='magenta', linestyle='--', linewidth=1.5)
plt.title(u'E)', fontweight='bold')
plt.xticks(np.arange(0, 18, 3))
plt.setp(ax6.get_xticklabels(), visible=False)
plt.ylabel(u'PDF', fontweight='bold')

ax7 = fig.add_subplot(3, 4, 7)
plt.plot(xcdf_tas_cru_eneb, cdf_tas_cru_eneb, color='black', linestyle='--', linewidth=1.5)
plt.plot(xcdf_tas_reg_eneb, cdf_tas_reg_eneb, color='green', linestyle='--', linewidth=1.5)
plt.plot(xcdf_tas_had_eneb, cdf_tas_had_eneb, color='magenta', linestyle='--', linewidth=1.5)
plt.title(u'H)', fontweight='bold')
plt.xticks(np.arange(20, 32, 3))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax7.get_xticklabels(), visible=False)
plt.ylabel(u'CDF', fontweight='bold')

ax8 = fig.add_subplot(3, 4, 8)
plt.plot(xpdf_tas_cru_eneb, pdf_tas_cru_eneb, color='black', linestyle='--', linewidth=1.5)
plt.plot(xpdf_tas_reg_eneb, pdf_tas_reg_eneb, color='green', linestyle='--', linewidth=1.5)
plt.plot(xpdf_tas_had_eneb, pdf_tas_had_eneb, color='magenta', linestyle='--', linewidth=1.5)
plt.title(u'K)', fontweight='bold')
plt.xticks(np.arange(20, 32, 3))
plt.setp(ax8.get_xticklabels(), visible=False)
plt.ylabel(u'PDF', fontweight='bold')

ax9 = fig.add_subplot(3, 4, 9)
plt.plot(xcdf_pre_cru_matopiba, cdf_pre_cru_matopiba, color='black', linestyle='--', linewidth=1.5)
plt.plot(xcdf_pre_reg_matopiba, cdf_pre_reg_matopiba, color='green', linestyle='--', linewidth=1.5)
plt.plot(xcdf_pre_had_matopiba, cdf_pre_had_matopiba, color='magenta', linestyle='--', linewidth=1.5)
plt.title(u'C)', fontweight='bold')
plt.xticks(np.arange(0, 18, 3))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.xlabel(u'Precipitation (mm d⁻¹)', fontweight='bold')

ax10 = fig.add_subplot(3, 4, 10)
plt.plot(xpdf_pre_cru_matopiba, pdf_pre_cru_matopiba, color='black', linestyle='--', linewidth=1.5)
plt.plot(xpdf_pre_reg_matopiba, pdf_pre_reg_matopiba, color='green', linestyle='--', linewidth=1.5)
plt.plot(xpdf_pre_had_matopiba, pdf_pre_had_matopiba, color='magenta', linestyle='--', linewidth=1.5)
plt.title(u'F)', fontweight='bold')
plt.xticks(np.arange(0, 18, 3))
plt.xlabel(u'Precipitation (mm d⁻¹)', fontweight='bold')

ax11 = fig.add_subplot(3, 4, 11)
plt.plot(xcdf_tas_cru_matopiba, cdf_tas_cru_matopiba, color='black', linestyle='--', linewidth=1.5)
plt.plot(xcdf_tas_reg_matopiba, cdf_tas_reg_matopiba, color='green', linestyle='--', linewidth=1.5)
plt.plot(xcdf_tas_had_matopiba, cdf_tas_had_matopiba, color='magenta', linestyle='--', linewidth=1.5)
plt.title(u'I)', fontweight='bold')
plt.xticks(np.arange(20, 32, 3))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.xlabel(u'Temperature (°C)', fontweight='bold')

ax12 = fig.add_subplot(3, 4, 12)
plt.plot(xpdf_tas_cru_matopiba, pdf_tas_cru_matopiba, color='black', linestyle='--', linewidth=1.5)
plt.plot(xpdf_tas_reg_matopiba, pdf_tas_reg_matopiba, color='green', linestyle='--', linewidth=1.5)
plt.plot(xpdf_tas_had_matopiba, pdf_tas_had_matopiba, color='magenta', linestyle='--', linewidth=1.5)
plt.title(u'L)', fontweight='bold')
plt.xticks(np.arange(20, 32, 3))
plt.xlabel(u'Temperature (°C)', fontweight='bold')

fig.tight_layout()
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.50, hspace=0.50)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_cdf_pdf_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()


plt.plot(sort_pre_reg_samz, cdf_pre_reg_samz, color='black', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_had_samz, cdf_pre_had_samz, color='blue', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_cru_samz, cdf_pre_cru_samz, color='red', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_udel_samz, cdf_pre_udel_samz, color='green', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_chirps_samz, cdf_pre_chirps_samz, color='yellow', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_era5_samz, cdf_pre_era5_samz, color='magenta', linestyle='--', linewidth=1.5)
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'A) SAMZ', fontweight='bold')
plt.xticks(np.arange(0, 18, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = fig.add_subplot(3, 2, 2)
plt.plot(sort_tas_reg_samz, cdf_tas_reg_samz, color='black', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_had_samz, cdf_tas_had_samz, color='blue', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_cru_samz, cdf_tas_cru_samz, color='red', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_udel_samz, cdf_tas_udel_samz, color='green', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_era5_samz, cdf_tas_era5_samz, color='magenta', linestyle='--', linewidth=1.5)
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'B) SAMZ', fontweight='bold')
plt.xticks(np.arange(20, 32, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax2.get_xticklabels(), visible=False)

ax3 = fig.add_subplot(3, 2, 3)
plt.plot(sort_pre_reg_eneb, cdf_pre_reg_eneb, color='black', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_had_eneb, cdf_pre_had_eneb, color='blue', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_cru_eneb, cdf_pre_cru_eneb, color='red', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_udel_eneb, cdf_pre_udel_eneb, color='green', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_chirps_eneb, cdf_pre_chirps_eneb, color='yellow', linestyle='--', linewidth=1.5)
plt.plot(sort_pre_era5_eneb, cdf_pre_era5_eneb, color='magenta', linestyle='--', linewidth=1.5)
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'C) ENEB', fontweight='bold')
plt.xticks(np.arange(0, 18, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax3.get_xticklabels(), visible=False)
plt.ylabel(u'CDF', fontweight='bold')

ax4 = fig.add_subplot(3, 2, 4)
plt.plot(sort_tas_reg_eneb, cdf_tas_reg_eneb, color='black', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_had_eneb, cdf_tas_had_eneb, color='blue', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_cru_eneb, cdf_tas_cru_eneb, color='red', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_udel_eneb, cdf_tas_udel_eneb, color='green', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_era5_eneb, cdf_tas_era5_eneb, color='magenta', linestyle='--', linewidth=1.5)
ax4.xaxis.grid(True, which='major', linestyle='--')
ax4.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'D) ENEB', fontweight='bold')
plt.xticks(np.arange(20, 32, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax4.get_xticklabels(), visible=False)
plt.ylabel(u'CDF', fontweight='bold')

ax5 = fig.add_subplot(3, 2, 5)
ax51 = plt.plot(sort_pre_reg_matopiba, cdf_pre_reg_matopiba, 
sort_pre_had_matopiba, cdf_pre_had_matopiba, 
sort_pre_cru_matopiba, cdf_pre_cru_matopiba, 
sort_pre_udel_matopiba, cdf_pre_udel_matopiba, 
sort_pre_chirps_matopiba, cdf_pre_chirps_matopiba, 
sort_pre_era5_matopiba, cdf_pre_era5_matopiba)
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'E) MATOPIBA', fontweight='bold')
plt.xticks(np.arange(0, 18, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.xlabel(u'Precipitation (mm d⁻¹)', fontweight='bold')

l1, l2, l3, l4, l5, l6 = ax51
plt.setp(l1, color='black', linestyle='--', linewidth=1.5)
plt.setp(l2, color='blue', linestyle='--', linewidth=1.5)
plt.setp(l3, color='red', linestyle='--', linewidth=1.5)
plt.setp(l4, color='green', linestyle='--', linewidth=1.5)
plt.setp(l5, color='yellow', linestyle='--', linewidth=1.5)
plt.setp(l6, color='magenta', linestyle='--', linewidth=1.5)

legend = ['Reg','Had','CRU','UDEL','CHIRPS','ERA5']
plt.legend(ax51, legend, loc='lower left', bbox_to_anchor=(-0.3, -0.8), shadow=True, ncol=6)

ax6 = fig.add_subplot(3, 2, 6)
plt.plot(sort_tas_reg_matopiba, cdf_tas_reg_matopiba, color='black', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_had_matopiba, cdf_tas_had_matopiba, color='blue', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_cru_matopiba, cdf_tas_cru_matopiba, color='red', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_udel_matopiba, cdf_tas_udel_matopiba, color='green', linestyle='--', linewidth=1.5)
plt.plot(sort_tas_era5_matopiba, cdf_tas_era5_matopiba, color='magenta', linestyle='--', linewidth=1.5)
ax6.xaxis.grid(True, which='major', linestyle='--')
ax6.yaxis.grid(True, which='major', linestyle='--')
plt.title(u'F) MATOPIBA', fontweight='bold')
plt.xticks(np.arange(20, 32, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.xlabel(u'Temperature (°C)', fontweight='bold')

fig.tight_layout()
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.50, hspace=0.50)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_cdf_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()

