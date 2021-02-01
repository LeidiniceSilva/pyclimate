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
	

def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/rcm'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	obs = []
	conf_int = []		
	for mon in range(1, 12 + 1):
		obs.append(np.nanmean(value[mon::12], axis=0))
		conf_int.append(np.percentile(value[mon::12], [2.5, 97.5]))
	print(conf_int)

	return value, obs
	
	
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
	
	return value, rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/hist'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	gcm = []
	for mon in range(1, 12 + 1):
		gcm.append(np.nanmean(value[mon::12], axis=0))
	
	return value, gcm

	              
# Import regcm exps model end obs database climatology
mon_pre_cru_samz, pre_cru_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
mon_pre_reg_samz, pre_reg_samz = import_rcm('pr', 'samz', 'hist', '1986-2005')
mon_pre_had_samz, pre_had_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')
mon_pre_cru_eneb, pre_cru_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
mon_pre_reg_eneb, pre_reg_eneb = import_rcm('pr', 'eneb', 'hist', '1986-2005')
mon_pre_had_eneb, pre_had_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')
mon_pre_cru_matopiba, pre_cru_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')
mon_pre_reg_matopiba, pre_reg_matopiba = import_rcm('pr', 'matopiba', 'hist', '1986-2005')
mon_pre_had_matopiba, pre_had_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

mon_tas_cru_samz, tas_cru_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
mon_tas_reg_samz, tas_reg_samz = import_rcm('tas', 'samz', 'hist', '1986-2005')
mon_tas_had_samz, tas_had_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')
mon_tas_cru_eneb, tas_cru_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
mon_tas_reg_eneb, tas_reg_eneb = import_rcm('tas', 'eneb', 'hist', '1986-2005')
mon_tas_had_eneb, tas_had_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')
mon_tas_cru_matopiba, tas_cru_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')
mon_tas_reg_matopiba, tas_reg_matopiba = import_rcm('tas', 'matopiba', 'hist', '1986-2005')
mon_tas_had_matopiba, tas_had_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

# Calculate CFD function
xcdf_pre_cru_samz, cdf_pre_cru_samz = cdf_function(mon_pre_cru_samz)
xcdf_pre_reg_samz, cdf_pre_reg_samz = cdf_function(mon_pre_reg_samz)
xcdf_pre_had_samz, cdf_pre_had_samz = cdf_function(mon_pre_had_samz)
xcdf_pre_cru_eneb, cdf_pre_cru_eneb = cdf_function(mon_pre_cru_eneb)
xcdf_pre_reg_eneb, cdf_pre_reg_eneb = cdf_function(mon_pre_reg_eneb)
xcdf_pre_had_eneb, cdf_pre_had_eneb = cdf_function(mon_pre_had_eneb)
xcdf_pre_cru_matopiba, cdf_pre_cru_matopiba = cdf_function(mon_pre_cru_matopiba)
xcdf_pre_reg_matopiba, cdf_pre_reg_matopiba = cdf_function(mon_pre_reg_matopiba)
xcdf_pre_had_matopiba, cdf_pre_had_matopiba = cdf_function(mon_pre_had_matopiba)

xcdf_tas_cru_samz, cdf_tas_cru_samz = cdf_function(mon_tas_cru_samz)
xcdf_tas_reg_samz, cdf_tas_reg_samz = cdf_function(np.nanmean(mon_tas_reg_samz, axis=0))
xcdf_tas_had_samz, cdf_tas_had_samz = cdf_function(mon_tas_had_samz)
xcdf_tas_cru_eneb, cdf_tas_cru_eneb = cdf_function(mon_tas_cru_eneb)
xcdf_tas_reg_eneb, cdf_tas_reg_eneb = cdf_function(np.nanmean(mon_tas_reg_eneb, axis=0))
xcdf_tas_had_eneb, cdf_tas_had_eneb = cdf_function(mon_tas_had_eneb)
xcdf_tas_cru_matopiba, cdf_tas_cru_matopiba = cdf_function(mon_tas_cru_matopiba)
xcdf_tas_reg_matopiba, cdf_tas_reg_matopiba = cdf_function(np.nanmean(mon_tas_reg_matopiba, axis=0))
xcdf_tas_had_matopiba, cdf_tas_had_matopiba = cdf_function(mon_tas_had_matopiba)

# Calculate PDF function
xpdf_pre_cru_samz, pdf_pre_cru_samz = pdf_function(mon_pre_cru_samz)
xpdf_pre_reg_samz, pdf_pre_reg_samz = pdf_function(mon_pre_reg_samz)
xpdf_pre_had_samz, pdf_pre_had_samz = pdf_function(mon_pre_had_samz)
xpdf_pre_cru_eneb, pdf_pre_cru_eneb = pdf_function(mon_pre_cru_eneb)
xpdf_pre_reg_eneb, pdf_pre_reg_eneb = pdf_function(mon_pre_reg_eneb)
xpdf_pre_had_eneb, pdf_pre_had_eneb = pdf_function(mon_pre_had_eneb)
xpdf_pre_cru_matopiba, pdf_pre_cru_matopiba = pdf_function(mon_pre_cru_matopiba)
xpdf_pre_reg_matopiba, pdf_pre_reg_matopiba = pdf_function(mon_pre_reg_matopiba)
xpdf_pre_had_matopiba, pdf_pre_had_matopiba = pdf_function(mon_pre_had_matopiba)

xpdf_tas_cru_samz, pdf_tas_cru_samz = pdf_function(mon_tas_cru_samz)
xpdf_tas_reg_samz, pdf_tas_reg_samz = pdf_function(np.nanmean(mon_tas_reg_samz, axis=0))
xpdf_tas_had_samz, pdf_tas_had_samz = pdf_function(mon_tas_had_samz)
xpdf_tas_cru_eneb, pdf_tas_cru_eneb = pdf_function(mon_tas_cru_eneb)
xpdf_tas_reg_eneb, pdf_tas_reg_eneb = pdf_function(np.nanmean(mon_tas_reg_eneb, axis=0))
xpdf_tas_had_eneb, pdf_tas_had_eneb = pdf_function(mon_tas_had_eneb)
xpdf_tas_cru_matopiba, pdf_tas_cru_matopiba = pdf_function(mon_tas_cru_matopiba)
xpdf_tas_reg_matopiba, pdf_tas_reg_matopiba = pdf_function(np.nanmean(mon_tas_reg_matopiba, axis=0))
xpdf_tas_had_matopiba, pdf_tas_had_matopiba = pdf_function(mon_tas_had_matopiba)

# Confidence interval 95%
pre_samz_ci1=[8.66463728, 9.4697933, 6.85738347, 3.80205981, 1.48250971, 0.92583204, 1.22552552, 2.63836498, 4.51578057, 5.77292271, 7.07533377, 8.21227317]
pre_samz_ci2=[10.75263259, 11.66272891, 8.50332589, 5.47224972, 2.52379874, 1.68387562, 2.02090873, 3.62090154, 5.93973876, 7.98845394, 10.45760779, 11.00553899]
pre_eneb_ci1=[1.00524572, 1.1608404, 1.77165509, 1.30258448, 1.71733001, 1.79458527, 0.79820421, 0.28923712, 0.29726263, 0.5088564 , 0.39716837, 0.84027369]
pre_eneb_ci2=[4.5356897, 7.152737, 7.64559294, 5.35725617, 5.83513105, 4.48715296, 2.64040059, 2.33545691, 1.56852434, 1.90245899, 3.78144468, 5.38589151]
pre_matopiba_ci1=[5.23713094, 5.05914688, 3.40499815, 1.34858458, 0.50445026, 0.22964761, 0.18256121, 0.42671296, 1.65213115, 2.92743726, 4.84032965, 4.16228163]
pre_matopiba_ci2=[10.35689692, 11.49636469, 7.45052783, 4.22976214, 1.72208759, 1.12680329, 0.80865735, 1.84274545, 4.14807294, 6.75314355, 10.0606626, 11.274613]

tas_samz_ci1=[25.69711699, 25.71469378, 25.97765923, 25.39793468, 24.91330719, 24.70398798, 25.73415108, 26.23648839, 26.75862598,26.37253828, 26.01073117, 25.61359653]
tas_samz_ci2=[26.92223835, 26.85330153, 27.19564362, 26.84376068 , 26.48522525, 26.68297215, 27.20409451, 27.39006267, 28.20131063, 27.48015785, 27.10853429, 27.11589193]
tas_eneb_ci1=[25.99318895, 25.83138466, 25.38425841, 24.71016569, 23.53299055, 23.11000266, 23.37063761, 24.47628813, 25.65004044, 26.16076336, 26.20762401, 26.22430735]
tas_eneb_ci2=[27.63031745, 27.23039775, 27.07775173, 25.82194371, 24.92942662, 24.34608192, 24.69219999, 25.40287352, 26.57566276, 27.12908754, 27.51666064, 27.58264961]
tas_matopiba_ci1=[25.55731535, 25.70387478, 25.79148459, 25.5457233 , 24.99058108, 24.97589211, 25.75307593, 26.66938272, 26.64489403, 26.10251489, 25.4382093 , 25.45791321]
tas_matopiba_ci2=[27.15601463, 26.68446865, 27.67121773, 27.08947453, 26.64261518, 26.57647147, 27.74292717, 28.66518302, 28.50709934, 27.67222061, 27.13844676, 26.91656485]

# Plot model end obs data climatology
fig = plt.figure(figsize=(10,5))
time = np.arange(0.5, 12 + 0.5)

ax10 = fig.add_subplot(3, 4, 1)
annual_cycle1 = ax10.plot(time, pre_cru_samz, time, pre_reg_samz, time, pre_had_samz, time, pre_samz_ci1, time, pre_samz_ci2)
plt.text(12, 13, u'A) SAMZ', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.ylim(0, 12)
plt.yticks(np.arange(0, 14, 2))
plt.setp(ax10.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5 = annual_cycle1
plt.setp(l1, linewidth=1.5, color='black', linestyle='--')
plt.setp(l2, linewidth=1.5, color='green', linestyle='--')
plt.setp(l3, linewidth=1.5, color='magenta', linestyle='--')
plt.setp(l4, linewidth=1.5, color='slategray', linestyle='--')
plt.setp(l5, linewidth=1.5, color='slategray', linestyle='--')
plt.fill_between(time, pre_samz_ci1, pre_samz_ci2, facecolor='slategray',
		 alpha=0.4, interpolate=True)		     
legend = ['CRU', 'Reg', 'Had', '95% CI']
plt.legend(annual_cycle1, legend, fontsize=5.5, loc=9, shadow=True, ncol=1)

ax1 = fig.add_subplot(3, 4, 2)
line1, = ax1.plot(xcdf_pre_cru_samz, cdf_pre_cru_samz, color='black', label='CRU', linestyle='--', linewidth=1.5)
line1, = ax1.plot(xcdf_pre_reg_samz, cdf_pre_reg_samz, color='green', label='Reg', linestyle='--', linewidth=1.5)
line1, = ax1.plot(xcdf_pre_had_samz, cdf_pre_had_samz, color='magenta', label='Had', linestyle='--', linewidth=1.5)
plt.xlim(0, 14)
plt.ylim(0, 1)
plt.xticks(np.arange(0, 16, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax1.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')

ax1_1 = ax1.twinx()
line1_1, = ax1_1.plot(xpdf_pre_cru_samz, pdf_pre_cru_samz, color='black', linestyle='--', linewidth=1.5)
line1_1, = ax1_1.plot(xpdf_pre_reg_samz, pdf_pre_reg_samz, color='green', linestyle='--', linewidth=1.5)
line1_1, = ax1_1.plot(xpdf_pre_had_samz, pdf_pre_had_samz, color='magenta', linestyle='--', linewidth=1.5)
plt.xlim(0, 14)
plt.ylim(0, 0.5)
plt.xticks(np.arange(0, 16, 2))
plt.yticks(np.arange(0, 0.6, 0.1))
plt.setp(ax1_1.get_xticklabels(), visible=False)

ax20 = fig.add_subplot(3, 4, 3)
annual_cycle2 = ax20.plot(time, tas_cru_samz, time, np.nanmean(tas_reg_samz, axis=1), time, tas_had_samz, time, tas_samz_ci1, time, tas_samz_ci2)
plt.text(12, 33, u'D) SAMZ', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.ylim(22, 32)
plt.yticks(np.arange(22, 34, 2))
plt.setp(ax20.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5 = annual_cycle2
plt.setp(l1, linewidth=1.5, color='black', linestyle='--')
plt.setp(l2, linewidth=1.5, color='green', linestyle='--')
plt.setp(l3, linewidth=1.5, color='magenta', linestyle='--')
plt.setp(l4, linewidth=1.5, color='slategray', linestyle='--')
plt.setp(l5, linewidth=1.5, color='slategray', linestyle='--')
plt.fill_between(time, tas_samz_ci1, tas_samz_ci2, facecolor='slategray',
		 alpha=0.4, interpolate=True)	

ax2 = fig.add_subplot(3, 4, 4)
line2, = ax2.plot(xcdf_tas_cru_samz, cdf_tas_cru_samz, color='black', linestyle='--', linewidth=1.5)
line2, = ax2.plot(xcdf_tas_reg_samz, cdf_tas_reg_samz, color='green', linestyle='--', linewidth=1.5)
line2, = ax2.plot(xcdf_tas_had_samz, cdf_tas_had_samz, color='magenta', linestyle='--', linewidth=1.5)
plt.xlim(20, 30)
plt.ylim(0, 1)
plt.xticks(np.arange(20, 32, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax2.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')

ax2_1 = ax2.twinx()
line2_1, = ax2_1.plot(xpdf_tas_cru_samz, pdf_tas_cru_samz, color='black', linestyle='--', linewidth=1.5)
line2_1, = ax2_1.plot(xpdf_tas_reg_samz, pdf_tas_reg_samz, color='green', linestyle='--', linewidth=1.5)
line2_1, = ax2_1.plot(xpdf_tas_had_samz, pdf_tas_had_samz, color='magenta', linestyle='--', linewidth=1.5)
plt.xlim(20, 30)
plt.xticks(np.arange(20, 32, 2))
plt.yticks(np.arange(0, 1.0, 0.2))
plt.setp(ax2_1.get_xticklabels(), visible=False)
		 
ax30 = fig.add_subplot(3, 4, 5)
annual_cycle3 = ax30.plot(time, pre_cru_eneb, time, pre_reg_eneb, time, pre_had_eneb, time, pre_eneb_ci1, time, pre_eneb_ci2)
plt.text(12, 13, u'B) ENEB', fontweight='bold')
plt.ylabel(u'Precipitation (mm d⁻¹)', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.ylim(0, 12)
plt.yticks(np.arange(0, 14, 2))
plt.setp(ax30.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5 = annual_cycle3
plt.setp(l1, linewidth=1.5, color='black', linestyle='--')
plt.setp(l2, linewidth=1.5, color='green', linestyle='--')
plt.setp(l3, linewidth=1.5, color='magenta', linestyle='--')
plt.setp(l4, linewidth=1.5, color='slategray', linestyle='--')
plt.setp(l5, linewidth=1.5, color='slategray', linestyle='--')
plt.fill_between(time, pre_eneb_ci1, pre_eneb_ci2, facecolor='slategray',
		 alpha=0.4, interpolate=True)

ax3 = fig.add_subplot(3, 4, 6)
line3, = ax3.plot(xcdf_pre_cru_eneb, cdf_pre_cru_eneb, color='black', linestyle='--', linewidth=1.5)
line3, = ax3.plot(xcdf_pre_reg_eneb, cdf_pre_reg_eneb, color='green', linestyle='--', linewidth=1.5)
line3, = ax3.plot(xcdf_pre_had_eneb, cdf_pre_had_eneb, color='magenta', linestyle='--', linewidth=1.5)
plt.xlim(0, 14)
plt.ylim(0, 1)
plt.xticks(np.arange(0, 16, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax3.get_xticklabels(), visible=False)
plt.ylabel(u'CDF', fontweight='bold')
plt.grid(True, which='major', linestyle='--')

ax3_1 = ax3.twinx()
line3_1, = ax3_1.plot(xpdf_pre_cru_eneb, pdf_pre_cru_eneb, color='black', linestyle='--', linewidth=1.5)
line3_1, = ax3_1.plot(xpdf_pre_reg_eneb, pdf_pre_reg_eneb, color='green', linestyle='--', linewidth=1.5)
line3_1, = ax3_1.plot(xpdf_pre_had_eneb, pdf_pre_had_eneb, color='magenta', linestyle='--', linewidth=1.5)
plt.xlim(0, 14)
plt.ylim(0, 0.5)
plt.xticks(np.arange(0, 16, 2))
plt.yticks(np.arange(0, 0.6, 0.1))
plt.setp(ax3_1.get_xticklabels(), visible=False)
plt.ylabel(u'PDF', fontweight='bold')
		 
ax40 = fig.add_subplot(3, 4, 7)
annual_cycle4 = ax40.plot(time, tas_cru_eneb, time, np.nanmean(tas_reg_eneb, axis=1), time, tas_had_eneb, time, tas_eneb_ci1, time, tas_eneb_ci2)
plt.text(12, 33, u'E) ENEB', fontweight='bold')
plt.ylabel(u'Temperature (°C)', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.ylim(22, 32)
plt.yticks(np.arange(22, 34, 2))
plt.setp(ax40.get_xticklabels(), visible=False)
plt.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5 = annual_cycle4
plt.setp(l1, linewidth=1.5, color='black', linestyle='--')
plt.setp(l2, linewidth=1.5, color='green', linestyle='--')
plt.setp(l3, linewidth=1.5, color='magenta', linestyle='--')
plt.setp(l4, linewidth=1.5, color='slategray', linestyle='--')
plt.setp(l5, linewidth=1.5, color='slategray', linestyle='--')
plt.fill_between(time, tas_eneb_ci1, tas_eneb_ci2, facecolor='slategray',
		 alpha=0.4, interpolate=True)

ax4 = fig.add_subplot(3, 4, 8)
line4, = ax4.plot(xcdf_tas_cru_eneb, cdf_tas_cru_eneb, color='black', linestyle='--', linewidth=1.5)
line4, = ax4.plot(xcdf_tas_reg_eneb, cdf_tas_reg_eneb, color='green', linestyle='--', linewidth=1.5)
line4, = ax4.plot(xcdf_tas_had_eneb, cdf_tas_had_eneb, color='magenta', linestyle='--', linewidth=1.5)
plt.xlim(20, 30)
plt.xticks(np.arange(20, 32, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.setp(ax4.get_xticklabels(), visible=False)
plt.ylabel(u'CDF', fontweight='bold')
plt.grid(True, which='major', linestyle='--')

ax4_1 = ax4.twinx()
line4_1, = ax4_1.plot(xpdf_tas_cru_eneb, pdf_tas_cru_eneb, color='black', linestyle='--', linewidth=1.5)
line4_1, = ax4_1.plot(xpdf_tas_reg_eneb, pdf_tas_reg_eneb, color='green', linestyle='--', linewidth=1.5)
line4_1, = ax4_1.plot(xpdf_tas_had_eneb, pdf_tas_had_eneb, color='magenta', linestyle='--', linewidth=1.5)
plt.xlim(20, 30)
plt.xticks(np.arange(20, 32, 2))
plt.yticks(np.arange(0, 1.0, 0.2))
plt.setp(ax4_1.get_xticklabels(), visible=False)
plt.ylabel(u'PDF', fontweight='bold')
		 
ax50 = fig.add_subplot(3, 4, 9)
annual_cycle5 = ax50.plot(time, pre_cru_matopiba, time, pre_reg_matopiba, time, pre_had_matopiba, time, pre_matopiba_ci1, time, pre_matopiba_ci2)
plt.text(12, 13, u'C) MATOPIBA', fontweight='bold')
plt.xlabel(u'Months', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.ylim(0, 12)
plt.yticks(np.arange(0, 14, 2))
plt.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5 = annual_cycle5
plt.setp(l1, linewidth=1.5, color='black', linestyle='--')
plt.setp(l2, linewidth=1.5, color='green', linestyle='--')
plt.setp(l3, linewidth=1.5, color='magenta', linestyle='--')
plt.setp(l4, linewidth=1.5, color='slategray', linestyle='--')
plt.setp(l5, linewidth=1.5, color='slategray', linestyle='--')
plt.fill_between(time, pre_matopiba_ci1, pre_matopiba_ci2, facecolor='slategray',
		 alpha=0.4, interpolate=True)

ax5 = fig.add_subplot(3, 4, 10)
line5, = ax5.plot(xcdf_pre_cru_matopiba, cdf_pre_cru_matopiba, color='black', linestyle='--', linewidth=1.5)
line5, = ax5.plot(xcdf_pre_reg_matopiba, cdf_pre_reg_matopiba, color='green', linestyle='--', linewidth=1.5)
line5, = ax5.plot(xcdf_pre_had_matopiba, cdf_pre_had_matopiba, color='magenta', linestyle='--', linewidth=1.5)
plt.xlim(0, 14)
plt.ylim(0, 1)
plt.xticks(np.arange(0, 16, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.xlabel(u'Precipitation (mm d⁻¹)', fontweight='bold')
plt.grid(True, which='major', linestyle='--')

ax5_1 = ax5.twinx()
line5_1, = ax5_1.plot(xpdf_pre_cru_matopiba, pdf_pre_cru_matopiba, color='black', linestyle='--', linewidth=1.5)
line5_1, = ax5_1.plot(xpdf_pre_reg_matopiba, pdf_pre_reg_matopiba, color='green', linestyle='--', linewidth=1.5)
line5_1, = ax5_1.plot(xpdf_pre_had_matopiba, pdf_pre_had_matopiba, color='magenta', linestyle='--', linewidth=1.5)
plt.xlim(0, 14)
plt.ylim(0, 0.5)
plt.xticks(np.arange(0, 16, 2))
plt.yticks(np.arange(0, 0.6, 0.1))
plt.xlabel(u'Precipitation (mm d⁻¹)', fontweight='bold')

ax60 = fig.add_subplot(3, 4, 11)
annual_cycle6 = ax60.plot(time, tas_cru_matopiba, time, np.nanmean(tas_reg_matopiba, axis=1), time, tas_had_matopiba, time, tas_matopiba_ci1, time, tas_matopiba_ci2)
plt.text(12, 33, u'F) MATOPIBA', fontweight='bold')
plt.xlabel(u'Months', fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
plt.ylim(22, 32)
plt.yticks(np.arange(22, 34, 2))
plt.grid(True, which='major', linestyle='--')
l1, l2, l3, l4, l5 = annual_cycle6
plt.setp(l1, linewidth=1.5, color='black', linestyle='--')
plt.setp(l2, linewidth=1.5, color='green', linestyle='--')
plt.setp(l3, linewidth=1.5, color='magenta', linestyle='--')
plt.setp(l4, linewidth=1.5, color='slategray', linestyle='--')
plt.setp(l5, linewidth=1.5, color='slategray', linestyle='--')
plt.fill_between(time, tas_matopiba_ci1, tas_matopiba_ci2, facecolor='slategray',
		 alpha=0.4, interpolate=True)

ax6 = fig.add_subplot(3, 4, 12)
line6, = ax6.plot(xcdf_tas_cru_matopiba, cdf_tas_cru_matopiba, color='black', linestyle='--', linewidth=1.5)
line6, = ax6.plot(xcdf_tas_reg_matopiba, cdf_tas_reg_matopiba, color='green', linestyle='--', linewidth=1.5)
line6, = ax6.plot(xcdf_tas_had_matopiba, cdf_tas_had_matopiba, color='magenta', linestyle='--', linewidth=1.5)
plt.xlim(20, 30)
plt.xticks(np.arange(20, 32, 2))
plt.yticks(np.arange(0, 1.2, 0.2))
plt.xlabel(u'Temperature (°C)', fontweight='bold')
plt.grid(True, which='major', linestyle='--')

ax6_1 = ax6.twinx()
line6_1, = ax6_1.plot(xpdf_tas_cru_matopiba, pdf_tas_cru_matopiba, color='black', linestyle='--', linewidth=1.5)
line6_1, = ax6_1.plot(xpdf_tas_reg_matopiba, pdf_tas_reg_matopiba, color='green', linestyle='--', linewidth=1.5)
line6_1, = ax6_1.plot(xpdf_tas_had_matopiba, pdf_tas_had_matopiba, color='magenta', linestyle='--', linewidth=1.5)
plt.xlim(20, 30)
plt.xticks(np.arange(20, 32, 2))
plt.yticks(np.arange(0, 1.0, 0.2))
plt.xlabel(u'Temperature (°C)', fontweight='bold')

#~ fig.tight_layout()
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95, wspace=0.80, hspace=0.40)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_annual_cycle_cdf_pdf_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()






