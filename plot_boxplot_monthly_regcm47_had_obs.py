# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot boxplot from Reg and Had models end obs database"

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

	
def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return obs
	
               
# Import regcm exps model end obs database climatology
p_reg = import_rcm('pr', 'amz_neb', 'hist', '1986-2005')
p_had = import_gcm('pr', 'amz_neb', 'hist', '1986-2005')
p_cru = import_obs('pre', 'amz_neb', 'cru_ts4.04', '1986-2005')
p_udel = import_obs('pre', 'amz_neb', 'udel_v301', '1986-2005')
p_chirps = import_obs('precip', 'amz_neb', 'chirps-v2.0', '1986-2005')
p_era5 = import_obs('mtpr', 'amz_neb', 'era5', '1986-2005')

t_reg = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')
t_had = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')
t_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', '1986-2005')
t_udel = import_obs('temp', 'amz_neb', 'udel_v301', '1986-2005')
t_era5 = import_obs('t2m', 'amz_neb', 'era5', '1986-2005')

# Plot model end obs data boxplot
fig = plt.figure()
time1 = np.arange(1,7)
time2 = np.arange(1,6)

ax1 = fig.add_subplot(3, 2, 1)
data = [p_reg, p_had, p_cru, p_udel, p_chirps, p_era5]
a1 = ax1.boxplot(data, patch_artist=True, notch=True, bootstrap=10000, vert=1)

# Change outline and fill color
for box in a1['boxes']:
    box.set( color='black', linewidth=1)
    box.set( facecolor = 'blue' )

# Change color and linewidth of the whiskers
for whisker in a1['whiskers']:
    whisker.set(color='black', linewidth=1)

# Change color and linewidth of the caps
for cap in a1['caps']:
    cap.set(color='black', linewidth=1)

# Change color and linewidth of the medians
for median in a1['medians']:
    median.set(color='black', linewidth=1)

# Change the style of fliers and their fill
for flier in a1['fliers']:
    flier.set(marker='o', color='black', alpha=1)
    
plt.title(u'A) SAMZ', fontweight='bold')
plt.xticks(time1, (u'Reg', u'Had', u'CRU', u'UDEL', 'CHIRPS', 'ERA5'))
plt.yticks(np.arange(0, 12, 2), fontsize=7)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.xaxis.grid(True, which='major', linestyle='--')
ax1.yaxis.grid(True, which='major', linestyle='--')

ax2 = fig.add_subplot(3, 2, 2)
data = [t_had, t_had, t_cru, t_udel, t_era5]
a2 = ax2.boxplot(data, patch_artist=True, notch=True, bootstrap=10000, vert=1)

# Change outline and fill color
for box in a2['boxes']:
    box.set( color='black', linewidth=1)
    box.set( facecolor = 'red' )

# Change color and linewidth of the whiskers
for whisker in a2['whiskers']:
    whisker.set(color='black', linewidth=1)

# Change color and linewidth of the caps
for cap in a2['caps']:
    cap.set(color='black', linewidth=1)

# Change color and linewidth of the medians
for median in a2['medians']:
    median.set(color='black', linewidth=1)

# Change the style of fliers and their fill
for flier in a2['fliers']:
    flier.set(marker='o', color='black', alpha=1)

plt.title(u'B) SAMZ', fontweight='bold')
plt.xticks(time2, (u'Reg', u'Had', u'CRU', u'UDEL', 'CHIRPS', 'ERA5'))
plt.yticks(np.arange(20, 30, 2), fontsize=7)
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.xaxis.grid(True, which='major', linestyle='--')
ax2.yaxis.grid(True, which='major', linestyle='--')

ax3 = fig.add_subplot(3, 2, 3)
data = [p_reg, p_had, p_cru, p_udel, p_chirps, p_era5]
a3 = ax3.boxplot(data, patch_artist=True, notch=True, bootstrap=10000, vert=1)

# Change outline and fill color
for box in a3['boxes']:
    box.set( color='black', linewidth=1)
    box.set( facecolor = 'blue' )

# Change color and linewidth of the whiskers
for whisker in a3['whiskers']:
    whisker.set(color='black', linewidth=1)

# Change color and linewidth of the caps
for cap in a3['caps']:
    cap.set(color='black', linewidth=1)

# Change color and linewidth of the medians
for median in a3['medians']:
    median.set(color='black', linewidth=1)

# Change the style of fliers and their fill
for flier in a3['fliers']:
    flier.set(marker='o', color='black', alpha=1)
    
plt.title(u'C) ENEB', fontweight='bold')
plt.ylabel(u'Precipitation (mm d⁻¹)', fontweight='bold')
plt.xticks(time1, (u'Reg', u'Had', u'CRU', u'UDEL', 'CHIRPS', 'ERA5'))
plt.yticks(np.arange(0, 12, 2), fontsize=7)
plt.setp(ax3.get_xticklabels(), visible=False)
ax3.xaxis.grid(True, which='major', linestyle='--')
ax3.yaxis.grid(True, which='major', linestyle='--')

ax4 = fig.add_subplot(3, 2, 4)
data = [t_had, t_had, t_cru, t_udel, t_era5]
a4 = ax4.boxplot(data, patch_artist=True, notch=True, bootstrap=10000, vert=1)

# Change outline and fill color
for box in a4['boxes']:
    box.set( color='black', linewidth=1)
    box.set( facecolor = 'red' )

# Change color and linewidth of the whiskers
for whisker in a4['whiskers']:
    whisker.set(color='black', linewidth=1)

# Change color and linewidth of the caps
for cap in a4['caps']:
    cap.set(color='black', linewidth=1)

# Change color and linewidth of the medians
for median in a4['medians']:
    median.set(color='black', linewidth=1)

# Change the style of fliers and their fill
for flier in a4['fliers']:
    flier.set(marker='o', color='black', alpha=1)
    
plt.title(u'D) ENEB', fontweight='bold')
plt.ylabel(u'Temperature (°C)', fontweight='bold')
plt.xticks(time2, (u'Reg', u'Had', u'CRU', u'UDEL', 'CHIRPS', 'ERA5'))
plt.yticks(np.arange(20, 30, 2), fontsize=7)
plt.setp(ax4.get_xticklabels(), visible=False)
ax4.xaxis.grid(True, which='major', linestyle='--')
ax4.yaxis.grid(True, which='major', linestyle='--')

ax5 = fig.add_subplot(3, 2, 5)
data = [p_reg, p_had, p_cru, p_udel, p_chirps, p_era5]
a5 = ax5.boxplot(data, patch_artist=True, notch=True, bootstrap=10000, vert=1)

# Change outline and fill color
for box in a5['boxes']:
    box.set( color='black', linewidth=1)
    box.set( facecolor = 'blue' )

# Change color and linewidth of the whiskers
for whisker in a5['whiskers']:
    whisker.set(color='black', linewidth=1)

# Change color and linewidth of the caps
for cap in a5['caps']:
    cap.set(color='black', linewidth=1)

# Change color and linewidth of the medians
for median in a5['medians']:
    median.set(color='black', linewidth=1)

# Change the style of fliers and their fill
for flier in a5['fliers']:
    flier.set(marker='o', color='black', alpha=1)

plt.title(u'E) MATOPIBA', fontweight='bold')
plt.xlabel('Dataset', fontweight='bold')
plt.xticks(time1, (u'Reg', u'Had', u'CRU', u'UDEL', 'CHIRPS', 'ERA5'), fontsize=7)
plt.yticks(np.arange(0, 12, 2), fontsize=7)
ax5.xaxis.grid(True, which='major', linestyle='--')
ax5.yaxis.grid(True, which='major', linestyle='--')

ax6 = fig.add_subplot(3, 2, 6)
data = [t_had, t_had, t_cru, t_udel, t_era5]
a6 = ax6.boxplot(data, patch_artist=True, notch=True, bootstrap=10000, vert=1)

# Change outline and fill color
for box in a6['boxes']:
    box.set( color='black', linewidth=1)
    box.set( facecolor = 'red' )

# Change color and linewidth of the whiskers
for whisker in a6['whiskers']:
    whisker.set(color='black', linewidth=1)

# Change color and linewidth of the caps
for cap in a6['caps']:
    cap.set(color='black', linewidth=1)

# Change color and linewidth of the medians
for median in a6['medians']:
    median.set(color='black', linewidth=1)

# Change the style of fliers and their fill
for flier in a6['fliers']:
    flier.set(marker='o', color='black', alpha=1)
    
plt.title(u'F) MATOPIBA', fontweight='bold')
plt.xlabel('Dataset', fontweight='bold')
plt.xticks(time2, (u'Reg', u'Had', u'CRU', u'UDEL', 'ERA5'), fontsize=7)
plt.yticks(np.arange(20, 30, 2), fontsize=7)
ax6.xaxis.grid(True, which='major', linestyle='--')
ax6.yaxis.grid(True, which='major', linestyle='--')

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.35, hspace=0.35)

# Path out to save bias figure
path_out = '/home/nice/Documents'
name_out = 'pyplt_boxplot_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()


