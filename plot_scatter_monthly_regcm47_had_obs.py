# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script scatter plot from Reg and Had models end obs database"

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

fig=plt.figure()

ax1 = plt.subplot(3, 5, 1)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax1.set_xticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax1.text(-2.0,17.3, '   Annual  ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
             
ax2 = plt.subplot(3, 5, 2)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax2.set_xticklabels([])
ax2.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax2.text(-2.0,17.3, '     DJF     ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
                
ax3 = plt.subplot(3, 5, 3)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax3.set_xticklabels([])
ax3.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax3.text(-2.0,17.3, '    MAM    ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})

ax4 = plt.subplot(3, 5, 4)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax4.set_xticklabels([])
ax4.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax4.text(-2.0,17.3, '      JJA     ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
   
ax5 = plt.subplot(3, 5, 5)
plt.scatter(3, 10, s=80, c='blue', marker='o', label='Had (RCP2.6)')
plt.scatter(3, 14, s=80, c='blue', marker='D', label='Reg_Had (RCP2.6)')
plt.scatter(-2, -10, s=80, c='red', marker='o', label='Had (RCP8.55)')
plt.scatter(1, -5, s=80, c='red', marker='D', label='Reg_Had (RCP8.55)')
ax5.set_xticklabels([])
ax5.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax5.text(-2.0,17.3, '     SON    ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
title = ax5.text(3.75,11., '      AMZ       ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
      
ax6 = plt.subplot(3, 5, 6)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax6.set_xticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
plt.ylabel(u'Precipitation change (%)', fontweight='bold')

ax7 = plt.subplot(3, 5, 7)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax7.set_xticklabels([])
ax7.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax8 = plt.subplot(3, 5, 8)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax8.set_xticklabels([])
ax8.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax9 = plt.subplot(3, 5, 9)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-1, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax9.set_xticklabels([])
ax9.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax10 = plt.subplot(3, 5, 10)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax10.set_xticklabels([])
ax10.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax10.text(3.75,11., '      NEB       ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
      
ax11 = plt.subplot(3, 5, 11)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax12 = plt.subplot(3, 5, 12)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax12.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax13 = plt.subplot(3, 5, 13)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax13.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
plt.xlabel(u'Temperature change (°C)', fontweight='bold')

ax14 = plt.subplot(3, 5, 14)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax14.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax15 = plt.subplot(3, 5, 15)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax15.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax15.text(3.75,11., ' MATOPIBA  ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.6, 'pad':4})
                
ax5.legend(bbox_to_anchor=(1.4, 1), loc=2, borderaxespad=0.5)

# Path out to save bias figure
path_out = '/home/nice/Documents'
name_out = 'pyplt_scatter_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()


