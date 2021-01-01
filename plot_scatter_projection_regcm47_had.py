# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script scatter plot from Reg and Had models"

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


def import_rcm(var, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_amz_neb_reg_had_{2}_mon_{3}_lonlat_seamask.nc'.format(path, var, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, rcm


def import_gcm(var, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)	
	arq  = '{0}/{1}_amz_neb_Amon_HadGEM2-ES_{2}_r1i1p1_mon_{3}_lonlat_seamask.nc'.format(path, var, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, gcm
	

# Import regcm and hadgem exp 
lat, lon, pre_reg_hist = import_rcm('pr', 'hist', '1986-2005')
lat, lon, pre_had_hist = import_gcm('pr', 'hist', '1986-2005')
lat, lon, pre_reg_r26 = import_rcm('pr', 'rcp26', '2080-2099')
lat, lon, pre_had_r26 = import_gcm('pr', 'rcp26', '2080-2099')
lat, lon, pre_reg_r85 = import_rcm('pr', 'rcp85', '2080-2099')
lat, lon, pre_had_r85 = import_gcm('pr', 'rcp85', '2080-2099')

lat, lon, tas_reg_hist = import_rcm('tas', 'hist', '1986-2005')
lat, lon, tas_had_hist = import_gcm('tas', 'hist', '1986-2005')
lat, lon, tas_reg_r26 = import_rcm('tas', 'rcp26', '2080-2099')
lat, lon, tas_had_r26 = import_gcm('tas', 'rcp26', '2080-2099')
lat, lon, tas_reg_r85 = import_rcm('tas', 'rcp85', '2080-2099')
lat, lon, tas_had_r85 = import_gcm('tas', 'rcp85', '2080-2099')

# Compute skill metrics
diff_pre_reg_r26_hist = pre_reg_r26 - pre_reg_hist
diff_pre_had_r26_hist = pre_had_r26 - pre_had_hist
diff_pre_reg_r85_hist = pre_reg_r85 - pre_reg_hist
diff_pre_had_r85_hist = pre_had_r85 - pre_had_hist

diff_tas_reg_r26_hist = tas_reg_r26 - tas_reg_hist
diff_tas_had_r26_hist = tas_had_r26 - tas_had_hist
diff_tas_reg_r85_hist = tas_reg_r85 - tas_reg_hist
diff_tas_had_r85_hist = tas_had_r85 - tas_had_hist

# Compute skill metrics
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
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax5.set_xticklabels([])
ax5.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax5.text(-2.0,17.3, '     SON    ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
title = ax5.text(3.75,11., '     SAMZ      ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
      
ax6 = plt.subplot(3, 5, 6)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax6.set_xticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
plt.ylabel(u'Precipitation (mm d⁻¹)', fontweight='bold')

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
title = ax10.text(3.75,11., '     ENEB      ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
      
ax11 = plt.subplot(3, 5, 11)
plt.scatter(3, 10, s=80, c='blue', marker='o', label='Reg (RCP26)')
plt.scatter(3, 14, s=80, c='blue', marker='D', label='Had (RCP26)')
plt.scatter(-2, -10, s=80, c='red', marker='o', label='Reg (RCP85)')
plt.scatter(1, -5, s=80, c='red', marker='D', label='Had (RCP85)')
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
ax11.legend(loc='lower left', bbox_to_anchor=(-0.7, -0.7), shadow=True, ncol=4)

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
plt.xlabel(u'Temperature (°C)', fontweight='bold')

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
                
# Path out to save bias figure
path_out = '/home/nice/Documents'
name_out = 'pyplt_scatter_reg_had_rcp-hist.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()


