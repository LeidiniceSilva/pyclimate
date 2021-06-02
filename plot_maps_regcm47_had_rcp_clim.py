# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot climatology maps from Reg and Had models end obs database"

import os
import conda
import netCDF4
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy.ma as ma

# mpl.use('Agg')

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from matplotlib import colors 
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
from os.path import expanduser
from scipy import stats


def import_rcm(var, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_amz_neb_reg_had_{2}_mon_{3}_lonlat_seamask.nc'.format(path, var, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = var[:][:,:,:]

	return lat, lon, rcm


def import_gcm(var, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)	
	arq  = '{0}/{1}_amz_neb_Amon_HadGEM2-ES_{2}_r1i1p1_mon_{3}_lonlat_seamask.nc'.format(path, var, exp, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm  = var[:][:,:,:]

	return lat, lon, gcm
	

def function_ttest(data1, data2):
	
	# Calculate the mean and standard error
	x1_bar, x2_bar = np.nanmean(data1, axis=0), np.nanmean(data2, axis=0)
	var_x1, var_x2= np.var(data1, ddof=1), np.var(data2, ddof=1)

	# Pooled sample variance
	n1, n2 = 240, 240
	pool_var = ( ((n1-1)*var_x1) + ((n2-1)*var_x2) ) / (n1+n2-2)

	# Standard error
	std_error = np.sqrt(pool_var * (1.0 / n1 + 1.0 / n2))

	# Calculate t statistics
	t = abs(x1_bar - x2_bar) / std_error

	# Calculate p value
	p_value = 1 - stats.t.cdf(x=t, df=238)
	
	return p_value

	
	
def basemap(lat, lon):
	
	aux_lon1 = []
	aux_lon2 = []
	for l in lon:
		if l <= 180:
			aux_lon1.append(l)
		else:
			aux_lon2.append(l-360)
		
	lon = np.array(aux_lon1[::-1] + aux_lon2[::-1])
	new_lat = lat
	new_lon = lon[::-1]
	
	map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-20., urcrnrlon=-15.,urcrnrlat=10., resolution='c')
	map.drawmapboundary(color='white')
	#~ map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4)
	#~ map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4)
	
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	lons, lats = np.meshgrid(new_lon, new_lat)

	xx, yy = map(lons,lats)
	
	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='black', linewidth=.5)
	map.readshapefile('{0}/lim_unid_fed/lim_unid_fed'.format(path), 'lim_unid_fed', drawbounds=True, color='black', linewidth=.5)
	
	return map, xx, yy


# Import regcm and hadgem exp
lat, lon, pre_reg_hist = import_rcm('pr', 'hist', '1986-2005')
lat, lon, pre_had_hist = import_gcm('pr', 'hist', '1986-2005')
lat, lon, pre_reg_rcp26 = import_rcm('pr', 'rcp26', '2080-2099')
lat, lon, pre_had_rcp26 = import_gcm('pr', 'rcp26', '2080-2099')
lat, lon, pre_reg_rcp85 = import_rcm('pr', 'rcp85', '2080-2099')
lat, lon, pre_had_rcp85 = import_gcm('pr', 'rcp85', '2080-2099')

lat, lon, tas_reg_hist = import_rcm('tas', 'hist', '1986-2005')
lat, lon, tas_had_hist = import_gcm('tas', 'hist', '1986-2005')
lat, lon, tas_reg_rcp26 = import_rcm('tas', 'rcp26', '2080-2099')
lat, lon, tas_had_rcp26 = import_gcm('tas', 'rcp26', '2080-2099')
lat, lon, tas_reg_rcp85 = import_rcm('tas', 'rcp85', '2080-2099')
lat, lon, tas_had_rcp85 = import_gcm('tas', 'rcp85', '2080-2099')

# Compute change from rcp and historical period
diff_pre_reg_rcp26_hist = np.nanmean(pre_reg_rcp26, axis=0) - np.nanmean(pre_reg_hist, axis=0)
diff_pre_had_rcp26_hist = np.nanmean(pre_had_rcp26, axis=0) - np.nanmean(pre_had_hist, axis=0)
diff_pre_reg_rcp85_hist = np.nanmean(pre_reg_rcp85, axis=0) - np.nanmean(pre_reg_hist, axis=0)
diff_pre_had_rcp85_hist = np.nanmean(pre_had_rcp85, axis=0) - np.nanmean(pre_had_hist, axis=0)

diff_tas_reg_rcp26_hist = np.nanmean(np.nanmean(tas_reg_rcp26, axis=0), axis=0) - np.nanmean(np.nanmean(tas_reg_hist, axis=0), axis=0)
diff_tas_had_rcp26_hist = np.nanmean(tas_had_rcp26, axis=0) - np.nanmean(tas_had_hist, axis=0)
diff_tas_reg_rcp85_hist = np.nanmean(np.nanmean(tas_reg_rcp85, axis=0), axis=0) - np.nanmean(np.nanmean(tas_reg_hist, axis=0), axis=0)
diff_tas_had_rcp85_hist = np.nanmean(tas_had_rcp85, axis=0) - np.nanmean(tas_had_hist, axis=0)

p_value_pre_reg_rcp26_hist = function_ttest(pre_reg_rcp26, pre_reg_hist)
p_value_pre_had_rcp26_hist = function_ttest(pre_had_rcp26, pre_had_hist)
p_value_pre_reg_rcp85_hist = function_ttest(pre_reg_rcp85, pre_reg_hist)
p_value_pre_had_rcp85_hist = function_ttest(pre_had_rcp85, pre_had_hist)

p_value_tas_reg_rcp26_hist = function_ttest(tas_reg_rcp26, tas_reg_hist)
p_value_tas_had_rcp26_hist = function_ttest(tas_had_rcp26, tas_had_hist)
p_value_tas_reg_rcp85_hist = function_ttest(tas_reg_rcp85, tas_reg_hist)
p_value_tas_had_rcp85_hist = function_ttest(tas_had_rcp85, tas_had_hist)

# Plot bias maps 
fig = plt.figure()
levs1 = [-2, -1.5, -1, -0.5, -0.1, 0.1, 0.5, 1, 1.5, 2]
levs2 = [0, 0.5, 1, 1.5, 2, 3, 4, 5, 6]

ax = fig.add_subplot(4, 2, 1)
plt.title(u'A)', loc='left', fontsize=6, fontweight='bold')
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, diff_pre_reg_rcp26_hist, latlon=True, levels=levs1, cmap=cm.BrBG, extend='both')
p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 2, 2)
plt.title(u'B)', loc='left', fontsize=6, fontweight='bold')
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, diff_pre_had_rcp26_hist, latlon=True, levels=levs1, cmap=plt.cm.BrBG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)                                
p_value_pre_had_rcp26_hist = ma.masked_where(p_value_pre_had_rcp26_hist >= 0.1, p_value_pre_had_rcp26_hist) 
map.contourf(xx, yy, p_value_pre_had_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 2, 3)
plt.title(u'C)', loc='left', fontsize=6, fontweight='bold')
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, diff_pre_reg_rcp85_hist, llatlon=True, levels=levs1, cmap=cm.BrBG, extend='both')
p_value_pre_reg_rcp85_hist = ma.masked_where(p_value_pre_reg_rcp85_hist >= 0.1, p_value_pre_reg_rcp85_hist) 
map.contourf(xx, yy, p_value_pre_reg_rcp85_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 2, 4)
plt.title(u'D)', loc='left', fontsize=6, fontweight='bold')
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, diff_pre_had_rcp85_hist, latlon=True, levels=levs1, cmap=cm.BrBG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
p_value_pre_had_rcp85_hist = ma.masked_where(p_value_pre_had_rcp85_hist >= 0.1, p_value_pre_had_rcp85_hist) 
map.contourf(xx, yy, p_value_pre_had_rcp85_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 2, 5)
plt.title(u'E)', loc='left', fontsize=6, fontweight='bold')
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, diff_tas_reg_rcp26_hist, latlon=True, levels=levs2, cmap=cm.YlOrRd, extend='both')
p_value_tas_reg_rcp26_hist = ma.masked_where(p_value_tas_reg_rcp26_hist >= 0.001, p_value_tas_reg_rcp26_hist) 
map.contourf(xx, yy, p_value_tas_reg_rcp26_hist[0,:,:], colors='none', hatches=["////"])

ax = fig.add_subplot(4, 2, 6)
plt.title(u'F)', loc='left', fontsize=6, fontweight='bold')
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, diff_tas_had_rcp26_hist, latlon=True, levels=levs2, cmap=plt.cm.YlOrRd, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)                                
p_value_tas_had_rcp26_hist = ma.masked_where(p_value_tas_had_rcp26_hist >= 0.001, p_value_tas_had_rcp26_hist) 
map.contourf(xx, yy, p_value_tas_had_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 2, 7)
plt.title(u'G)', loc='left', fontsize=6, fontweight='bold')
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, diff_tas_reg_rcp85_hist, llatlon=True, levels=levs2, cmap=cm.YlOrRd, extend='both')
p_value_tas_reg_rcp85_hist = ma.masked_where(p_value_tas_reg_rcp85_hist >= 0.001, p_value_tas_reg_rcp85_hist) 
map.contourf(xx, yy, p_value_tas_reg_rcp85_hist[0,:,:], colors='none', hatches=["////"])

ax = fig.add_subplot(4, 2, 8)
plt.title(u'H)', loc='left', fontsize=6, fontweight='bold')
map, xx, yy = basemap(lat, lon)
map.contourf(xx, yy, diff_tas_had_rcp85_hist, latlon=True, levels=levs2, cmap=cm.YlOrRd, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)
p_value_tas_had_rcp85_hist = ma.masked_where(p_value_tas_had_rcp85_hist >= 0.001, p_value_tas_had_rcp85_hist) 
map.contourf(xx, yy, p_value_tas_had_rcp85_hist, colors='none', hatches=["////"])

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_diff_reg_had_rcp-hist.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')

plt.show()
exit()


