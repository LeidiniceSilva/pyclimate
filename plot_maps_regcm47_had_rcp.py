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

# mpl.use('Agg')

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from matplotlib import colors as c
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
from os.path import expanduser
    
    
def import_rcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_rcm = value[2:240:3,:,:]
	djf_rcm = np.nanmean(season_rcm[3:80:4], axis=0)
	mam_rcm = np.nanmean(season_rcm[0:80:4], axis=0)
	jja_rcm = np.nanmean(season_rcm[1:80:4], axis=0)
	son_rcm = np.nanmean(season_rcm[2:80:4], axis=0)
	annual_rcm = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, djf_rcm, mam_rcm, jja_rcm, son_rcm, annual_rcm


def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/{0}'.format(exp)
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_gcm = value[2:240:3,:,:]
	djf_gcm = np.nanmean(season_gcm[3:80:4], axis=0)
	mam_gcm = np.nanmean(season_gcm[0:80:4], axis=0)
	jja_gcm = np.nanmean(season_gcm[1:80:4], axis=0)
	son_gcm = np.nanmean(season_gcm[2:80:4], axis=0)
	annual_gcm = np.nanmean(value[0:240:12,:,:], axis=0)

	return lat, lon, djf_gcm, mam_gcm, jja_gcm, son_gcm, annual_gcm
	

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
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	
	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='gray', linewidth=.5)
	map.readshapefile('{0}/lim_unid_fed/lim_unid_fed'.format(path), 'lim_unid_fed', drawbounds=True, color='black', linewidth=.5)
	
	return map, xx, yy
	

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
	
	
# Import regcm and hadgem exp
lat, lon, pre_djf_rcm_hist, pre_mam_rcm_hist, pre_jja_rcm_hist, pre_son_rcm_hist, pre_annual_rcm_hist = import_rcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, pre_djf_gcm_hist, pre_mam_gcm_hist, pre_jja_gcm_hist, pre_son_gcm_hist, pre_annual_gcm_hist = import_gcm('pr', 'amz_neb', 'hist', '1986-2005')
lat, lon, pre_djf_rcm_rcp26, pre_mam_rcm_rcp26, pre_jja_rcm_rcp26, pre_son_rcm_rcp26, pre_annual_rcm_rcp26 = import_rcm('pr', 'amz_neb', 'rcp26', '2080-2099')
lat, lon, pre_djf_gcm_rcp26, pre_mam_gcm_rcp26, pre_jja_gcm_rcp26, pre_son_gcm_rcp26, pre_annual_gcm_rcp26 = import_gcm('pr', 'amz_neb', 'rcp26', '2080-2099')
lat, lon, pre_djf_rcm_rcp85, pre_mam_rcm_rcp85, pre_jja_rcm_rcp85, pre_son_rcm_rcp85, pre_annual_rcm_rcp85 = import_rcm('pr', 'amz_neb', 'rcp85', '2080-2099')
lat, lon, pre_djf_gcm_rcp85, pre_mam_gcm_rcp85, pre_jja_gcm_rcp85, pre_son_gcm_rcp85, pre_annual_gcm_rcp85 = import_gcm('pr', 'amz_neb', 'rcp85', '2080-2099')

lat, lon, tas_djf_rcm_hist, tas_mam_rcm_hist, tas_jja_rcm_hist, tas_son_rcm_hist, tas_annual_rcm_hist = import_rcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, tas_djf_gcm_hist, tas_mam_gcm_hist, tas_jja_gcm_hist, tas_son_gcm_hist, tas_annual_gcm_hist = import_gcm('tas', 'amz_neb', 'hist', '1986-2005')
lat, lon, tas_djf_rcm_rcp26, tas_mam_rcm_rcp26, tas_jja_rcm_rcp26, tas_son_rcm_rcp26, tas_annual_rcm_rcp26 = import_rcm('tas', 'amz_neb', 'rcp26', '2080-2099')
lat, lon, tas_djf_gcm_rcp26, tas_mam_gcm_rcp26, tas_jja_gcm_rcp26, tas_son_gcm_rcp26, tas_annual_gcm_rcp26 = import_gcm('tas', 'amz_neb', 'rcp26', '2080-2099')
lat, lon, tas_djf_rcm_rcp85, tas_mam_rcm_rcp85, tas_jja_rcm_rcp85, tas_son_rcm_rcp85, tas_annual_rcm_rcp85 = import_rcm('tas', 'amz_neb', 'rcp85', '2080-2099')
lat, lon, tas_djf_gcm_rcp85, tas_mam_gcm_rcp85, tas_jja_gcm_rcp85, tas_son_gcm_rcp85, tas_annual_gcm_rcp85 = import_gcm('tas', 'amz_neb', 'rcp85', '2080-2099')

# Compute change from rcp and historical period
pre_djf_rcm_rcp26_hist = pre_djf_rcm_rcp26 - pre_djf_rcm_hist
pre_mam_rcm_rcp26_hist = pre_mam_rcm_rcp26 - pre_mam_rcm_hist
pre_jja_rcm_rcp26_hist = pre_jja_rcm_rcp26 - pre_jja_rcm_hist
pre_son_rcm_rcp26_hist = pre_son_rcm_rcp26 - pre_son_rcm_hist
pre_annual_rcm_rcp26_hist = pre_annual_rcm_rcp26 - pre_annual_rcm_hist
pre_djf_rcm_rcp85_hist = pre_djf_rcm_rcp85 - pre_djf_rcm_hist
pre_mam_rcm_rcp85_hist = pre_mam_rcm_rcp85 - pre_mam_rcm_hist
pre_jja_rcm_rcp85_hist = pre_jja_rcm_rcp85 - pre_jja_rcm_hist
pre_son_rcm_rcp85_hist = pre_son_rcm_rcp85 - pre_son_rcm_hist
pre_annual_rcm_rcp85_hist = pre_annual_rcm_rcp85 - pre_annual_rcm_hist
pre_djf_gcm_rcp26_hist = pre_djf_gcm_rcp26 - pre_djf_gcm_hist
pre_mam_gcm_rcp26_hist = pre_mam_gcm_rcp26 - pre_mam_gcm_hist
pre_jja_gcm_rcp26_hist = pre_jja_gcm_rcp26 - pre_jja_gcm_hist
pre_son_gcm_rcp26_hist = pre_son_gcm_rcp26 - pre_son_gcm_hist
pre_annual_gcm_rcp26_hist = pre_annual_gcm_rcp26 - pre_annual_gcm_hist
pre_djf_gcm_rcp85_hist = pre_djf_gcm_rcp85 - pre_djf_gcm_hist
pre_mam_gcm_rcp85_hist = pre_mam_gcm_rcp85 - pre_mam_gcm_hist
pre_jja_gcm_rcp85_hist = pre_jja_gcm_rcp85 - pre_jja_gcm_hist
pre_son_gcm_rcp85_hist = pre_son_gcm_rcp85 - pre_son_gcm_hist
pre_annual_gcm_rcp85_hist = pre_annual_gcm_rcp85 - pre_annual_gcm_hist

tas_djf_rcm_rcp26_hist = np.nanmean(tas_djf_rcm_rcp26, axis=0) - np.nanmean(tas_djf_rcm_hist, axis=0)
tas_mam_rcm_rcp26_hist = np.nanmean(tas_mam_rcm_rcp26, axis=0) - np.nanmean(tas_mam_rcm_hist, axis=0)
tas_jja_rcm_rcp26_hist = np.nanmean(tas_jja_rcm_rcp26, axis=0) - np.nanmean(tas_jja_rcm_hist, axis=0)
tas_son_rcm_rcp26_hist = np.nanmean(tas_son_rcm_rcp26, axis=0) - np.nanmean(tas_son_rcm_hist, axis=0)
tas_annual_rcm_rcp26_hist = np.nanmean(tas_annual_rcm_rcp26, axis=0) - np.nanmean(tas_annual_rcm_hist, axis=0)
tas_djf_rcm_rcp85_hist = np.nanmean(tas_djf_rcm_rcp85, axis=0) - np.nanmean(tas_djf_rcm_hist, axis=0)
tas_mam_rcm_rcp85_hist = np.nanmean(tas_mam_rcm_rcp85, axis=0) - np.nanmean(tas_mam_rcm_hist, axis=0)
tas_jja_rcm_rcp85_hist = np.nanmean(tas_jja_rcm_rcp85, axis=0) - np.nanmean(tas_jja_rcm_hist, axis=0)
tas_son_rcm_rcp85_hist = np.nanmean(tas_son_rcm_rcp85, axis=0) - np.nanmean(tas_son_rcm_hist, axis=0)
tas_annual_rcm_rcp85_hist = np.nanmean(tas_annual_rcm_rcp85, axis=0) - np.nanmean(tas_annual_rcm_hist, axis=0)
tas_djf_gcm_rcp26_hist = tas_djf_gcm_rcp26 - tas_djf_gcm_hist
tas_mam_gcm_rcp26_hist = tas_mam_gcm_rcp26 - tas_mam_gcm_hist
tas_jja_gcm_rcp26_hist = tas_jja_gcm_rcp26 - tas_jja_gcm_hist
tas_son_gcm_rcp26_hist = tas_son_gcm_rcp26 - tas_son_gcm_hist
tas_annual_gcm_rcp26_hist = tas_annual_gcm_rcp26 - tas_annual_gcm_hist
tas_djf_gcm_rcp85_hist = tas_djf_gcm_rcp85 - tas_djf_gcm_hist
tas_mam_gcm_rcp85_hist = tas_mam_gcm_rcp85 - tas_mam_gcm_hist
tas_jja_gcm_rcp85_hist = tas_jja_gcm_rcp85 - tas_jja_gcm_hist
tas_son_gcm_rcp85_hist = tas_son_gcm_rcp85 - tas_son_gcm_hist
tas_annual_gcm_rcp85_hist = tas_annual_gcm_rcp85 - tas_annual_gcm_hist

# Compute ttest
#~ p_value_pre_djf_rcm_rcp26_hist = function_ttest(pre_reg_rcp26, pre_reg_hist)
#~ p_value_pre_had_rcp26_hist = function_ttest(pre_had_rcp26, pre_had_hist)
#~ p_value_pre_reg_rcp85_hist = function_ttest(pre_reg_rcp85, pre_reg_hist)
#~ p_value_pre_had_rcp85_hist = function_ttest(pre_had_rcp85, pre_had_hist)

#~ p_value_tas_reg_rcp26_hist = function_ttest(tas_reg_rcp26, tas_reg_hist)
#~ p_value_tas_had_rcp26_hist = function_ttest(tas_had_rcp26, tas_had_hist)
#~ p_value_tas_reg_rcp85_hist = function_ttest(tas_reg_rcp85, tas_reg_hist)
#~ p_value_tas_had_rcp85_hist = function_ttest(tas_had_rcp85, tas_had_hist)

# Plot bias maps 
fig = plt.figure(figsize=(8, 4))
levs1 = [0.5, 1, 1.5, 2.5, 3, 3.5, 4]
levs2 = [1, 1.5, 2, 4, 6, 8, 10]

ax = fig.add_subplot(4, 5, 1)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_djf_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist_djf = ma.masked_where(p_value_pre_reg_rcp26_hist_djf >= 0.3, p_value_pre_reg_rcp26_hist_djf) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist_djf, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 2)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_mam_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist_mam = ma.masked_where(p_value_pre_reg_rcp26_hist_mam >= 0.3, p_value_pre_reg_rcp26_hist_mam) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist_mam, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 3)
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_jja_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist_jja = ma.masked_where(p_value_pre_reg_rcp26_hist_jja >= 0.3, p_value_pre_reg_rcp26_hist_jja) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist_jja, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 4)
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_son_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 5)
plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_annual_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 6)
plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_djf_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 7)
plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_mam_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 8)
plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_jja_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 9)
plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_son_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 10)
plt.title(u'J)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_annual_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.YlOrRd, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 11)
plt.title(u'K)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_djf_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 12)
plt.title(u'L)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_mam_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 13)
plt.title(u'M)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_jja_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 14)
plt.title(u'N)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_son_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 15)
plt.title(u'O)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_annual_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 16)
plt.title(u'P)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_djf_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 17)
plt.title(u'Q)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_mam_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 18)
plt.title(u'R)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_jja_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 19)
plt.title(u'S)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_son_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 20)
plt.title(u'T)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, tas_annual_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])
	
# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_diff_reg_had_rcp-hist_tas.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()

levs1 = [-2, -1, -0.5, 0.5, 1, 2]
levs2 = [-4, -2, -1, 1, 2, 4]

ax = fig.add_subplot(4, 5, 1)
plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_djf_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist_djf = ma.masked_where(p_value_pre_reg_rcp26_hist_djf >= 0.3, p_value_pre_reg_rcp26_hist_djf) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist_djf, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 2)
plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_mam_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist_mam = ma.masked_where(p_value_pre_reg_rcp26_hist_mam >= 0.3, p_value_pre_reg_rcp26_hist_mam) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist_mam, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 3)
plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_jja_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist_jja = ma.masked_where(p_value_pre_reg_rcp26_hist_jja >= 0.3, p_value_pre_reg_rcp26_hist_jja) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist_jja, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 4)
plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_son_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 5)
plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_annual_rcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 6)
plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_djf_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 7)
plt.title(u'G)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_mam_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 8)
plt.title(u'H)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_jja_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 9)
plt.title(u'I)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_son_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 10)
plt.title(u'J)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_annual_gcm_rcp26_hist, levels=levs1, latlon=True, cmap=cm.BrBG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 11)
plt.title(u'K)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_djf_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 12)
plt.title(u'L)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_mam_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 13)
plt.title(u'M)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_jja_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 14)
plt.title(u'N)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_son_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 15)
plt.title(u'O)', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_annual_rcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 16)
plt.title(u'P)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_djf_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 17)
plt.title(u'Q)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_mam_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 18)
plt.title(u'R)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_jja_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 19)
plt.title(u'S)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_son_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

ax = fig.add_subplot(4, 5, 20)
plt.title(u'T)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_diff = map.contourf(xx, yy, pre_annual_gcm_rcp85_hist, levels=levs2, latlon=True, cmap=cm.BrBG, extend='both')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6)     
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
#~ p_value_pre_reg_rcp26_hist = ma.masked_where(p_value_pre_reg_rcp26_hist >= 0.3, p_value_pre_reg_rcp26_hist) 
#~ map.contourf(xx, yy, p_value_pre_reg_rcp26_hist, colors='none', hatches=["////"])

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_diff_reg_had_rcp-hist_pre.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()







