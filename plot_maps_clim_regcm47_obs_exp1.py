# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot seasonal climatology maps from regcm47 and hadgem models and obs database"

import os
import conda
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib as mpl 
import matplotlib.cm as cm
import scipy.stats as stats
import matplotlib.pyplot as plt

# mpl.use('Agg')

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from os.path import expanduser
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm


def import_obs(var, area, dataset, period, dt):
	
	path = '/home/nice/Documents/dataset/obs/reg_exp1'
	arq  = '{0}/{1}_{2}_{3}_obs_{4}_{5}_lonlat.nc'.format(path, var, area, dataset, period, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	mean_obs = np.nanmean(value, axis=0)

	return lat, lon, mean_obs
	
	
def import_rcm(var, area, exp, period, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp1/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_{4}_{5}_lonlat_seamask.nc'.format(path, var, area, exp, period, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	mean_rcm = np.nanmean(value, axis=0)

	return lat, lon, mean_rcm
	

def import_gcm(var, area, exp, period, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp1/{0}'.format(exp)
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_{4}_{5}_lonlat_seamask.nc'.format(path, var, area, exp, period, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	mean_gcm = np.nanmean(value, axis=0)

	return lat, lon, mean_gcm
	
	
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
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=.8)

	return map, xx, yy


# Import models and obs database 
lat, lon, pre_djf_cru = import_obs('pre', 'amz_neb', 'cru_ts4.04', 'djf', '1986-2005')	   
lat, lon, pre_mam_cru = import_obs('pre', 'amz_neb', 'cru_ts4.04', 'mam', '1986-2005')	   
lat, lon, pre_jja_cru = import_obs('pre', 'amz_neb', 'cru_ts4.04', 'jja', '1986-2005')	   
lat, lon, pre_son_cru = import_obs('pre', 'amz_neb', 'cru_ts4.04', 'son', '1986-2005')	   
lat, lon, pre_ann_cru = import_obs('pre', 'amz_neb', 'cru_ts4.04', 'ann', '1986-2005')	   

lat, lon, pre_djf_rcm = import_rcm('pr', 'amz_neb', 'hist', 'djf', '1986-2005')	   
lat, lon, pre_mam_rcm = import_rcm('pr', 'amz_neb', 'hist', 'mam', '1986-2005')	   
lat, lon, pre_jja_rcm = import_rcm('pr', 'amz_neb', 'hist', 'jja', '1986-2005')	   
lat, lon, pre_son_rcm = import_rcm('pr', 'amz_neb', 'hist', 'son', '1986-2005')	   
lat, lon, pre_ann_rcm = import_rcm('pr', 'amz_neb', 'hist', 'ann', '1986-2005')	   

lat, lon, pre_djf_gcm = import_gcm('pr', 'amz_neb', 'hist', 'djf', '1986-2005')	   
lat, lon, pre_mam_gcm = import_gcm('pr', 'amz_neb', 'hist', 'mam', '1986-2005')	   
lat, lon, pre_jja_gcm = import_gcm('pr', 'amz_neb', 'hist', 'jja', '1986-2005')	   
lat, lon, pre_son_gcm = import_gcm('pr', 'amz_neb', 'hist', 'son', '1986-2005')	   
lat, lon, pre_ann_gcm = import_gcm('pr', 'amz_neb', 'hist', 'ann', '1986-2005')	   

lat, lon, tas_djf_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', 'djf', '1986-2005')	   
lat, lon, tas_mam_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', 'mam', '1986-2005')	   
lat, lon, tas_jja_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', 'jja', '1986-2005')	   
lat, lon, tas_son_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', 'son', '1986-2005')	   
lat, lon, tas_ann_cru = import_obs('tmp', 'amz_neb', 'cru_ts4.04', 'ann', '1986-2005')	   

lat, lon, tas_djf_rcm = import_rcm('tas', 'amz_neb', 'hist', 'djf', '1986-2005')	   
lat, lon, tas_mam_rcm = import_rcm('tas', 'amz_neb', 'hist', 'mam', '1986-2005')	   
lat, lon, tas_jja_rcm = import_rcm('tas', 'amz_neb', 'hist', 'jja', '1986-2005')	   
lat, lon, tas_son_rcm = import_rcm('tas', 'amz_neb', 'hist', 'son', '1986-2005')	   
lat, lon, tas_ann_rcm = import_rcm('tas', 'amz_neb', 'hist', 'ann', '1986-2005')	   

lat, lon, tas_djf_gcm = import_gcm('tas', 'amz_neb', 'hist', 'djf', '1986-2005')	   
lat, lon, tas_mam_gcm = import_gcm('tas', 'amz_neb', 'hist', 'mam', '1986-2005')	   
lat, lon, tas_jja_gcm = import_gcm('tas', 'amz_neb', 'hist', 'jja', '1986-2005')	   
lat, lon, tas_son_gcm = import_gcm('tas', 'amz_neb', 'hist', 'son', '1986-2005')	   
lat, lon, tas_ann_gcm = import_gcm('tas', 'amz_neb', 'hist', 'ann', '1986-2005')	   

# Plot models and obs database 
fig = plt.figure(figsize=(7,7))
levs1 = [0, 3, 6, 9, 12, 15]
levs2 = [19, 21, 24, 27, 30, 33]

ax = fig.add_subplot(5, 3, 1)
plt.title(u'A) CRU DJF', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_djf_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 2)
plt.title(u'B) RegCM4.7 DJF', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_djf_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 3)
plt.title(u'C) HadGEM2-ES DJF', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_djf_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu, extend='max') 
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Precipitation \n (mm d⁻¹)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 4)
plt.title(u'D) CRU MAM', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_mam_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 5)
plt.title(u'E) RegCM4.7 MAM', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_mam_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu)  
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(5, 3, 6)
plt.title(u'F) HadGEM2-ES MAM', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_mam_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu, extend='max')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Precipitation \n (mm d⁻¹)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 7)
plt.title(u'G) CRU JJA', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_jja_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 8)
plt.title(u'H) RegCM4.7 JJA', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_jja_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
map.drawmeridians(np.arange(-85.,-15.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,10.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 9)
plt.title(u'I) HadGEM2-ES JJA', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_jja_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu, extend='max') 
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Precipitation \n (mm d⁻¹)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(5, 3, 10)
plt.title(u'J) CRU SON', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_son_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 11)
plt.title(u'K) RegCM4.7 SON', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_son_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu)
map.drawmeridians(np.arange(-85.,-15.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,10.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 12)
plt.title(u'L) HadGEM2-ES SON', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_son_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu, extend='max')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Precipitation \n (mm d⁻¹)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 13)
plt.title(u'M) CRU ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_ann_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 14)
plt.title(u'N) RegCM4.7 ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_ann_rcm, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(5, 3, 15)
plt.title(u'O) HadGEM2-ES ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, pre_ann_gcm, levels=levs1, latlon=True, cmap=cm.YlGnBu, extend='max')  
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Precipitation \n (mm d⁻¹)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_pre_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

ax = fig.add_subplot(5, 3, 1)
plt.title(u'A) CRU DJF', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_djf_cru, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 2)
plt.title(u'B) RegCM4.7 DJF', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_djf_rcm[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 3)
plt.title(u'C) HadGEM2-ES DJF', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_djf_gcm, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='max') 
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Temperature \n (°C)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 4)
plt.title(u'D) CRU MAM', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_mam_cru, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(5, 3, 5)
plt.title(u'E) RegCM4.7 MAM', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_mam_rcm[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd)  
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(5, 3, 6)
plt.title(u'F) HadGEM2-ES MAM', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_mam_gcm, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='max')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Temperature \n (°C)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 7)
plt.title(u'G) CRU JJA', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_jja_cru, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 8)
plt.title(u'H) RegCM4.7 JJA', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_jja_rcm[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-15.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,10.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 9)
plt.title(u'I) HadGEM2-ES JJA', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_jja_gcm, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='max') 
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Temperature \n (°C)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(5, 3, 10)
plt.title(u'J) CRU SON', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_son_cru, levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 11)
plt.title(u'K) RegCM4.7 SON', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_son_rcm[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd)
map.drawmeridians(np.arange(-85.,-15.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,10.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 12)
plt.title(u'L) HadGEM2-ES SON', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_son_gcm, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='max')
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Temperature \n (°C)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 13)
plt.title(u'M) CRU ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_ann_cru, levels=levs2, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

ax = fig.add_subplot(5, 3, 14)
plt.title(u'N) RegCM4.7 ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_ann_rcm[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd) 
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')
	
ax = fig.add_subplot(5, 3, 15)
plt.title(u'O) HadGEM2-ES ANN', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plot_maps_mean = map.contourf(xx, yy, tas_ann_gcm, levels=levs2, latlon=True, cmap=cm.YlOrRd, extend='max')  
cbar = map.colorbar()
cbar.ax.tick_params(labelsize=6) 
cbar.set_label('Temperature \n (°C)', size=6, labelpad=5, fontweight='bold', rotation=90)
map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[0,0,0,0], linewidth=0.4, color='black')

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_tas_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()









