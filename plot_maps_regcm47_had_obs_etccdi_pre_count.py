# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/01/2021"
__description__ = "This script plot maps from extremes index"

import os
import conda
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import warnings ; warnings.filterwarnings("ignore")
import matplotlib.cm as cm

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from matplotlib.path import Path
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import PathPatch


def import_obs(var, area, dataset, freq, dt):
	
	path = '/home/nice/Documents/dataset/obs/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, var, area, dataset, freq, dt)	

	dict_var = {u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}

	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs  = np.nanmean(var[:][:,:,:], axis=0)
	
	return lat, lon, obs
	
	
def import_rcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/rcm/eca'	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}

	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, rcm


def import_gcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/gcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}

	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, gcm


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
	
	map = Basemap(projection='cyl', llcrnrlat=-16, urcrnrlat=-2, llcrnrlon=-52, urcrnrlon=-41, resolution=None, suppress_ticks=True, lon_0=0, celestial=False)
	map.drawmeridians(np.arange(-52.,-41.,4.), size=6, labels=[0,0,0,1], linewidth=0.5)
	map.drawparallels(np.arange(-16.,-2.,4.), size=6, labels=[1,0,0,0], linewidth=0.5)
 
	lons, lats = np.meshgrid(new_lon, new_lat)

	xx, yy = map(lons,lats)

	# Import shapefile from word and matopiba 
	map.readshapefile('/home/nice/Documents/github_projects/shp/shp_world/world', 'world', drawbounds=True, color='white')
	map.readshapefile('/home/nice/Documents/github_projects/shp/shp_matopiba/matopiba', 'matopiba', drawbounds=True, color='black')
	x0, x1 = plt.xlim()
	y0, y1 = plt.ylim()
	map_edges = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
	polys = [map_edges]
	map.readshapefile('/home/nice/Documents/github_projects/shp/shp_matopiba/matopiba', 'matopiba2', drawbounds=False)
	polys = polys + getattr(map, 'matopiba2')
	codes = [[Path.MOVETO] + [Path.LINETO for p in p[1:]] for p in polys] # creating a PathPatch
	polys_lin = [v for p in polys for v in p]
	codes_lin = [cdg for cdgs in codes for cdg in cdgs]
	path  = Path(polys_lin, codes_lin)
	patch = PathPatch(path, facecolor='white', lw=0)
	plt.gca().add_patch(patch)
	
	return map, xx, yy
	
	
# Import regcm exp and cru databases 	
lat, lon, obs_cdd = import_obs('eca_cdd', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
lat, lon, rcm_cdd = import_rcm('eca_cdd', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_cdd = import_gcm('eca_cdd', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_cwd = import_obs('eca_cwd', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
lat, lon, rcm_cwd = import_rcm('eca_cwd', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_cwd = import_gcm('eca_cwd', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_r10mm = import_obs('eca_r10mm', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
lat, lon, rcm_r10mm = import_rcm('eca_r10mm', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r10mm = import_gcm('eca_r10mm', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_r20mm = import_obs('eca_r20mm', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
lat, lon, rcm_r20mm = import_rcm('eca_r20mm', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_r20mm = import_gcm('eca_r20mm', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# Plot maps with the function
fig = plt.figure(figsize=(4, 5.5))
levs1 = [10, 20, 40, 60, 80, 100, 120, 150, 175, 200]
levs2 = [5, 10, 20, 30, 40, 50, 60, 70, 80, 100]
levs3 = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80]
levs4 = [0, 1, 5, 10, 15, 20, 25, 30, 35, 40]

ax = fig.add_subplot(4, 3, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'A) CDD Obs (days)', loc='left', fontsize=6, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_cdd, levels=levs1, latlon=True, cmap=cm.YlGnBu)
	
ax = fig.add_subplot(4, 3, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'B) CDD Reg (days)', loc='left', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_cdd, levels=levs1, latlon=True, cmap=cm.YlGnBu)

ax = fig.add_subplot(4, 3, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'C) CDD Had (days)', loc='left', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_cdd, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(4, 3, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'D) CWD Obs (days)', loc='left', fontsize=6, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_cwd, levels=levs2, latlon=True, cmap=cm.YlGnBu)

ax = fig.add_subplot(4, 3, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'E) CWD Reg (days)', loc='left', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_cwd, levels=levs2, latlon=True, cmap=cm.YlGnBu) 

ax = fig.add_subplot(4, 3, 6)
map, xx, yy = basemap(lat, lon)
plt.title(u'F) CWD Had (days)', loc='left', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_cwd, levels=levs2, latlon=True, cmap=cm.YlGnBu)
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(4, 3, 7)
map, xx, yy = basemap(lat, lon)
plt.title(u'G) R10mm Obs (days)', loc='left', fontsize=6, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_r10mm, levels=levs3, latlon=True, cmap=cm.YlGnBu) 

ax = fig.add_subplot(4, 3, 8)
map, xx, yy = basemap(lat, lon)
plt.title(u'H) R10mm Reg (days)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_r10mm, levels=levs3, latlon=True, cmap=cm.YlGnBu)

ax = fig.add_subplot(4, 3, 9)
map, xx, yy = basemap(lat, lon)
plt.title(u'I) R10mm Had (days)', loc='left', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_r10mm, levels=levs3, latlon=True, cmap=cm.YlGnBu)
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(4, 3, 10) 
map, xx, yy = basemap(lat, lon)
plt.title(u'J) R20mm Obs (days)', loc='left', fontsize=6, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_r20mm, levels=levs4, latlon=True, cmap=cm.YlGnBu)

ax = fig.add_subplot(4, 3, 11)
map, xx, yy = basemap(lat, lon)
plt.title(u'K) R20mm Reg (days)', loc='left', fontsize=6, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_r20mm, levels=levs4, latlon=True, cmap=cm.YlGnBu) 

ax = fig.add_subplot(4, 3, 12)
map, xx, yy = basemap(lat, lon)
plt.title(u'L) R20mm Had (days)', loc='left', fontsize=6, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_r20mm, levels=levs4, latlon=True, cmap=cm.YlGnBu) 
cbar = map.colorbar(ticks=levs4, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

plt.subplots_adjust(left=0.10, bottom=0.10, right=0.99, top=0.99, wspace=0.30, hspace=0.30)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_etccdi_pre_count_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()	
	
	
	
	
	
