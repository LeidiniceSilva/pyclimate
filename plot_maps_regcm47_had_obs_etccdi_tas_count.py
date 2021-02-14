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

	dict_var = {u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}

	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs  = np.nanmean(var[:][:,:,:], axis=0)
	
	return lat, lon, obs
	
	
def import_rcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/rcm/eca'	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	
	
	dict_var = {u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, rcm


def import_gcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/gcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}
	
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
lat, lon, obs_su = import_obs('eca_su', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
lat, lon, rcm_su = import_rcm('eca_su', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_su = import_gcm('eca_su', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tr = import_obs('eca_tr', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
lat, lon, rcm_tr = import_rcm('eca_tr', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tr = import_gcm('eca_tr', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tx10p = import_obs('eca_tx10p', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
lat, lon, rcm_tx10p = import_rcm('eca_tx10p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tx10p = import_gcm('eca_tx10p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tx90p = import_obs('eca_tx90p', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
lat, lon, rcm_tx90p = import_rcm('eca_tx90p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tx90p = import_gcm('eca_tx90p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tn10p = import_obs('eca_tn10p', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
lat, lon, rcm_tn10p = import_rcm('eca_tn10p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tn10p = import_gcm('eca_tn10p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

lat, lon, obs_tn90p = import_obs('eca_tn90p', 'amz_neb', 'xavier_obs', 'yr', '1986-2005')   
lat, lon, rcm_tn90p = import_rcm('eca_tn90p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
lat, lon, gcm_tn90p = import_gcm('eca_tn90p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

# Plot maps with the function
fig = plt.figure(figsize=(4, 8))
levs1 = [338, 341, 344, 347, 350, 353, 356, 359, 362, 365]
levs2 = [65, 95, 125, 155, 185, 215, 245, 275, 305, 335]
levs3 = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]

ax = fig.add_subplot(6, 3, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'A) SU Obs (days)', loc='left', fontsize=6, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_su, levels=levs1, latlon=True, cmap=cm.YlOrRd)
	
ax = fig.add_subplot(6, 3, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'B) SU Reg (days)', loc='left', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_su[0,:,:], levels=levs1, latlon=True, cmap=cm.YlOrRd)

ax = fig.add_subplot(6, 3, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'C) SU Had (days)', loc='left', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_su, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(6, 3, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'D) TR Obs (days)', loc='left', fontsize=6, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_tr, levels=levs2, latlon=True, cmap=cm.YlOrRd)

ax = fig.add_subplot(6, 3, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'E) TR Reg (days)', loc='left', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_tr[0,:,:], levels=levs2, latlon=True, cmap=cm.YlOrRd) 

ax = fig.add_subplot(6, 3, 6)
map, xx, yy = basemap(lat, lon)
plt.title(u'F) TR Had (mm)', loc='left', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_tr, levels=levs2, latlon=True, cmap=cm.YlOrRd)
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(6, 3, 7)
map, xx, yy = basemap(lat, lon)
plt.title(u'G) TX10p Obs (%)', loc='left', fontsize=6, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_tx10p, levels=levs3, latlon=True, cmap=cm.YlOrRd) 

ax = fig.add_subplot(6, 3, 8)
map, xx, yy = basemap(lat, lon)
plt.title(u'H) TX10p Reg (%)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_tx10p[0,:,:], levels=levs3, latlon=True, cmap=cm.YlOrRd)

ax = fig.add_subplot(6, 3, 9)
map, xx, yy = basemap(lat, lon)
plt.title(u'I) TX10p Had (%)', loc='left', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_tx10p, levels=levs3, latlon=True, cmap=cm.YlOrRd)
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(6, 3, 10) 
map, xx, yy = basemap(lat, lon)
plt.title(u'J) TX90p Obs (%)', loc='left', fontsize=6, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_tx90p, levels=levs3, latlon=True, cmap=cm.YlOrRd)

ax = fig.add_subplot(6, 3, 11)
map, xx, yy = basemap(lat, lon)
plt.title(u'K) TX90p Reg (%)', loc='left', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_tx90p[0,:,:], levels=levs3, latlon=True, cmap=cm.YlOrRd) 

ax = fig.add_subplot(6, 3, 12)
map, xx, yy = basemap(lat, lon)
plt.title(u'L) TX90p Had (%)', loc='left', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_tx90p, levels=levs3, latlon=True, cmap=cm.YlOrRd) 
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(6, 3, 13)
map, xx, yy = basemap(lat, lon)
plt.title(u'M) TN10p Obs (%)', loc='left', fontsize=6, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_tn10p, levels=levs3, latlon=True, cmap=cm.YlOrRd) 

ax = fig.add_subplot(6, 3, 14)
map, xx, yy = basemap(lat, lon)
plt.title(u'N) TN10p Reg (%)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_tn10p[0,:,:], levels=levs3, latlon=True, cmap=cm.YlOrRd)

ax = fig.add_subplot(6, 3, 15)
map, xx, yy = basemap(lat, lon)
plt.title(u'O) TN10p Had (%)', loc='left', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_tn10p, levels=levs3, latlon=True, cmap=cm.YlOrRd)
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(6, 3, 16) 
map, xx, yy = basemap(lat, lon)
plt.title(u'P) TN90p Obs (%)', loc='left', fontsize=6, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
plt.ylabel(u'Latitude', fontsize=6, fontweight='bold', labelpad=15)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_tn90p, levels=levs3, latlon=True, cmap=cm.YlOrRd)

ax = fig.add_subplot(6, 3, 17)
map, xx, yy = basemap(lat, lon)
plt.title(u'Q) TN90p Reg (%)', loc='left', fontsize=6, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_tn90p[0,:,:], levels=levs3, latlon=True, cmap=cm.YlOrRd) 

ax = fig.add_subplot(6, 3, 18)
map, xx, yy = basemap(lat, lon)
plt.title(u'R) TN90p Had (%)', loc='left', fontsize=6, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, fontweight='bold', labelpad=10)
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_tn90p, levels=levs3, latlon=True, cmap=cm.YlOrRd) 
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

plt.subplots_adjust(left=0.10, bottom=0.10, right=0.99, top=0.99, wspace=0.30, hspace=0.30)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_etccdi_tas_count_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()	
	
	
	
	
	
