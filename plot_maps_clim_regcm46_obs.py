# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot climatology maps from regcm46 and obs database"

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

from matplotlib import colors as c
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
from os.path import expanduser


def import_obs(area, obs, time):
	
	param = 'precip' 
	date  = '2001-2010'

	path  = '/home/nice/Documents/dataset/obs/reg_pbl'
	arq   = '{0}/{1}_{2}_{3}_{4}_{5}.nc'.format(path, param, area, obs, time, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean_obs = np.nanmean(value, axis=0)

	return lat, lon, mean_obs
	
	
def import_sim(area, exp, time):
	
	param = 'pr' 
	date  = '2001-2010'

	path  = '/home/nice/Documents/dataset/rcm/reg_pbl'
	arq   = '{0}/{1}_{2}_{3}_{4}_{5}.nc'.format(path, param, area, exp, time, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean_sim = np.nanmean(value, axis=0)

	return lat, lon, mean_sim
	

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
	
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	lons, lats = np.meshgrid(new_lon, new_lat)

	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=.8)

	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
def plot_maps_clim(djf_obs, djf_exp1, djf_exp2, jja_obs, jja_exp1, jja_exp2, ann_obs, ann_exp1, ann_exp2):
		
	fig = plt.figure(figsize=(8,4))
	levs = [0, 1, 2, 4, 6, 8, 10, 12, 14, 16]

	ax = fig.add_subplot(331)
	plt.title(u'A) GPCP DJF', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, djf_obs, levels=levs, latlon=True, cmap=cm.Blues) 
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')

	ax = fig.add_subplot(332)
	plt.title(u'B) Reg_H DJF', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, djf_exp1, levels=levs, latlon=True, cmap=cm.Blues)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')

	ax = fig.add_subplot(333)
	plt.title(u'C) Reg_UW DJF', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, djf_exp2, levels=levs, latlon=True, cmap=cm.Blues)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
		
	ax = fig.add_subplot(334)
	plt.title(u'D) GPCP JJA', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, jja_obs, levels=levs, latlon=True, cmap=cm.Blues) 
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')
		
	ax = fig.add_subplot(335)
	plt.title(u'E) Reg_H JJA', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, jja_exp1, levels=levs, latlon=True, cmap=cm.Blues)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')

	ax = fig.add_subplot(336)
	plt.title(u'F) Reg_UW JJA', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, jja_exp2, levels=levs, latlon=True, cmap=cm.Blues)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
		
	ax = fig.add_subplot(337)
	plt.title(u'G) GPCP ANN', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, ann_obs, levels=levs, latlon=True, cmap=cm.Blues)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')

	ax = fig.add_subplot(338)
	plt.title(u'H) Reg_H ANN', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, ann_exp1, levels=levs, latlon=True, cmap=cm.Blues)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
			
	ax = fig.add_subplot(339)
	plt.title(u'I) Reg_UW ANN', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, ann_exp2, levels=levs, latlon=True, cmap=cm.Blues, extend='max')
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')

	cb_ax = fig.add_axes([0.92, 0.2, 0.016, 0.6])
	bounds=[0, 1, 2, 4, 6, 8, 10, 12, 14, 16]
	cbar = fig.colorbar(plt_maps_clim, cax=cb_ax, orientation='vertical', boundaries=bounds, shrink=0.5, pad=0.5)
	cbar.set_label(u'Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
	cbar.ax.tick_params(labelsize=8)  

	return plt_maps_clim


# Import regcm exps and obs database 
lat, lon, djf_obs = import_obs(u'amz_neb', u'gpcp_v2.3_obs', 'djf')
lat, lon, djf_exp1 = import_sim(u'amz_neb', u'regcm_exp1', 'djf')
lat, lon, djf_exp2 = import_sim(u'amz_neb', u'regcm_exp2', 'djf')

lat, lon, jja_obs = import_obs(u'amz_neb', u'gpcp_v2.3_obs', 'jja')
lat, lon, jja_exp1 = import_sim(u'amz_neb', u'regcm_exp1', 'jja')
lat, lon, jja_exp2 = import_sim(u'amz_neb', u'regcm_exp2', 'jja')

lat, lon, ann_obs = import_obs(u'amz_neb', u'gpcp_v2.3_obs', 'ann')
lat, lon, ann_exp1 = import_sim(u'amz_neb', u'regcm_exp1', 'ann')
lat, lon, ann_exp2 = import_sim(u'amz_neb', u'regcm_exp2', 'ann')

# Plot maps with the function
plt_map = plot_maps_clim(djf_obs, djf_exp1, djf_exp2, jja_obs, jja_exp1, jja_exp2, ann_obs, ann_exp1, ann_exp2)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()



