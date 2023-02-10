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

from matplotlib import colors as c
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
from os.path import expanduser


def import_obs(area, obs, time):
	
	param = 'tcc' # tcc or msshf
	date  = '2001-2005'

	path  = '/home/nice/Documentos/dataset/obs/reg_pbl'
	arq   = '{0}/{1}_{2}_{3}_{4}_{5}.nc'.format(path, param, area, obs, time, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean_obs = np.nanmean(value, axis=0)

	return lat, lon, mean_obs
	
	
def import_sim(area, exp, time):
	
	param = 'clt' # clt or hfss
	date  = '2001-2005'

	path  = '/home/nice/Documentos/dataset/rcm/reg_pbl'
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

	path = '/home/nice/Documentos/github_projects/shp'
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=.8)

	xx, yy = map(lons,lats)
	
	return map, xx, yy


def plot_maps_clim(djf_obs, bias_djf_exp1, bias_djf_exp2, jja_obs, bias_jja_exp1, bias_jja_exp2, ann_obs, bias_ann_exp1, bias_ann_exp2):
		
	fig = plt.figure(figsize=(8,4))
	
	var = 'tcc' # tcc or msshf 
	
	if var == 'tcc':
		levs1 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
		color1 = cm.Greys
		text1 = u'Cloud area fraction (%)'
		levs2 = [-50, -40, -30, -20, -10, 10, 20, 30, 40, 50]
		color2 = cm.RdGy
		text2 = u'Bias of cloud area fraction (%)'
	else:
		levs1 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]
		color1 = cm.turbo
		text1 = u'Sensible heat flux (W m⁻²)'
		levs2 = [-50, -40, -30, -20, -10, 10, 20, 30, 40, 50]
		color2 = cm.Spectral_r
		text2 = u'Bias of sensible heat flux (W m⁻²)'
		
	ax = fig.add_subplot(331)
	plt.title(u'A) ERA5 DJF', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plot_maps_clim = map.contourf(xx, yy, djf_obs, levels=levs1, latlon=True, cmap=color1) 
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')

	cb_ax = fig.add_axes([0.92, 0.2, 0.016, 0.6])
	bounds = levs1
	cbar = fig.colorbar(plot_maps_clim, cax=cb_ax, orientation='vertical', boundaries=bounds, shrink=0.5, pad=0.5)
	cbar.set_label('{0}'.format(text1), fontsize=8, fontweight='bold')
	cbar.ax.tick_params(labelsize=8)  
	
	ax = fig.add_subplot(332)
	plt.title(u'B) Reg_Holtslag - ERA5 DJF', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	map.contourf(xx, yy, bias_djf_exp1, levels=levs2, latlon=True, cmap=color2)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')

	ax = fig.add_subplot(333)
	plt.title(u'C) Reg_UW-PBL - ERA5 DJF', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	map.contourf(xx, yy, bias_djf_exp2, levels=levs2, latlon=True, cmap=color2)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
		
	ax = fig.add_subplot(334)
	plt.title(u'D) ERA5 JJA', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, jja_obs, levels=levs1, latlon=True, cmap=color1) 
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')
		
	ax = fig.add_subplot(335)
	plt.title(u'E) Reg_Holtslag - ERA5 JJA', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plot_maps_clim = map.contourf(xx, yy, bias_jja_exp1, levels=levs2, latlon=True, cmap=color2)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')

	ax = fig.add_subplot(336)
	plt.title(u'F) Reg_UW-PBL - ERA5 JJA', loc='left', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plot_maps_clim = map.contourf(xx, yy, bias_jja_exp2, levels=levs2, latlon=True, cmap=color2)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
		
	ax = fig.add_subplot(337)
	plt.title(u'G) ERA5 ANN', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	map.contourf(xx, yy, ann_obs, levels=levs1, latlon=True, cmap=color1)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')

	ax = fig.add_subplot(338)
	plt.title(u'H) Reg_Holtslag - ERA5 ANN', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_ann_exp1, levels=levs2, latlon=True, cmap=color2, extend='both')
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')

	cb_ax = fig.add_axes([1, 0.2, 0.016, 0.6])
	bounds = levs2
	cbar = fig.colorbar(plt_maps_bias, cax=cb_ax, orientation='vertical', boundaries=bounds, shrink=0.5, pad=0.5)
	cbar.set_label('{0}'.format(text2), fontsize=8, fontweight='bold')
	cbar.ax.tick_params(labelsize=8)  
	
	ax = fig.add_subplot(339)
	plt.title(u'I) Reg_UW-PBL - ERA5 ANN', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	map.contourf(xx, yy, bias_ann_exp2, levels=levs2, latlon=True, cmap=color2)
	map.drawmeridians(np.arange(-85.,-5.,20.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-20.,15.,10.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')

	return plt_maps_clim


# Import regcm exps and obs database 
lat, lon, djf_obs = import_obs(u'amz_neb', u'era5', 'djf')
lat, lon, djf_exp1 = import_sim(u'amz_neb', u'regcm_exp1', 'djf')
lat, lon, djf_exp2 = import_sim(u'amz_neb', u'regcm_exp2', 'djf')

lat, lon, jja_obs = import_obs(u'amz_neb', u'era5', 'jja')
lat, lon, jja_exp1 = import_sim(u'amz_neb', u'regcm_exp1', 'jja')
lat, lon, jja_exp2 = import_sim(u'amz_neb', u'regcm_exp2', 'jja')

lat, lon, ann_obs = import_obs(u'amz_neb', u'era5', 'ann')
lat, lon, ann_exp1 = import_sim(u'amz_neb', u'regcm_exp1', 'ann')
lat, lon, ann_exp2 = import_sim(u'amz_neb', u'regcm_exp2', 'ann')

var = 'tcc' # tcc or msshf 

if var == 'tcc':
	djf_obs = djf_obs*100
	djf_exp1 = djf_exp1*100
	djf_exp2 = djf_exp2*100
	jja_obs = jja_obs*100
	jja_exp1 = jja_exp1*100
	jja_exp2 = jja_exp2*100
	ann_obs = ann_obs*100
	ann_exp1 = ann_exp1*100
	ann_exp2 = ann_exp2*100
else:
	djf_obs = djf_obs*(-1)
	djf_exp1 = djf_exp1
	djf_exp2 = djf_exp2
	jja_obs = jja_obs*(-1)
	jja_exp1 = jja_exp1
	jja_exp2 = jja_exp2
	ann_obs = ann_obs*(-1)
	ann_exp1 = ann_exp1
	ann_exp2 = ann_exp2

bias_djf_exp1 = djf_exp1 - djf_obs
bias_djf_exp2 = djf_exp2 - djf_obs

bias_jja_exp1 = jja_exp1 - jja_obs
bias_jja_exp2 = jja_exp2 - jja_obs

bias_ann_exp1 = ann_exp1 - ann_obs
bias_ann_exp2 = ann_exp2 - ann_obs

# Plot maps with the function
plt_map = plot_maps_clim(djf_obs, bias_djf_exp1, bias_djf_exp2, jja_obs, bias_jja_exp1, bias_jja_exp2, ann_obs, bias_ann_exp1, bias_ann_exp2)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_bias_{0}_regcm_pbl_obs_2001-2005.png'.format(var)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
