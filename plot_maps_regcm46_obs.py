# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/26/2018"
__description__ = "This script plot climatology maps from Rec_EXP models end OBS basedata"

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


def import_sim(param, exp):
	
	path = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq  = '{0}/{1}_amz_neb_{2}_mon_2001-2010.nc'.format(path, param, exp)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	exp  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, exp


def import_obs(param, obs):
	
	path = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq  = '{0}/{1}_amz_neb_{2}_mon_2001-2010.nc'.format(path, param, obs)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs  = np.nanmean(var[:][:,:,:], axis=0)
	
	return lat, lon, obs
	

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
	map.drawmeridians(np.arange(-85.,-5.,10.), size=6, labels=[0,0,0,1], linewidth=0.4)
	map.drawparallels(np.arange(-20.,15.,5.), size=6, labels=[1,0,0,0], linewidth=0.4)
	map.drawcoastlines(linewidth=1, color='k')
	map.drawcountries(linewidth=1, color='k')
	
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	lons, lats = np.meshgrid(new_lon, new_lat)

	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
def plot_maps_clim(exp1, exp2, obs1):
		
	fig = plt.figure()

	levs = [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

	# Plot firt maps reg_exp1 model 
	ax = fig.add_subplot(321)
	plt.title(u'A) Reg_Exp1 DJF/2001-2010', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, exp1[:,:], levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot firt maps reg_exp2 model 
	ax = fig.add_subplot(322)
	plt.title(u'B) Reg_Exp1 JJA/2001-2010', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, exp2[:,:], levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot thirth maps cru obs 
	ax = fig.add_subplot(323)
	plt.title(u'C) Reg_Exp2 DJF/2001-2010', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, obs1[:,:], levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot firt maps reg_exp1 model 
	ax = fig.add_subplot(324)
	plt.title(u'D) Reg_Exp2 JJA/2001-2010', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, exp1[:,:], levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot firt maps reg_exp2 model 
	ax = fig.add_subplot(325)
	plt.title(u'E) GPCP DJF/2001-2010', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, exp2[:,:], levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot thirth maps cru obs 
	ax = fig.add_subplot(326)
	plt.title(u'F) GPCP JJA/2001-2010', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, obs1[:,:], levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	return plt_maps_clim
	
	
def plot_maps_bias(bias_exp1, bias_exp2):
		
	fig = plt.figure(figsize=(8,4))

	levs = [-6, -5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5, 6]

	# Plot firt maps reg_exp1 model 
	ax = fig.add_subplot(221)
	plt.title(u'A) Reg_Exp1 - GPCP DJF/2001-2010', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp1[:,:], levels=levs, latlon=True, cmap=cm.RdBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot firt maps reg_exp2 model 
	ax = fig.add_subplot(222)
	plt.title(u'B) Reg_Exp1 - GPCP JJA/2001-2010', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp2[:,:], levels=levs, latlon=True, cmap=cm.RdBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 

	# Plot firt maps reg_exp1 model 
	ax = fig.add_subplot(223)
	plt.title(u'C) Reg_Exp2 - GPCP DJF/2001-2010', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp1[:,:], levels=levs, latlon=True, cmap=cm.RdBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot firt maps reg_exp2 model 
	ax = fig.add_subplot(224)
	plt.title(u'D) Reg_Exp2 - GPCP JJA/2001-2010', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp2[:,:], levels=levs, latlon=True, cmap=cm.RdBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	return plt_maps_bias


# Import regcm exp and cru databases 	   
lat, lon, exp1 = import_sim('pr', 'regcm_exp1')
lat, lon, exp2 = import_sim('pr', 'regcm_exp2')
lat, lon, obs1 = import_obs('precip', 'gpcp_v2.2_obs')

#~ # Plot maps with the function
plt_map = plot_maps_clim(exp1, exp2, obs1)

# Compute and plot bias from regcm exp and cru database
#~ bias_exp1 = exp1 - obs1
#~ bias_exp2 = exp2 - obs1
#~ plt_map = plot_maps_bias(bias_exp1, bias_exp2)

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/regcm_pbl/results'
name_out = 'pyplt_maps_clim_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()



