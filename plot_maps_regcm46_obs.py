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


def import_sim(exp, season):
	
	path = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq  = '{0}/pr_amz_neb_{1}_{2}_2001-2010_clim.nc'.format(path, exp, season)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables['pr'][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	reg  = var[:][:,:,:]

	return lat, lon, reg


def import_obs(obs, season):
	
	path = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq  = '{0}/precip_amz_neb_{1}_{2}_2001-2010_clim.nc'.format(path, obs, season)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables['precip'][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gpcp  = var[:][:,:,:]
	
	return lat, lon, gpcp
	

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
	
	
def plot_maps_clim(exp1_djf, exp1_jja, exp2_djf, exp2_jja, obs_djf, obs_jja):
		
	fig = plt.figure()

	levs = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16]

	# Plot firt maps reg_exp1 model 
	ax = fig.add_subplot(321)
	plt.title(u'A) Reg_Exp1 DJF/2001-2010', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, exp1_djf[0,:,:], levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot firt maps reg_exp2 model 
	ax = fig.add_subplot(322)
	plt.title(u'B) Reg_Exp1 JJA/2001-2010', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, exp1_jja[0,:,:], levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot thirth maps cru obs 
	ax = fig.add_subplot(323)
	plt.title(u'C) Reg_Exp2 DJF/2001-2010', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, exp2_djf[0,:,:], levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot firt maps reg_exp1 model 
	ax = fig.add_subplot(324)
	plt.title(u'D) Reg_Exp2 JJA/2001-2010', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_map = map.contourf(xx, yy, exp2_jja[0,:,:], levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot firt maps reg_exp2 model 
	ax = fig.add_subplot(325)
	plt.title(u'E) GPCP DJF/2001-2010', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, obs_djf[0,:,:], levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot thirth maps cru obs 
	ax = fig.add_subplot(326)
	plt.title(u'F) GPCP JJA/2001-2010', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_clim = map.contourf(xx, yy, obs_jja[0,:,:], levels=levs, latlon=True, cmap=cm.YlGnBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	return plt_maps_clim
	
	
def plot_maps_bias(bias_exp1_djf, bias_exp1_jja, bias_exp2_djf, bias_exp2_jja):
		
	fig = plt.figure(figsize=(8,4))

	levs = [-8, -6, -5, -4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4, 5, 6, 8]

	# Plot firt maps reg_exp1 model 
	ax = fig.add_subplot(221)
	plt.title(u'A) Reg_Exp1 - GPCP DJF/2001-2010', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp1_djf[0,:,:], levels=levs, latlon=True, cmap=cm.RdBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot firt maps reg_exp2 model 
	ax = fig.add_subplot(222)
	plt.title(u'B) Reg_Exp1 - GPCP JJA/2001-2010', fontsize=8, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp1_jja[0,:,:], levels=levs, latlon=True, cmap=cm.RdBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 

	# Plot firt maps reg_exp1 model 
	ax = fig.add_subplot(223)
	plt.title(u'C) Reg_Exp2 - GPCP DJF/2001-2010', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp2_djf[0,:,:], levels=levs, latlon=True, cmap=cm.RdBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	# Plot firt maps reg_exp2 model 
	ax = fig.add_subplot(224)
	plt.title(u'D) Reg_Exp2 - GPCP JJA/2001-2010', fontsize=8, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	plt_maps_bias = map.contourf(xx, yy, bias_exp2_jja[0,:,:], levels=levs, latlon=True, cmap=cm.RdBu)
	cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	return plt_maps_bias


# Import regcm exp and cru databases 	   
lat, lon, exp1_djf = import_sim('regcm_exp1', 'djf')
lat, lon, exp1_jja = import_sim('regcm_exp1', 'jja')

lat, lon, exp2_djf = import_sim('regcm_exp1', 'djf')
lat, lon, exp2_jja = import_sim('regcm_exp1', 'jja')

lat, lon, obs_djf = import_obs('gpcp_v2.2_obs', 'djf')
lat, lon, obs_jja = import_obs('gpcp_v2.2_obs', 'jja')

# Compute and plot bias from regcm exp and cru database
bias_exp1_djf = exp1_djf - obs_djf
bias_exp1_jja = exp1_jja - obs_jja

bias_exp2_djf = exp2_djf - obs_djf
bias_exp2_jja = exp2_jja - obs_jja

# Plot maps with the function
#~ plt_map = plot_maps_clim(exp1_djf, exp1_jja, exp2_djf, exp2_jja, obs_djf, obs_jja)
plt_map = plot_maps_bias(bias_exp1_djf, bias_exp1_jja, bias_exp2_djf, bias_exp2_jja)

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/regcm_pbl/results'
name_out = 'pyplt_maps_bias_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()



