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
	
	path  = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq   = '{0}/{1}_amz_neb_{2}_mon_2001-2010.nc'.format(path, param, exp)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	exp = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, exp


def import_obs(param, obs):
	
	path  = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq   = '{0}/{1}_amz_neb_{2}_mon_2001-2010.nc'.format(path, param, obs)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	obs   = np.nanmean(var[:][:,:,:], axis=0)
	
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
	
	
def colormap():
	
	var_colormap = cm.YlGnBu  # If precipitation
	levs = [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
	
	return var_colormap, levs
	
	
def function_plot(exp1, exp2, obs1):
		
	fig = plt.figure()

	# Plot firt maps reg_exp1 model 
	ax = fig.add_subplot(311)
	plt.title('Precipitação (mm/d) 2001-2010', fontsize=8, fontweight='bold')
	plt.text(-37, -18, u'A) Reg_Exp1', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	var_colormap, levs = colormap()
	plt_map = map.contourf(xx, yy, exp1[:,:], levels=levs, latlon=True, cmap=var_colormap)
	map.colorbar(ticks=levs, drawedges=True)

	# Plot firt maps reg_exp2 model 
	ax = fig.add_subplot(312)
	plt.text(-37, -18, u'B) Reg_Exp2', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	var_colormap, levs = colormap()
	plt_map = map.contourf(xx, yy, exp2[:,:], levels=levs, latlon=True, cmap=var_colormap)
	map.colorbar(ticks=levs, drawedges=True)

	# Plot thirth maps cru obs 
	ax = fig.add_subplot(313)
	plt.text(-37, -18, u'C) GPCP', fontsize=8, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=8, labelpad=16, fontweight='bold')
	map, xx, yy = basemap(lat, lon)
	var_colormap, levs = colormap()
	plt_map = map.contourf(xx, yy, obs1[:,:], levels=levs, latlon=True, cmap=var_colormap)
	map.colorbar(ticks=levs, drawedges=True, ax=ax)
	
	return plt_map

# Import regcm exp and cru databases 	   
lat, lon, exp1 = import_sim('pr', 'regcm_exp1')
lat, lon, exp2 = import_sim('pr', 'regcm_exp2')
lat, lon, obs1  = import_obs('precip', 'gpcp_v2.2_obs')

# Plot maps with the function
plt_map = function_plot(exp1, exp2, obs1)

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/regcm_pbl/results'
name_out = 'pyplt_maps_pr_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()



