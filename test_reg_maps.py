# -*- coding: utf-8 -*-

import os
import conda
import netCDF4
import numpy as np
import matplotlib as mpl  
import matplotlib.pyplot as plt

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset


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
	
	bmap = Basemap(projection='cyl', llcrnrlat=-20, urcrnrlat=10, llcrnrlon=-85, urcrnrlon=-15, resolution='c')
	bmap.drawparallels(np.arange(-20.,15.,5.), labels=[1, 0, 0, 0], color='gray', linewidth=0.5, dashes=[6, 6], fontsize=8, zorder=1)
	bmap.drawmeridians(np.arange(-85.,-5.,10.), labels=[0, 0, 0, 1], color='gray', linewidth=0.5, dashes=[6, 6], fontsize=8, zorder=1)
	
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
		
	lons, lats = np.meshgrid(new_lon, new_lat)
	x, y = bmap(lons, lats)
	
	bmap.readshapefile('/home/nice/Documents/github_projects/shp/shp_world/world', 'world', drawbounds=True, color='gray', linewidth=1.)
	bmap.readshapefile('/home/nice/Documents/github_projects/shp/lim_unid_fed/lim_unid_fed', 'lim_unid_fed', drawbounds=True, color='black', linewidth=1.)

	return bmap, x, y
	
	
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

# Plot map RegCM and HadGEM2-ES
fig = plt.figure(figsize=(10, 8))
levs = [-8, -6, -4, -2, -1, 1, 2, 4, 6, 8]

# Plot firt maps reg_exp1 model 
ax = fig.add_subplot(221)
bmap, x, y = basemap(lat, lon)
plt_maps_bias = bmap.contourf(x, y, bias_exp1_djf[0,:,:], levels=levs, latlon=True, cmap='BrBG')

# Plot second maps reg_exp1 model 
ax = fig.add_subplot(222)
bmap, x, y = basemap(lat, lon)
plt_maps_bias = bmap.contourf(x, y, bias_exp2_djf[0,:,:], levels=levs, latlon=True, cmap='BrBG')
cbar = bmap.colorbar(ticks=levs, drawedges=True, ax=ax, extend='both', shrink=0.8)
cbar.ax.tick_params(labelsize=6) 
cbar.set_label(u'Precipitation (mm d⁻¹)', rotation=90, fontsize=8, fontweight='bold')
cbar.ax.yaxis.set_label_position('right')

# Plot thirth maps reg_exp1 model 
ax = fig.add_subplot(223)
bmap, x, y = basemap(lat, lon)
plt_maps_bias = bmap.contourf(x, y, bias_exp1_jja[0,:,:], levels=levs, latlon=True, cmap='BrBG')

# Plot fourth maps reg_exp1 model 
ax = fig.add_subplot(224)
bmap, x, y = basemap(lat, lon)
plt_maps_bias = bmap.contourf(x, y, bias_exp2_jja[0,:,:], levels=levs, latlon=True, cmap='BrBG')
cbar = bmap.colorbar(ticks=levs, drawedges=True, ax=ax, extend='both', shrink=0.8)
cbar.ax.tick_params(labelsize=6) 
cbar.set_label(u'Precipitation (mm d⁻¹)', rotation=90, fontsize=8, fontweight='bold')
cbar.ax.yaxis.set_label_position('right')

plt.show()
exit()





