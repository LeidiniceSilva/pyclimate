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
	exp = var[:,:,:]

	return lat, lon, exp


def import_obs(param, obs):
	
	path  = '/home/nice/Documents/ufrn/papers/regcm_pbl/datas'
	arq   = '{0}/{1}_amz_neb_{2}_mon_2001-2010.nc'.format(path, param, obs)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	obs = var[:,:,:]
	
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
	new_lat = lat[::-1]
	new_lon = lon[::-1]
	
	map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-20., urcrnrlon=-15.,urcrnrlat=10., resolution='c')
	map.drawmeridians(np.arange(-85.,-5.,10.), labels=[0,0,0,1], linewidth=0.2)
	map.drawparallels(np.arange(-20.,15.,5.), labels=[1,0,0,0], linewidth=0.2)
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
	
	param = 'tas'
	
	if param == 'pr':
		precip_colors = ['#FDFDFD', '#04E9E7', '#019FF4', '#0300F4',
		'#02FD02', '#01C501', '#008E00', '#FDF802', '#E5BC00', '#FD9500',
		'#FD0000', '#D40000', '#BF0000', '#F800FD', '#9854C6']
		#~ precip_colors = ['#FFFFFF', '#78F573','#37D23C','#0FA00F',
		#~ '#96D2FA','#50A5F5','#1464D2','#FFFAAA','#FDD802','#FFC03C',
		#~ '#FF6000','#E11400','#BC0000','#F800FD','#9854C6']
		
		mpl.colors.ListedColormap(precip_colors)
		levs = [0, 0.5, 1, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10., 15.0, 20.0]
		norm = mpl.colors.BoundaryNorm(levs, 16)

	else:
		levs = [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]
		var_colormap = cm.coolwarm
		norm = mpl.colors.BoundaryNorm(levs, 16)	
	
	return var_colormap, levs, norm
	
	
def function_plot(exp1, exp2, obs1):
		
	fig = plt.figure(figsize=(24,34))
	
	# Plot firt maps reg_exp1 model 
	ax = fig.add_subplot(311)
	plt.title('Temperatura 2m Reg_Exp1 ($^\circ$C)', fontsize=30, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=30, labelpad=30, fontweight='bold')

	map, xx, yy = basemap(lat, lon)
	var_colormap, levs, norm = colormap()
	plt_map = map.contourf(xx, yy, exp1[2,:,:], levels=levs, latlon=True, norm=norm, cmap=var_colormap)
	map.colorbar(ticks=levs, drawedges=True)

	# Plot firt maps reg_exp2 model 
	ax = fig.add_subplot(312)
	plt.title(u'Temperatura 2m Reg_Exp2 ($^\circ$C)', fontsize=30, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=30, labelpad=30, fontweight='bold')

	map, xx, yy = basemap(lat, lon)
	precip_colormap, levs, norm = colormap()
	plt_map = map.contourf(xx, yy, exp2[2,:,:], levels=levs, latlon=True, norm=norm, cmap=var_colormap)
	map.colorbar(ticks=levs, drawedges=True)

	# Plot thirth maps cru obs 
	ax = fig.add_subplot(313)
	plt.title(u'Temperatura 2m CRU ($^\circ$C)', fontsize=30, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=30, labelpad=30, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=30, labelpad=30, fontweight='bold')

	map, xx, yy = basemap(lat, lon)
	precip_colormap, levs, norm = colormap()
	plt_map = map.contourf(xx, yy, obs1[2,:,:], levels=levs, latlon=True, norm=norm, cmap=var_colormap)
	map.colorbar(ticks=levs, drawedges=True)
	
	return plt_map


# Import regcm exp and cru databases 	   
lat, lon, exp1 = import_sim('tas', 'regcm_exp1')
lat, lon, exp2 = import_sim('tas', 'regcm_exp2')
lat, lon, obs1  = import_obs('tmp', 'cru_ts4.02_obs')

# Plot maps with the function
plt_map = function_plot(exp1, exp2, obs1)

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/regcm_pbl/results'
name_out = 'pyplt_maps_tas_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()



