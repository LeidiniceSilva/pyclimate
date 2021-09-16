# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/26/2018"
__description__ = "This script plot vertically integrated moisture flux"

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


def import_data(var):
	
	path  = '/home/nice/Downloads'
	arq   = '{0}/{1}_lonlat.nc'.format(path, var)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[var][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][0,:,:]

	return lat, lon, value
	

def import_data_lev1(var):
	
	path  = '/home/nice/Downloads'
	arq   = '{0}/{1}_lonlat.nc'.format(path, var)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[var][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	lev   = data.variables['kz'][:]
	value = var[:][0,17,:,:]

	return lat, lon, value


def import_data_lev2(var):
	
	path  = '/home/nice/Downloads'
	arq   = '{0}/{1}_lonlat.nc'.format(path, var)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[var][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	lev   = data.variables['kz'][:]
	value = var[:][0,4,:,:]

	return lat, lon, value
		

def comp_vimf(ua, va, hus):
	
	p1 = ua * hus
	p2 = va * hus
	p3 = p1 + p2
	
	return p3
	
		
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
	
	xin = np.linspace(map.xmin,map.xmax,10) 
	yin = np.linspace(map.ymin,map.ymax,10) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	lons, lats = np.meshgrid(new_lon, new_lat)

	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='black', linewidth=1.)
	
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
print('Import data')	
# Import regcm exps and obs database 
lat, lon, ps = import_data(u'ps')
#~ lat, lon, ua_lev1 = import_data_lev1(u'ua')
#~ lat, lon, va_lev1 = import_data_lev1(u'va')
#~ lat, lon, hus_lev1 = import_data_lev1(u'hus')

#~ lat, lon, ua_lev2 = import_data_lev2(u'ua')
#~ lat, lon, va_lev2 = import_data_lev2(u'va')
#~ lat, lon, hus_lev2 = import_data_lev2(u'hus')

#~ print('Compute VIMF')	
#~ # Compute VIMF
#~ vimf_lev1 = comp_vimf(ua_lev1, va_lev1, hus_lev1)
#~ vimf_lev2 = comp_vimf(ua_lev2, va_lev2, hus_lev2)
#~ vimf = 75 * (vimf_lev1 + vimf_lev2)
#~ ua_vimf = 75 * (vimf_lev1)
#~ va_vimf = 75 * (vimf_lev1)

#~ print('Compute wind speed')
#~ # Compute wind speed
#~ windspeed = (ua_lev1 ** 2 + va_lev1 ** 2) ** 0.5

print('Plot maps')	 
# Plot maps                                      
fig = plt.figure()
levels = np.arange(575.1194, 1017.6015, 100)

ax = fig.add_subplot(111)
map, xx, yy = basemap(lat, lon) 
ps = map.contour(xx, yy, ps, levels=levels, linewidths=1, colors='gray')
plt.clabel(ps, levels[1::2], inline=1, fmt='%1.1f', fontsize=12)
plt.title(u'A) Surface pressure (hPa)', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=25, fontsize=8, fontweight='bold')
map.drawmeridians(np.arange(-85.,-5.,10.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-20.,15.,5.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')

print('Save figure')	
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_sp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()


           
#~ print('Plot maps')	 
#~ # Plot maps                                      
#~ fig = plt.figure()

#~ ax = fig.add_subplot(111)
#~ map, xx, yy = basemap(lat, lon) 
#~ map.streamplot(xx, yy, ua_lev1, va_lev1, density=0.5, color='black', linewidths=0.5)
#~ plt.title(u'A) Streamline of Wind (m s⁻¹) 1000hPa', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', labelpad=15, fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', labelpad=25, fontsize=8, fontweight='bold')
#~ map.drawmeridians(np.arange(-85.,-5.,10.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
#~ map.drawparallels(np.arange(-20.,15.,5.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')

#~ print('Save figure')	
#~ # Path out to save figure
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_streamline_uv.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
#~ plt.show()
#~ exit()

#~ print('Plot maps')	 
#~ # Plot maps                          
#~ fig = plt.figure()
#~ levs = [0, 1, 2, 3, 4, 5, 6, 7, 8]

#~ ax = fig.add_subplot(111)
#~ map, xx, yy = basemap(lat, lon) 
#~ map.contourf(xx, yy, windspeed, levels=levs, latlon=True, cmap=cm.Greens, extend='max')
#~ q=map.quiver(xx[::8,::8], yy[::8,::8], ua_lev1[::8,::8], va_lev1[::8,::8]) 
#~ plt.title(u'A) Wind (m s⁻¹) 1000hPa', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', labelpad=15, fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', labelpad=25, fontsize=8, fontweight='bold')
#~ cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=8) 
#~ map.drawmeridians(np.arange(-85.,-5.,10.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
#~ map.drawparallels(np.arange(-20.,15.,5.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')
#~ qk = plt.quiverkey(q, 0.85, -0.15, 20, '8 m s⁻¹', labelpos='N')

#~ print('Save figure')	
#~ # Path out to save figure
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_uv.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
#~ plt.show()
#~ exit()

#~ print('Plot maps')	 
#~ # Plot maps          
#~ fig = plt.figure()
#~ levs = [0, 5, 10, 15, 20, 25, 30]

#~ ax = fig.add_subplot(111)
#~ map, xx, yy = basemap(lat, lon) 
#~ map.contourf(xx, yy, hus_lev1, levels=levs, latlon=True, cmap=cm.Greens, extend='max')
#~ map.quiver(xx[::8,::8], yy[::8,::8], ua_lev1[::8,::8], va_lev1[::8,::8]) 
#~ plt.title(u'A) Specific Humidity (g kg⁻¹) Wind (––> 8 m s⁻¹) 1000hPa', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', labelpad=15, fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', labelpad=25, fontsize=8, fontweight='bold')
#~ cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=8) 
#~ map.drawmeridians(np.arange(-85.,-5.,10.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
#~ map.drawparallels(np.arange(-20.,15.,5.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')

#~ print('Save figure')	
#~ # Path out to save figure
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_sh_uv.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
#~ plt.show()
#~ exit()

#~ print('Plot maps')	
#~ # Plot maps 
#~ fig = plt.figure()
#~ levs = [-16, -12, -8, -4, 4, 8, 12, 16]

#~ ax = fig.add_subplot(111)
#~ map, xx, yy = basemap(lat, lon)
#~ map.contourf(xx, yy, vimf, levels=levs, latlon=True, cmap=cm.PiYG, extend='both')
#~ map.quiver(xx[::10,::10], yy[::10,::10], ua_vimf[::10,::10], va_vimf[::10,::10]) 
#~ plt.title(u'A) Vertically Integrated Moisture Flux (1000hPa – 250hPa)', loc='left', fontsize=8, fontweight='bold')
#~ plt.xlabel(u'Longitude', labelpad=15, fontsize=8, fontweight='bold')
#~ plt.ylabel(u'Latitude', labelpad=25, fontsize=8, fontweight='bold')
#~ cbar = map.colorbar(ticks=levs, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=8) 
#~ map.drawmeridians(np.arange(-85.,-5.,10.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
#~ map.drawparallels(np.arange(-20.,15.,5.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')

#~ print('Save figure')	
#~ # Path out to save figure
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_wimf.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
#~ plt.show()
#~ exit()
