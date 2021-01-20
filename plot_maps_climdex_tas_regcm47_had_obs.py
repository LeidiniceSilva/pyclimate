# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/01/2021"
__description__ = "This script plot maps from extremes index"

import os
import conda
import numpy as np
import matplotlib.pyplot as plt
import warnings ; warnings.filterwarnings("ignore")

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from matplotlib.path import Path
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import PathPatch
from matplotlib.colors import ListedColormap, BoundaryNorm


def import_obs(var, area, freq,  model, exp, dt):
	
	path = '/home/nice/Documents/hadgem2-es_rclimdex/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_r1i1p1_{6}.nc'.format(path, var, area, freq,  model, exp, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs  = np.nanmean(var[:][:,:,:], axis=0)
	
	return lat, lon, obs
	
	
def import_rcm(var, area, freq,  model, exp, dt):
	
	path = '/home/nice/Documents/hadgem2-es_rclimdex/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_r1i1p1_{6}.nc'.format(path, var, area, freq,  model, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, rcm


def import_gcm(var, area, freq,  model, exp, dt):
	
	path = '/home/nice/Documents/hadgem2-es_rclimdex/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_r1i1p1_{6}.nc'.format(path, var, area, freq,  model, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
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
	
	#~ xin = np.linspace(map.xmin,map.xmax,20) 
	#~ yin = np.linspace(map.ymin,map.ymax,20) 
	#~ lons = np.arange(-85.,-5.,0.25) 
	#~ lats = np.arange(-20.,15.,-0.25) 
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
lat, lon, obs_prcptot = import_obs('prcptotETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')   
lat, lon, rcm_prcptot = import_rcm('prcptotETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
lat, lon, gcm_prcptot = import_gcm('prcptotETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

lat, lon, obs_r95p = import_obs('r95pETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')   
lat, lon, rcm_r95p = import_rcm('r95pETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
lat, lon, gcm_r95p = import_gcm('r95pETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

lat, lon, obs_r99p = import_obs('r99pETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')   
lat, lon, rcm_r99p = import_rcm('r99pETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
lat, lon, gcm_r99p = import_gcm('r99pETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

lat, lon, obs_rx1day = import_obs('rx1dayETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')   
lat, lon, rcm_rx1day = import_rcm('rx1dayETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
lat, lon, gcm_rx1day = import_gcm('rx1dayETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

lat, lon, obs_rx5day = import_obs('rx5dayETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')   
lat, lon, rcm_rx5day = import_rcm('rx5dayETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
lat, lon, gcm_rx5day = import_gcm('rx5dayETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

lat, lon, obs_sdii = import_obs('sdiiETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')   
lat, lon, rcm_sdii = import_rcm('sdiiETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
lat, lon, gcm_sdii = import_gcm('sdiiETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

# Plot maps with the function
fig = plt.figure(figsize=(4, 8))
levs1 = [50, 100, 150, 200, 300, 400, 500, 1000, 1500, 2000]
levs2 = [10, 20, 30, 40, 50, 60, 80, 100, 150, 200]
levs3 = [1, 2, 4, 6, 8, 10, 15, 16, 18, 20]

ax = fig.add_subplot(6, 3, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'A) prcptot CRU (mm)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_prcptot, levels=levs1, latlon=True, cmap=cm.YlGnBu)
	
ax = fig.add_subplot(6, 3, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'B) prcptot Reg (mm)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_prcptot, levels=levs1, latlon=True, cmap=cm.YlGnBu)

ax = fig.add_subplot(6, 3, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'C) prcptot Had (mm)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_prcptot, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(6, 3, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'D) r95p CRU (mm)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_r95p, levels=levs1, latlon=True, cmap=cm.YlGnBu)

ax = fig.add_subplot(6, 3, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'E) r95p Reg (mm)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_r95p, levels=levs1, latlon=True, cmap=cm.YlGnBu) 

ax = fig.add_subplot(6, 3, 6)
map, xx, yy = basemap(lat, lon)
plt.title(u'F) r95p Had (mm)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_r95p, levels=levs1, latlon=True, cmap=cm.YlGnBu)
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(6, 3, 7)
map, xx, yy = basemap(lat, lon)
plt.title(u'G) r99p CRU (mm)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_r99p, levels=levs2, latlon=True, cmap=cm.YlGnBu) 

ax = fig.add_subplot(6, 3, 8)
map, xx, yy = basemap(lat, lon)
plt.title(u'H) r99p Reg (mm)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_r99p, levels=levs2, latlon=True, cmap=cm.YlGnBu)

ax = fig.add_subplot(6, 3, 9)
map, xx, yy = basemap(lat, lon)
plt.title(u'I) r99p Had (mm)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_r99p, levels=levs2, latlon=True, cmap=cm.YlGnBu)
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(6, 3, 10) 
map, xx, yy = basemap(lat, lon)
plt.title(u'J) rx1day CRU (mm d⁻¹)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_rx1day, levels=levs2, latlon=True, cmap=cm.YlGnBu)

ax = fig.add_subplot(6, 3, 11)
map, xx, yy = basemap(lat, lon)
plt.title(u'K) rx1day Reg (mm d⁻¹)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_rx1day, levels=levs2, latlon=True, cmap=cm.YlGnBu) 

ax = fig.add_subplot(6, 3, 12)
map, xx, yy = basemap(lat, lon)
plt.title(u'L) rx1day Had (mm d⁻¹)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_rx1day, levels=levs2, latlon=True, cmap=cm.YlGnBu) 
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(6, 3, 13) 
map, xx, yy = basemap(lat, lon)
plt.title(u'J) rx5day CRU (mm d⁻¹)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_rx5day, levels=levs2, latlon=True, cmap=cm.YlGnBu)

ax = fig.add_subplot(6, 3, 14)
map, xx, yy = basemap(lat, lon)
plt.title(u'K) rx5day Reg (mm d⁻¹)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_rx5day, levels=levs2, latlon=True, cmap=cm.YlGnBu) 

ax = fig.add_subplot(6, 3, 15)
map, xx, yy = basemap(lat, lon)
plt.title(u'L) rx5day Had (mm d⁻¹)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_rx5day, levels=levs2, latlon=True, cmap=cm.YlGnBu) 
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(6, 3, 16) 
map, xx, yy = basemap(lat, lon)
plt.title(u'M) sdii CRU (mm d⁻¹)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, obs_sdii, levels=levs3, latlon=True, cmap=cm.YlGnBu)

ax = fig.add_subplot(6, 3, 17)
map, xx, yy = basemap(lat, lon)
plt.title(u'N) sdii Reg (mm d⁻¹)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, rcm_sdii, levels=levs3, latlon=True, cmap=cm.YlGnBu) 

ax = fig.add_subplot(6, 3, 18)
map, xx, yy = basemap(lat, lon)
plt.title(u'O) sdii Had (mm d⁻¹)', fontsize=6, fontweight='bold')
plt.text(-51, -5, u'\u25B2 \nN ', fontsize=6)
map.contourf(xx, yy, gcm_sdii, levels=levs3, latlon=True, cmap=cm.YlGnBu) 
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

#~ fig.tight_layout()
plt.subplots_adjust(left=0.10, bottom=0.10, right=0.99, top=0.99, wspace=0.30, hspace=0.30)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_climdex_pre_reg_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()	
	
	
	
	
	
