# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/21/2021"
__description__ = "This script plot vertical across section"

import os
import conda
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.patches as mpatches

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from netCDF4 import Dataset as nc
from matplotlib.colors import Normalize
from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from wrf import to_np, CoordPair, vertcross, latlon_coords, interpline, interp2dxy, xy


def find_ind_latlon2d(lats, lat_pt, lons, lon_pt):
	
	a = abs(lats - lat_pt) + abs(lons - lon_pt)
	i,j = np.unravel_index(a.argmin(),a.shape)

	return i, j
  
 
def map_RegCMdomain(lat_start, lat_end, lon_start, lon_end, lon0, lat0):

    m = Basemap(ax=ax, llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end, urcrnrlat=lat_end,
                resolution='c', area_thresh=1000., projection='cyl', lon_0=lon0, lat_0=lat0, lat_ts=0)

    m.drawparallels(range(-20, 10, 10), labels=[1,0,0,0], fontsize=8, dashes=[1, 2],linewidth=1, color='black', zorder=12)
    m.drawmeridians(range(-85, -15, 10),labels=[0,0,0,1], fontsize=8, dashes=[1, 2],linewidth=1, color='black', zorder=12)
    
    path = '/home/nice/Documents/github_projects/shp'
    m.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='black', linewidth=.5)
    m.readshapefile('{0}/lim_unid_fed/lim_unid_fed'.format(path), 'lim_unid_fed', drawbounds=True, color='black', linewidth=.5)
    
    return m
	 
	 
# Open and read models and obs database 	         
dirIN  = '/home/nice/Downloads'
infMOD = os.path.join(dirIN,'atm_RegCM4_HadG_historical_yr_1986-2005.nc')

inf_mod   = nc(infMOD, mode='r')
lat       = inf_mod.variables['xlat'][:]
lon       = inf_mod.variables['xlon'][:]
PS  	  = inf_mod.variables['ps'][0,:,:]
plev      = inf_mod.variables['plev'][:]
RH        = inf_mod.variables['rh'][0,:,:,:] 
inf_mod.close()

# Get number of pressure level
nlev = len(plev)

# Latitude and longitude of projection origin 
lat0     = -6
lon0     = -50

# Boundaries of the domain map
lat_start = -20
lat_end   = 10
lon_start = -85
lon_end   = -15

# Start and end coordinates for the cross section
start_lat = -20.0; start_lon = -85.0
end_lat = 10.; end_lon = -15.0   

# Compute the cross section 
# 1. Find the starting and ending point of the transect
start_i, start_j = find_ind_latlon2d(lat, start_lat, lon, start_lon)
end_i, end_j     = find_ind_latlon2d(lat, end_lat, lon, end_lon)
start = (start_j, start_i)
end   = (end_j, end_i)

# 2. Compute the vertical cross-section interpolation
xy   	   = xy(PS, start_point=start, end_point=end)
rh_cross   = interp2dxy(RH, xy, meta=False)

# 3. Interpolate LAT-LON along the cross section and build LAT-LON pairs 
lat_interp   = interpline(lat, start_point=CoordPair(x=start_j, y=start_i), end_point=CoordPair(x=end_j, y=end_i), latlon=True, meta=False)
lon_interp   = interpline(lon, start_point=CoordPair(x=start_j, y=start_i), end_point=CoordPair(x=end_j, y=end_i), latlon=False, meta=False)
latlon_pairs = np.dstack((lat_interp, lon_interp))
ps_interp    = interpline(PS, start_point=CoordPair(x=start_j, y=start_i), end_point=CoordPair(x=end_j, y=end_i), latlon=False, meta=False)

# Plot models and obs database 
fig = plt.figure(1, figsize=(10,8))
clevs = [5., 10., 20., 30., 40., 50., 60., 70., 80., 90.]
norm  = colors.BoundaryNorm(boundaries=clevs, ncolors=256)
cmap  = "BrBG" 

ax = fig.add_subplot(3,2,1)
m = map_RegCMdomain(lat_start, lat_end, lon_start, lon_end, lon0, lat0)
lon_S, lat_S = m(lon, lat) 
plot_RH = m.contourf(lon_S, lat_S, RH[2,:,:], norm=norm, levels=clevs, ax=ax, zorder=1, cmap=cmap, extend='both')
plt.title('A)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel('Latitude', fontsize=8, fontweight='bold', labelpad=25)
plt.xlabel('Longitude', fontsize=8, fontweight='bold', labelpad=15)
cbar = fig.colorbar(plot_RH, ax=ax, orientation='horizontal', aspect=10)

nxy    = np.shape(rh_cross)[1]
x      = np.arange(0,nxy,1)
X,Z    = np.meshgrid(x, plev)
x_axis = np.arange(0,nxy,1)
y1     = np.empty(nxy,'float')
y1.fill(1000.)
y2     = np.minimum(ps_interp[:]/100.,y1)

ax = fig.add_subplot(3,2,2)   
cross_RH = ax.contourf(x_axis, plev, to_np(rh_cross[:,:]), levels=clevs, cmap=cmap, extend='both')               
plt.title('B)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel('Vertical levels (hPa)', fontsize=8, fontweight='bold', labelpad=10)

# Add the orography
YY = np.arange(0, plev.shape[0], 2)
XX = np.arange(0, x_axis.shape[0], 2)
points = np.meshgrid(YY, XX)
xx, yy = np.meshgrid(YY, XX)
x, y   = m(xx, yy) 
ax.fill_between(x_axis, y2, 1000.0, facecolor='black')
plt.gca().invert_yaxis()
plt.colorbar(cross_RH)

# Set the x-ticks and y-ticks 
coord_pairs = latlon_pairs[0,:,:]
x_ticks     = np.arange(coord_pairs.shape[0])
x_labels    = coord_pairs
ax.set_xticks(x_ticks[::5])
ax.set_xticklabels(x_labels[::5], rotation=70, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
vert_vals = plev[::-1]
ax.set_yticks(vert_vals)
ax.set_yticklabels(vert_vals[:], fontsize=8)

ax = fig.add_subplot(3,2,3)
m = map_RegCMdomain(lat_start, lat_end, lon_start, lon_end, lon0, lat0)
lon_S, lat_S = m(lon, lat) 
plot_RH = m.contourf(lon_S, lat_S, RH[2,:,:], norm=norm, levels=clevs, ax=ax, zorder=1, cmap=cmap, extend='both')
plt.title('C)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel('Latitude', fontsize=8, fontweight='bold', labelpad=25)
plt.xlabel('Longitude', fontsize=8, fontweight='bold', labelpad=15)
cbar = fig.colorbar(plot_RH, ax=ax, orientation='horizontal')

nxy    = np.shape(rh_cross)[1]
x      = np.arange(0,nxy,1)
X,Z    = np.meshgrid(x, plev)
x_axis = np.arange(0,nxy,1)
y1     = np.empty(nxy,'float')
y1.fill(1000.)
y2     = np.minimum(ps_interp[:]/100.,y1)

ax = fig.add_subplot(3,2,4)   
cross_RH = ax.contourf(x_axis, plev, to_np(rh_cross[:,:]), levels=clevs, cmap=cmap, extend='both')               
plt.title('D)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel('Vertical levels (hPa)', fontsize=8, fontweight='bold', labelpad=10)

# Add the orography
YY = np.arange(0, plev.shape[0], 2)
XX = np.arange(0, x_axis.shape[0], 2)
points = np.meshgrid(YY, XX)
xx, yy = np.meshgrid(YY, XX)
x, y   = m(xx, yy) 
ax.fill_between(x_axis, y2, 1000.0, facecolor='black')
plt.gca().invert_yaxis()
plt.colorbar(cross_RH)

# Set the x-ticks and y-ticks
coord_pairs = latlon_pairs[0,:,:]
x_ticks     = np.arange(coord_pairs.shape[0])
x_labels    = coord_pairs
ax.set_xticks(x_ticks[::5])
ax.set_xticklabels(x_labels[::5], rotation=70, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
vert_vals = plev[::-1]
ax.set_yticks(vert_vals)
ax.set_yticklabels(vert_vals[:], fontsize=8)

ax = fig.add_subplot(3,2,5)
m = map_RegCMdomain(lat_start, lat_end, lon_start, lon_end, lon0, lat0)
lon_S, lat_S = m(lon, lat) 
plot_RH = m.contourf(lon_S, lat_S, RH[2,:,:], norm=norm, levels=clevs, ax=ax, zorder=1, cmap=cmap, extend='both')
plt.title('E)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel('Latitude', fontsize=8, fontweight='bold', labelpad=25)
plt.xlabel('Longitude', fontsize=8, fontweight='bold', labelpad=15)
cbar = fig.colorbar(plot_RH, ax=ax, orientation='horizontal')

nxy    = np.shape(rh_cross)[1]
x      = np.arange(0,nxy,1)
X,Z    = np.meshgrid(x, plev)
x_axis = np.arange(0,nxy,1)
y1     = np.empty(nxy,'float')
y1.fill(1000.)
y2     = np.minimum(ps_interp[:]/100.,y1)

ax = fig.add_subplot(3,2,6)   
cross_RH = ax.contourf(x_axis, plev, to_np(rh_cross[:,:]), levels=clevs, cmap=cmap, extend='both')               
plt.title('F)', loc='left', fontsize=8, fontweight='bold')
plt.ylabel('Vertical levels (hPa)', fontsize=8, fontweight='bold', labelpad=10)

# Add the orography
YY = np.arange(0, plev.shape[0], 2)
XX = np.arange(0, x_axis.shape[0], 2)
points = np.meshgrid(YY, XX)
xx, yy = np.meshgrid(YY, XX)
x, y   = m(xx, yy) 
ax.fill_between(x_axis, y2, 1000.0, facecolor='black')
plt.gca().invert_yaxis()
plt.colorbar(cross_RH)

# Set the x-ticks and y-ticks
coord_pairs = latlon_pairs[0,:,:]
x_ticks     = np.arange(coord_pairs.shape[0])
x_labels    = coord_pairs
ax.set_xticks(x_ticks[::5])
ax.set_xticklabels(x_labels[::5], rotation=70, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
vert_vals = plev[::-1]
ax.set_yticks(vert_vals)
ax.set_yticklabels(vert_vals[:], fontsize=8)

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_vertical_cross_reg_exp2.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()
