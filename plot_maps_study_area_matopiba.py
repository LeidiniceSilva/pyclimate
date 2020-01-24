# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "11/07/2019"
__description__ = "This script plot maps from Rclimdex"

import os
import conda
import numpy as np
import matplotlib as mpl  #; mpl.use('Agg')
import matplotlib.pyplot as plt
import warnings ; warnings.filterwarnings("ignore")

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from matplotlib.path import Path
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import PathPatch

# Import rclimdex database 
obse = 'pr_cmap.nc'
obse_d   = Dataset(obse)
obse_v   = obse_d.variables['pr'][0,...]
lat = obse_d.variables['Y'][:]
lon = obse_d.variables['X'][:]

y1 = -18
y2 =   0
x1 = -53
x2 = -40

# Preparing lat lon
lat = lat - np.mean(np.diff(lat))/2.
lon = lon - np.mean(np.diff(lon))/2.
lons, lats = np.meshgrid(lon, lat)

# Plot map from study area
fig = plt.figure(figsize=(10, 8))

parallels = np.arange( -90,  91, 2)
meridians = np.arange(-180, 180, 2)
bmap = Basemap(projection='cyl', llcrnrlat=y1, urcrnrlat=y2, llcrnrlon=x1, urcrnrlon=x2, resolution=None, suppress_ticks=True, lon_0=0, celestial=False)
bmap.drawparallels(parallels, labels=[1, 0, 0, 0], color='gray', linewidth=0.01, dashes=[4, 6], fontsize=10, zorder=1)
bmap.drawmeridians(meridians, labels=[0, 0, 0, 1], color='gray', linewidth=0.01, dashes=[4, 6], fontsize=10, zorder=1)
x, y = bmap(lons, lats)

# Import shapefile from word and matopiba 
bmap.readshapefile('./world', 'world', drawbounds=True)
bmap.readshapefile('./matopiba', 'matopiba', drawbounds=True, color='black', linewidth=2.5)

# Masking map with shape
x0, x1 = plt.xlim()
y0, y1 = plt.ylim()

map_edges = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
polys = [map_edges]

bmap.readshapefile('./matopiba', 'matopiba2', drawbounds=False)
polys = polys + getattr(bmap, 'matopiba2')
codes = [[Path.MOVETO] + [Path.LINETO for p in p[1:]] for p in polys] # creating a PathPatch

polys_lin = [v for p in polys for v in p]
codes_lin = [cdg for cdgs in codes for cdg in cdgs]
path  = Path(polys_lin, codes_lin)
patch = PathPatch(path, facecolor='white', lw=0)

plt.gca().add_patch(patch) # masking the data

# Bar colors and levels
colors  = ('#ffffff', '#b2b2f2', '#8484f2', '#4363d8', '#4143f2', '#0006f2') # blue
pale_lev= (0, 300, 500, 600, 700, 800)

pa_over = colors[ -1]
mycolor = colors[:-1]
my_cmap = ListedColormap(mycolor)
my_cmap.set_over(pa_over)

norm = BoundaryNorm(pale_lev, ncolors=my_cmap.N, clip=False)

cs0 = plt.pcolormesh(x, y, obse_v*0, cmap=my_cmap, norm=norm)
#~ plt.colorbar(cs0, orientation='vertical', extend='max', spacing='uniform', ticks=pale_lev, extendfrac='auto', cmap=my_cmap)
plt.title('Área de Estudo (MATOPIBA)', fontsize=15)
plt.xlabel('Longitude', fontsize=15, labelpad=30)
plt.ylabel('Latitude', fontsize=15, labelpad=30)

plt.text(-42, -17, u'\u25B2 \nN ', fontsize=15)
plt.text(-45.5, -9, u'MA', color='black', fontsize=15)
plt.text(-49, -10, u'TO', color='black', fontsize=15)
plt.text(-45, -6.5, u'PI', color='black', fontsize=15)
plt.text(-45, -13, u'BA', color='black', fontsize=15)

plt.text(-45.56, -9.70, u'Alto Parnaiba', color='blue', fontsize=10)
plt.text(-46.02, -7.00, u'Balsas', color='blue', fontsize=10)
plt.text(-44.88, -5.62, u'Barra do Corda', color='blue', fontsize=10)
plt.text(-50.34, -10.95, u'Porto Nacional', color='blue', fontsize=10)
plt.text(-49.3, -8.80, u'Pedro Afonso', color='blue', fontsize=10)
plt.text(-48.22, -11.80, u'Peixe', color='blue', fontsize=10)
plt.text(-49.25, -12.60, u'Taguatinga', color='blue', fontsize=10)
plt.text(-42.67, -6.95, u'Floriano', color='blue', fontsize=10)
plt.text(-44.07, -8.80, u'Bom Jesus', color='blue', fontsize=10)
plt.text(-42.81, -4.80, u'Teresina', color='blue', fontsize=10)
plt.text(-45.00, -11.85, u'Barreiras', color='blue', fontsize=10)
plt.text(-43.00, -13.20, u'Correntina', color='blue', fontsize=10)

bmap.plot(-45.56, -9.06, 'ko', mfc='red', markersize=5)
bmap.plot(-46.02, -7.32, 'ko', mfc='red', markersize=5)
bmap.plot(-45.13, -5.30, 'ko', mfc='red', markersize=5)
bmap.plot(-48.25, -10.43, 'ko', mfc='red', markersize=5)
bmap.plot(-48.18, -9.06, 'ko', mfc='red', markersize=5)
bmap.plot(-48.22, -12.00, 'ko', mfc='red', markersize=5)
bmap.plot(-46.25, -12.24, 'ko', mfc='red', markersize=5)
bmap.plot(-43.01, -6.46, 'ko', mfc='red', markersize=5)
bmap.plot(-44.07, -9.06, 'ko', mfc='red', markersize=5)
bmap.plot(-42.81, -5.08, 'ko', mfc='red', markersize=5)
bmap.plot(-45.00, -12.09, 'ko', mfc='red', markersize=5)
bmap.plot(-43.37, -13.20, 'ko', mfc='red', markersize=5)

# Path out to save figure
path_out = '/home/nice'
name_out = 'pyplt_maps_study_area_matopiba.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()



