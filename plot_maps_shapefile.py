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

print
print('# Importing data')
print

obse = 'pr_cmap.nc'
obse_d   = Dataset(obse)
obse_v   = obse_d.variables['pr'][0,...]
lat = obse_d.variables['Y'][:]
lon = obse_d.variables['X'][:]

#~ print 'input -->', obse
#~ print 'shape -->', obse_v.shape

#~ print 'min -->', np.nanmin(obse_v)
#~ print 'max -->', np.nanmax(obse_v)

y1 = -18
y2 =   0
x1 = -53
x2 = -40

    
# Preparing lat lon
lat = lat - np.mean(np.diff(lat))/2.
lon = lon - np.mean(np.diff(lon))/2.
lons, lats = np.meshgrid(lon, lat)

print
print('# Stating figure')

fig = plt.figure(figsize=(10, 8))

plt.title('A) Rx5day', fontsize=15)

parallels = np.arange( -90,  91, 2)
meridians = np.arange(-180, 180, 2)

bmap = Basemap(projection='cyl', llcrnrlat=y1, urcrnrlat=y2, llcrnrlon=x1, urcrnrlon=x2, resolution=None, suppress_ticks=True, lon_0=0, celestial=False)
bmap.drawparallels(parallels, labels=[1, 0, 0, 0], color='gray', linewidth=0.01, dashes=[4, 6], fontsize=8, zorder=1)
bmap.drawmeridians(meridians, labels=[0, 0, 0, 1], color='gray', linewidth=0.01, dashes=[4, 6], fontsize=8, zorder=1)

x, y = bmap(lons, lats)

print
print('# Drawing shapes contour')
print

bmap.readshapefile('./world', 'world', drawbounds=True)

bmap.readshapefile('./matopiba', 'matopiba', drawbounds=True, color='black', linewidth=2.5)


print
print('# Masking map with shape')
print

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


print
print('# Creating map')
print

cs0 = plt.pcolormesh(x, y, obse_v, cmap=my_cmap, norm=norm)

plt.colorbar(cs0, orientation='vertical', extend='max', spacing='uniform', ticks=pale_lev, extendfrac='auto', cmap=my_cmap)

print
print('# Saving figure')
print

plt.savefig('matopiba_test.png.png', bbox_inches='tight', dpi=200)
plt.close()

if os.path.exists('matopiba_test.png'):
    print('done --> matopiba.png')




