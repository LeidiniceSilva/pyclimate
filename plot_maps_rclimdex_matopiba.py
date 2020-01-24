# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "24/01/2020"
__description__ = "This script plot trend from extremes index"

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


def basemap():
	
	lat = np.arange( -90,  91, 2)
	lon = np.arange(-180, 180, 2)
	
	lat = lat - np.mean(np.diff(lat))/2.
	lon = lon - np.mean(np.diff(lon))/2.

	bmap = Basemap(projection='cyl', llcrnrlat=-18, urcrnrlat=0, llcrnrlon=-53, urcrnrlon=-40, resolution=None, suppress_ticks=True, lon_0=0, celestial=False)
	bmap.drawparallels(np.arange( -90,  91, 2), labels=[1, 0, 0, 0], color='gray', linewidth=0.01, dashes=[4, 6], fontsize=8, zorder=1)
	bmap.drawmeridians(np.arange(-180, 180, 4), labels=[0, 0, 0, 1], color='gray', linewidth=0.01, dashes=[4, 6], fontsize=8, zorder=1)
	
	lons, lats = np.meshgrid(lon, lat)
	x, y = bmap(lons, lats)
	
	return bmap, x, y
	

def shipefile():
	
	bmap, x, y = basemap()

	# Import shapefile from word and matopiba 
	bmap.readshapefile('/home/nice/Documents/github_projects/shp/shp_world/world', 'world', drawbounds=True, color='white')
	bmap.readshapefile('/home/nice/Documents/github_projects/shp/shp_matopiba/matopiba', 'matopiba', drawbounds=True, color='black', linewidth=2)

	# Masking map with shape
	x0, x1 = plt.xlim()
	y0, y1 = plt.ylim()

	map_edges = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
	polys = [map_edges]

	bmap.readshapefile('/home/nice/Documents/github_projects/shp/shp_matopiba/matopiba', 'matopiba2', drawbounds=False)
	polys = polys + getattr(bmap, 'matopiba2')
	codes = [[Path.MOVETO] + [Path.LINETO for p in p[1:]] for p in polys] # creating a PathPatch

	polys_lin = [v for p in polys for v in p]
	codes_lin = [cdg for cdgs in codes for cdg in cdgs]
	path  = Path(polys_lin, codes_lin)
	patch = PathPatch(path, facecolor='white', lw=0)
	
	return patch


lons = [-36.70, -44.61, -40.79, -41.86, -39.00, -47.48, -45.93, -37.26, -43.35, -40.46, -43.71, -37.04]
lats = [-9.44, -13.33, -14.88, -11.3, -4.28, -5.53, -9.1, -7.01, -4.86, -9.36, -8.41, -10.95]
trend_cdd = [0.872, 1.190, -0.654, 0.444, 3.16, -0.789, -1.198, 0.718, 0.68, -1.592, 0.596, 0.137]
trend_cwd = [-3.388, 0.296, 0.098, -0.428, -1.242, -0.226, 0.854, 0.059, -0.004, -0.002, 1.0212, 0.0273]
trend_r95p = [12,14, -2.03, -1.515, 3.4922, 4.2, 9.478, 6.717, -4.784, 6.522, 2.04, 1.021, 14.152]
trend_txx = [0.077, 0.062, 0.087, 0.11, -0.121, 0.066, 0.12, 0.112, 0.169, 0.126, 0.138, -0.068]

#~ # Define size from marker
#~ def get_marker_color(trend):
    #~ if trend > 0.0:
        #~ return ('r^')
    #~ elif trend < 0.0:
        #~ return ('bv')
    #~ else:
        #~ return ('ro')

#~ for n in range(len(lats)):
    #~ x,y = map(lons, lats)
    #~ marker_string = get_marker_color(trend_cwd[n])
    #~ map.plot(x[n], y[n], marker_string, markersize=trend_cwd[n]*1)
    #~ print(x[n], y[n], marker_string, trend_cwd[n]*1)
    
# Plot maps with trend from extremes index - Rainfall
fig = plt.figure()

# Plot firt map 
ax = fig.add_subplot(231)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'A) Pr95p', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)
plt.ylabel(u'Latitude', fontsize=8, labelpad=25)

bmap.plot(-45.56, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-46.02, -7.32, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-45.13, -5.30, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.25, -10.43, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.18, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.22, -12.00, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-46.25, -12.24, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-43.01, -6.46, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-44.07, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-42.81, -5.08, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-45.00, -12.09, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-43.37, -13.20, 'b^', mfc='cornflowerblue', markersize=5)

# Plot second map 
ax = fig.add_subplot(232)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'B) Pr99p', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)

bmap.plot(-45.56, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-46.02, -7.32, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-45.13, -5.30, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.25, -10.43, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.18, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.22, -12.00, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-46.25, -12.24, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-43.01, -6.46, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-44.07, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-42.81, -5.08, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-45.00, -12.09, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-43.37, -13.20, 'b^', mfc='cornflowerblue', markersize=5)

# Plot third map 
ax = fig.add_subplot(233)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'C) Rx1dia', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)

bmap.plot(-45.56, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-46.02, -7.32, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-45.13, -5.30, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.25, -10.43, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.18, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.22, -12.00, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-46.25, -12.24, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-43.01, -6.46, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-44.07, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-42.81, -5.08, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-45.00, -12.09, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-43.37, -13.20, 'b^', mfc='cornflowerblue', markersize=5)

# Plot fourth map 
ax = fig.add_subplot(234)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'D) Rx5dia', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)
plt.xlabel(u'Longitude', fontsize=8, labelpad=20)
plt.ylabel(u'Latitude', fontsize=8, labelpad=25)

bmap.plot(-45.56, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-46.02, -7.32, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-45.13, -5.30, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.25, -10.43, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.18, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.22, -12.00, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-46.25, -12.24, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-43.01, -6.46, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-44.07, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-42.81, -5.08, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-45.00, -12.09, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-43.37, -13.20, 'b^', mfc='cornflowerblue', markersize=5)

# Plot fifth map 
ax = fig.add_subplot(235)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'E) SDII', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)
plt.xlabel(u'Longitude', fontsize=8, labelpad=20)

bmap.plot(-45.56, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-46.02, -7.32, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-45.13, -5.30, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.25, -10.43, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.18, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.22, -12.00, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-46.25, -12.24, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-43.01, -6.46, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-44.07, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-42.81, -5.08, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-45.00, -12.09, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-43.37, -13.20, 'b^', mfc='cornflowerblue', markersize=5)

# Plot sixth map 
ax = fig.add_subplot(236)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'F) PRECTOT', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)
plt.xlabel(u'Longitude', fontsize=8, labelpad=20)

bmap.plot(-45.56, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-46.02, -7.32, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-45.13, -5.30, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.25, -10.43, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.18, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-48.22, -12.00, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-46.25, -12.24, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-43.01, -6.46, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-44.07, -9.06, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-42.81, -5.08, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-45.00, -12.09, 'b^', mfc='cornflowerblue', markersize=5)
bmap.plot(-43.37, -13.20, 'b^', mfc='cornflowerblue', markersize=5)

# Path out to save figure
path_out = '/home/nice'
name_out = 'pyplt_maps_trend_rclimdex_xavier_matopiba.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()











