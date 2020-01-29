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


lons = [-45.56, -46.02, -45.13, -48.25, -48.18, -48.22, -46.25, -43.01, -44.07, -42.81, -45.00, -43.37]
lats = [-9.06, -7.32, -5.30, -10.43, -9.06, -12.00, -12.24, -6.46, -9.06, -5.08, -12.09, -13.20]

trend_rx1day = [1.009, 0.323, 0.054, -0.062, 0.434, 0.464, 1.309, 0.577, -0.196, -0.218, 0.575, -0.237, -0.034]
trend_rx5day = [1.441, 1.752, 0.757, 0.252, 1.368, 0.963, 1.977, 0.853, 0.846, -0.239, -2.188, 0.183]
trend_sdii = [0.084, 0.105, 0.105, 0.099, 0.079, 0.084, 0.039, 0.041, 0.049, 0.004, -0.03, -0.024]
trend_pr95p = [10.449, 11.013, 1.862, 5.502, 1.509, 5.14, 5.067, -2.594, 2.255, -2.176, -5.465, 1.233]
trend_pr99p = [4.958, 3.658, 1.734, -1.113, 1.978, 1.367, 1.551, 2.949, 0.918, 2.354, -6.231, 0.335]
trend_prectot = [3.958, 5.893, 5.632, 2.789, -0.025, 8.274, -1.751, -3.34, 1.918, -3.67, -10.198, -3.365]

trend_txx = [0.072, 0.077, 0.042, 0.062, 0.163, 0.06, 0.213, 0.044, 0.046, 0.047, 0.032, 0.146]
trend_txn = [0.045, 0.101, 0.09, 0.045, 0.084, 0.074, -0.003, 0.079, 0.033, 0.035, 0.095, -0.025]
trend_tnx = [-0.004, 0.05, 0.095, 0.002, 0.021, 0.022, 0.178, 0.036, 0.038, 0.043, 0.054, 0.068]
trend_tnn = [-0.072, 0.004, 0.089, -1.199, 0.016, -0.02, 0.019, 0.032, 0.001, 0.047, -0.074, 0.013]
trend_tx90p = [0.574, 0.818, 0.595, 0.299, 1.493, 0.138, 1.007, 0.659, 0.304, 0.321, 0.178, 0.534]
trend_tn90p = [-0.206, 0.679, -0.504, -0.214, -0.065, -0.156, 1.288, 0.701, 0.407, 0.705, 0.394, 1.002]

# Plot maps with trend from extremes index - Temperatuure
fig = plt.figure()

# Plot firt map 
ax = fig.add_subplot(231)
bmap, x, y = basemap()
patch = shipefile()

#~ def get_marker_color(trend_tn90p):
    #~ if trend_tn90p > 0.0:
        #~ return ('r^')
    #~ elif trend_tn90p < 0.0:
        #~ return ('rv')
    #~ else:
        #~ return ('ro')

#~ for n in range(len(lats)):
    #~ marker_string = get_marker_color(trend_tn90p[n])
    #~ print(marker_string, trend_tn90p[n]*20)
#~ exit()

plt.text(-52.5, -2, u'A) Txx', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)
plt.ylabel(u'Latitude', fontsize=8, labelpad=25)

bmap.plot(-45.56, -9.06, 'k^', mfc='lightcoral', markersize=2.44)
bmap.plot(-46.02, -7.32, 'k^', mfc='lightcoral', markersize=2.54)
bmap.plot(-45.13, -5.30, 'r^', mfc='lightcoral', markersize=1.85)
bmap.plot(-48.25, -10.43, 'k^', mfc='lightcoral', markersize=2.24)
bmap.plot(-48.18, -9.06, 'k^', mfc='lightcoral', markersize=3.26)
bmap.plot(-48.22, -12.00, 'r^', mfc='lightcoral', markersize=2.2)
bmap.plot(-46.25, -12.24, 'k^', mfc='lightcoral', markersize=4.26)
bmap.plot(-43.01, -6.46, 'k^', mfc='lightcoral', markersize=1.87)
bmap.plot(-44.07, -9.06, 'k^', mfc='lightcoral', markersize=1.91)
bmap.plot(-42.81, -5.08, 'k^', mfc='lightcoral', markersize=1.94)
bmap.plot(-45.00, -12.09, 'r^', mfc='lightcoral', markersize=1.64)
bmap.plot(-43.37, -13.20, 'k^', mfc='lightcoral', markersize=3.92)

# Plot second map 
ax = fig.add_subplot(232)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'B) Txn', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)

bmap.plot(-45.56, -9.06, 'r^', mfc='lightcoral', markersize=1.89)
bmap.plot(-46.02, -7.32, 'k^', mfc='lightcoral', markersize=3.02)
bmap.plot(-45.13, -5.30, 'k^', mfc='lightcoral', markersize=2.79)
bmap.plot(-48.25, -10.43, 'r^', mfc='lightcoral', markersize=1.89)
bmap.plot(-48.18, -9.06, 'r^', mfc='lightcoral', markersize=2.68)
bmap.plot(-48.22, -12.00, 'r^', mfc='lightcoral', markersize=2.48)
bmap.plot(-46.25, -12.24, 'rv', mfc='lightcoral', markersize=1.5)
bmap.plot(-43.01, -6.46, 'k^', mfc='lightcoral', markersize=2.58)
bmap.plot(-44.07, -9.06, 'r^', mfc='lightcoral', markersize=1.66)
bmap.plot(-42.81, -5.08, 'r^', mfc='lightcoral', markersize=1.70)
bmap.plot(-45.00, -12.09, 'k^', mfc='lightcoral', markersize=2.9)
bmap.plot(-43.37, -13.20, 'rv', mfc='lightcoral', markersize=1.5)

# Plot third map 
ax = fig.add_subplot(233)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'C) Tnx', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)

bmap.plot(-45.56, -9.06, 'rv', mfc='lightcoral', markersize=1.8)
bmap.plot(-46.02, -7.32, 'k^', mfc='lightcoral', markersize=2.0)
bmap.plot(-45.13, -5.30, 'k^', mfc='lightcoral', markersize=3.9)
bmap.plot(-48.25, -10.43, 'r^', mfc='lightcoral', markersize=1.5)
bmap.plot(-48.18, -9.06, 'r^', mfc='lightcoral', markersize=1.52)
bmap.plot(-48.22, -12.00, 'r^', mfc='lightcoral', markersize=1.53)
bmap.plot(-46.25, -12.24, 'k^', mfc='lightcoral', markersize=3.55)
bmap.plot(-43.01, -6.46, 'k^', mfc='lightcoral', markersize=1.72)
bmap.plot(-44.07, -9.06, 'r^', mfc='lightcoral', markersize=1.76)
bmap.plot(-42.81, -5.08, 'k^', mfc='lightcoral', markersize=1.85)
bmap.plot(-45.00, -12.09, 'k^', mfc='lightcoral', markersize=2.08)
bmap.plot(-43.37, -13.20, 'k^', mfc='lightcoral', markersize=2.36)
# Plot fourth map 
ax = fig.add_subplot(234)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'D) Tnn', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)
plt.xlabel(u'Longitude', fontsize=8, labelpad=20)
plt.ylabel(u'Latitude', fontsize=8, labelpad=25)

bmap.plot(-45.56, -9.06, 'rv', mfc='lightcoral', markersize=2.44)
bmap.plot(-46.02, -7.32, 'r^', mfc='lightcoral', markersize=1.8)
bmap.plot(-45.13, -5.30, 'k^', mfc='lightcoral', markersize=2.77)
bmap.plot(-48.25, -10.43, 'kv', mfc='lightcoral', markersize=12.0)
bmap.plot(-48.18, -9.06, 'r^', mfc='lightcoral', markersize=1.32)
bmap.plot(-48.22, -12.00, 'rv', mfc='lightcoral', markersize=1.4)
bmap.plot(-46.25, -12.24, 'r^', mfc='lightcoral', markersize=1.38)
bmap.plot(-43.01, -6.46, 'r^', mfc='lightcoral', markersize=1.64)
bmap.plot(-44.07, -9.06, 'r^', mfc='lightcoral', markersize=1.2)
bmap.plot(-42.81, -5.08, 'r^', mfc='lightcoral', markersize=1.94)
bmap.plot(-45.00, -12.09, 'kv', mfc='lightcoral', markersize=1.48)
bmap.plot(-43.37, -13.20, 'r^', mfc='lightcoral', markersize=1.3)

# Plot fifth map 
ax = fig.add_subplot(235)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'E) Tx90p', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)
plt.xlabel(u'Longitude', fontsize=8, labelpad=20)

bmap.plot(-45.56, -9.06, 'r^', mfc='lightcoral', markersize=6.447)
bmap.plot(-46.02, -7.32, 'k^', mfc='lightcoral', markersize=8.364)
bmap.plot(-45.13, -5.30, 'k^', mfc='lightcoral', markersize=6.85)
bmap.plot(-48.25, -10.43, 'r^', mfc='lightcoral', markersize=3.44)
bmap.plot(-48.18, -9.06, 'k^', mfc='lightcoral', markersize=14.76)
bmap.plot(-48.22, -12.00, 'r^', mfc='lightcoral', markersize=1.69)
bmap.plot(-46.25, -12.24, 'k^', mfc='lightcoral', markersize=10.01)
bmap.plot(-43.01, -6.46, 'k^', mfc='lightcoral', markersize=7.09)
bmap.plot(-44.07, -9.06, 'r^', mfc='lightcoral', markersize=3.04)
bmap.plot(-42.81, -5.08, 'r^', mfc='lightcoral', markersize=3.21)
bmap.plot(-45.00, -12.09, 'r^', mfc='lightcoral', markersize=1.77)
bmap.plot(-43.37, -13.20, 'r^', mfc='lightcoral', markersize=5.68)

# Plot sixth map 
ax = fig.add_subplot(236)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'F) Tn90p', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)
plt.xlabel(u'Longitude', fontsize=8, labelpad=20)

bmap.plot(-45.56, -9.06, 'kv', mfc='lightcoral', markersize=2.12)
bmap.plot(-46.02, -7.32, 'k^', mfc='lightcoral', markersize=7.50)
bmap.plot(-45.13, -5.30, 'rv', mfc='lightcoral', markersize=5.05)
bmap.plot(-48.25, -10.43, 'rv', mfc='lightcoral', markersize=2.28)
bmap.plot(-48.18, -9.06, 'rv', mfc='lightcoral', markersize=0.65)
bmap.plot(-48.22, -12.00, 'rv', mfc='lightcoral', markersize=1.6)
bmap.plot(-46.25, -12.24, 'k^', mfc='lightcoral', markersize=12.76)
bmap.plot(-43.01, -6.46, 'k^', mfc='lightcoral', markersize=7.02)
bmap.plot(-44.07, -9.06, 'k^', mfc='lightcoral', markersize=4.13)
bmap.plot(-42.81, -5.08, 'k^', mfc='lightcoral', markersize=7.1)
bmap.plot(-45.00, -12.09, 'k^', mfc='lightcoral', markersize=3.99)
bmap.plot(-43.37, -13.20, 'k^', mfc='lightcoral', markersize=10.0)

# Path out to save figure
path_out = '/home/nice'
name_out = 'pyplt_maps_trend_temp_rclimdex_xavier_matopiba.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()


# Plot maps with trend from extremes index - Rainfall
fig = plt.figure()

# Plot firt map 
ax = fig.add_subplot(231)
bmap, x, y = basemap()
patch = shipefile()

#~ def get_marker_color(trend_prectot):
    #~ if trend_prectot > 0.0:
        #~ return ('b^')
    #~ elif trend_prectot < 0.0:
        #~ return ('bv')
    #~ else:
        #~ return ('ro')

#~ for n in range(len(lats)):
    #~ marker_string = get_marker_color(trend_prectot[n])
    #~ print(marker_string, trend_prectot[n]*1)
#~ exit()

plt.text(-52.5, -2, u'A) Rx1dia', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)
plt.ylabel(u'Latitude', fontsize=8, labelpad=25)

bmap.plot(-45.56, -9.06, 'k^', mfc='cornflowerblue', markersize=10.09)
bmap.plot(-46.02, -7.32, 'b^', mfc='cornflowerblue', markersize=3.23)
bmap.plot(-45.13, -5.30, 'b^', mfc='cornflowerblue', markersize=0.54)
bmap.plot(-48.25, -10.43, 'bv', mfc='cornflowerblue', markersize=0.62)
bmap.plot(-48.18, -9.06, 'b^', mfc='cornflowerblue', markersize=4.34)
bmap.plot(-48.22, -12.00, 'b^', mfc='cornflowerblue', markersize=4.64)
bmap.plot(-46.25, -12.24, 'k^', mfc='cornflowerblue', markersize=13.09)
bmap.plot(-43.01, -6.46, 'b^', mfc='cornflowerblue', markersize=5.77)
bmap.plot(-44.07, -9.06, 'bv', mfc='cornflowerblue', markersize=1.96)
bmap.plot(-42.81, -5.08, 'bv', mfc='cornflowerblue', markersize=2.18)
bmap.plot(-45.00, -12.09, 'b^', mfc='cornflowerblue', markersize=5.75)
bmap.plot(-43.37, -13.20, 'bv', mfc='slategrey', markersize=2.37)

# Plot second map 
ax = fig.add_subplot(232)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'B) Rx5dia', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)

bmap.plot(-45.56, -9.06, 'b^', mfc='cornflowerblue', markersize=7.41)
bmap.plot(-46.02, -7.32, 'k^', mfc='cornflowerblue', markersize=9.52)
bmap.plot(-45.13, -5.30, 'b^', mfc='cornflowerblue', markersize=4.52)
bmap.plot(-48.25, -10.43, 'b^', mfc='cornflowerblue', markersize=1.26)
bmap.plot(-48.18, -9.06, 'b^', mfc='cornflowerblue', markersize=7.34)
bmap.plot(-48.22, -12.00, 'b^', mfc='cornflowerblue', markersize=4.81)
bmap.plot(-46.25, -12.24, 'b^', mfc='cornflowerblue', markersize=9.99)
bmap.plot(-43.01, -6.46, 'b^', mfc='cornflowerblue', markersize=4.53)
bmap.plot(-44.07, -9.06, 'b^', mfc='cornflowerblue', markersize=4.45)
bmap.plot(-42.81, -5.08, 'bv', mfc='cornflowerblue', markersize=1.20)
bmap.plot(-45.00, -12.09, 'kv', mfc='cornflowerblue', markersize=10.44)
bmap.plot(-43.37, -13.20, 'b^', mfc='cornflowerblue', markersize=0.93)

# Plot third map 
ax = fig.add_subplot(233)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'C) SDII', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)

bmap.plot(-45.56, -9.06, 'k^', mfc='cornflowerblue', markersize=1.84)
bmap.plot(-46.02, -7.32, 'k^', mfc='cornflowerblue', markersize=2.05)
bmap.plot(-45.13, -5.30, 'k^', mfc='cornflowerblue', markersize=2.05)
bmap.plot(-48.25, -10.43, 'b^', mfc='cornflowerblue', markersize=1.99)
bmap.plot(-48.18, -9.06, 'b^', mfc='cornflowerblue', markersize=1.79)
bmap.plot(-48.22, -12.00, 'b^', mfc='cornflowerblue', markersize=1.84)
bmap.plot(-46.25, -12.24, 'b^', mfc='cornflowerblue', markersize=1.39)
bmap.plot(-43.01, -6.46, 'b^', mfc='cornflowerblue', markersize=1.41)
bmap.plot(-44.07, -9.06, 'b^', mfc='cornflowerblue', markersize=1.49)
bmap.plot(-42.81, -5.08, 'b^', mfc='cornflowerblue', markersize=1.04)
bmap.plot(-45.00, -12.09, 'bv', mfc='cornflowerblue', markersize=1.3)
bmap.plot(-43.37, -13.20, 'bv', mfc='cornflowerblue', markersize=1.24)

# Plot fourth map 
ax = fig.add_subplot(234)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'D) Pr95p', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)
plt.xlabel(u'Longitude', fontsize=8, labelpad=20)
plt.ylabel(u'Latitude', fontsize=8, labelpad=25)

bmap.plot(-45.56, -9.06, 'k^', mfc='cornflowerblue', markersize=10.45)
bmap.plot(-46.02, -7.32, 'k^', mfc='cornflowerblue', markersize=11.01)
bmap.plot(-45.13, -5.30, 'b^', mfc='cornflowerblue', markersize=1.862)
bmap.plot(-48.25, -10.43, 'b^', mfc='cornflowerblue', markersize=5.502)
bmap.plot(-48.18, -9.06, 'b^', mfc='cornflowerblue', markersize=1.509)
bmap.plot(-48.22, -12.00, 'b^', mfc='cornflowerblue', markersize=5.14)
bmap.plot(-46.25, -12.24, 'b^', mfc='cornflowerblue', markersize=5.067)
bmap.plot(-43.01, -6.46, 'bv', mfc='cornflowerblue', markersize=2.594)
bmap.plot(-44.07, -9.06, 'b^', mfc='cornflowerblue', markersize=2.255)
bmap.plot(-42.81, -5.08, 'bv', mfc='cornflowerblue', markersize=2.176)
bmap.plot(-45.00, -12.09, 'bv', mfc='cornflowerblue', markersize=5.465)
bmap.plot(-43.37, -13.20, 'b^', mfc='cornflowerblue', markersize=1.233)

# Plot fifth map 
ax = fig.add_subplot(235)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'E) Pr99p', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)
plt.xlabel(u'Longitude', fontsize=8, labelpad=20)

bmap.plot(-45.56, -9.06, 'b^', mfc='cornflowerblue', markersize=4.958)
bmap.plot(-46.02, -7.32, 'b^', mfc='cornflowerblue', markersize=3.658)
bmap.plot(-45.13, -5.30, 'b^', mfc='cornflowerblue', markersize=1.734)
bmap.plot(-48.25, -10.43, 'bv', mfc='cornflowerblue', markersize=1.113)
bmap.plot(-48.18, -9.06, 'b^', mfc='cornflowerblue', markersize=1.978)
bmap.plot(-48.22, -12.00, 'b^', mfc='cornflowerblue', markersize=1.367)
bmap.plot(-46.25, -12.24, 'b^', mfc='cornflowerblue', markersize=1.551)
bmap.plot(-43.01, -6.46, 'bv', mfc='cornflowerblue', markersize=2.949)
bmap.plot(-44.07, -9.06, 'b^', mfc='cornflowerblue', markersize=0.918)
bmap.plot(-42.81, -5.08, 'b^', mfc='cornflowerblue', markersize=2.354)
bmap.plot(-45.00, -12.09, 'bv', mfc='cornflowerblue', markersize=6.231)
bmap.plot(-43.37, -13.20, 'b^', mfc='cornflowerblue', markersize=0.335)

# Plot sixth map 
ax = fig.add_subplot(236)
bmap, x, y = basemap()
patch = shipefile()

plt.text(-52.5, -2, u'F) PRECTOT', fontsize=10)
plt.text(-42, -17, u'\u25B2 \nN ', fontsize=8)
plt.xlabel(u'Longitude', fontsize=8, labelpad=20)

bmap.plot(-45.56, -9.06, 'b^', mfc='cornflowerblue', markersize=3.958)
bmap.plot(-46.02, -7.32, 'b^', mfc='cornflowerblue', markersize=5.893)
bmap.plot(-45.13, -5.30, 'b^', mfc='cornflowerblue', markersize=5.632)
bmap.plot(-48.25, -10.43, 'b^', mfc='cornflowerblue', markersize=2.789)
bmap.plot(-48.18, -9.06, 'bv', mfc='cornflowerblue', markersize=0.025)
bmap.plot(-48.22, -12.00, 'b^', mfc='cornflowerblue', markersize=8.274)
bmap.plot(-46.25, -12.24, 'bv', mfc='cornflowerblue', markersize=1.751)
bmap.plot(-43.01, -6.46, 'bv', mfc='cornflowerblue', markersize=3.34)
bmap.plot(-44.07, -9.06, 'b^', mfc='cornflowerblue', markersize=1.918)
bmap.plot(-42.81, -5.08, 'bv', mfc='cornflowerblue', markersize=3.67)
bmap.plot(-45.00, -12.09, 'bv', mfc='cornflowerblue', markersize=10.198)
bmap.plot(-43.37, -13.20, 'bv', mfc='cornflowerblue', markersize=3.365)

# Path out to save figure
path_out = '/home/nice'
name_out = 'pyplt_maps_trend_prec_rclimdex_xavier_matopiba.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()







