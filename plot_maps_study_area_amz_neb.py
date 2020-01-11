# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/10/2020"
__description__ = "This script plot study area AMZ na NEB domain"

import os
import conda
import cartopy
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# Create a large figure:
fig = plt.figure()
ax = fig.add_subplot(111)

# First plot
map1 = Basemap(projection='cyl', lat_0=0, lon_0=0)
map1.drawmeridians(np.arange(-180.,181.,60.), labels=[0,0,0,1], linewidth=0.4)
map1.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,0], linewidth=0.4)
map1.drawcoastlines(linewidth=1, color='k')
map1.drawcountries(linewidth=1, color='k')

plt.title(u'Área de Estudo (AMZ_NEB)', fontsize=15)
plt.xlabel(u'Longitude', fontsize=15, labelpad=20)
plt.ylabel(u'Latitude', fontsize=15, labelpad=30)

lons = np.arange(-180,180,0.25) 
lats  = np.arange(90,-90,-0.25) 
x, y = map1(lons, lats)

axins = zoomed_inset_axes(ax, 4, loc="center", bbox_to_anchor=(0,0))
axins.set_xlim(-85, -15)
axins.set_ylim(-20, 10)

# Second plot
map2 = Basemap(llcrnrlon=-85, llcrnrlat=-20, urcrnrlon=-15, urcrnrlat=10, ax=axins)
map2.etopo()
map2.drawcoastlines(linewidth=1, color='k')
mark_inset(ax, axins, loc1=4, loc2=2)

plt.text(-22, 3.5, u'\u25B2 \nN ', ha='center', fontsize=10, family='Arial', rotation = 0)


# Path out to save figure
path_out = '/home/nice'
name_out = 'pyplt_maps_study_area_amz_neb.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()


# Create a large figure:
fig = plt.figure()

# Add an axes set and draw coastlines:
ax1 = plt.axes([0.01, 0.49, 0.8, 0.5], projection=ccrs.PlateCarree())
ax1.set_global()
ax1.coastlines()
ax1.set_title(u'Área de Estudo (AMZ_NEB)', fontsize=15)

# Draw the rectangular extent of the second plot on the first:
ax1.add_patch(mpatches.Rectangle(xy=[-85, -20], width=70, height=30, facecolor='none',
                                    edgecolor='blue', transform=ccrs.PlateCarree()))
                                    
# Add a second axes set (overlaps first) and draw coastlines:
ax2 = plt.axes([0.35, 0.15, 0.5, 0.3], projection=ccrs.PlateCarree())
ax2.set_extent([-85, -15, -20, 10], crs=ccrs.PlateCarree())
ax2.coastlines()
ax2.set_xlabel('Longitude')
ax2.set_ylabel('Latitude')
ax2.set_xticks(np.linspace(-85, -15, 5), crs=cartopy.crs.PlateCarree())
ax2.set_yticks(np.linspace(-20, 10, 5), crs=cartopy.crs.PlateCarree())
ax2.text(-22, 3.5, u'\u25B2 \nN ', ha='center', fontsize=10, family='Arial', rotation = 0)

# Path out to save figure
path_out = '/home/nice'
name_out = 'pyplt_maps_study_area_amz_neb.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()
