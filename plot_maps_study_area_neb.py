# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "11/22/2019"
__description__ = "This script plot study area NEB domain"

import os
import cartopy
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import matplotlib.patches as mpatches

from matplotlib.axes import Axes
from cartopy.feature import NaturalEarthFeature, BORDERS
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.geoaxes import GeoAxes

GeoAxes._pcolormesh_patched = Axes.pcolormesh

# Open xavier basedata
path_var = '/home/nice/Documents/ufrn/papers/wmrn/data/xavier/'
prec = xr.open_mfdataset(path_var + 'prec_daily_UT_Brazil_v2.2_1986-2005.nc').prec

# set states limits
states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none', name='admin_1_states_provinces_shp')

# Plot maps xavier basedata

# set up the plot
fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()))
prec.plot(transform=ccrs.PlateCarree(), cmap='Spectral', cbar_kwargs={'label': 'Precipitação (mm)'})

ax.axes.axis('tight')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_xticks(np.linspace(-75, -33, 5), crs=cartopy.crs.PlateCarree())
ax.set_yticks(np.linspace(-34, 6, 5), crs=cartopy.crs.PlateCarree())

ax.add_patch(mpatches.Rectangle(xy=[-48, -18], width=13, height=18, facecolor='none',
                                    edgecolor='black', transform=ccrs.PlateCarree()))
                                    
ax.text(-35, -33,u'\u25B2 \nN ', ha='center', fontsize=10, family='Arial', rotation = 0)
                          
ax.set_title(u'Climatogia de Precipitação (1986-2005)')
ax.coastlines()
ax.add_feature(states, edgecolor='black', facecolor='none')
ax.add_feature(BORDERS)

plt.tight_layout()
plt.draw()

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/wmrn/results'
name_out = 'pyplt_maps_study_area.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400)

plt.show()
exit()
