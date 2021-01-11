import os
import conda
#~ import numpy as np
#~ import matplotlib as mpl  #; mpl.use('Agg')
#~ import matplotlib.pyplot as plt
#~ import warnings ; warnings.filterwarnings("ignore")

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cartopy.feature import NaturalEarthFeature, LAND, COASTLINE
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def brazil_states(projection=ccrs.PlateCarree()):
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw=dict(projection=projection))
    ax.set_extent([-82, -32, -45, 10])
    ax.stock_img()
    ax.add_feature(LAND)
    ax.add_feature(COASTLINE)
    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return fig, ax


#~ fig, ax = brazil_states()
#~ states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                             #~ name='admin_1_states_provinces_shp')

#~ _ = ax.add_feature(states, edgecolor='gray')

    
dx, dy = 0.5, 0.5

# generate 2 2d grids for the x & y bounds
y, x = np.mgrid[slice(-20, 15 + dy, dy),
                slice(-90, -30 + dx, dx)]

z = np.sin(x) + np.cos(5 + y*x) * np.cos(x) + x*y/100
n, m = z.shape
z = z - np.exp(np.random.rand(n,m))

#Setup the map
m = Basemap(projection='merc', llcrnrlat=-20, urcrnrlat=15,\
            llcrnrlon=-90, urcrnrlon=-30, resolution='l')
m.drawcoastlines()
m.drawcountries()
m.drawstates()


xx, yy = m(x,y)
cs = m.contourf(xx,yy,z,cmap=plt.cm.Spectral_r)
cbar = plt.colorbar(cs, orientation='horizontal', shrink=0.5)

##adding hatches:
cs = m.contourf(xx, yy, z, levels=[np.min(z),5,np.max(z)], colors='none',
                  hatches=[None,'.',],
                  extend='lower')

plt.show()
