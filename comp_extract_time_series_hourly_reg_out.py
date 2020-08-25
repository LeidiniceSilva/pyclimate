# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/24/2020"
__description__ = "This script compute time serie from RegCM output"

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import numpy as np
import xarray as xr
import pandas as pd
import time

# Set latitude and longitude
lat = [-6.57]
lon = [-37.25]

# Stations 
list_city = ['serra_negra']

# Variable name
var_names = ['rsnl', 'rsns']

# Set correct path of the netcdf files
path_var = '/home/nice/Documents/choluck/rttm/'

# Function to read the netcdf files
def rawData(var2get_xr, var_name2get):
    return var2get_xr[var_name2get].sel(lon=xr.DataArray(lon, dims='z'),
                                          lat=xr.DataArray(lat, dims='z'),
                                          method='nearest').values

# Getting data from NetCDF files
for n, var_name2get in enumerate(var_names):
    var2get_yr = xr.open_mfdataset(path_var + var_name2get + '_rttm_hourly_2014_lonlat.nc')
    var2get_xr = (xr.combine_by_coords([var2get_yr]))

    if n == 0:
        var_ar = rawData(var2get_xr, var_name2get)
        n_lines = var_ar.shape[0]
        time = var2get_xr.time.values
    elif var_name2get != 'prec':
        var_ar = np.c_[var_ar, rawData(var2get_xr, var_name2get)]
    else:
        prec = rawData(var2get_xr, var_name2get)
        prec_size_var_ar = np.zeros((var_ar.shape[0], len(lon))) * np.nan
        prec_size_var_ar[:prec.shape[0], :] = prec
        var_ar = np.c_[var_ar, prec_size_var_ar]

# Saving file.csv
for n in range(len(lat)):
    print('File {} of {}'.format(n+1, len(lat)))
    name_file = 'lat{}_lon{}_{}_rttm_reg_2014.csv'.format(lat[n], lon[n], list_city[n])
    if ~np.isnan(var_ar[0, n]):
        file = var_ar[:, n::len(lon)]
        pd.DataFrame(file, index=time, columns=var_names).to_csv(name_file, float_format='%.1f')
        
        
      
