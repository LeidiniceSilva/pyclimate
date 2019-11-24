# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "11/22/2019"
__description__ = "This script compute time series from Xavier basedata"

import numpy as np
import xarray as xr
import pandas as pd
import time

# Set latitude and longitude
lat = [-6.43, -4.28, -7.01, -9.44]
lon = [-36.58, -39.0, -37.26, -36.7]

# Capitals (AL - BA - CE - MA - PB - PE - PI - RN - SE)
list_city = ['cruzeta', 'guaramiranga', 'patos', 'alagoas']

# Variables names
var_names = ['Tmax', 'Tmin', 'prec']

# Set correct path of the netcdf files
path_var = '/home/nice/Downloads/'

# Function to read the netcdf files
def rawData(var2get_xr, var_name2get):
    return var2get_xr[var_name2get].sel(longitude=xr.DataArray(lon, dims='z'),
                                          latitude=xr.DataArray(lat, dims='z'),
                                          method='nearest').values

# Getting data from NetCDF files
for n, var_name2get in enumerate(var_names):
    var2get_yr = xr.open_mfdataset(path_var + var_name2get + '_daily_UT_Brazil_v2.1_19860101_20051231.nc')
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
    name_file = 'lat{}_lon{}_{}.csv'.format(lat[n], lon[n], list_city[n])
    if ~np.isnan(var_ar[0, n]):
        file = var_ar[:, n::len(lon)]
        pd.DataFrame(file, index=time, columns=var_names).to_csv(name_file, float_format='%.1f')

