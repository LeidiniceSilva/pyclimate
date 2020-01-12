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
#~ lat = [-9.06, -7.32, -5.30, -10.43, -8.96, -12.00]
#~ lon = [-45.56, -46.02, -45.13, -48.25, -48.18, -48.22]
lat = [-12.24, -6.46, -9.06, -5.08, -12.09, -13.20]
lon = [-46.25, -43.01, -44.07, -42.81, -45.00, -43.37]

# MATOPIBA stations 
#~ list_city = ['alto_parnaiba', 'balsas', 'barra_corda', 'porto_nacional', 'pedro_afonso', 'peixe']
list_city = ['taguatinga', 'floriano', 'bom_jesus', 'teresina', 'barreiras', 'correntina']

# Variables names
var_names = ['Tmax', 'Tmin', 'prec']

# Set correct path of the netcdf files
path_var = '/home/nice/Documents/'

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

