# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "01/08/2019"
__description__ = "This script plot boxplot and time series from CMIP5"


import netCDF4
import numpy as np

from pylab import *
from netCDF4 import Dataset


model = 'bcc-csm1-1'
var   = 'pr'
	
mod_path = '/home/nice/Downloads/{0}'.format(model.upper())
arq1     = '{0}/{1}_Amon_{2}_historical_r1i1p1_185001-201212.nc'.format(mod_path, var, model)	
data1    = netCDF4.Dataset(arq1)
print data1
exit()

var1     = data1.variables['pr'][:]
lat      = data1.variables['lat'][:]
lon      = data1.variables['lon'][:]
mod      = var1[:][:,:,:]
mod_ini1 = np.nanmean(mod, axis=1)
mod_end1 = np.nanmean(mod_ini1, axis=1)

print mod_end1
print

exit()

