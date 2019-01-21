# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "01/21/2019"
__description__ = "This script to extract time series from CMIP5 models"


import netCDF4
import numpy as np

from pylab import *
from netCDF4 import Dataset

var   = 'pr'
model = 'BCC-CSM1-1'
exp   = 'historical_r1i1p1'
data  = '198001-200512'

mod_path = '/vol3/disco1/nice/cmip_data/cmip5_hist/BCC-CSM1.1'
arq1     = '{0}/{1}_amz_neb_Amon_{2}_{3}_mon_{4}.nc'.format(mod_path, var, model, exp, data)	
data1    = netCDF4.Dataset(arq1)

var1     = data1.variables['pr'][:]
lat      = data1.variables['lat'][:]
lon      = data1.variables['lon'][:]
mod      = var1[:][:,:,:]
mod_ini1 = np.nanmean(mod, axis=1)
mod_end1 = np.nanmean(mod_ini1, axis=1)

print len(mod_end1)
print
exit()

print(mod_end1, file=open('test.txt', 'w'))



