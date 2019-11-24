# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/13/2019"
__description__ = "This script to extract time series from Xavier database"


import netCDF4
import numpy as np
# mpl.use('Agg')

var=u'Tmin'

arq   = '{0}_sf_2016.nc'.format(var)
data  = netCDF4.Dataset(arq)
var   = data.variables[var][:] 
lat   = data.variables['latitude'][:]
lon   = data.variables['longitude'][:]
value = var[:][:,:,:]

ts_day = np.nanmean(np.nanmean(value[:,:,:], axis=1), axis=1)
data_day  = ts_day.reshape(366, 1)

data_day=[]
for i in ts_day:
	print(round(i,2))
	data_day.append(i)

path_out  = '/home/nice/'
file_out = '{0}/{1}_sf_2016.txt'.format(path_out, var)	
np.savetxt(file_out, (data_day), fmt="%s")

exit()
		


