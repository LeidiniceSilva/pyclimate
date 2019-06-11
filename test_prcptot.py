# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilvae@gmail.com"
__date__        = "06/03/2019"
__description__ = "This script plot Rclimdex based in HadGEM2-ES model CMIP5"

import netCDF4
import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from matplotlib.font_manager import FontProperties

path = '/home/nice/Documentos/data_file/cmip_data/cmip5/hadgem2-es_rclimdex/historical'
arq  = '{0}/prcptotETCCDI_yr_HadGEM2-ES_historical_r1i1p1_1859-2005.nc'.format(path)

data = netCDF4.Dataset(arq)
var  = data.variables['prcptotETCCDI'][:] 
lat  = data.variables['lat'][:]
lon  = data.variables['lon'][:]
idx  = np.nanmean(np.nanmean(var[:,:,:], axis=1), axis=1) 

fig, ax = plt.subplots(figsize=(12, 8))

time = np.arange(0.5, 30 + 0.5)

x = [i for i in range(0, len(idx[116:146]))]
x = np.reshape(x, (len(x), 1))
y = idx[116:146]
model = LinearRegression()
model.fit(x, y)
trend = model.predict(x)
z = np.polyfit(y.flatten(), trend.flatten(), 1)

r2 = r2_score(y, trend)

plt_idx = plt.plot(time, y, marker='o', lw=1.5, color='darkblue')
plt_trd = plt.plot(time, trend, lw=1.5, color='red')
plt.title('Annual Total Precipitation in Wet Days - MATOPIBA \n R2={0} y={1}x+{2}'.format(round(r2,2),round(z[0],2),round(z[1],2)), fontsize = 20) 
plt.xlabel('Year', fontsize = 15)
plt.ylabel('PRCPTOT (mm)', fontsize = 15)

plt.xticks(time, np.arange(1976, 2006), rotation=90, fontsize=15)
plt.yticks(np.arange(824, 860, 2), fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=15, length=6, width=1, pad=4)
plt.legend([u'PRCPTOT', u'Trend'], loc=4, shadow=True, ncol=1, prop=FontProperties(size=15))
plt.grid()

plt.savefig('test.png', bbox_inches='tight')
plt.show()
exit()






