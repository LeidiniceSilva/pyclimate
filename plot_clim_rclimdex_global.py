# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "06/03/2019"
__description__ = "This script plot Rclimdex based in HadGEM2-ES model CMIP5"

import os
import netCDF4
import statistics
import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from matplotlib.font_manager import FontProperties


def import_rclimdex(variable):
	
	path = '/home/nice/Documentos/data_file/cmip_data/cmip5/hadgem2-es_rclimdex/historical'
	arq  = '{0}/{1}_yr_HadGEM2-ES_historical_r1i1p1_1859-2005.nc'.format(path, variable)

	data = netCDF4.Dataset(arq)
	var  = data.variables[variable][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	idx  = np.nanmean(np.nanmean(var[:,:,:], axis=1), axis=1) 

	return idx
	

# Import hadgem rclimdex cmip5 	 
# r1mm  
r1mm = import_rclimdex('r1mmETCCDI')
r1mm_x = [i for i in range(0, len(r1mm[116:146]))]
r1mm_x = np.reshape(r1mm_x, (len(r1mm_x), 1))
r1mm_y = r1mm[116:146]
model = LinearRegression()
model.fit(r1mm_x, r1mm_y)
r1mm_trend = model.predict(r1mm_x)
r1mm_z = np.polyfit(r1mm_y.flatten(), r1mm_trend.flatten(), 1)
r1mm_r2 = r2_score(r1mm_y, r1mm_trend)
r1mm_median = statistics.median(r1mm_y)

# r10mm  
r10mm = import_rclimdex('r10mmETCCDI')
r10mm_x = [i for i in range(0, len(r10mm[116:146]))]
r10mm_x = np.reshape(r10mm_x, (len(r10mm_x), 1))
r10mm_y = r10mm[116:146]
model = LinearRegression()
model.fit(r10mm_x, r10mm_y)
r10mm_trend = model.predict(r10mm_x)
r10mm_z = np.polyfit(r10mm_y.flatten(), r10mm_trend.flatten(), 1)
r10mm_r2 = r2_score(r10mm_y, r10mm_trend)
r10mm_median = statistics.median(r10mm_y)

# r20mm  
r20mm = import_rclimdex('r20mmETCCDI')
r20mm_x = [i for i in range(0, len(r20mm[116:146]))]
r20mm_x = np.reshape(r20mm_x, (len(r20mm_x), 1))
r20mm_y = r20mm[116:146]
model = LinearRegression()
model.fit(r20mm_x, r20mm_y)
r20mm_trend = model.predict(r20mm_x)
r20mm_z = np.polyfit(r20mm_y.flatten(), r20mm_trend.flatten(), 1)
r20mm_r2 = r2_score(r20mm_y, r20mm_trend)
r20mm_median = statistics.median(r20mm_y)

# Plot maps hadgem2-es model 
fig = plt.figure()
time = np.arange(1, 30 + 1)
objects = [u'1975', u'1980', u'1985', u'1990', u'1995', u'2000', u'2005', u'2010']

ax1 = fig.add_subplot(311)
ax1.plot(time, r1mm_y, marker='o', lw=1.5, color='darkblue')
ax1.plot(time, r1mm_trend, lw=1.5, color='red')
ax1.axhline(r1mm_median, lw=1.5, linestyle='dashed', color='gray', alpha=2)
plt.title('A) R1mm', loc='left', fontsize=10, fontweight='bold')
plt.yticks(np.arange(126, 129.5, 0.5))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.text(0.5, 128.5, u'R2={0} y={1}x+{2}'.format(round(r1mm_r2,2),round(r1mm_z[0],2),round(r1mm_z[1],2)), fontweight='bold')

ax2 = fig.add_subplot(312)
ax2.plot(time, r10mm_y, marker='o', lw=1.5, color='darkblue')
ax2.plot(time, r10mm_trend, lw=1.5, color='red')
ax2.axhline(r10mm_median, lw=1.5, linestyle='dashed', color='gray', alpha=2)
ax2.set_ylabel('Number of days', fontweight='bold')
plt.title('B) R10mm', loc='left', fontsize=10, fontweight='bold')
plt.yticks(np.arange(22, 23.5, 0.25))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.text(0.5, 23, u'R2={0} y={1}x+{2}'.format(round(r10mm_r2,2),round(r10mm_z[0],2),round(r10mm_z[1],2)), fontweight='bold')

ax3 = fig.add_subplot(313)
ax3.plot(time, r20mm_y, marker='o', lw=1.5, color='darkblue')
ax3.plot(time, r20mm_trend, lw=1.5, color='red')
ax3.axhline(r20mm_median, lw=1.5, linestyle='dashed', color='gray', alpha=2)
ax3.set_xlabel('Years', fontweight='bold')
plt.title('C) R20mm', loc='left', fontsize=10, fontweight='bold')
plt.xticks(np.arange(0.5, 32, 5), objects)
plt.yticks(np.arange(7, 8.5, 0.251))
plt.text(0.5, 8, u'R2={0} y={1}x+{2}'.format(round(r20mm_r2,2),round(r20mm_z[0],2),round(r20mm_z[1],2)), fontweight='bold')

path_out = '/home/nice/Documentos/ufrn/PhD_project/results/rclimdex'
name_out = 'pyplt_clim_prcptotETCCDI_rclimdex_hadgem_1976-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')

plt.show()
exit()






