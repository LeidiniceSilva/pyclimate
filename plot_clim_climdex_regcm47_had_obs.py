# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/01/2021"
__description__ = "This script plot trends from extremes index"

import os
import netCDF4
import statistics
import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from matplotlib.font_manager import FontProperties


def import_obs(var, area, freq,  model, exp, dt):
	
	path = '/home/nice/Documents/hadgem2-es_rclimdex/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_r1i1p1_{6}.nc'.format(path, var, area, freq,  model, exp, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs  = np.nanmean(np.nanmean(var[:-1,:,:], axis=1), axis=1) 
	
	return obs
	
	
def import_rcm(var, area, freq,  model, exp, dt):
	
	path = '/home/nice/Documents/hadgem2-es_rclimdex/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_r1i1p1_{6}.nc'.format(path, var, area, freq,  model, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = np.nanmean(np.nanmean(var[:-1,:,:], axis=1), axis=1) 

	return rcm


def import_gcm(var, area, freq,  model, exp, dt):
	
	path = '/home/nice/Documents/hadgem2-es_rclimdex/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_r1i1p1_{6}.nc'.format(path, var, area, freq,  model, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	gcm  = np.nanmean(np.nanmean(var[:-1,:,:], axis=1), axis=1) 

	return gcm
	

# Import regcm exp and cru databases 	
obs_prcptot = import_obs('prcptotETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')   
rcm_prcptot = import_rcm('prcptotETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
gcm_prcptot = import_gcm('prcptotETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

obs_r95p = import_obs('r95pETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')   
rcm_r95p = import_rcm('r95pETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
gcm_r95p = import_gcm('r95pETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

obs_r99p = import_obs('r99pETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')   
rcm_r99p = import_rcm('r99pETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
gcm_r99p = import_gcm('r99pETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

obs_rx1day = import_obs('rx1dayETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')   
rcm_rx1day = import_rcm('rx1dayETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
gcm_rx1day = import_gcm('rx1dayETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

obs_rx5day = import_obs('rx5dayETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')   
rcm_rx5day = import_rcm('rx5dayETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
gcm_rx5day = import_gcm('rx5dayETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')

obs_sdii = import_obs('sdiiETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')   
rcm_sdii = import_rcm('sdiiETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')
gcm_sdii = import_gcm('sdiiETCCDI', 'amz_neb', 'yr', 'HadGEM2-ES', 'historical', '1986-2005')


# Import hadgem rclimdex cmip5 	 
# prcptot  
prcptot_x = [i for i in range(0, len(obs_prcptot))]
prcptot_x = np.reshape(prcptot_x, (len(prcptot_x), 1))
prcptot_y = obs_prcptot
model = LinearRegression()
model.fit(prcptot_x, prcptot_y)
prcptot_trend = model.predict(prcptot_x)
prcptot_z = np.polyfit(prcptot_y.flatten(), prcptot_trend.flatten(), 1)
prcptot_r2 = r2_score(prcptot_y, prcptot_trend)
prcptot_median = statistics.median(prcptot_y)

# r95p  
r95p_x = [i for i in range(0, len(obs_r95p))]
r95p_x = np.reshape(r95p_x, (len(r95p_x), 1))
r95p_y = obs_r95p
model = LinearRegression()
model.fit(r95p_x, r95p_y)
r95p_trend = model.predict(r95p_x)
r95p_z = np.polyfit(r95p_y.flatten(), r95p_trend.flatten(), 1)
r95p_r2 = r2_score(r95p_y, r95p_trend)
r95p_median = statistics.median(r95p_y)

# r99p  
r99p_x = [i for i in range(0, len(obs_r99p))]
r99p_x = np.reshape(r99p_x, (len(r99p_x), 1))
r99p_y = obs_r99p
model = LinearRegression()
model.fit(r99p_x, r99p_y)
r99p_trend = model.predict(r99p_x)
r99p_z = np.polyfit(r99p_y.flatten(), r99p_trend.flatten(), 1)
r99p_r2 = r2_score(r99p_y, r99p_trend)
r99p_median = statistics.median(r99p_y)

#~ # r99p  
#~ rx1day_x = [i for i in range(0, len(obs_rx1day))]
#~ rx1day_x = np.reshape(rx1day_x, (len(rx1day_x), 1))
#~ rx1day_y = obs_rx1day
#~ model = LinearRegression()
#~ model.fit(rx1day_x, rx1dayp_y)
#~ rx1day_trend = model.predict(rx1day_x)
#~ rx1day_z = np.polyfit(rx1day_y.flatten(), rx1day_trend.flatten(), 1)
#~ rx1day_r2 = r2_score(rx1day_y, rx1day_trend)
#~ rx1day_median = statistics.median(rx1day_y)

#~ # r99p  
#~ rx5day_x = [i for i in range(0, len(obs_r99p))]
#~ r99p_x = np.reshape(r99p_x, (len(r99p_x), 1))
#~ r99p_y = obs_r99p
#~ model = LinearRegression()
#~ model.fit(r99p_x, r99p_y)
#~ r99p_trend = model.predict(r99p_x)
#~ r99p_z = np.polyfit(r99p_y.flatten(), r99p_trend.flatten(), 1)
#~ r99p_r2 = r2_score(r99p_y, r99p_trend)
#~ r99p_median = statistics.median(r99p_y)

#~ # r99p  
#~ r99p_x = [i for i in range(0, len(obs_r99p))]
#~ r99p_x = np.reshape(r99p_x, (len(r99p_x), 1))
#~ r99p_y = obs_r99p
#~ model = LinearRegression()
#~ model.fit(r99p_x, r99p_y)
#~ r99p_trend = model.predict(r99p_x)
#~ r99p_z = np.polyfit(r99p_y.flatten(), r99p_trend.flatten(), 1)
#~ r99p_r2 = r2_score(r99p_y, r99p_trend)
#~ r99p_median = statistics.median(r99p_y)

# Plot maps hadgem2-es model 
fig = plt.figure()
time = np.arange(1, 19 + 1)
objects = [u'1986',u'1987',u'1988',u'1989',u'1990',u'1991',u'1992',u'1993',u'1994',u'1995',u'1996',u'1997',u'1998',u'1999',u'2000',u'2001',u'2002',u'2003',u'2004',u'2005']

ax1 = fig.add_subplot(611)
ax1.plot(time, prcptot_y, marker='o', lw=1.5, color='darkblue')
ax1.plot(time, prcptot_trend, lw=1.5, color='red')
ax1.axhline(prcptot_median, lw=1.5, linestyle='dashed', color='gray', alpha=2)
plt.title('A) prcptot', loc='left', fontsize=10, fontweight='bold')
#~ plt.yticks(np.arange(126, 129.5, 0.5))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.text(0.5, 128.5, u'R2={0} y={1}x+{2}'.format(round(prcptot_r2,2),round(prcptot_z[0],2),round(prcptot_z[1],2)), fontweight='bold')

ax2 = fig.add_subplot(612)
ax2.plot(time, r95p_y, marker='o', lw=1.5, color='darkblue')
ax2.plot(time, r95p_trend, lw=1.5, color='red')
ax2.axhline(r95p_median, lw=1.5, linestyle='dashed', color='gray', alpha=2)
ax2.set_ylabel('Number of days', fontweight='bold')
plt.title('B) r95p', loc='left', fontsize=10, fontweight='bold')
#~ plt.yticks(np.arange(22, 23.5, 0.25))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.text(0.5, 23, u'R2={0} y={1}x+{2}'.format(round(r95p_r2,2),round(r95p_z[0],2),round(r95p_z[1],2)), fontweight='bold')

ax3 = fig.add_subplot(613)
ax3.plot(time, r99p_y, marker='o', lw=1.5, color='darkblue')
ax3.plot(time, r99p_trend, lw=1.5, color='red')
ax3.axhline(r99p_median, lw=1.5, linestyle='dashed', color='gray', alpha=2)
ax3.set_xlabel('Years', fontweight='bold')
plt.title('C) r99p', loc='left', fontsize=10, fontweight='bold')
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.xticks(np.arange(0.5, 32, 5), objects)
#~ plt.yticks(np.arange(7, 8.5, 0.251))
#~ plt.text(0.5, 8, u'R2={0} y={1}x+{2}'.format(round(r99p_r2,2),round(r99p_z[0],2),round(r99p_z[1],2)), fontweight='bold')

ax4 = fig.add_subplot(614)
ax4.plot(time, prcptot_y, marker='o', lw=1.5, color='darkblue')
ax4.plot(time, prcptot_trend, lw=1.5, color='red')
ax4.axhline(prcptot_median, lw=1.5, linestyle='dashed', color='gray', alpha=2)
plt.title('D) prcptot', loc='left', fontsize=10, fontweight='bold')
#~ plt.yticks(np.arange(126, 129.5, 0.5))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.text(0.5, 128.5, u'R2={0} y={1}x+{2}'.format(round(prcptot_r2,2),round(prcptot_z[0],2),round(prcptot_z[1],2)), fontweight='bold')

ax5 = fig.add_subplot(615)
ax5.plot(time, r95p_y, marker='o', lw=1.5, color='darkblue')
ax5.plot(time, r95p_trend, lw=1.5, color='red')
ax5.axhline(r95p_median, lw=1.5, linestyle='dashed', color='gray', alpha=2)
ax5.set_ylabel('Number of days', fontweight='bold')
plt.title('E) r95p', loc='left', fontsize=10, fontweight='bold')
#~ plt.yticks(np.arange(22, 23.5, 0.25))
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#~ plt.text(0.5, 23, u'R2={0} y={1}x+{2}'.format(round(r95p_r2,2),round(r95p_z[0],2),round(r95p_z[1],2)), fontweight='bold')

ax6 = fig.add_subplot(616)
ax6.plot(time, r99p_y, marker='o', lw=1.5, color='darkblue')
ax6.plot(time, r99p_trend, lw=1.5, color='red')
ax6.axhline(r99p_median, lw=1.5, linestyle='dashed', color='gray', alpha=2)
ax6.set_xlabel('Years', fontweight='bold')
plt.title('F) r99p', loc='left', fontsize=10, fontweight='bold')
plt.xticks(np.arange(1, 19, 1), objects)
#~ plt.yticks(np.arange(7, 8.5, 0.251))
#~ plt.text(0.5, 8, u'R2={0} y={1}x+{2}'.format(round(r99p_r2,2),round(r99p_z[0],2),round(r99p_z[1],2)), fontweight='bold')

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.10, hspace=0.60)

path_out = '/home/nice/Downloads'
name_out = 'pyplt_climdex_pre_regcm47_had_obs_1986-2005.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()






