# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/08/2019"
__description__ = "This script plot taylor diagram from Reg and Had models end obs database to ETCCDI"

import os
import netCDF4
import numpy as np
import numpy.ma as ma
import scipy.stats as st
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as FA
import mpl_toolkits.axisartist.grid_finder as GF

from matplotlib.projections import PolarAxes


class TaylorDiagram(object):
    """
    Taylor diagram.
    Plot model standard deviation and correlation to reference (data)
    sample in a single-quadrant polar plot, with r=stddev and
    theta=arccos(correlation).
    """

    def __init__(self, refstd, fig=None, rect=313, label='_', marker='', color='', srange=(0., 3.), extend=False):
        """
        Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using `mpl_toolkits.axisartist.floating_axes`.
        Parameters:
        * refstd: reference standard deviation to be compared to
        * fig: input Figure or None
        * rect: subplot definition
        * label: reference label
        * srange: stddev axis extension, in units of *refstd*
        * extend: extend diagram to negative correlations
        """

        self.refstd = refstd    
               
        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = np.array([0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.98, 1])
        
        if extend:
            # Diagram extended to negative correlations
            self.tmax = np.pi
            rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
        else:
            # Diagram limited to positive correlations
            self.tmax = np.pi/2
            
        tlocs = np.arccos(rlocs)        # Conversion to polar angles
        gl1 = GF.FixedLocator(tlocs)    # Positions
        tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        # Standard deviation axis extent (in units of reference stddev)
        self.smin = srange[0] * self.refstd
        self.smax = srange[1] * self.refstd
        

        ghelper = FA.GridHelperCurveLinear(tr,
                                           extremes=(0,self.tmax, # 1st quadrant
                                                     self.smin,self.smax),
                                           grid_locator1=gl1,
                                           tick_formatter1=tf1,
                                           )
        
        if fig is None:
            fig = plt.figure()
        
        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")   # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text(u'Correlation')

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].toggle(ticklabels=True, label=True)
        ax.axis["left"].major_ticklabels.set_rotation(-90)
        ax.axis["left"].label.set_pad(10) 
        ax.axis["left"].label.set_text(u'  Standard Deviation')
        		
        ax.axis["right"].set_axis_direction("top")    # "Y-axis"
        ax.axis["right"].toggle(ticklabels=True, label=True)
        ax.axis["right"].major_ticklabels.set_rotation(90)
        ax.axis["right"].major_ticklabels.set_pad(12)
        ax.axis["bottom"].set_visible(False)      # Unused

        ax.grid(color='k', axis='x', linestyle='--', linewidth=1)
        
        #~ if self.smin:
            #~ ax.axis["bottom"].toggle(ticklabels=True, label=True)
        #~ else:

        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        # Add reference point and stddev contour
        l, = self.ax.plot([0], self.refstd, 'k*', ls='', ms=8, label=label)
        t = np.linspace(0, self.tmax)
        r = np.zeros_like(t) + self.refstd
        self.ax.plot(t, r, 'k--', label='_')

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [l]


    def add_sample(self, stddev, corrcoef, *args, **kwargs):
        """
        Add sample (*stddev*, *corrcoeff*) to the Taylor
        diagram. *args* and *kwargs* are directly propagated to the
        'Figure.plot' command.
        """
	
        l, = self.ax.plot(np.arccos(corrcoef), stddev, *args, **kwargs)  # (theta, radius)
        self.samplePoints.append(l)
   
        return l
        

    def add_grid(self, *args, **kwargs):
        """Add a grid."""
        
        self._ax.grid(*args, **kwargs)
        

    def add_contours(self, levels=5, **kwargs):
        """
        Add constant centered RMS difference contours, defined by *levels*.
        """
        
        rs, ts = np.meshgrid(np.linspace(self.smin, self.smax), np.linspace(0, self.tmax))
        rms = np.sqrt(self.refstd**2 + rs**2 - 2*self.refstd*rs*np.cos(ts))
        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)
        
        return contours


def import_obs(var, area, dataset, freq, dt):
	
	path = '/home/nice/Documents/dataset/obs/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_lonlat.nc'.format(path, var, area, dataset, freq, dt)	

	dict_var = {u'eca_prectot': u'precip', 
	u'eca_r95p': u'precip',
	u'eca_r99p': u'precip', 
	u'eca_rx1day': u'highest_one_day_precipitation_amount_per_time_period',
	u'eca_rx5day': u'highest_five_day_precipitation_amount_per_time_period',
	u'eca_sdii': u'simple_daily_intensitiy_index_per_time_period',
	u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}

	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	annual_obs = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return annual_obs
	
	
def import_rcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/rcm/eca'	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_prectot': u'pr', 
	u'eca_r95p': u'pr',
	u'eca_r99p': u'pr', 
	u'eca_rx1day': u'highest_one_day_precipitation_amount_per_time_period',
	u'eca_rx5day': u'highest_five_day_precipitation_amount_per_time_period',
	u'eca_sdii': u'simple_daily_intensitiy_index_per_time_period',
	u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	annual_rcm = np.nanmean(np.nanmean(var[:][0:20,:,:], axis=1), axis=1)
	
	return annual_rcm


def import_gcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/gcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_prectot': u'pr', 
	u'eca_r95p': u'pr',
	u'eca_r99p': u'pr', 
	u'eca_rx1day': u'highest_one_day_precipitation_amount_per_time_period',
	u'eca_rx5day': u'highest_five_day_precipitation_amount_per_time_period',
	u'eca_sdii': u'simple_daily_intensitiy_index_per_time_period',
	u'eca_cdd': u'consecutive_dry_days_index_per_time_period', 
	u'eca_cwd': u'consecutive_wet_days_index_per_time_period',
	u'eca_r10mm': u'heavy_precipitation_days_index_per_time_period', 
	u'eca_r20mm': u'very_heavy_precipitation_days_index_per_time_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	annual_gcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	return annual_gcm


if __name__=='__main__':
	
	# Import regcm exp and cru databases 
	obs_prcptot = import_obs('eca_prectot', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_prcptot = import_rcm('eca_prectot', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_prcptot = import_gcm('eca_prectot', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_r95p = import_obs('eca_r95p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r95p = import_rcm('eca_r95p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r95p = import_gcm('eca_r95p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_r99p = import_obs('eca_r99p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r99p = import_rcm('eca_r99p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r99p = import_gcm('eca_r99p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_rx1day = import_obs('eca_rx1day', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_rx1day = import_rcm('eca_rx1day', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_rx1day = import_gcm('eca_rx1day', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_rx5day = import_obs('eca_rx5day', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_rx5day = import_rcm('eca_rx5day', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_rx5day = import_gcm('eca_rx5day', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_sdii = import_obs('eca_sdii', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_sdii = import_rcm('eca_sdii', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_sdii = import_gcm('eca_sdii', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_cdd = import_obs('eca_cdd', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_cdd = import_rcm('eca_cdd', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_cdd = import_gcm('eca_cdd', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_cwd = import_obs('eca_cwd', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_cwd = import_rcm('eca_cwd', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_cwd = import_gcm('eca_cwd', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_r10mm = import_obs('eca_r10mm', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r10mm = import_rcm('eca_r10mm', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r10mm = import_gcm('eca_r10mm', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_r20mm = import_obs('eca_r20mm', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r20mm = import_rcm('eca_r20mm', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r20mm = import_gcm('eca_r20mm', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	# Reference database standard desviation
	stdrefs = dict(SAMZ = 1,
				 ENEB = 1,
				 MATOPIBA = 1)       

	text1 = dict(SAMZ='A)',
				 ENEB='B)',
				 MATOPIBA='C)')
				 
	# Compute stddev and correlation coefficient of models
	# Sample std, rho: Be sure to check order and that correct numbers are placed!
	samples = dict(SAMZ=[[obs_prcptot.std(ddof=1), np.corrcoef(obs_prcptot, rcm_prcptot)[0,1], 'PRCPTOT-Reg', 'o', 'b'],
                          [obs_r95p.std(ddof=1), np.corrcoef(obs_r95p, rcm_r95p)[0,1], 'R95p-Reg', 's', 'b'],
                          [obs_r99p.std(ddof=1), np.corrcoef(obs_r99p, rcm_r99p)[0,1], 'R99p-Reg', 'v', 'b'],
                          [obs_rx1day.std(ddof=1), np.corrcoef(obs_rx1day, rcm_rx1day)[0,1], 'Rx1day-Reg', '^', 'b'],
                          [obs_rx5day.std(ddof=1), np.corrcoef(obs_rx5day, rcm_rx5day)[0,1], 'Rx5day-Reg', 'd', 'b'],
                          [obs_sdii.std(ddof=1), np.corrcoef(obs_sdii, rcm_sdii)[0,1], 'SDII-Reg', 'h', 'b'],
                          [obs_cdd.std(ddof=1), np.corrcoef(obs_cdd, rcm_cdd)[0,1], 'CDD-Reg', '*', 'b'],
                          [obs_cwd.std(ddof=1), np.corrcoef(obs_cwd, rcm_cwd)[0,1], 'CWD-Reg', 'p', 'b'],
                          [obs_r10mm.std(ddof=1), np.corrcoef(obs_r10mm, rcm_r10mm)[0,1], 'R10mm-Reg', '>', 'b'],
                          [obs_r20mm.std(ddof=1), np.corrcoef(obs_r20mm, rcm_r20mm)[0,1], 'R20mm-Reg', '<', 'b'],
                          [obs_prcptot.std(ddof=1), np.corrcoef(obs_prcptot, gcm_prcptot)[0,1], 'PRCPTOT-Had', 'o', 'r'],
                          [obs_r95p.std(ddof=1), np.corrcoef(obs_r95p, gcm_r95p)[0,1], 'R95p-Had', 's', 'r'],
                          [obs_r99p.std(ddof=1), np.corrcoef(obs_r99p, gcm_r99p)[0,1], 'R99p-Had', 'v', 'r'],
                          [obs_rx1day.std(ddof=1), np.corrcoef(obs_rx1day, gcm_rx1day)[0,1], 'Rx1day-Had', '^', 'r'],
                          [obs_rx5day.std(ddof=1), np.corrcoef(obs_rx5day, gcm_rx5day)[0,1], 'Rx5day-Had', 'd', 'r'],
                          [obs_sdii.std(ddof=1), np.corrcoef(obs_sdii, gcm_sdii)[0,1], 'SDII-Had', 'h', 'r'],
                          [obs_cdd.std(ddof=1), np.corrcoef(obs_cdd, gcm_cdd)[0,1], 'CDD-Had', '*', 'r'],
                          [obs_cwd.std(ddof=1), np.corrcoef(obs_cwd, gcm_cwd)[0,1], 'CWD-Had', 'p', 'r'],
                          [obs_r10mm.std(ddof=1), np.corrcoef(obs_r10mm, gcm_r10mm)[0,1], 'R10mm-Had', '>', 'r'],
                          [obs_r20mm.std(ddof=1), np.corrcoef(obs_r20mm, gcm_r20mm)[0,1], 'R20mm-Had', '<', 'r']],
                   ENEB=[[obs_prcptot.std(ddof=1), np.corrcoef(obs_prcptot, rcm_prcptot)[0,1], 'PRCPTOT-Reg', 'o', 'b'],
                          [obs_r95p.std(ddof=1), np.corrcoef(obs_r95p, rcm_r95p)[0,1], 'R95p-Reg', 's', 'b'],
                          [obs_r99p.std(ddof=1), np.corrcoef(obs_r99p, rcm_r99p)[0,1], 'R99p-Reg', 'v', 'b'],
                          [obs_rx1day.std(ddof=1), np.corrcoef(obs_rx1day, rcm_rx1day)[0,1], 'Rx1day-Reg', '^', 'b'],
                          [obs_rx5day.std(ddof=1), np.corrcoef(obs_rx5day, rcm_rx5day)[0,1], 'Rx5day-Reg', 'd', 'b'],
                          [obs_sdii.std(ddof=1), np.corrcoef(obs_sdii, rcm_sdii)[0,1], 'SDII-Reg', 'h', 'b'],
                          [obs_cdd.std(ddof=1), np.corrcoef(obs_cdd, rcm_cdd)[0,1], 'CDD-Reg', '*', 'b'],
                          [obs_cwd.std(ddof=1), np.corrcoef(obs_cwd, rcm_cwd)[0,1], 'CWD-Reg', 'p', 'b'],
                          [obs_r10mm.std(ddof=1), np.corrcoef(obs_r10mm, rcm_r10mm)[0,1], 'R10mm-Reg', '>', 'b'],
                          [obs_r20mm.std(ddof=1), np.corrcoef(obs_r20mm, rcm_r20mm)[0,1], 'R20mm-Reg', '<', 'b'],
                          [obs_prcptot.std(ddof=1), np.corrcoef(obs_prcptot, gcm_prcptot)[0,1], 'PRCPTOT-Had', 'o', 'r'],
                          [obs_r95p.std(ddof=1), np.corrcoef(obs_r95p, gcm_r95p)[0,1], 'R95p-Had', 's', 'r'],
                          [obs_r99p.std(ddof=1), np.corrcoef(obs_r99p, gcm_r99p)[0,1], 'R99p-Had', 'v', 'r'],
                          [obs_rx1day.std(ddof=1), np.corrcoef(obs_rx1day, gcm_rx1day)[0,1], 'Rx1day-Had', '^', 'r'],
                          [obs_rx5day.std(ddof=1), np.corrcoef(obs_rx5day, gcm_rx5day)[0,1], 'Rx5day-Had', 'd', 'r'],
                          [obs_sdii.std(ddof=1), np.corrcoef(obs_sdii, gcm_sdii)[0,1], 'SDII-Had', 'h', 'r'],
                          [obs_cdd.std(ddof=1), np.corrcoef(obs_cdd, gcm_cdd)[0,1], 'CDD-Had', '*', 'r'],
                          [obs_cwd.std(ddof=1), np.corrcoef(obs_cwd, gcm_cwd)[0,1], 'CWD-Had', 'p', 'r'],
                          [obs_r10mm.std(ddof=1), np.corrcoef(obs_r10mm, gcm_r10mm)[0,1], 'R10mm-Had', '>', 'r'],
                          [obs_r20mm.std(ddof=1), np.corrcoef(obs_r20mm, gcm_r20mm)[0,1], 'R20mm-Had', '<', 'r']], 
               MATOPIBA=[[obs_prcptot.std(ddof=1), np.corrcoef(obs_prcptot, rcm_prcptot)[0,1], 'PRCPTOT-Reg', 'o', 'b'],
                          [obs_r95p.std(ddof=1), np.corrcoef(obs_r95p, rcm_r95p)[0,1], 'R95p-Reg', 's', 'b'],
                          [obs_r99p.std(ddof=1), np.corrcoef(obs_r99p, rcm_r99p)[0,1], 'R99p-Reg', 'v', 'b'],
                          [obs_rx1day.std(ddof=1), np.corrcoef(obs_rx1day, rcm_rx1day)[0,1], 'Rx1day-Reg', '^', 'b'],
                          [obs_rx5day.std(ddof=1), np.corrcoef(obs_rx5day, rcm_rx5day)[0,1], 'Rx5day-Reg', 'd', 'b'],
                          [obs_sdii.std(ddof=1), np.corrcoef(obs_sdii, rcm_sdii)[0,1], 'SDII-Reg', 'h', 'b'],
                          [obs_cdd.std(ddof=1), np.corrcoef(obs_cdd, rcm_cdd)[0,1], 'CDD-Reg', '*', 'b'],
                          [obs_cwd.std(ddof=1), np.corrcoef(obs_cwd, rcm_cwd)[0,1], 'CWD-Reg', 'p', 'b'],
                          [obs_r10mm.std(ddof=1), np.corrcoef(obs_r10mm, rcm_r10mm)[0,1], 'R10mm-Reg', '>', 'b'],
                          [obs_r20mm.std(ddof=1), np.corrcoef(obs_r20mm, rcm_r20mm)[0,1], 'R20mm-Reg', '<', 'b'],
                          [obs_prcptot.std(ddof=1), np.corrcoef(obs_prcptot, gcm_prcptot)[0,1], 'PRCPTOT-Had', 'o', 'r'],
                          [obs_r95p.std(ddof=1), np.corrcoef(obs_r95p, gcm_r95p)[0,1], 'R95p-Had', 's', 'r'],
                          [obs_r99p.std(ddof=1), np.corrcoef(obs_r99p, gcm_r99p)[0,1], 'R99p-Had', 'v', 'r'],
                          [obs_rx1day.std(ddof=1), np.corrcoef(obs_rx1day, gcm_rx1day)[0,1], 'Rx1day-Had', '^', 'r'],
                          [obs_rx5day.std(ddof=1), np.corrcoef(obs_rx5day, gcm_rx5day)[0,1], 'Rx5day-Had', 'd', 'r'],
                          [obs_sdii.std(ddof=1), np.corrcoef(obs_sdii, gcm_sdii)[0,1], 'SDII-Had', 'h', 'r'],
                          [obs_cdd.std(ddof=1), np.corrcoef(obs_cdd, gcm_cdd)[0,1], 'CDD-Had', '*', 'r'],
                          [obs_cwd.std(ddof=1), np.corrcoef(obs_cwd, gcm_cwd)[0,1], 'CWD-Had', 'p', 'r'],
                          [obs_r10mm.std(ddof=1), np.corrcoef(obs_r10mm, gcm_r10mm)[0,1], 'R10mm-Had', '>', 'r'],
                          [obs_r20mm.std(ddof=1), np.corrcoef(obs_r20mm, gcm_r20mm)[0,1], 'R20mm-Had', '<', 'r']])

	# Plot Taylor Diagram
	rects = dict(SAMZ=311,
				 ENEB=312,
				 MATOPIBA=313)

	# Plot model end obs data taylor diagram 			 
	fig = plt.figure(figsize=(6, 6.5))
	
	for var in ['SAMZ', 'ENEB', 'MATOPIBA']:

		dia = TaylorDiagram(stdrefs[var], fig=fig, rect=rects[var], label='Reference', srange=(0., 3.), extend=True)
		dia.samplePoints[0].set_color('black')
		
		# Add samples to Taylor diagram

		for i, (stddev,corrcoef,name,mark,cor) in enumerate(samples[var]):
			dia.add_sample(stddev, corrcoef,
						   label=name, marker=mark, color='black', mfc=cor, ms=8, ls='')		
			plt.text(-3., 2.8, text1[var], fontweight='bold')

		# Add RMS contours, and label them
		contours = dia.add_contours(levels=5, colors='0.5')
		plt.clabel(contours, inline=1, fontsize=8, fmt='%.1f')

	# Add a figure legend
	fig.legend(dia.samplePoints, 
			   [ p.get_label() for p in dia.samplePoints ], 
			   prop=dict(size=10), ncol=1, numpoints=1, loc=7)
	
	plt.subplots_adjust(left=0.10, bottom=0.10, right=0.90, top=0.90, wspace=0.65, hspace=0.65)

	# Path out to save figure
	path_out = '/home/nice/Downloads'
	name_out = 'pyplt_taylor_diagram_etccdi_pre_reg_had_obs_1986-2005.png'
	if not os.path.exists(path_out):
		create_path(path_out)
	plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
	plt.show()
	exit()
