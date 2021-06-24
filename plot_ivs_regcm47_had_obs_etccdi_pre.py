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

    def __init__(self, refstd, fig=None, rect=313, label='_', marker='', color='', srange=(0., 9), extend=False):
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
        rlocs = np.array([0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.98, 1])
        
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
        #~ ax.axis["left"].toggle(ticklabels=True, label=True)
        #~ ax.axis["left"].major_ticklabels.set_rotation(-90)
        #~ ax.axis["left"].label.set_pad(10) 
        ax.axis["left"].label.set_text(u'  Standard Deviation')
        		
        ax.axis["right"].set_axis_direction("top")    # "Y-axis"
        ax.axis["right"].toggle(ticklabels=True, label=True)
        #~ ax.axis["right"].major_ticklabels.set_rotation(90)
        #~ ax.axis["right"].major_ticklabels.set_pad(12)
        ax.axis["right"].major_ticklabels.set_axis_direction(
			"bottom" if extend else "left")

        if self.smin:
            ax.axis["bottom"].toggle(ticklabels=False, label=False)
        else:
            ax.axis["bottom"].set_visible(False)      # Unused

        ax.grid(color='k', axis='x', linestyle='--', linewidth=1)
        
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

	dict_var = {u'eca_prcptot': u'precip', 
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

	dict_var = {u'eca_prcptot': u'pr', 
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

	dict_var = {u'eca_prcptot': u'pr', 
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
	obs_prcptot_samz = import_obs('eca_prcptot', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_prcptot_samz = import_rcm('eca_prcptot', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_prcptot_samz = import_gcm('eca_prcptot', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_r95p_samz = import_obs('eca_r95p', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r95p_samz = import_rcm('eca_r95p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r95p_samz = import_gcm('eca_r95p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_r99p_samz = import_obs('eca_r99p', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r99p_samz = import_rcm('eca_r99p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r99p_samz = import_gcm('eca_r99p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_rx1day_samz = import_obs('eca_rx1day', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_rx1day_samz = import_rcm('eca_rx1day', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_rx1day_samz = import_gcm('eca_rx1day', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_rx5day_samz = import_obs('eca_rx5day', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_rx5day_samz = import_rcm('eca_rx5day', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_rx5day_samz = import_gcm('eca_rx5day', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_sdii_samz = import_obs('eca_sdii', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_sdii_samz = import_rcm('eca_sdii', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_sdii_samz = import_gcm('eca_sdii', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_cdd_samz = import_obs('eca_cdd', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_cdd_samz = import_rcm('eca_cdd', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_cdd_samz = import_gcm('eca_cdd', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_cwd_samz = import_obs('eca_cwd', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_cwd_samz = import_rcm('eca_cwd', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_cwd_samz = import_gcm('eca_cwd', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_r10mm_samz = import_obs('eca_r10mm', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r10mm_samz = import_rcm('eca_r10mm', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r10mm_samz = import_gcm('eca_r10mm', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_r20mm_samz = import_obs('eca_r20mm', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r20mm_samz = import_rcm('eca_r20mm', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r20mm_samz = import_gcm('eca_r20mm', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_prcptot_eneb = import_obs('eca_prcptot', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_prcptot_eneb = import_rcm('eca_prcptot', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_prcptot_eneb = import_gcm('eca_prcptot', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_r95p_eneb = import_obs('eca_r95p', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r95p_eneb = import_rcm('eca_r95p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r95p_eneb = import_gcm('eca_r95p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_r99p_eneb = import_obs('eca_r99p', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r99p_eneb = import_rcm('eca_r99p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r99p_eneb = import_gcm('eca_r99p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_rx1day_eneb = import_obs('eca_rx1day', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_rx1day_eneb = import_rcm('eca_rx1day', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_rx1day_eneb = import_gcm('eca_rx1day', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_rx5day_eneb = import_obs('eca_rx5day', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_rx5day_eneb = import_rcm('eca_rx5day', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_rx5day_eneb = import_gcm('eca_rx5day', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_sdii_eneb = import_obs('eca_sdii', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_sdii_eneb = import_rcm('eca_sdii', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_sdii_eneb = import_gcm('eca_sdii', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_cdd_eneb = import_obs('eca_cdd', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_cdd_eneb = import_rcm('eca_cdd', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_cdd_eneb = import_gcm('eca_cdd', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_cwd_eneb = import_obs('eca_cwd', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_cwd_eneb = import_rcm('eca_cwd', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_cwd_eneb = import_gcm('eca_cwd', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_r10mm_eneb = import_obs('eca_r10mm', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r10mm_eneb = import_rcm('eca_r10mm', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r10mm_eneb = import_gcm('eca_r10mm', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_r20mm_eneb = import_obs('eca_r20mm', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r20mm_eneb = import_rcm('eca_r20mm', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r20mm_eneb = import_gcm('eca_r20mm', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	
	obs_prcptot_matopiba = import_obs('eca_prcptot', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_prcptot_matopiba = import_rcm('eca_prcptot', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_prcptot_matopiba = import_gcm('eca_prcptot', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_r95p_matopiba = import_obs('eca_r95p', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r95p_matopiba = import_rcm('eca_r95p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r95p_matopiba = import_gcm('eca_r95p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_r99p_matopiba = import_obs('eca_r99p', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r99p_matopiba = import_rcm('eca_r99p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r99p_matopiba = import_gcm('eca_r99p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_rx1day_matopiba = import_obs('eca_rx1day', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_rx1day_matopiba = import_rcm('eca_rx1day', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_rx1day_matopiba = import_gcm('eca_rx1day', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_rx5day_matopiba = import_obs('eca_rx5day', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_rx5day_matopiba = import_rcm('eca_rx5day', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_rx5day_matopiba = import_gcm('eca_rx5day', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_sdii_matopiba = import_obs('eca_sdii', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_sdii_matopiba = import_rcm('eca_sdii', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_sdii_matopiba = import_gcm('eca_sdii', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_cdd_matopiba = import_obs('eca_cdd', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_cdd_matopiba = import_rcm('eca_cdd', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_cdd_matopiba = import_gcm('eca_cdd', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_cwd_matopiba = import_obs('eca_cwd', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_cwd_matopiba = import_rcm('eca_cwd', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_cwd_matopiba = import_gcm('eca_cwd', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_r10mm_matopiba = import_obs('eca_r10mm', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r10mm_matopiba = import_rcm('eca_r10mm', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r10mm_matopiba = import_gcm('eca_r10mm', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_r20mm_matopiba = import_obs('eca_r20mm', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_r20mm_matopiba = import_rcm('eca_r20mm', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_r20mm_matopiba = import_gcm('eca_r20mm', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	
	# Reference database standard desviation
	stdrefs = dict(SAMZ = 1,
				 ENEB = 1,
				 MATOPIBA = 1)       

	text1 = dict(SAMZ='A)',
				 ENEB='B)',
				 MATOPIBA='C)')
				 
	# Compute stddev and correlation coefficient of models
	# Sample std, rho: Be sure to check order and that correct numbers are placed!
	samples = dict(SAMZ=[[obs_prcptot_samz.std(ddof=1), np.corrcoef(obs_prcptot_samz, rcm_prcptot_samz)[0,1], 'PRCPTOT-Reg', 'o', 'b'],
                          [obs_r95p_samz.std(ddof=1), np.corrcoef(obs_r95p_samz, rcm_r95p_samz)[0,1], 'R95p-Reg', 's', 'b'],
                          [obs_r99p_samz.std(ddof=1), np.corrcoef(obs_r99p_samz, rcm_r99p_samz)[0,1], 'R99p-Reg', 'v', 'b'],
                          [obs_rx1day_samz.std(ddof=1), np.corrcoef(obs_rx1day_samz, rcm_rx1day_samz)[0,1], 'Rx1day-Reg', '^', 'b'],
                          [obs_rx5day_samz.std(ddof=1), np.corrcoef(obs_rx5day_samz, rcm_rx5day_samz)[0,1], 'Rx5day-Reg', 'd', 'b'],
                          [obs_sdii_samz.std(ddof=1), np.corrcoef(obs_sdii_samz, rcm_sdii_samz)[0,1], 'SDII-Reg', 'h', 'b'],
                          [obs_cdd_samz.std(ddof=1), np.corrcoef(obs_cdd_samz, rcm_cdd_samz)[0,1], 'CDD-Reg', '*', 'b'],
                          [obs_cwd_samz.std(ddof=1), np.corrcoef(obs_cwd_samz, rcm_cwd_samz)[0,1], 'CWD-Reg', 'p', 'b'],
                          [obs_r10mm_samz.std(ddof=1), np.corrcoef(obs_r10mm_samz, rcm_r10mm_samz)[0,1], 'R10mm-Reg', '>', 'b'],
                          [obs_r20mm_samz.std(ddof=1), np.corrcoef(obs_r20mm_samz, rcm_r20mm_samz)[0,1], 'R20mm-Reg', '<', 'b'],
                          [obs_prcptot_samz.std(ddof=1), np.corrcoef(obs_prcptot_samz, gcm_prcptot_samz)[0,1], 'PRCPTOT-Had', 'o', 'r'],
                          [obs_r95p_samz.std(ddof=1), np.corrcoef(obs_r95p_samz, gcm_r95p_samz)[0,1], 'R95p-Had', 's', 'r'],
                          [obs_r99p_samz.std(ddof=1), np.corrcoef(obs_r99p_samz, gcm_r99p_samz)[0,1], 'R99p-Had', 'v', 'r'],
                          [obs_rx1day_samz.std(ddof=1), np.corrcoef(obs_rx1day_samz, gcm_rx1day_samz)[0,1], 'Rx1day-Had', '^', 'r'],
                          [obs_rx5day_samz.std(ddof=1), np.corrcoef(obs_rx5day_samz, gcm_rx5day_samz)[0,1], 'Rx5day-Had', 'd', 'r'],
                          [obs_sdii_samz.std(ddof=1), np.corrcoef(obs_sdii_samz, gcm_sdii_samz)[0,1], 'SDII-Had', 'h', 'r'],
                          [obs_cdd_samz.std(ddof=1), np.corrcoef(obs_cdd_samz, gcm_cdd_samz)[0,1], 'CDD-Had', '*', 'r'],
                          [obs_cwd_samz.std(ddof=1), np.corrcoef(obs_cwd_samz, gcm_cwd_samz)[0,1], 'CWD-Had', 'p', 'r'],
                          [obs_r10mm_samz.std(ddof=1), np.corrcoef(obs_r10mm_samz, gcm_r10mm_samz)[0,1], 'R10mm-Had', '>', 'r'],
                          [obs_r20mm_samz.std(ddof=1), np.corrcoef(obs_r20mm_samz, gcm_r20mm_samz)[0,1], 'R20mm-Had', '<', 'r']],
                   ENEB=[[obs_prcptot_eneb.std(ddof=1), np.corrcoef(obs_prcptot_eneb, rcm_prcptot_eneb)[0,1], 'PRCPTOT-Reg', 'o', 'b'],
                          [obs_r95p_eneb.std(ddof=1), np.corrcoef(obs_r95p_eneb, rcm_r95p_eneb)[0,1], 'R95p-Reg', 's', 'b'],
                          [obs_r99p_eneb.std(ddof=1), np.corrcoef(obs_r99p_eneb, rcm_r99p_eneb)[0,1], 'R99p-Reg', 'v', 'b'],
                          [obs_rx1day_eneb.std(ddof=1), np.corrcoef(obs_rx1day_eneb, rcm_rx1day_eneb)[0,1], 'Rx1day-Reg', '^', 'b'],
                          [obs_rx5day_eneb.std(ddof=1), np.corrcoef(obs_rx5day_eneb, rcm_rx5day_eneb)[0,1], 'Rx5day-Reg', 'd', 'b'],
                          [obs_sdii_eneb.std(ddof=1), np.corrcoef(obs_sdii_eneb, rcm_sdii_eneb)[0,1], 'SDII-Reg', 'h', 'b'],
                          [obs_cdd_eneb.std(ddof=1), np.corrcoef(obs_cdd_eneb, rcm_cdd_eneb)[0,1], 'CDD-Reg', '*', 'b'],
                          [obs_cwd_eneb.std(ddof=1), np.corrcoef(obs_cwd_eneb, rcm_cwd_eneb)[0,1], 'CWD-Reg', 'p', 'b'],
                          [obs_r10mm_eneb.std(ddof=1), np.corrcoef(obs_r10mm_eneb, rcm_r10mm_eneb)[0,1], 'R10mm-Reg', '>', 'b'],
                          [obs_r20mm_eneb.std(ddof=1), np.corrcoef(obs_r20mm_eneb, rcm_r20mm_eneb)[0,1], 'R20mm-Reg', '<', 'b'],
                          [obs_prcptot_eneb.std(ddof=1), np.corrcoef(obs_prcptot_eneb, gcm_prcptot_eneb)[0,1], 'PRCPTOT-Had', 'o', 'r'],
                          [obs_r95p_eneb.std(ddof=1), np.corrcoef(obs_r95p_eneb, gcm_r95p_eneb)[0,1], 'R95p-Had', 's', 'r'],
                          [obs_r99p_eneb.std(ddof=1), np.corrcoef(obs_r99p_eneb, gcm_r99p_eneb)[0,1], 'R99p-Had', 'v', 'r'],
                          [obs_rx1day_eneb.std(ddof=1), np.corrcoef(obs_rx1day_eneb, gcm_rx1day_eneb)[0,1], 'Rx1day-Had', '^', 'r'],
                          [obs_rx5day_eneb.std(ddof=1), np.corrcoef(obs_rx5day_eneb, gcm_rx5day_eneb)[0,1], 'Rx5day-Had', 'd', 'r'],
                          [obs_sdii_eneb.std(ddof=1), np.corrcoef(obs_sdii_eneb, gcm_sdii_eneb)[0,1], 'SDII-Had', 'h', 'r'],
                          [obs_cdd_eneb.std(ddof=1), np.corrcoef(obs_cdd_eneb, gcm_cdd_eneb)[0,1], 'CDD-Had', '*', 'r'],
                          [obs_cwd_eneb.std(ddof=1), np.corrcoef(obs_cwd_eneb, gcm_cwd_eneb)[0,1], 'CWD-Had', 'p', 'r'],
                          [obs_r10mm_eneb.std(ddof=1), np.corrcoef(obs_r10mm_eneb, gcm_r10mm_eneb)[0,1], 'R10mm-Had', '>', 'r'],
                          [obs_r20mm_eneb.std(ddof=1), np.corrcoef(obs_r20mm_eneb, gcm_r20mm_eneb)[0,1], 'R20mm-Had', '<', 'r']], 
               MATOPIBA=[[obs_prcptot_matopiba.std(ddof=1), np.corrcoef(obs_prcptot_matopiba, rcm_prcptot_matopiba)[0,1], 'PRCPTOT-Reg', 'o', 'b'],
                          [obs_r95p_matopiba.std(ddof=1), np.corrcoef(obs_r95p_matopiba, rcm_r95p_matopiba)[0,1], 'R95p-Reg', 's', 'b'],
                          [obs_r99p_matopiba.std(ddof=1), np.corrcoef(obs_r99p_matopiba, rcm_r99p_matopiba)[0,1], 'R99p-Reg', 'v', 'b'],
                          [obs_rx1day_matopiba.std(ddof=1), np.corrcoef(obs_rx1day_matopiba, rcm_rx1day_matopiba)[0,1], 'Rx1day-Reg', '^', 'b'],
                          [obs_rx5day_matopiba.std(ddof=1), np.corrcoef(obs_rx5day_matopiba, rcm_rx5day_matopiba)[0,1], 'Rx5day-Reg', 'd', 'b'],
                          [obs_sdii_matopiba.std(ddof=1), np.corrcoef(obs_sdii_matopiba, rcm_sdii_matopiba)[0,1], 'SDII-Reg', 'h', 'b'],
                          [obs_cdd_matopiba.std(ddof=1), np.corrcoef(obs_cdd_matopiba, rcm_cdd_matopiba)[0,1], 'CDD-Reg', '*', 'b'],
                          [obs_cwd_matopiba.std(ddof=1), np.corrcoef(obs_cwd_matopiba, rcm_cwd_matopiba)[0,1], 'CWD-Reg', 'p', 'b'],
                          [obs_r10mm_matopiba.std(ddof=1), np.corrcoef(obs_r10mm_matopiba, rcm_r10mm_matopiba)[0,1], 'R10mm-Reg', '>', 'b'],
                          [obs_r20mm_matopiba.std(ddof=1), np.corrcoef(obs_r20mm_matopiba, rcm_r20mm_matopiba)[0,1], 'R20mm-Reg', '<', 'b'],
                          [obs_prcptot_matopiba.std(ddof=1), np.corrcoef(obs_prcptot_matopiba, gcm_prcptot_matopiba)[0,1], 'PRCPTOT-Had', 'o', 'r'],
                          [obs_r95p_matopiba.std(ddof=1), np.corrcoef(obs_r95p_matopiba, gcm_r95p_matopiba)[0,1], 'R95p-Had', 's', 'r'],
                          [obs_r99p_matopiba.std(ddof=1), np.corrcoef(obs_r99p_matopiba, gcm_r99p_matopiba)[0,1], 'R99p-Had', 'v', 'r'],
                          [obs_rx1day_matopiba.std(ddof=1), np.corrcoef(obs_rx1day_matopiba, gcm_rx1day_matopiba)[0,1], 'Rx1day-Had', '^', 'r'],
                          [obs_rx5day_matopiba.std(ddof=1), np.corrcoef(obs_rx5day_matopiba, gcm_rx5day_matopiba)[0,1], 'Rx5day-Had', 'd', 'r'],
                          [obs_sdii_matopiba.std(ddof=1), np.corrcoef(obs_sdii_matopiba, gcm_sdii_matopiba)[0,1], 'SDII-Had', 'h', 'r'],
                          [obs_cdd_matopiba.std(ddof=1), np.corrcoef(obs_cdd_matopiba, gcm_cdd_matopiba)[0,1], 'CDD-Had', '*', 'r'],
                          [obs_cwd_matopiba.std(ddof=1), np.corrcoef(obs_cwd_matopiba, gcm_cwd_matopiba)[0,1], 'CWD-Had', 'p', 'r'],
                          [obs_r10mm_matopiba.std(ddof=1), np.corrcoef(obs_r10mm_matopiba, gcm_r10mm_matopiba)[0,1], 'R10mm-Had', '>', 'r'],
                          [obs_r20mm_matopiba.std(ddof=1), np.corrcoef(obs_r20mm_matopiba, gcm_r20mm_matopiba)[0,1], 'R20mm-Had', '<', 'r']])

	# Plot Taylor Diagram
	rects = dict(SAMZ=311,
				 ENEB=312,
				 MATOPIBA=313)

	# Plot model end obs data taylor diagram 			 
	fig = plt.figure(figsize=(6, 6.5))
	
	for var in ['SAMZ', 'ENEB', 'MATOPIBA']:

		dia = TaylorDiagram(stdrefs[var], fig=fig, rect=rects[var], label='Reference', srange=(0., 9), extend=True)
		dia.samplePoints[0].set_color('black')
		
		# Add samples to Taylor diagram

		for i, (stddev,corrcoef,name,mark,cor) in enumerate(samples[var]):
			dia.add_sample(stddev, corrcoef,
						   label=name, marker=mark, color='black', mfc=cor, ms=8, ls='')		
			plt.text(-9., 8.8, text1[var], fontweight='bold')

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
	plt.close('all')
	plt.cla()
	exit()	

