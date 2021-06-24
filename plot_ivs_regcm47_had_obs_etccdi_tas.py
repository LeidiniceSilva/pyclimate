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

    def __init__(self, refstd, fig=None, rect=313, label='_', marker='', color='', srange=(0., 9.), extend=False):
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
        rlocs = np.array([0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.98, 1.])
        
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

	dict_var = {u'eca_txx': u'tmax', 
	u'eca_txn': u'tmax',
	u'eca_tnx': u'tmin', 
	u'eca_tnn': u'tmin',
	u'eca_dtr': u'tmax',
	u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}

	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	annual_obs = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	return annual_obs
	
	
def import_rcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/rcm/eca'	
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_txx': u'tasmax', 
	u'eca_txn': u'tasmax',
	u'eca_tnx': u'tasmin', 
	u'eca_tnn': u'tasmin',
	u'eca_dtr': u'tasmax',
	u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	annual_rcm = np.nanmean(np.nanmean(var[:][0:20,:,:], axis=1), axis=1)
	
	return annual_rcm


def import_gcm(var, area, model, exp, freq, dt):
	
	path = '/home/nice/Documents/dataset/gcm/eca'
	arq  = '{0}/{1}_{2}_{3}_{4}_{5}_{6}_lonlat.nc'.format(path, var, area, model, exp, freq, dt)	

	dict_var = {u'eca_txx': u'tasmax', 
	u'eca_txn': u'tasmax',
	u'eca_tnx': u'tasmin', 
	u'eca_tnn': u'tasmin',
	u'eca_dtr': u'tasmax',
	u'eca_su': u'summer_days_index_per_time_period', 
	u'eca_tr': u'tropical_nights_index_per_time_period',
	u'eca_tx10p': u'very_cold_days_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tx90p': u'very_warm_days_percent_wrt_90th_percentile_of_reference_period', 
	u'eca_tn10p': u'cold_nights_percent_wrt_10th_percentile_of_reference_period',
	u'eca_tn90p': u'warm_nights_percent_wrt_90th_percentile_of_reference_period'}
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	annual_gcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	return annual_gcm


if __name__=='__main__':
	
	# Import regcm exp and cru databases 
	obs_txx_samz = import_obs('eca_txx', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_txx_samz = import_rcm('eca_txx', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_txx_samz = import_gcm('eca_txx', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_txn_samz = import_obs('eca_txn', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_txn_samz = import_rcm('eca_txn', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_txn_samz = import_gcm('eca_txn', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tnx_samz = import_obs('eca_tnx', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tnx_samz = import_rcm('eca_tnx', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tnx_samz = import_gcm('eca_tnx', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tnn_samz = import_obs('eca_tnn', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tnn_samz = import_rcm('eca_tnn', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tnn_samz = import_gcm('eca_tnn', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_dtr_samz = import_obs('eca_dtr', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_dtr_samz = import_rcm('eca_dtr', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_dtr_samz = import_gcm('eca_dtr', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_su_samz = import_obs('eca_su', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_su_samz = import_rcm('eca_su', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_su_samz = import_gcm('eca_su', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tr_samz = import_obs('eca_tr', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tr_samz = import_rcm('eca_tr', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tr_samz = import_gcm('eca_tr', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tx10p_samz = import_obs('eca_tx10p', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tx10p_samz = import_rcm('eca_tx10p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tx10p_samz = import_gcm('eca_tx10p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tx90p_samz = import_obs('eca_tx90p', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tx90p_samz = import_rcm('eca_tx90p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tx90p_samz = import_gcm('eca_tx90p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tn10p_samz = import_obs('eca_tn10p', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tn10p_samz = import_rcm('eca_tn10p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tn10p_samz = import_gcm('eca_tn10p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tn90p_samz = import_obs('eca_tn90p', 'samz', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tn90p_samz = import_rcm('eca_tn90p', 'samz', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tn90p_samz = import_gcm('eca_tn90p', 'samz', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_txx_eneb = import_obs('eca_txx', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_txx_eneb = import_rcm('eca_txx', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_txx_eneb = import_gcm('eca_txx', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_txn_eneb = import_obs('eca_txn', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_txn_eneb = import_rcm('eca_txn', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_txn_eneb = import_gcm('eca_txn', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tnx_eneb = import_obs('eca_tnx', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tnx_eneb = import_rcm('eca_tnx', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tnx_eneb = import_gcm('eca_tnx', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tnn_eneb = import_obs('eca_tnn', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tnn_eneb = import_rcm('eca_tnn', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tnn_eneb = import_gcm('eca_tnn', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_dtr_eneb = import_obs('eca_dtr', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_dtr_eneb = import_rcm('eca_dtr', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_dtr_eneb = import_gcm('eca_dtr', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_su_eneb = import_obs('eca_su', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_su_eneb = import_rcm('eca_su', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_su_eneb = import_gcm('eca_su', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tr_eneb = import_obs('eca_tr', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tr_eneb = import_rcm('eca_tr', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tr_eneb = import_gcm('eca_tr', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tx10p_eneb = import_obs('eca_tx10p', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tx10p_eneb = import_rcm('eca_tx10p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tx10p_eneb = import_gcm('eca_tx10p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tx90p_eneb = import_obs('eca_tx90p', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tx90p_eneb = import_rcm('eca_tx90p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tx90p_eneb = import_gcm('eca_tx90p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tn10p_eneb = import_obs('eca_tn10p', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tn10p_eneb = import_rcm('eca_tn10p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tn10p_eneb = import_gcm('eca_tn10p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tn90p_eneb = import_obs('eca_tn90p', 'eneb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tn90p_eneb = import_rcm('eca_tn90p', 'eneb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tn90p_eneb = import_gcm('eca_tn90p', 'eneb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_txx_matopiba = import_obs('eca_txx', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_txx_matopiba = import_rcm('eca_txx', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_txx_matopiba = import_gcm('eca_txx', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_txn_matopiba = import_obs('eca_txn', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_txn_matopiba = import_rcm('eca_txn', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_txn_matopiba = import_gcm('eca_txn', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tnx_matopiba = import_obs('eca_tnx', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tnx_matopiba = import_rcm('eca_tnx', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tnx_matopiba = import_gcm('eca_tnx', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tnn_matopiba = import_obs('eca_tnn', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tnn_matopiba = import_rcm('eca_tnn', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tnn_matopiba = import_gcm('eca_tnn', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_dtr_matopiba = import_obs('eca_dtr', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_dtr_matopiba = import_rcm('eca_dtr', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_dtr_matopiba = import_gcm('eca_dtr', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_su_matopiba = import_obs('eca_su', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_su_matopiba = import_rcm('eca_su', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_su_matopiba = import_gcm('eca_su', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tr_matopiba = import_obs('eca_tr', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tr_matopiba = import_rcm('eca_tr', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tr_matopiba = import_gcm('eca_tr', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tx10p_matopiba = import_obs('eca_tx10p', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tx10p_matopiba = import_rcm('eca_tx10p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tx10p_matopiba = import_gcm('eca_tx10p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tx90p_matopiba = import_obs('eca_tx90p', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tx90p_matopiba = import_rcm('eca_tx90p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tx90p_matopiba = import_gcm('eca_tx90p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tn10p_matopiba = import_obs('eca_tn10p', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tn10p_matopiba = import_rcm('eca_tn10p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tn10p_matopiba = import_gcm('eca_tn10p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	obs_tn90p_matopiba = import_obs('eca_tn90p', 'matopiba', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tn90p_matopiba = import_rcm('eca_tn90p', 'matopiba', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tn90p_matopiba = import_gcm('eca_tn90p', 'matopiba', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')
	
	# Reference database standard desviation
	stdrefs = dict(SAMZ = 1,
				 ENEB = 1,
				 MATOPIBA = 1)       

	text1 = dict(SAMZ='A)',
				 ENEB='B)',
				 MATOPIBA='C)')
				 
	# Compute stddev and correlation coefficient of models
	# Sample std, rho: Be sure to check order and that correct numbers are placed!
	samples = dict(SAMZ=[[obs_txx_samz.std(ddof=1), np.corrcoef(obs_txx_samz, np.nanmean(rcm_txx_samz, axis=1))[0,1], 'TXx-Reg', 'o', 'b'],
                          [obs_txn_samz.std(ddof=1), np.corrcoef(obs_txn_samz, np.nanmean(rcm_txn_samz, axis=1))[0,1], 'TXn-Reg', 's', 'b'],
                          [obs_tnx_samz.std(ddof=1), np.corrcoef(obs_tnx_samz, np.nanmean(rcm_tnx_samz, axis=1))[0,1], 'TNx-Reg', 'v', 'b'],
                          [obs_tnn_samz.std(ddof=1), np.corrcoef(obs_tnn_samz, np.nanmean(rcm_tnn_samz, axis=1))[0,1], 'TNn-Reg', '^', 'b'],
                          [obs_dtr_samz.std(ddof=1), np.corrcoef(obs_dtr_samz, np.nanmean(rcm_dtr_samz, axis=1))[0,1], 'DTR-Reg', 'd', 'b'],
                          [obs_su_samz.std(ddof=1), np.corrcoef(obs_su_samz, np.nanmean(rcm_su_samz, axis=1))[0,1], 'SU-Reg', 'h', 'b'],
                          [obs_tr_samz.std(ddof=1), np.corrcoef(obs_tr_samz, np.nanmean(rcm_tr_samz, axis=1))[0,1], 'TR-Reg', '*', 'b'],
                          [obs_tx10p_samz.std(ddof=1), np.corrcoef(obs_tx10p_samz, np.nanmean(rcm_tx10p_samz, axis=1))[0,1], 'TX10p-Reg', 'p', 'b'],
                          [obs_tx90p_samz.std(ddof=1), np.corrcoef(obs_tx90p_samz, np.nanmean(rcm_tx90p_samz, axis=1))[0,1], 'TX90p-Reg', '>', 'b'],
                          [obs_tn10p_samz.std(ddof=1), np.corrcoef(obs_tn10p_samz, np.nanmean(rcm_tn10p_samz, axis=1))[0,1], 'TN10p-Reg', '<', 'b'],
                          [obs_tn90p_samz.std(ddof=1), np.corrcoef(obs_tn90p_samz, np.nanmean(rcm_tn90p_samz, axis=1))[0,1], 'TN90p-Had', 'H', 'b'],
                          [obs_txx_samz.std(ddof=1), np.corrcoef(obs_txx_samz, gcm_txx_samz)[0,1], 'TXx-Reg', 'o', 'r'],
                          [obs_txn_samz.std(ddof=1), np.corrcoef(obs_txn_samz, gcm_txn_samz)[0,1], 'TXn-Reg', 's', 'r'],
                          [obs_tnx_samz.std(ddof=1), np.corrcoef(obs_tnx_samz, gcm_tnx_samz)[0,1], 'TNx-Reg', 'v', 'r'],
                          [obs_tnn_samz.std(ddof=1), np.corrcoef(obs_tnn_samz, gcm_tnn_samz)[0,1], 'TNn-Reg', '^', 'r'],
                          [obs_dtr_samz.std(ddof=1), np.corrcoef(obs_dtr_samz, gcm_dtr_samz)[0,1], 'DTR-Reg', 'd', 'r'],
                          [obs_su_samz.std(ddof=1), np.corrcoef(obs_su_samz, gcm_su_samz)[0,1], 'SU-Reg', 'h', 'r'],
                          [obs_tr_samz.std(ddof=1), np.corrcoef(obs_tr_samz, gcm_tr_samz)[0,1], 'TR-Reg', '*', 'r'],
                          [obs_tx10p_samz.std(ddof=1), np.corrcoef(obs_tx10p_samz, gcm_tx10p_samz)[0,1], 'TX10p-Reg', 'p', 'r'],
                          [obs_tx90p_samz.std(ddof=1), np.corrcoef(obs_tx90p_samz, gcm_tx90p_samz)[0,1], 'TX90p-Reg', '>', 'r'],
                          [obs_tn10p_samz.std(ddof=1), np.corrcoef(obs_tn10p_samz, gcm_tn10p_samz)[0,1], 'TN10p-Reg', '<', 'r'],
                          [obs_tn90p_samz.std(ddof=1), np.corrcoef(obs_tn90p_samz, gcm_tn90p_samz)[0,1], 'TN90p-Had', 'H', 'r']],
                   ENEB=[[obs_txx_eneb.std(ddof=1), np.corrcoef(obs_txx_eneb, np.nanmean(rcm_txx_eneb, axis=1))[0,1], 'TXx-Reg', 'o', 'b'],
                          [obs_txn_eneb.std(ddof=1), np.corrcoef(obs_txn_eneb, np.nanmean(rcm_txn_eneb, axis=1))[0,1], 'TXn-Reg', 's', 'b'],
                          [obs_tnx_eneb.std(ddof=1), np.corrcoef(obs_tnx_eneb, np.nanmean(rcm_tnx_eneb, axis=1))[0,1], 'TNx-Reg', 'v', 'b'],
                          [obs_tnn_eneb.std(ddof=1), np.corrcoef(obs_tnn_eneb, np.nanmean(rcm_tnn_eneb, axis=1))[0,1], 'TNn-Reg', '^', 'b'],
                          [obs_dtr_eneb.std(ddof=1), np.corrcoef(obs_dtr_eneb, np.nanmean(rcm_dtr_eneb, axis=1))[0,1], 'DTR-Reg', 'd', 'b'],
                          [obs_su_eneb.std(ddof=1), np.corrcoef(obs_su_eneb, np.nanmean(rcm_su_eneb, axis=1))[0,1], 'SU-Reg', 'h', 'b'],
                          [obs_tr_eneb.std(ddof=1), np.corrcoef(obs_tr_eneb, np.nanmean(rcm_tr_eneb, axis=1))[0,1], 'TR-Reg', '*', 'b'],
                          [obs_tx10p_eneb.std(ddof=1), np.corrcoef(obs_tx10p_eneb, np.nanmean(rcm_tx10p_eneb, axis=1))[0,1], 'TX10p-Reg', 'p', 'b'],
                          [obs_tx90p_eneb.std(ddof=1), np.corrcoef(obs_tx90p_eneb, np.nanmean(rcm_tx90p_eneb, axis=1))[0,1], 'TX90p-Reg', '>', 'b'],
                          [obs_tn10p_eneb.std(ddof=1), np.corrcoef(obs_tn10p_eneb, np.nanmean(rcm_tn10p_eneb, axis=1))[0,1], 'TN10p-Reg', '<', 'b'],
                          [obs_tn90p_eneb.std(ddof=1), np.corrcoef(obs_tn90p_eneb, np.nanmean(rcm_tn90p_eneb, axis=1))[0,1], 'TN90p-Had', 'H', 'b'],
                          [obs_txx_eneb.std(ddof=1), np.corrcoef(obs_txx_eneb, gcm_txx_eneb)[0,1], 'TXx-Reg', 'o', 'r'],
                          [obs_txn_eneb.std(ddof=1), np.corrcoef(obs_txn_eneb, gcm_txn_eneb)[0,1], 'TXn-Reg', 's', 'r'],
                          [obs_tnx_eneb.std(ddof=1), np.corrcoef(obs_tnx_eneb, gcm_tnx_eneb)[0,1], 'TNx-Reg', 'v', 'r'],
                          [obs_tnn_eneb.std(ddof=1), np.corrcoef(obs_tnn_eneb, gcm_tnn_eneb)[0,1], 'TNn-Reg', '^', 'r'],
                          [obs_dtr_eneb.std(ddof=1), np.corrcoef(obs_dtr_eneb, gcm_dtr_eneb)[0,1], 'DTR-Reg', 'd', 'r'],
                          [obs_su_eneb.std(ddof=1), np.corrcoef(obs_su_eneb, gcm_su_eneb)[0,1], 'SU-Reg', 'h', 'r'],
                          [obs_tr_eneb.std(ddof=1), np.corrcoef(obs_tr_eneb, gcm_tr_eneb)[0,1], 'TR-Reg', '*', 'r'],
                          [obs_tx10p_eneb.std(ddof=1), np.corrcoef(obs_tx10p_eneb, gcm_tx10p_eneb)[0,1], 'TX10p-Reg', 'p', 'r'],
                          [obs_tx90p_eneb.std(ddof=1), np.corrcoef(obs_tx90p_eneb, gcm_tx90p_eneb)[0,1], 'TX90p-Reg', '>', 'r'],
                          [obs_tn10p_eneb.std(ddof=1), np.corrcoef(obs_tn10p_eneb, gcm_tn10p_eneb)[0,1], 'TN10p-Reg', '<', 'r'],
                          [obs_tn90p_eneb.std(ddof=1), np.corrcoef(obs_tn90p_eneb, gcm_tn90p_eneb)[0,1], 'TN90p-Had', 'H', 'r']], 
               MATOPIBA=[[obs_txx_matopiba.std(ddof=1), np.corrcoef(obs_txx_matopiba, np.nanmean(rcm_txx_matopiba, axis=1))[0,1], 'TXx-Reg', 'o', 'b'],
                          [obs_txn_matopiba.std(ddof=1), np.corrcoef(obs_txn_matopiba, np.nanmean(rcm_txn_matopiba, axis=1))[0,1], 'TXn-Reg', 's', 'b'],
                          [obs_tnx_matopiba.std(ddof=1), np.corrcoef(obs_tnx_matopiba, np.nanmean(rcm_tnx_matopiba, axis=1))[0,1], 'TNx-Reg', 'v', 'b'],
                          [obs_tnn_matopiba.std(ddof=1), np.corrcoef(obs_tnn_matopiba, np.nanmean(rcm_tnn_matopiba, axis=1))[0,1], 'TNn-Reg', '^', 'b'],
                          [obs_dtr_matopiba.std(ddof=1), np.corrcoef(obs_dtr_matopiba, np.nanmean(rcm_dtr_matopiba, axis=1))[0,1], 'DTR-Reg', 'd', 'b'],
                          [obs_su_matopiba.std(ddof=1), np.corrcoef(obs_su_matopiba, np.nanmean(rcm_su_matopiba, axis=1))[0,1], 'SU-Reg', 'h', 'b'],
                          [obs_tr_matopiba.std(ddof=1), np.corrcoef(obs_tr_matopiba, np.nanmean(rcm_tr_matopiba, axis=1))[0,1], 'TR-Reg', '*', 'b'],
                          [obs_tx10p_matopiba.std(ddof=1), np.corrcoef(obs_tx10p_matopiba, np.nanmean(rcm_tx10p_matopiba, axis=1))[0,1], 'TX10p-Reg', 'p', 'b'],
                          [obs_tx90p_matopiba.std(ddof=1), np.corrcoef(obs_tx90p_matopiba, np.nanmean(rcm_tx90p_matopiba, axis=1))[0,1], 'TX90p-Reg', '>', 'b'],
                          [obs_tn10p_matopiba.std(ddof=1), np.corrcoef(obs_tn10p_matopiba, np.nanmean(rcm_tn10p_matopiba, axis=1))[0,1], 'TN10p-Reg', '<', 'b'],
                          [obs_tn90p_matopiba.std(ddof=1), np.corrcoef(obs_tn90p_matopiba, np.nanmean(rcm_tn90p_matopiba, axis=1))[0,1], 'TN90p-Had', 'H', 'b'],
                          [obs_txx_matopiba.std(ddof=1), np.corrcoef(obs_txx_matopiba, gcm_txx_matopiba)[0,1], 'TXx-Reg', 'o', 'r'],
                          [obs_txn_matopiba.std(ddof=1), np.corrcoef(obs_txn_matopiba, gcm_txn_matopiba)[0,1], 'TXn-Reg', 's', 'r'],
                          [obs_tnx_matopiba.std(ddof=1), np.corrcoef(obs_tnx_matopiba, gcm_tnx_matopiba)[0,1], 'TNx-Reg', 'v', 'r'],
                          [obs_tnn_matopiba.std(ddof=1), np.corrcoef(obs_tnn_matopiba, gcm_tnn_matopiba)[0,1], 'TNn-Reg', '^', 'r'],
                          [obs_dtr_matopiba.std(ddof=1), np.corrcoef(obs_dtr_matopiba, gcm_dtr_matopiba)[0,1], 'DTR-Reg', 'd', 'r'],
                          [obs_su_matopiba.std(ddof=1), np.corrcoef(obs_su_matopiba, gcm_su_matopiba)[0,1], 'SU-Reg', 'h', 'r'],
                          [obs_tr_matopiba.std(ddof=1), np.corrcoef(obs_tr_matopiba, gcm_tr_matopiba)[0,1], 'TR-Reg', '*', 'r'],
                          [obs_tx10p_matopiba.std(ddof=1), np.corrcoef(obs_tx10p_matopiba, gcm_tx10p_matopiba)[0,1], 'TX10p-Reg', 'p', 'r'],
                          [obs_tx90p_matopiba.std(ddof=1), np.corrcoef(obs_tx90p_matopiba, gcm_tx90p_matopiba)[0,1], 'TX90p-Reg', '>', 'r'],
                          [obs_tn10p_matopiba.std(ddof=1), np.corrcoef(obs_tn10p_matopiba, gcm_tn10p_matopiba)[0,1], 'TN10p-Reg', '<', 'r'],
                          [obs_tn90p_matopiba.std(ddof=1), np.corrcoef(obs_tn90p_matopiba, gcm_tn90p_matopiba)[0,1], 'TN90p-Had', 'H', 'r']])

	# Plot Taylor Diagram
	rects = dict(SAMZ=311,
				 ENEB=312,
				 MATOPIBA=313)

	# Plot model end obs data taylor diagram 			 
	fig = plt.figure(figsize=(6, 6.5))
	
	for var in ['SAMZ', 'ENEB', 'MATOPIBA']:

		dia = TaylorDiagram(stdrefs[var], fig=fig, rect=rects[var], label='Reference', srange=(0., 9.), extend=True)
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
	name_out = 'pyplt_taylor_diagram_etccdi_tas_reg_had_obs_1986-2005.png'
	if not os.path.exists(path_out):
		create_path(path_out)
	plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
	plt.close('all')
	plt.cla()
	exit()	
