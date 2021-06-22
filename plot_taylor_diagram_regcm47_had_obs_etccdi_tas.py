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
	obs_txx = import_obs('eca_txx', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_txx = import_rcm('eca_txx', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_txx = import_gcm('eca_txx', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_txn = import_obs('eca_txn', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_txn = import_rcm('eca_txn', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_txn = import_gcm('eca_txn', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_tnx = import_obs('eca_tnx', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tnx = import_rcm('eca_tnx', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tnx = import_gcm('eca_tnx', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_tnn = import_obs('eca_tnn', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tnn = import_rcm('eca_tnn', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tnn = import_gcm('eca_tnn', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_dtr = import_obs('eca_dtr', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_dtr = import_rcm('eca_dtr', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_dtr = import_gcm('eca_dtr', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_su = import_obs('eca_su', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_su = import_rcm('eca_su', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_su = import_gcm('eca_su', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_tr = import_obs('eca_tr', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tr = import_rcm('eca_tr', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tr = import_gcm('eca_tr', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_tx10p = import_obs('eca_tx10p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tx10p = import_rcm('eca_tx10p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tx10p = import_gcm('eca_tx10p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_tx90p = import_obs('eca_tx90p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tx90p = import_rcm('eca_tx90p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tx90p = import_gcm('eca_tx90p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_tn10p = import_obs('eca_tn10p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tn10p = import_rcm('eca_tn10p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tn10p = import_gcm('eca_tn10p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	obs_tn90p = import_obs('eca_tn90p', 'amz_neb', 'cpc_obs', 'yr', '1986-2005')   
	rcm_tn90p = import_rcm('eca_tn90p', 'amz_neb', 'RegCM47_had', 'historical', 'yr', '1986-2005')
	gcm_tn90p = import_gcm('eca_tn90p', 'amz_neb', 'HadGEM2-ES', 'historical', 'yr', '1986-2005')

	# Reference database standard desviation
	stdrefs = dict(SAMZ = 1,
				 ENEB = 1,
				 MATOPIBA = 1)       

	text1 = dict(SAMZ='A)',
				 ENEB='B)',
				 MATOPIBA='C)')
				 
	# Compute stddev and correlation coefficient of models
	# Sample std, rho: Be sure to check order and that correct numbers are placed!
	samples = dict(SAMZ=[[obs_txx.std(ddof=1), np.corrcoef(obs_txx, np.nanmean(rcm_txx, axis=1))[0,1], 'TXx-Reg', 'o', 'b'],
                          [obs_txn.std(ddof=1), np.corrcoef(obs_txn, np.nanmean(rcm_txn, axis=1))[0,1], 'TXn-Reg', 's', 'b'],
                          [obs_tnx.std(ddof=1), np.corrcoef(obs_tnx, np.nanmean(rcm_tnx, axis=1))[0,1], 'TNx-Reg', 'v', 'b'],
                          [obs_tnn.std(ddof=1), np.corrcoef(obs_tnn, np.nanmean(rcm_tnn, axis=1))[0,1], 'TNn-Reg', '^', 'b'],
                          [obs_dtr.std(ddof=1), np.corrcoef(obs_dtr, np.nanmean(rcm_dtr, axis=1))[0,1], 'DTR-Reg', 'd', 'b'],
                          [obs_su.std(ddof=1), np.corrcoef(obs_su, np.nanmean(rcm_su, axis=1))[0,1], 'SU-Reg', 'h', 'b'],
                          [obs_tr.std(ddof=1), np.corrcoef(obs_tr, np.nanmean(rcm_tr, axis=1))[0,1], 'TR-Reg', '*', 'b'],
                          [obs_tx10p.std(ddof=1), np.corrcoef(obs_tx10p, np.nanmean(rcm_tx10p, axis=1))[0,1], 'TX10p-Reg', 'p', 'b'],
                          [obs_tx90p.std(ddof=1), np.corrcoef(obs_tx90p, np.nanmean(rcm_tx90p, axis=1))[0,1], 'TX90p-Reg', '>', 'b'],
                          [obs_tn10p.std(ddof=1), np.corrcoef(obs_tn10p, np.nanmean(rcm_tn10p, axis=1))[0,1], 'TN10p-Reg', '<', 'b'],
                          [obs_tn90p.std(ddof=1), np.corrcoef(obs_tn90p, np.nanmean(rcm_tn90p, axis=1))[0,1], 'TN90p-Had', 'H', 'b'],
                          [obs_txx.std(ddof=1), np.corrcoef(obs_txx, gcm_txx)[0,1], 'TXx-Reg', 'o', 'r'],
                          [obs_txn.std(ddof=1), np.corrcoef(obs_txn, gcm_txn)[0,1], 'TXn-Reg', 's', 'r'],
                          [obs_tnx.std(ddof=1), np.corrcoef(obs_tnx, gcm_tnx)[0,1], 'TNx-Reg', 'v', 'r'],
                          [obs_tnn.std(ddof=1), np.corrcoef(obs_tnn, gcm_tnn)[0,1], 'TNn-Reg', '^', 'r'],
                          [obs_dtr.std(ddof=1), np.corrcoef(obs_dtr, gcm_dtr)[0,1], 'DTR-Reg', 'd', 'r'],
                          [obs_su.std(ddof=1), np.corrcoef(obs_su, gcm_su)[0,1], 'SU-Reg', 'h', 'r'],
                          [obs_tr.std(ddof=1), np.corrcoef(obs_tr, gcm_tr)[0,1], 'TR-Reg', '*', 'r'],
                          [obs_tx10p.std(ddof=1), np.corrcoef(obs_tx10p, gcm_tx10p)[0,1], 'TX10p-Reg', 'p', 'r'],
                          [obs_tx90p.std(ddof=1), np.corrcoef(obs_tx90p, gcm_tx90p)[0,1], 'TX90p-Reg', '>', 'r'],
                          [obs_tn10p.std(ddof=1), np.corrcoef(obs_tn10p, gcm_tn10p)[0,1], 'TN10p-Reg', '<', 'r'],
                          [obs_tn90p.std(ddof=1), np.corrcoef(obs_tn90p, gcm_tn90p)[0,1], 'TN90p-Had', 'H', 'r']],
                   ENEB=[[obs_txx.std(ddof=1), np.corrcoef(obs_txx, np.nanmean(rcm_txx, axis=1))[0,1], 'TXx-Reg', 'o', 'b'],
                          [obs_txn.std(ddof=1), np.corrcoef(obs_txn, np.nanmean(rcm_txn, axis=1))[0,1], 'TXn-Reg', 's', 'b'],
                          [obs_tnx.std(ddof=1), np.corrcoef(obs_tnx, np.nanmean(rcm_tnx, axis=1))[0,1], 'TNx-Reg', 'v', 'b'],
                          [obs_tnn.std(ddof=1), np.corrcoef(obs_tnn, np.nanmean(rcm_tnn, axis=1))[0,1], 'TNn-Reg', '^', 'b'],
                          [obs_dtr.std(ddof=1), np.corrcoef(obs_dtr, np.nanmean(rcm_dtr, axis=1))[0,1], 'DTR-Reg', 'd', 'b'],
                          [obs_su.std(ddof=1), np.corrcoef(obs_su, np.nanmean(rcm_su, axis=1))[0,1], 'SU-Reg', 'h', 'b'],
                          [obs_tr.std(ddof=1), np.corrcoef(obs_tr, np.nanmean(rcm_tr, axis=1))[0,1], 'TR-Reg', '*', 'b'],
                          [obs_tx10p.std(ddof=1), np.corrcoef(obs_tx10p, np.nanmean(rcm_tx10p, axis=1))[0,1], 'TX10p-Reg', 'p', 'b'],
                          [obs_tx90p.std(ddof=1), np.corrcoef(obs_tx90p, np.nanmean(rcm_tx90p, axis=1))[0,1], 'TX90p-Reg', '>', 'b'],
                          [obs_tn10p.std(ddof=1), np.corrcoef(obs_tn10p, np.nanmean(rcm_tn10p, axis=1))[0,1], 'TN10p-Reg', '<', 'b'],
                          [obs_tn90p.std(ddof=1), np.corrcoef(obs_tn90p, np.nanmean(rcm_tn90p, axis=1))[0,1], 'TN90p-Had', 'H', 'b'],
                          [obs_txx.std(ddof=1), np.corrcoef(obs_txx, gcm_txx)[0,1], 'TXx-Reg', 'o', 'r'],
                          [obs_txn.std(ddof=1), np.corrcoef(obs_txn, gcm_txn)[0,1], 'TXn-Reg', 's', 'r'],
                          [obs_tnx.std(ddof=1), np.corrcoef(obs_tnx, gcm_tnx)[0,1], 'TNx-Reg', 'v', 'r'],
                          [obs_tnn.std(ddof=1), np.corrcoef(obs_tnn, gcm_tnn)[0,1], 'TNn-Reg', '^', 'r'],
                          [obs_dtr.std(ddof=1), np.corrcoef(obs_dtr, gcm_dtr)[0,1], 'DTR-Reg', 'd', 'r'],
                          [obs_su.std(ddof=1), np.corrcoef(obs_su, gcm_su)[0,1], 'SU-Reg', 'h', 'r'],
                          [obs_tr.std(ddof=1), np.corrcoef(obs_tr, gcm_tr)[0,1], 'TR-Reg', '*', 'r'],
                          [obs_tx10p.std(ddof=1), np.corrcoef(obs_tx10p, gcm_tx10p)[0,1], 'TX10p-Reg', 'p', 'r'],
                          [obs_tx90p.std(ddof=1), np.corrcoef(obs_tx90p, gcm_tx90p)[0,1], 'TX90p-Reg', '>', 'r'],
                          [obs_tn10p.std(ddof=1), np.corrcoef(obs_tn10p, gcm_tn10p)[0,1], 'TN10p-Reg', '<', 'r'],
                          [obs_tn90p.std(ddof=1), np.corrcoef(obs_tn90p, gcm_tn90p)[0,1], 'TN90p-Had', 'H', 'r']], 
               MATOPIBA=[[obs_txx.std(ddof=1), np.corrcoef(obs_txx, np.nanmean(rcm_txx, axis=1))[0,1], 'TXx-Reg', 'o', 'b'],
                          [obs_txn.std(ddof=1), np.corrcoef(obs_txn, np.nanmean(rcm_txn, axis=1))[0,1], 'TXn-Reg', 's', 'b'],
                          [obs_tnx.std(ddof=1), np.corrcoef(obs_tnx, np.nanmean(rcm_tnx, axis=1))[0,1], 'TNx-Reg', 'v', 'b'],
                          [obs_tnn.std(ddof=1), np.corrcoef(obs_tnn, np.nanmean(rcm_tnn, axis=1))[0,1], 'TNn-Reg', '^', 'b'],
                          [obs_dtr.std(ddof=1), np.corrcoef(obs_dtr, np.nanmean(rcm_dtr, axis=1))[0,1], 'DTR-Reg', 'd', 'b'],
                          [obs_su.std(ddof=1), np.corrcoef(obs_su, np.nanmean(rcm_su, axis=1))[0,1], 'SU-Reg', 'h', 'b'],
                          [obs_tr.std(ddof=1), np.corrcoef(obs_tr, np.nanmean(rcm_tr, axis=1))[0,1], 'TR-Reg', '*', 'b'],
                          [obs_tx10p.std(ddof=1), np.corrcoef(obs_tx10p, np.nanmean(rcm_tx10p, axis=1))[0,1], 'TX10p-Reg', 'p', 'b'],
                          [obs_tx90p.std(ddof=1), np.corrcoef(obs_tx90p, np.nanmean(rcm_tx90p, axis=1))[0,1], 'TX90p-Reg', '>', 'b'],
                          [obs_tn10p.std(ddof=1), np.corrcoef(obs_tn10p, np.nanmean(rcm_tn10p, axis=1))[0,1], 'TN10p-Reg', '<', 'b'],
                          [obs_tn90p.std(ddof=1), np.corrcoef(obs_tn90p, np.nanmean(rcm_tn90p, axis=1))[0,1], 'TN90p-Had', 'H', 'b'],
                          [obs_txx.std(ddof=1), np.corrcoef(obs_txx, gcm_txx)[0,1], 'TXx-Reg', 'o', 'r'],
                          [obs_txn.std(ddof=1), np.corrcoef(obs_txn, gcm_txn)[0,1], 'TXn-Reg', 's', 'r'],
                          [obs_tnx.std(ddof=1), np.corrcoef(obs_tnx, gcm_tnx)[0,1], 'TNx-Reg', 'v', 'r'],
                          [obs_tnn.std(ddof=1), np.corrcoef(obs_tnn, gcm_tnn)[0,1], 'TNn-Reg', '^', 'r'],
                          [obs_dtr.std(ddof=1), np.corrcoef(obs_dtr, gcm_dtr)[0,1], 'DTR-Reg', 'd', 'r'],
                          [obs_su.std(ddof=1), np.corrcoef(obs_su, gcm_su)[0,1], 'SU-Reg', 'h', 'r'],
                          [obs_tr.std(ddof=1), np.corrcoef(obs_tr, gcm_tr)[0,1], 'TR-Reg', '*', 'r'],
                          [obs_tx10p.std(ddof=1), np.corrcoef(obs_tx10p, gcm_tx10p)[0,1], 'TX10p-Reg', 'p', 'r'],
                          [obs_tx90p.std(ddof=1), np.corrcoef(obs_tx90p, gcm_tx90p)[0,1], 'TX90p-Reg', '>', 'r'],
                          [obs_tn10p.std(ddof=1), np.corrcoef(obs_tn10p, gcm_tn10p)[0,1], 'TN10p-Reg', '<', 'r'],
                          [obs_tn90p.std(ddof=1), np.corrcoef(obs_tn90p, gcm_tn90p)[0,1], 'TN90p-Had', 'H', 'r']])

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
	name_out = 'pyplt_taylor_diagram_etccdi_tas_reg_had_obs_1986-2005.png'
	if not os.path.exists(path_out):
		create_path(path_out)
	plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
	plt.show()
	exit()
