# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot taylor diagram from regcm47 and hadgem models and obs database"

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

    def __init__(self, refstd, fig=None, rect=336, label='_', marker='', color='', srange=(0., 1.5), extend=False):
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
        #~ plt.rcParams.update({'font.size':8})
        plt.rcParams.update({'axes.titlesize': 'small'})

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")   # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text(u'r')       

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].toggle(ticklabels=True, label=True)
        ax.axis["left"].major_ticklabels.set_rotation(-90)
        ax.axis["left"].label.set_pad(10) 
        ax.axis["left"].label.set_text(u'SD')
        		
        ax.axis["right"].set_axis_direction("top")    # "Y-axis"
        ax.axis["right"].toggle(ticklabels=True, label=True)
        ax.axis["right"].major_ticklabels.set_rotation(90)
        ax.axis["right"].major_ticklabels.set_pad(12)
        ax.axis["bottom"].set_visible(False)      # Unused

        #~ ax.grid(color='k', axis='x', linestyle='--', linewidth=1)
        
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


def import_obs(var, area, dataset, period, dt):
	
	path = '/home/nice/Documents/dataset/obs/reg_exp1'
	arq  = '{0}/{1}_{2}_{3}_obs_{4}_{5}_lonlat.nc'.format(path, var, area, dataset, period, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	mean_obs = np.nanmean(np.nanmean(value, axis=1), axis=1) 
	
	return mean_obs
	
	
def import_rcm(var, area, exp, period, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp1/hist'
	arq  = '{0}/{1}_{2}_reg_had_{3}_{4}_{5}_lonlat_seamask.nc'.format(path, var, area, exp, period, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	mean_rcm = np.nanmean(np.nanmean(value, axis=1), axis=1) 

	return mean_rcm
	

def import_gcm(var, area, exp, period, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp1/hist'
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_{4}_{5}_lonlat_seamask.nc'.format(path, var, area, exp, period, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	mean_gcm = np.nanmean(np.nanmean(value, axis=1), axis=1) 

	return mean_gcm
	
	
if __name__=='__main__':
	
	# Import models and obs database 
	# Precipitation
	pre_djf_obs_samz = import_obs('pre', 'samz', 'cru_ts4.04', 'djf', '1986-2005')	   
	pre_mam_obs_samz = import_obs('pre', 'samz', 'cru_ts4.04', 'mam', '1986-2005')	   
	pre_jja_obs_samz = import_obs('pre', 'samz', 'cru_ts4.04', 'jja', '1986-2005')	   
	pre_son_obs_samz = import_obs('pre', 'samz', 'cru_ts4.04', 'son', '1986-2005')	   
	pre_ann_obs_samz = import_obs('pre', 'samz', 'cru_ts4.04', 'ann', '1986-2005')	   

	pre_djf_obs_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', 'djf', '1986-2005')	   
	pre_mam_obs_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', 'mam', '1986-2005')	   
	pre_jja_obs_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', 'jja', '1986-2005')	   
	pre_son_obs_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', 'son', '1986-2005')	   
	pre_ann_obs_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', 'ann', '1986-2005')	   

	pre_djf_obs_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', 'djf', '1986-2005')	   
	pre_mam_obs_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', 'mam', '1986-2005')	   
	pre_jja_obs_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', 'jja', '1986-2005')	   
	pre_son_obs_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', 'son', '1986-2005')	   
	pre_ann_obs_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', 'ann', '1986-2005')	

	pre_djf_rcm_samz = import_rcm('pr', 'samz', 'hist', 'djf', '1986-2005')	   
	pre_mam_rcm_samz = import_rcm('pr', 'samz', 'hist', 'mam', '1986-2005')	   
	pre_jja_rcm_samz = import_rcm('pr', 'samz', 'hist', 'jja', '1986-2005')	   
	pre_son_rcm_samz = import_rcm('pr', 'samz', 'hist', 'son', '1986-2005')	   
	pre_ann_rcm_samz = import_rcm('pr', 'samz', 'hist', 'ann', '1986-2005')	   

	pre_djf_rcm_eneb = import_rcm('pr', 'eneb', 'hist', 'djf', '1986-2005')	   
	pre_mam_rcm_eneb = import_rcm('pr', 'eneb', 'hist', 'mam', '1986-2005')	   
	pre_jja_rcm_eneb = import_rcm('pr', 'eneb', 'hist', 'jja', '1986-2005')	   
	pre_son_rcm_eneb = import_rcm('pr', 'eneb', 'hist', 'son', '1986-2005')	   
	pre_ann_rcm_eneb = import_rcm('pr', 'eneb', 'hist', 'ann', '1986-2005')	   

	pre_djf_rcm_matopiba = import_rcm('pr', 'matopiba', 'hist', 'djf', '1986-2005')	   
	pre_mam_rcm_matopiba = import_rcm('pr', 'matopiba', 'hist', 'mam', '1986-2005')	   
	pre_jja_rcm_matopiba = import_rcm('pr', 'matopiba', 'hist', 'jja', '1986-2005')	   
	pre_son_rcm_matopiba = import_rcm('pr', 'matopiba', 'hist', 'son', '1986-2005')	   
	pre_ann_rcm_matopiba = import_rcm('pr', 'matopiba', 'hist', 'ann', '1986-2005')	

	pre_djf_gcm_samz = import_gcm('pr', 'samz', 'hist', 'djf', '1986-2005')	   
	pre_mam_gcm_samz = import_gcm('pr', 'samz', 'hist', 'mam', '1986-2005')	   
	pre_jja_gcm_samz = import_gcm('pr', 'samz', 'hist', 'jja', '1986-2005')	   
	pre_son_gcm_samz = import_gcm('pr', 'samz', 'hist', 'son', '1986-2005')	   
	pre_ann_gcm_samz = import_gcm('pr', 'samz', 'hist', 'ann', '1986-2005')	   

	pre_djf_gcm_eneb = import_gcm('pr', 'eneb', 'hist', 'djf', '1986-2005')	   
	pre_mam_gcm_eneb = import_gcm('pr', 'eneb', 'hist', 'mam', '1986-2005')	   
	pre_jja_gcm_eneb = import_gcm('pr', 'eneb', 'hist', 'jja', '1986-2005')	   
	pre_son_gcm_eneb = import_gcm('pr', 'eneb', 'hist', 'son', '1986-2005')	   
	pre_ann_gcm_eneb = import_gcm('pr', 'eneb', 'hist', 'ann', '1986-2005')	   

	pre_djf_gcm_matopiba = import_gcm('pr', 'matopiba', 'hist', 'djf', '1986-2005')	   
	pre_mam_gcm_matopiba = import_gcm('pr', 'matopiba', 'hist', 'mam', '1986-2005')	   
	pre_jja_gcm_matopiba = import_gcm('pr', 'matopiba', 'hist', 'jja', '1986-2005')	   
	pre_son_gcm_matopiba = import_gcm('pr', 'matopiba', 'hist', 'son', '1986-2005')	   
	pre_ann_gcm_matopiba = import_gcm('pr', 'matopiba', 'hist', 'ann', '1986-2005')	

	# Temperature
	tas_djf_obs_samz = import_obs('tmp', 'samz', 'cru_ts4.04', 'djf', '1986-2005')	   
	tas_mam_obs_samz = import_obs('tmp', 'samz', 'cru_ts4.04', 'mam', '1986-2005')	   
	tas_jja_obs_samz = import_obs('tmp', 'samz', 'cru_ts4.04', 'jja', '1986-2005')	   
	tas_son_obs_samz = import_obs('tmp', 'samz', 'cru_ts4.04', 'son', '1986-2005')	   
	tas_ann_obs_samz = import_obs('tmp', 'samz', 'cru_ts4.04', 'ann', '1986-2005')	   

	tas_djf_obs_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', 'djf', '1986-2005')	   
	tas_mam_obs_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', 'mam', '1986-2005')	   
	tas_jja_obs_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', 'jja', '1986-2005')	   
	tas_son_obs_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', 'son', '1986-2005')	   
	tas_ann_obs_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', 'ann', '1986-2005')	   

	tas_djf_obs_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', 'djf', '1986-2005')	   
	tas_mam_obs_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', 'mam', '1986-2005')	   
	tas_jja_obs_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', 'jja', '1986-2005')	   
	tas_son_obs_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', 'son', '1986-2005')	   
	tas_ann_obs_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', 'ann', '1986-2005')	

	tas_djf_rcm_samz = import_rcm('tas', 'samz', 'hist', 'djf', '1986-2005')	   
	tas_mam_rcm_samz = import_rcm('tas', 'samz', 'hist', 'mam', '1986-2005')	   
	tas_jja_rcm_samz = import_rcm('tas', 'samz', 'hist', 'jja', '1986-2005')	   
	tas_son_rcm_samz = import_rcm('tas', 'samz', 'hist', 'son', '1986-2005')	   
	tas_ann_rcm_samz = import_rcm('tas', 'samz', 'hist', 'ann', '1986-2005')	   

	tas_djf_rcm_eneb = import_rcm('tas', 'eneb', 'hist', 'djf', '1986-2005')	   
	tas_mam_rcm_eneb = import_rcm('tas', 'eneb', 'hist', 'mam', '1986-2005')	   
	tas_jja_rcm_eneb = import_rcm('tas', 'eneb', 'hist', 'jja', '1986-2005')	   
	tas_son_rcm_eneb = import_rcm('tas', 'eneb', 'hist', 'son', '1986-2005')	   
	tas_ann_rcm_eneb = import_rcm('tas', 'eneb', 'hist', 'ann', '1986-2005')	   

	tas_djf_rcm_matopiba = import_rcm('tas', 'matopiba', 'hist', 'djf', '1986-2005')	   
	tas_mam_rcm_matopiba = import_rcm('tas', 'matopiba', 'hist', 'mam', '1986-2005')	   
	tas_jja_rcm_matopiba = import_rcm('tas', 'matopiba', 'hist', 'jja', '1986-2005')	   
	tas_son_rcm_matopiba = import_rcm('tas', 'matopiba', 'hist', 'son', '1986-2005')	   
	tas_ann_rcm_matopiba = import_rcm('tas', 'matopiba', 'hist', 'ann', '1986-2005')	

	tas_djf_gcm_samz = import_gcm('tas', 'samz', 'hist', 'djf', '1986-2005')	   
	tas_mam_gcm_samz = import_gcm('tas', 'samz', 'hist', 'mam', '1986-2005')	   
	tas_jja_gcm_samz = import_gcm('tas', 'samz', 'hist', 'jja', '1986-2005')	   
	tas_son_gcm_samz = import_gcm('tas', 'samz', 'hist', 'son', '1986-2005')	   
	tas_ann_gcm_samz = import_gcm('tas', 'samz', 'hist', 'ann', '1986-2005')	   

	tas_djf_gcm_eneb = import_gcm('tas', 'eneb', 'hist', 'djf', '1986-2005')	   
	tas_mam_gcm_eneb = import_gcm('tas', 'eneb', 'hist', 'mam', '1986-2005')	   
	tas_jja_gcm_eneb = import_gcm('tas', 'eneb', 'hist', 'jja', '1986-2005')	   
	tas_son_gcm_eneb = import_gcm('tas', 'eneb', 'hist', 'son', '1986-2005')	   
	tas_ann_gcm_eneb = import_gcm('tas', 'eneb', 'hist', 'ann', '1986-2005')	   

	tas_djf_gcm_matopiba = import_gcm('tas', 'matopiba', 'hist', 'djf', '1986-2005')	   
	tas_mam_gcm_matopiba = import_gcm('tas', 'matopiba', 'hist', 'mam', '1986-2005')	   
	tas_jja_gcm_matopiba = import_gcm('tas', 'matopiba', 'hist', 'jja', '1986-2005')	   
	tas_son_gcm_matopiba = import_gcm('tas', 'matopiba', 'hist', 'son', '1986-2005')	   
	tas_ann_gcm_matopiba = import_gcm('tas', 'matopiba', 'hist', 'ann', '1986-2005')	
		
	# Reference database standard desviation
	stdrefs = dict(SAMZ1=1,
				 SAMZ2=1,
				 ENEB1=1,
				 ENEB2=1,
				 MATOPIBA1=1,
				 MATOPIBA2=1)       

	text1 = dict(SAMZ1='A) SAMZ',
				 SAMZ2='D) SAMZ',
				 ENEB1='B) ENEB',
				 ENEB2='E) ENEB',
				 MATOPIBA1='C) MATOPIBA',
				 MATOPIBA2='F) MATOPIBA')
				 
	# Compute stddev and correlation coefficient of models
	samples = dict(SAMZ1=[[pre_djf_rcm_samz.std(ddof=1), np.corrcoef(pre_djf_obs_samz, pre_djf_rcm_samz)[0,1], "$1$", 'black'],
                          [pre_mam_rcm_samz.std(ddof=1), np.corrcoef(pre_mam_obs_samz, pre_mam_rcm_samz)[0,1], "$2$", 'black'],
                          [pre_jja_rcm_samz.std(ddof=1), np.corrcoef(pre_jja_obs_samz, pre_jja_rcm_samz)[0,1], "$3$", 'black'],
                          [pre_son_rcm_samz.std(ddof=1), np.corrcoef(pre_son_obs_samz, pre_son_rcm_samz)[0,1], "$4$", 'black'],
                          [pre_ann_rcm_samz.std(ddof=1), np.corrcoef(pre_ann_obs_samz, pre_ann_rcm_samz)[0,1], "$5$", 'black'],
                          [pre_djf_gcm_samz.std(ddof=1), np.corrcoef(pre_djf_obs_samz, pre_djf_gcm_samz)[0,1], "$1$", 'dimgray'],
                          [pre_mam_gcm_samz.std(ddof=1), np.corrcoef(pre_mam_obs_samz, pre_mam_gcm_samz)[0,1], "$2$", 'dimgray'],
                          [pre_jja_gcm_samz.std(ddof=1), np.corrcoef(pre_jja_obs_samz, pre_jja_gcm_samz)[0,1], "$3$", 'dimgray'],
                          [pre_son_gcm_samz.std(ddof=1), np.corrcoef(pre_son_obs_samz, pre_son_gcm_samz)[0,1], "$4$", 'dimgray'],
                          [pre_ann_gcm_samz.std(ddof=1), np.corrcoef(pre_ann_obs_samz, pre_ann_gcm_samz)[0,1], "$5$", 'dimgray']],
                   SAMZ2=[[tas_djf_rcm_samz.std(ddof=1), np.corrcoef(tas_djf_obs_samz, np.nanmean(tas_djf_rcm_samz, axis=1))[0,1], "$1$", 'black'],
                          [tas_mam_rcm_samz.std(ddof=1), np.corrcoef(tas_mam_obs_samz, np.nanmean(tas_mam_rcm_samz, axis=1))[0,1], "$2$", 'black'],
                          [tas_jja_rcm_samz.std(ddof=1), np.corrcoef(tas_jja_obs_samz, np.nanmean(tas_jja_rcm_samz, axis=1))[0,1], "$3$", 'black'],
                          [tas_son_rcm_samz.std(ddof=1), np.corrcoef(tas_son_obs_samz, np.nanmean(tas_son_rcm_samz, axis=1))[0,1], "$4$", 'black'],
                          [tas_ann_rcm_samz.std(ddof=1), np.corrcoef(tas_ann_obs_samz, np.nanmean(tas_ann_rcm_samz, axis=1))[0,1], "$5$", 'black'],
                          [tas_djf_gcm_samz.std(ddof=1), np.corrcoef(tas_djf_obs_samz, tas_djf_gcm_samz)[0,1], "$1$", 'dimgray'],
                          [tas_mam_gcm_samz.std(ddof=1), np.corrcoef(tas_mam_obs_samz, tas_mam_gcm_samz)[0,1], "$2$", 'dimgray'],
                          [tas_jja_gcm_samz.std(ddof=1), np.corrcoef(tas_jja_obs_samz, tas_jja_gcm_samz)[0,1], "$3$", 'dimgray'],
                          [tas_son_gcm_samz.std(ddof=1), np.corrcoef(tas_son_obs_samz, tas_son_gcm_samz)[0,1], "$4$", 'dimgray'],
                          [tas_ann_gcm_samz.std(ddof=1), np.corrcoef(tas_ann_obs_samz, tas_ann_gcm_samz)[0,1], "$5$", 'dimgray']],
                   ENEB1=[[pre_djf_rcm_eneb.std(ddof=1), np.corrcoef(pre_djf_obs_eneb, pre_djf_rcm_eneb)[0,1], "$1$", 'black'],
                          [pre_mam_rcm_eneb.std(ddof=1), np.corrcoef(pre_mam_obs_eneb, pre_mam_rcm_eneb)[0,1], "$2$", 'black'],
                          [pre_jja_rcm_eneb.std(ddof=1), np.corrcoef(pre_jja_obs_eneb, pre_jja_rcm_eneb)[0,1], "$3$", 'black'],
                          [pre_son_rcm_eneb.std(ddof=1), np.corrcoef(pre_son_obs_eneb, pre_son_rcm_eneb)[0,1], "$4$", 'black'],
                          [pre_ann_rcm_eneb.std(ddof=1), np.corrcoef(pre_ann_obs_eneb, pre_ann_rcm_eneb)[0,1], "$5$", 'black'],
                          [pre_djf_gcm_eneb.std(ddof=1), np.corrcoef(pre_djf_obs_eneb, pre_djf_gcm_eneb)[0,1], "$1$", 'dimgray'],
                          [pre_mam_gcm_eneb.std(ddof=1), np.corrcoef(pre_mam_obs_eneb, pre_mam_gcm_eneb)[0,1], "$2$", 'dimgray'],
                          [pre_jja_gcm_eneb.std(ddof=1), np.corrcoef(pre_jja_obs_eneb, pre_jja_gcm_eneb)[0,1], "$3$", 'dimgray'],
                          [pre_son_gcm_eneb.std(ddof=1), np.corrcoef(pre_son_obs_eneb, pre_son_gcm_eneb)[0,1], "$4$", 'dimgray'],
                          [pre_ann_gcm_eneb.std(ddof=1), np.corrcoef(pre_ann_obs_eneb, pre_ann_gcm_eneb)[0,1], "$5$", 'dimgray']],                  
                   ENEB2=[[tas_djf_rcm_eneb.std(ddof=1), np.corrcoef(tas_djf_obs_eneb, np.nanmean(tas_djf_rcm_eneb, axis=1))[0,1], "$1$", 'black'],
                          [tas_mam_rcm_eneb.std(ddof=1), np.corrcoef(tas_mam_obs_eneb, np.nanmean(tas_mam_rcm_eneb, axis=1))[0,1], "$2$", 'black'],
                          [tas_jja_rcm_eneb.std(ddof=1), np.corrcoef(tas_jja_obs_eneb, np.nanmean(tas_jja_rcm_eneb, axis=1))[0,1], "$3$", 'black'],
                          [tas_son_rcm_eneb.std(ddof=1), np.corrcoef(tas_son_obs_eneb, np.nanmean(tas_son_rcm_eneb, axis=1))[0,1], "$4$", 'black'],
                          [tas_ann_rcm_eneb.std(ddof=1), np.corrcoef(tas_ann_obs_eneb, np.nanmean(tas_ann_rcm_eneb, axis=1))[0,1], "$5$", 'black'],
                          [tas_djf_gcm_eneb.std(ddof=1), np.corrcoef(tas_djf_obs_eneb, tas_djf_gcm_eneb)[0,1], "$1$", 'dimgray'],
                          [tas_mam_gcm_eneb.std(ddof=1), np.corrcoef(tas_mam_obs_eneb, tas_mam_gcm_eneb)[0,1], "$2$", 'dimgray'],
                          [tas_jja_gcm_eneb.std(ddof=1), np.corrcoef(tas_jja_obs_eneb, tas_jja_gcm_eneb)[0,1], "$3$", 'dimgray'],
                          [tas_son_gcm_eneb.std(ddof=1), np.corrcoef(tas_son_obs_eneb, tas_son_gcm_eneb)[0,1], "$4$", 'dimgray'],
                          [tas_ann_gcm_eneb.std(ddof=1), np.corrcoef(tas_ann_obs_eneb, tas_ann_gcm_eneb)[0,1], "$5$", 'dimgray']],
               MATOPIBA1=[[pre_djf_rcm_matopiba.std(ddof=1), np.corrcoef(pre_djf_obs_matopiba, pre_djf_rcm_matopiba)[0,1], "$1$", 'black'],
                          [pre_mam_rcm_matopiba.std(ddof=1), np.corrcoef(pre_mam_obs_matopiba, pre_mam_rcm_matopiba)[0,1], "$2$", 'black'],
                          [pre_jja_rcm_matopiba.std(ddof=1), np.corrcoef(pre_jja_obs_matopiba, pre_jja_rcm_matopiba)[0,1], "$3$", 'black'],
                          [pre_son_rcm_matopiba.std(ddof=1), np.corrcoef(pre_son_obs_matopiba, pre_son_rcm_matopiba)[0,1], "$4$", 'black'],
                          [pre_ann_rcm_matopiba.std(ddof=1), np.corrcoef(pre_ann_obs_matopiba, pre_ann_rcm_matopiba)[0,1], "$5$", 'black'],
                          [pre_djf_gcm_matopiba.std(ddof=1), np.corrcoef(pre_djf_obs_matopiba, pre_djf_gcm_matopiba)[0,1], "$1$", 'dimgray'],
                          [pre_mam_gcm_matopiba.std(ddof=1), np.corrcoef(pre_mam_obs_matopiba, pre_mam_gcm_matopiba)[0,1], "$2$", 'dimgray'],
                          [pre_jja_gcm_matopiba.std(ddof=1), np.corrcoef(pre_jja_obs_matopiba, pre_jja_gcm_matopiba)[0,1], "$3$", 'dimgray'],
                          [pre_son_gcm_matopiba.std(ddof=1), np.corrcoef(pre_son_obs_matopiba, pre_son_gcm_matopiba)[0,1], "$4$", 'dimgray'],
                          [pre_ann_gcm_matopiba.std(ddof=1), np.corrcoef(pre_ann_obs_matopiba, pre_ann_gcm_matopiba)[0,1], "$5$", 'dimgray']],                  
               MATOPIBA2=[[tas_djf_rcm_matopiba.std(ddof=1), np.corrcoef(tas_djf_obs_matopiba, np.nanmean(tas_djf_rcm_matopiba, axis=1))[0,1], "$1$", 'black'],
                          [tas_mam_rcm_matopiba.std(ddof=1), np.corrcoef(tas_mam_obs_matopiba, np.nanmean(tas_mam_rcm_matopiba, axis=1))[0,1], "$2$", 'black'],
                          [tas_jja_rcm_matopiba.std(ddof=1), np.corrcoef(tas_jja_obs_matopiba, np.nanmean(tas_jja_rcm_matopiba, axis=1))[0,1], "$3$", 'black'],
                          [tas_son_rcm_matopiba.std(ddof=1), np.corrcoef(tas_son_obs_matopiba, np.nanmean(tas_son_rcm_matopiba, axis=1))[0,1], "$4$", 'black'],
                          [tas_ann_rcm_matopiba.std(ddof=1), np.corrcoef(tas_ann_obs_matopiba, np.nanmean(tas_ann_rcm_matopiba, axis=1))[0,1], "$5$", 'black'],
                          [tas_djf_gcm_matopiba.std(ddof=1), np.corrcoef(tas_djf_obs_matopiba, tas_djf_gcm_matopiba)[0,1], "$1$", 'dimgray'],
                          [tas_mam_gcm_matopiba.std(ddof=1), np.corrcoef(tas_mam_obs_matopiba, tas_mam_gcm_matopiba)[0,1], "$2$", 'dimgray'],
                          [tas_jja_gcm_matopiba.std(ddof=1), np.corrcoef(tas_jja_obs_matopiba, tas_jja_gcm_matopiba)[0,1], "$3$", 'dimgray'],
                          [tas_son_gcm_matopiba.std(ddof=1), np.corrcoef(tas_son_obs_matopiba, tas_son_gcm_matopiba)[0,1], "$4$", 'dimgray'],
                          [tas_ann_gcm_matopiba.std(ddof=1), np.corrcoef(tas_ann_obs_matopiba, tas_ann_gcm_matopiba)[0,1], "$5$", 'dimgray']])

	x95 = [0.01, 8.5] # For Tair, this is for 95th level (r = 0.195)
	y95 = [0.0, 1.5]
	x99 = [0.01, 0.9] # For Tair, this is for 99th level (r = 0.254)
	y99 = [0.0, 1.5]

	rects = dict(SAMZ1=321,
				 SAMZ2=322,
				 ENEB1=323,
				 ENEB2=324,
				 MATOPIBA1=325,
				 MATOPIBA2=326)

	# Plot models and obs database taylor diagram
	fig = plt.figure(figsize=(6, 7))
	
	for var in ['SAMZ1', 'SAMZ2', 'ENEB1', 'ENEB2', 'MATOPIBA1', 'MATOPIBA2']:

		dia = TaylorDiagram(stdrefs[var], fig=fig, rect=rects[var], label=u'Referência', srange=(0., 1.5), extend=True)
		dia.samplePoints[0].set_color('black')
		dia.ax.plot(x95,y95,color='black')
		dia.ax.plot(x99,y99,color='black')

		# Add samples to Taylor diagram
		for i, (stddev,corrcoef,name,cor) in enumerate(samples[var]):
			dia.add_sample(stddev, corrcoef,
						   marker=name, color=cor, markersize=8, ls='')			   
			plt.text(-1.7, 1.9, text1[var], fontweight='bold', fontsize=9)
			#~ plt.text(-1.4, 0.4, 'RMSE', fontweight='bold', color='0.6', fontsize=9)			

		# Add RMS contours, and label them
		contours = dia.add_contours(levels=5, colors='0.5')
		plt.clabel(contours, inline=1, fontsize=8, fmt='%.1f')

	plt.text(-2.3, 4.6, u'1. DJF', fontsize=7.5)
	plt.text(-2.3, 4.4, u'2. MAM', fontsize=7.5)
	plt.text(-2.3, 4.2, u'3. JJA', fontsize=7.5)
	plt.text(-2.3, 4.0, u'4. SON', fontsize=7.5)
	plt.text(-2.3, 3.8, u'5. ANN', fontsize=7.5)

	plt.text(-2.3, 1.6, u'RegCM4.7', color='black', fontsize=7.5, fontweight='bold')
	plt.text(-2.3, 1.4, u'HadGEM2-ES', color='dimgray', fontsize=7.5, fontweight='bold')

	#~ # Add a figure legend
	#~ fig.legend(dia.samplePoints, 
			   #~ [ p.get_label() for p in dia.samplePoints ], 
			   #~ prop=dict(size=8), ncol=6, numpoints=1, loc='lower center')
		
	# Path out to save figure
	path_out = '/home/nice/Downloads'
	name_out = 'pyplt_taylor_diagram_reg_had_obs_1986-2005.png'
	if not os.path.exists(path_out):
		create_path(path_out)
	plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
	plt.show()
	exit()
