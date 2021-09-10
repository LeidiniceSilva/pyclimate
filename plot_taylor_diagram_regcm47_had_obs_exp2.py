# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "08/05/2021"
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

    def __init__(self, refstd, fig=None, rect=336, label='_', marker='', color='', srange=(0., 3.), extend=False):
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
        ax.axis["left"].label.set_text(u' Standard Deviation')
        		
        ax.axis["right"].set_axis_direction("top")    # "Y-axis"
        ax.axis["right"].toggle(ticklabels=True, label=True)
        ax.axis["right"].major_ticklabels.set_rotation(90)
        ax.axis["right"].major_ticklabels.set_pad(12)
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


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/rcm_exp1'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	ann_obs = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1) 
	season_obs = np.nanmean(np.nanmean(value[2:239:3,:,:], axis=1), axis=1)
	djf_obs = season_obs[3:80:4]
	jja_obs = season_obs[1:80:4]

	return ann_obs, djf_obs, jja_obs
	
	
def import_rcm_exp1(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp1/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	ann_rcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1),axis=1) 
	season_rcm = np.nanmean(np.nanmean(value[2:239:3,:,:], axis=1), axis=1)
	djf_rcm = season_rcm[3:80:4]
	jja_rcm = season_rcm[1:80:4]

	return ann_rcm, djf_rcm, jja_rcm


def import_rcm_exp2(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/rcm/rcm_exp2/{0}'.format(exp)	
	arq  = '{0}/{1}_{2}_reg_had_{3}_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	ann_rcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1),axis=1) 
	season_rcm = np.nanmean(np.nanmean(value[2:239:3,:,:], axis=1), axis=1)
	djf_rcm = season_rcm[3:80:4]
	jja_rcm = season_rcm[1:80:4]

	return ann_rcm, djf_rcm, jja_rcm
	
	
def import_gcm(var, area, exp, dt):
	
	path = '/home/nice/Documents/dataset/gcm/rcm_exp1/hist'	
	arq  = '{0}/{1}_{2}_Amon_HadGEM2-ES_{3}_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, exp, dt)	
	
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	ann_gcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1) 
	season_gcm = np.nanmean(np.nanmean(value[2:239:3,:,:], axis=1), axis=1)
	djf_gcm = season_gcm[3:80:4]
	jja_gcm = season_gcm[1:80:4]

	return ann_gcm, djf_gcm, jja_gcm


if __name__=='__main__':
	
	# Import models and obs database 
	# Precipitation
	pre_ann_cru_samz, pre_djf_cru_samz, pre_jja_cru_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')
	pre_ann_rcm_exp1_samz, pre_djf_rcm_exp1_samz, pre_jja_rcm_exp1_samz = import_rcm_exp1('pr', 'samz', 'hist', '1986-2005')
	pre_ann_rcm_exp2_samz, pre_djf_rcm_exp2_samz, pre_jja_rcm_exp2_samz = import_rcm_exp2('pr', 'samz', 'hist', '1986-2005')
	pre_ann_gcm_samz, pre_djf_gcm_samz, pre_jja_gcm_samz = import_gcm('pr', 'samz', 'hist', '1986-2005')

	pre_ann_cru_eneb, pre_djf_cru_eneb, pre_jja_cru_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')
	pre_ann_rcm_exp1_eneb, pre_djf_rcm_exp1_eneb, pre_jja_rcm_exp1_eneb = import_rcm_exp1('pr', 'eneb', 'hist', '1986-2005')
	pre_ann_rcm_exp2_eneb, pre_djf_rcm_exp2_eneb, pre_jja_rcm_exp2_eneb = import_rcm_exp2('pr', 'eneb', 'hist', '1986-2005')
	pre_ann_gcm_eneb, pre_djf_gcm_eneb, pre_jja_gcm_eneb = import_gcm('pr', 'eneb', 'hist', '1986-2005')

	pre_ann_cru_matopiba, pre_djf_cru_matopiba, pre_jja_cru_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')
	pre_ann_rcm_exp1_matopiba, pre_djf_rcm_exp1_matopiba, pre_jja_rcm_exp1_matopiba = import_rcm_exp1('pr', 'matopiba', 'hist', '1986-2005')
	pre_ann_rcm_exp2_matopiba, pre_djf_rcm_exp2_matopiba, pre_jja_rcm_exp2_matopiba = import_rcm_exp2('pr', 'matopiba', 'hist', '1986-2005')
	pre_ann_gcm_matopiba, pre_djf_gcm_matopiba, pre_jja_gcm_matopiba = import_gcm('pr', 'matopiba', 'hist', '1986-2005')

	# Temperature
	tas_ann_cru_samz, tas_djf_cru_samz, tas_jja_cru_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')
	tas_ann_rcm_exp1_samz, tas_djf_rcm_exp1_samz, tas_jja_rcm_exp1_samz = import_rcm_exp1('tas', 'samz', 'hist', '1986-2005')
	tas_ann_rcm_exp2_samz, tas_djf_rcm_exp2_samz, tas_jja_rcm_exp2_samz = import_rcm_exp2('tas', 'samz', 'hist', '1986-2005')
	tas_ann_gcm_samz, tas_djf_gcm_samz, tas_jja_gcm_samz = import_gcm('tas', 'samz', 'hist', '1986-2005')

	tas_ann_cru_eneb, tas_djf_cru_eneb, tas_jja_cru_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')
	tas_ann_rcm_exp1_eneb, tas_djf_rcm_exp1_eneb, tas_jja_rcm_exp1_eneb = import_rcm_exp1('tas', 'eneb', 'hist', '1986-2005')
	tas_ann_rcm_exp2_eneb, tas_djf_rcm_exp2_eneb, tas_jja_rcm_exp2_eneb = import_rcm_exp2('tas', 'eneb', 'hist', '1986-2005')
	tas_ann_gcm_eneb, tas_djf_gcm_eneb, tas_jja_gcm_eneb = import_gcm('tas', 'eneb', 'hist', '1986-2005')

	tas_ann_cru_matopiba, tas_djf_cru_matopiba, tas_jja_cru_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')
	tas_ann_rcm_exp1_matopiba, tas_djf_rcm_exp1_matopiba, tas_jja_rcm_exp1_matopiba = import_rcm_exp1('tas', 'matopiba', 'hist', '1986-2005')
	tas_ann_rcm_exp2_matopiba, tas_djf_rcm_exp2_matopiba, tas_jja_rcm_exp2_matopiba = import_rcm_exp2('tas', 'matopiba', 'hist', '1986-2005')
	tas_ann_gcm_matopiba, tas_djf_gcm_matopiba, tas_jja_gcm_matopiba = import_gcm('tas', 'matopiba', 'hist', '1986-2005')

	# Reference database standard desviation
	stdrefs = dict(SAMZ1=1,
				 SAMZ2=1,
				 ENEB1=1,
				 ENEB2=1,
				 MATOPIBA1=1,
				 MATOPIBA2=1)       

	text1 = dict(SAMZ1='A)',
				 SAMZ2='D)',
				 ENEB1='B)',
				 ENEB2='E)',
				 MATOPIBA1='C)',
				 MATOPIBA2='F)')
				 
	# Compute stddev and correlation coefficient of models
	samples = dict(SAMZ1=[[pre_djf_cru_samz.std(ddof=1), np.corrcoef(pre_djf_cru_samz, pre_djf_rcm_exp1_samz)[0,1], 'RegCM4.7_EXP1', 'o', 'blue'],
                          [pre_jja_cru_samz.std(ddof=1), np.corrcoef(pre_jja_cru_samz, pre_jja_rcm_exp1_samz)[0,1], 'RegCM4.7_EXP1', '^', 'blue'],
                          [pre_ann_cru_samz.std(ddof=1), np.corrcoef(pre_ann_cru_samz, pre_ann_rcm_exp1_samz)[0,1], 'RegCM4.7_EXP1', 's', 'blue'],
                          [pre_djf_cru_samz.std(ddof=1), np.corrcoef(pre_djf_cru_samz, pre_djf_rcm_exp2_samz)[0,1], 'RegCM4.7_EXP2', 'o', 'red'],
                          [pre_jja_cru_samz.std(ddof=1), np.corrcoef(pre_jja_cru_samz, pre_jja_rcm_exp2_samz)[0,1], 'RegCM4.7_EXP2', '^', 'red'],
                          [pre_ann_cru_samz.std(ddof=1), np.corrcoef(pre_ann_cru_samz, pre_ann_rcm_exp2_samz)[0,1], 'RegCM4.7_EXP2', 's', 'red'],
                          [pre_djf_cru_samz.std(ddof=1), np.corrcoef(pre_djf_cru_samz, pre_djf_gcm_samz)[0,1], 'HadGEM2-ES', 'o', 'gray'],
                          [pre_jja_cru_samz.std(ddof=1), np.corrcoef(pre_jja_cru_samz, pre_jja_gcm_samz)[0,1], 'HadGEM2-ES', '^', 'gray'],
                          [pre_ann_cru_samz.std(ddof=1), np.corrcoef(pre_ann_cru_samz, pre_ann_gcm_samz)[0,1], 'HadGEM2-ES', 's', 'gray']],
                   SAMZ2=[[tas_djf_cru_samz.std(ddof=1), np.corrcoef(tas_djf_cru_samz, np.nanmean(tas_djf_rcm_exp1_samz, axis=1))[0,1], 'RegCM4.7_EXP1', 'o', 'blue'],
                          [tas_jja_cru_samz.std(ddof=1), np.corrcoef(tas_jja_cru_samz, np.nanmean(tas_jja_rcm_exp1_samz, axis=1))[0,1], 'RegCM4.7_EXP1', '^', 'blue'],
                          [tas_ann_cru_samz.std(ddof=1), np.corrcoef(tas_ann_cru_samz, np.nanmean(tas_ann_rcm_exp1_samz, axis=1))[0,1], 'RegCM4.7_EXP1', 's', 'blue'],
                          [tas_djf_cru_samz.std(ddof=1), np.corrcoef(tas_djf_cru_samz, np.nanmean(tas_djf_rcm_exp2_samz, axis=1))[0,1], 'RegCM4.7_EXP2', 'o', 'red'],
                          [tas_jja_cru_samz.std(ddof=1), np.corrcoef(tas_jja_cru_samz, np.nanmean(tas_jja_rcm_exp2_samz, axis=1))[0,1], 'RegCM4.7_EXP2', '^', 'red'],
                          [tas_ann_cru_samz.std(ddof=1), np.corrcoef(tas_ann_cru_samz, np.nanmean(tas_ann_rcm_exp2_samz, axis=1))[0,1], 'RegCM4.7_EXP2', 's', 'red'],
                          [tas_djf_cru_samz.std(ddof=1), np.corrcoef(tas_djf_cru_samz, tas_djf_gcm_samz)[0,1], 'HadGEM2-ES', 'o', 'gray'],
                          [tas_jja_cru_samz.std(ddof=1), np.corrcoef(tas_jja_cru_samz, tas_jja_gcm_samz)[0,1], 'HadGEM2-ES', '^', 'gray'],
                          [tas_ann_cru_samz.std(ddof=1), np.corrcoef(tas_ann_cru_samz, tas_ann_gcm_samz)[0,1], 'HadGEM2-ES', 's', 'gray']],
                   ENEB1=[[pre_djf_cru_eneb.std(ddof=1), np.corrcoef(pre_djf_cru_eneb, pre_djf_rcm_exp1_eneb)[0,1], 'RegCM4.7_EXP1', 'o', 'blue'],
                          [pre_jja_cru_eneb.std(ddof=1), np.corrcoef(pre_jja_cru_eneb, pre_jja_rcm_exp1_eneb)[0,1], 'RegCM4.7_EXP1', '^', 'blue'],
                          [pre_ann_cru_eneb.std(ddof=1), np.corrcoef(pre_ann_cru_eneb, pre_ann_rcm_exp1_eneb)[0,1], 'RegCM4.7_EXP1', 's', 'blue'],
                          [pre_djf_cru_eneb.std(ddof=1), np.corrcoef(pre_djf_cru_eneb, pre_djf_rcm_exp2_eneb)[0,1], 'RegCM4.7_EXP2', 'o', 'red'],
                          [pre_jja_cru_eneb.std(ddof=1), np.corrcoef(pre_jja_cru_eneb, pre_jja_rcm_exp2_eneb)[0,1], 'RegCM4.7_EXP2', '^', 'red'],
                          [pre_ann_cru_eneb.std(ddof=1), np.corrcoef(pre_ann_cru_eneb, pre_ann_rcm_exp2_eneb)[0,1], 'RegCM4.7_EXP2', 's', 'red'],
                          [pre_djf_cru_eneb.std(ddof=1), np.corrcoef(pre_djf_cru_eneb, pre_djf_gcm_eneb)[0,1], 'HadGEM2-ES', 'o', 'gray'],
                          [pre_jja_cru_eneb.std(ddof=1), np.corrcoef(pre_jja_cru_eneb, pre_jja_gcm_eneb)[0,1], 'HadGEM2-ES', '^', 'gray'],
                          [pre_ann_cru_eneb.std(ddof=1), np.corrcoef(pre_ann_cru_eneb, pre_ann_gcm_eneb)[0,1], 'HadGEM2-ES', 's', 'gray']],                  
                   ENEB2=[[tas_djf_cru_eneb.std(ddof=1), np.corrcoef(tas_djf_cru_eneb, np.nanmean(tas_djf_rcm_exp1_eneb, axis=1))[0,1], 'RegCM4.7_EXP1', 'o', 'blue'],
                          [tas_jja_cru_eneb.std(ddof=1), np.corrcoef(tas_jja_cru_eneb, np.nanmean(tas_jja_rcm_exp1_eneb, axis=1))[0,1], 'RegCM4.7_EXP1', '^', 'blue'],
                          [tas_ann_cru_eneb.std(ddof=1), np.corrcoef(tas_ann_cru_eneb, np.nanmean(tas_ann_rcm_exp1_eneb, axis=1))[0,1], 'RegCM4.7_EXP1', 's', 'blue'],
                          [tas_djf_cru_eneb.std(ddof=1), np.corrcoef(tas_djf_cru_eneb, np.nanmean(tas_djf_rcm_exp2_eneb, axis=1))[0,1], 'RegCM4.7_EXP2', 'o', 'red'],
                          [tas_jja_cru_eneb.std(ddof=1), np.corrcoef(tas_jja_cru_eneb, np.nanmean(tas_jja_rcm_exp2_eneb, axis=1))[0,1], 'RegCM4.7_EXP2', '^', 'red'],
                          [tas_ann_cru_eneb.std(ddof=1), np.corrcoef(tas_ann_cru_eneb, np.nanmean(tas_ann_rcm_exp2_eneb, axis=1))[0,1], 'RegCM4.7_EXP2', 's', 'red'],
                          [tas_djf_cru_eneb.std(ddof=1), np.corrcoef(tas_djf_cru_eneb, tas_djf_gcm_eneb)[0,1], 'HadGEM2-ES', 'o', 'gray'],
                          [tas_jja_cru_eneb.std(ddof=1), np.corrcoef(tas_jja_cru_eneb, tas_jja_gcm_eneb)[0,1], 'HadGEM2-ES', '^', 'gray'],
                          [tas_ann_cru_eneb.std(ddof=1), np.corrcoef(tas_ann_cru_eneb, tas_ann_gcm_eneb)[0,1], 'HadGEM2-ES', 's', 'gray']],
               MATOPIBA1=[[pre_djf_cru_matopiba.std(ddof=1), np.corrcoef(pre_djf_cru_matopiba, pre_djf_rcm_exp1_matopiba)[0,1], 'RegCM4.7_EXP1', 'o', 'blue'],
                          [pre_jja_cru_matopiba.std(ddof=1), np.corrcoef(pre_jja_cru_matopiba, pre_jja_rcm_exp1_matopiba)[0,1], 'RegCM4.7_EXP1', '^', 'blue'],
                          [pre_ann_cru_matopiba.std(ddof=1), np.corrcoef(pre_ann_cru_matopiba, pre_ann_rcm_exp1_matopiba)[0,1], 'RegCM4.7_EXP1', 's', 'blue'],
                          [pre_djf_cru_matopiba.std(ddof=1), np.corrcoef(pre_djf_cru_matopiba, pre_djf_rcm_exp2_matopiba)[0,1], 'RegCM4.7_EXP2', 'o', 'red'],
                          [pre_jja_cru_matopiba.std(ddof=1), np.corrcoef(pre_jja_cru_matopiba, pre_jja_rcm_exp2_matopiba)[0,1], 'RegCM4.7_EXP2', '^', 'red'],
                          [pre_ann_cru_matopiba.std(ddof=1), np.corrcoef(pre_ann_cru_matopiba, pre_ann_rcm_exp2_matopiba)[0,1], 'RegCM4.7_EXP2', 's', 'red'],
                          [pre_djf_cru_matopiba.std(ddof=1), np.corrcoef(pre_djf_cru_matopiba, pre_djf_gcm_matopiba)[0,1], 'HadGEM2-ES', 'o', 'gray'],
                          [pre_jja_cru_matopiba.std(ddof=1), np.corrcoef(pre_jja_cru_matopiba, pre_jja_gcm_matopiba)[0,1], 'HadGEM2-ES', '^', 'gray'],
                          [pre_ann_cru_matopiba.std(ddof=1), np.corrcoef(pre_ann_cru_matopiba, pre_ann_gcm_matopiba)[0,1], 'HadGEM2-ES', 's', 'gray']],                  
               MATOPIBA2=[[tas_djf_cru_matopiba.std(ddof=1), np.corrcoef(tas_djf_cru_matopiba, np.nanmean(tas_djf_rcm_exp1_matopiba, axis=1))[0,1], 'DJF', 'o', 'blue'],
                          [tas_jja_cru_matopiba.std(ddof=1), np.corrcoef(tas_jja_cru_matopiba, np.nanmean(tas_jja_rcm_exp1_matopiba, axis=1))[0,1], 'JJA', '^', 'blue'],
                          [tas_ann_cru_matopiba.std(ddof=1), np.corrcoef(tas_ann_cru_matopiba, np.nanmean(tas_ann_rcm_exp1_matopiba, axis=1))[0,1], 'ANN', 's', 'blue'],
                          [tas_djf_cru_matopiba.std(ddof=1), np.corrcoef(tas_djf_cru_matopiba, np.nanmean(tas_djf_rcm_exp2_matopiba, axis=1))[0,1], 'DJF', 'o', 'red'],
                          [tas_jja_cru_matopiba.std(ddof=1), np.corrcoef(tas_jja_cru_matopiba, np.nanmean(tas_jja_rcm_exp2_matopiba, axis=1))[0,1], 'JJA', '^', 'red'],
                          [tas_ann_cru_matopiba.std(ddof=1), np.corrcoef(tas_ann_cru_matopiba, np.nanmean(tas_ann_rcm_exp2_matopiba, axis=1))[0,1], 'ANN', 's', 'red'],
                          [tas_djf_cru_matopiba.std(ddof=1), np.corrcoef(tas_djf_cru_matopiba, tas_djf_gcm_matopiba)[0,1], 'DJF', 'o', 'gray'],
                          [tas_jja_cru_matopiba.std(ddof=1), np.corrcoef(tas_jja_cru_matopiba, tas_jja_gcm_matopiba)[0,1], 'JJA', '^', 'gray'],
                          [tas_ann_cru_matopiba.std(ddof=1), np.corrcoef(tas_ann_cru_matopiba, tas_ann_gcm_matopiba)[0,1], 'ANN', 's', 'gray']])

	rects = dict(SAMZ1=321,
				 SAMZ2=322,
				 ENEB1=323,
				 ENEB2=324,
				 MATOPIBA1=325,
				 MATOPIBA2=326)

	# Plot models and obs database taylor diagram
	fig = plt.figure(figsize=(6, 6.5))
	
	for var in ['SAMZ1', 'SAMZ2', 'ENEB1', 'ENEB2', 'MATOPIBA1', 'MATOPIBA2']:

		dia = TaylorDiagram(stdrefs[var], fig=fig, rect=rects[var], label='CRU', srange=(0., 3.), extend=False)
		dia.samplePoints[0].set_color('r')
		
		# Add samples to Taylor diagram
		for i, (stddev,corrcoef,name,mark,cor) in enumerate(samples[var]):
			dia.add_sample(stddev, corrcoef,
						   label=name, marker=mark, mfc=cor, color='black', ms=5, ls='')			   
			plt.text(-0.7, 3.2, text1[var], fontweight='bold')

		# Add RMS contours, and label them
		contours = dia.add_contours(levels=5, colors='0.5')
		plt.clabel(contours, inline=1, fontsize=8, fmt='%.1f')

	plt.text(-3.4, 10.5, 'RegCM4.7_EXP1', color='blue', fontsize=8, fontweight='bold')
	plt.text(-3.4, 10, 'RegCM4.7_EXP2', color='red', fontsize=8, fontweight='bold')
	plt.text(-3.4, 9.5, 'HadGEM2-ES', color='gray', fontsize=8, fontweight='bold')
		
	# Add a figure legend
	fig.legend(dia.samplePoints, 
			   [ p.get_label() for p in dia.samplePoints ], 
			   prop=dict(size=10), ncol=1, numpoints=1, loc='center', frameon=False)

	plt.subplots_adjust(left=0.10, bottom=0.10, right=0.90, top=0.90, wspace=0.50, hspace=0.60)
 
	# Path out to save figure
	path_out = '/home/nice/Downloads'
	name_out = 'pyplt_taylor_diagram_reg_exp2.png'
	if not os.path.exists(path_out):
		create_path(path_out)
	plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
	plt.show()
	exit()
