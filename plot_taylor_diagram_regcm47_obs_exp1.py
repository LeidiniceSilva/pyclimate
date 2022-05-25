# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/08/2019"
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
        ax.axis["top"].label.set_text(u'                    Correlação')

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].toggle(ticklabels=True, label=True)
        ax.axis["left"].major_ticklabels.set_rotation(-90)
        ax.axis["left"].label.set_pad(10) 
        ax.axis["left"].label.set_text(u'  Desvio padrão')
        		
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


def import_obs(var, area, dataset, dt):
	
	path = '/home/nice/Documents/dataset/obs/reg_exp1'
	arq  = '{0}/{1}_{2}_{3}_obs_mon_{4}_lonlat.nc'.format(path, var, area, dataset, dt)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	
	season_obs = value[2:239:3,:,:]
	djf_obs = np.nanmean(np.nanmean(season_obs[3:80:4], axis=1), axis=1)
	mam_obs = np.nanmean(np.nanmean(season_obs[0:80:4], axis=1), axis=1)
	jja_obs = np.nanmean(np.nanmean(season_obs[1:80:4], axis=1), axis=1)
	son_obs = np.nanmean(np.nanmean(season_obs[2:80:4], axis=1), axis=1)
	annual_obs = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1)

	return djf_obs, mam_obs, jja_obs, son_obs, annual_obs
	

def import_rcm(var, area, mod, dt):
	
	path = '/home/nice/Documents/dataset/rcm/reg_exp1/hist'
	arq  = '{0}/{1}_{2}_{3}_hist_mon_{4}_lonlat_seamask.nc'.format(path, var, area, mod, dt)	

	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]
	
	season_rcm = value[2:239:3,:,:]
	djf_rcm = np.nanmean(np.nanmean(season_rcm[3:80:4], axis=1), axis=1)
	mam_rcm = np.nanmean(np.nanmean(season_rcm[0:80:4], axis=1), axis=1)
	jja_rcm = np.nanmean(np.nanmean(season_rcm[1:80:4], axis=1), axis=1)
	son_rcm = np.nanmean(np.nanmean(season_rcm[2:80:4], axis=1), axis=1)
	annual_rcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1)

	return djf_rcm, mam_rcm, jja_rcm, son_rcm, annual_rcm


def import_gcm(var, area, mod, dt):
	
	path = '/home/nice/Documents/dataset/gcm/reg_exp1/hist'	
	arq  = '{0}/{1}_{2}_Amon_{3}_hist_r1i1p1_mon_{4}_lonlat_seamask.nc'.format(path, var, area, mod, dt)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value  = var[:][:,:,:]

	season_gcm = value[2:240:3,:,:]
	djf_gcm = np.nanmean(np.nanmean(season_gcm[3:80:4], axis=1), axis=1)
	mam_gcm = np.nanmean(np.nanmean(season_gcm[0:80:4], axis=1), axis=1)
	jja_gcm = np.nanmean(np.nanmean(season_gcm[1:80:4], axis=1), axis=1)
	son_gcm = np.nanmean(np.nanmean(season_gcm[2:80:4], axis=1), axis=1)
	annual_gcm = np.nanmean(np.nanmean(value[0:240:12,:,:], axis=1), axis=1)

	return djf_gcm, mam_gcm, jja_gcm, son_gcm, annual_gcm


if __name__=='__main__':
	
	# Import models and obs database 
	# Precipitation
	pre_djf_cru_samz, pre_mam_cru_samz, pre_jja_cru_samz, pre_son_cru_samz, pre_annual_cru_samz = import_obs('pre', 'samz', 'cru_ts4.04', '1986-2005')	   
	pre_djf_cru_eneb, pre_mam_cru_eneb, pre_jja_cru_eneb, pre_son_cru_eneb, pre_annual_cru_eneb = import_obs('pre', 'eneb', 'cru_ts4.04', '1986-2005')	   
	pre_djf_cru_matopiba, pre_mam_cru_matopiba, pre_jja_cru_matopiba, pre_son_cru_matopiba, pre_annual_cru_matopiba = import_obs('pre', 'matopiba', 'cru_ts4.04', '1986-2005')	   

	pre_djf_rcm_samz, pre_mam_rcm_samz, pre_jja_rcm_samz, pre_son_rcm_samz, pre_annual_rcm_samz = import_rcm('pr', 'samz', 'reg_had', '1986-2005')	   
	pre_djf_rcm_eneb, pre_mam_rcm_eneb, pre_jja_rcm_eneb, pre_son_rcm_eneb, pre_annual_rcm_eneb = import_rcm('pr', 'eneb', 'reg_had', '1986-2005')	   
	pre_djf_rcm_matopiba, pre_mam_rcm_matopiba, pre_jja_rcm_matopiba, pre_son_rcm_matopiba, pre_annual_rcm_matopiba = import_rcm('pr', 'matopiba', 'reg_had', '1986-2005')	   

	pre_djf_gcm_samz, pre_mam_gcm_samz, pre_jja_gcm_samz, pre_son_gcm_samz, pre_annual_gcm_samz = import_gcm('pr', 'samz', 'HadGEM2-ES', '1986-2005')	   
	pre_djf_gcm_eneb, pre_mam_gcm_eneb, pre_jja_gcm_eneb, pre_son_gcm_eneb, pre_annual_gcm_eneb = import_gcm('pr', 'eneb', 'HadGEM2-ES', '1986-2005')	   
	pre_djf_gcm_matopiba, pre_mam_gcm_matopiba, pre_jja_gcm_matopiba, pre_son_gcm_matopiba, pre_annual_gcm_matopiba = import_gcm('pr', 'matopiba', 'HadGEM2-ES', '1986-2005')	   

	# Temperature
	tas_djf_cru_samz, tas_mam_cru_samz, tas_jja_cru_samz, tas_son_cru_samz, tas_annual_cru_samz = import_obs('tmp', 'samz', 'cru_ts4.04', '1986-2005')	   
	tas_djf_cru_eneb, tas_mam_cru_eneb, tas_jja_cru_eneb, tas_son_cru_eneb, tas_annual_cru_eneb = import_obs('tmp', 'eneb', 'cru_ts4.04', '1986-2005')	   
	tas_djf_cru_matopiba, tas_mam_cru_matopiba, tas_jja_cru_matopiba, tas_son_cru_matopiba, tas_annual_cru_matopiba = import_obs('tmp', 'matopiba', 'cru_ts4.04', '1986-2005')	   

	tas_djf_rcm_samz, tas_mam_rcm_samz, tas_jja_rcm_samz, tas_son_rcm_samz, tas_annual_rcm_samz = import_rcm('tas', 'samz', 'reg_had', '1986-2005')	   
	tas_djf_rcm_eneb, tas_mam_rcm_eneb, tas_jja_rcm_eneb, tas_son_rcm_eneb, tas_annual_rcm_eneb = import_rcm('tas', 'eneb', 'reg_had', '1986-2005')	   
	tas_djf_rcm_matopiba, tas_mam_rcm_matopiba, tas_jja_rcm_matopiba, tas_son_rcm_matopiba, tas_annual_rcm_matopiba = import_rcm('tas', 'matopiba', 'reg_had', '1986-2005')	   

	tas_djf_gcm_samz, tas_mam_gcm_samz, tas_jja_gcm_samz, tas_son_gcm_samz, tas_annual_gcm_samz = import_gcm('tas', 'samz', 'HadGEM2-ES', '1986-2005')	   
	tas_djf_gcm_eneb, tas_mam_gcm_eneb, tas_jja_gcm_eneb, tas_son_gcm_eneb, tas_annual_gcm_eneb = import_gcm('tas', 'eneb', 'HadGEM2-ES', '1986-2005')	   
	tas_djf_gcm_matopiba, tas_mam_gcm_matopiba, tas_jja_gcm_matopiba, tas_son_gcm_matopiba, tas_annual_gcm_matopiba = import_gcm('tas', 'matopiba', 'HadGEM2-ES', '1986-2005')	   
		
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
	samples = dict(SAMZ1=[[pre_djf_cru_samz.std(ddof=1), np.corrcoef(pre_djf_cru_samz, pre_djf_rcm_samz)[0,1], 'DJF', 'o', 'g'],
                          [pre_mam_cru_samz.std(ddof=1), np.corrcoef(pre_mam_cru_samz, pre_mam_rcm_samz)[0,1], 'MAM', 'o', 'b'],
                          [pre_jja_cru_samz.std(ddof=1), np.corrcoef(pre_jja_cru_samz, pre_jja_rcm_samz)[0,1], 'JJA', 'o', 'y'],
                          [pre_son_cru_samz.std(ddof=1), np.corrcoef(pre_son_cru_samz, pre_son_rcm_samz)[0,1], 'SON', 'o', 'c'],
                          [pre_annual_cru_samz.std(ddof=1), np.corrcoef(pre_annual_cru_samz, pre_annual_rcm_samz)[0,1], 'ANN', 'o', 'r'],
                          [pre_djf_cru_samz.std(ddof=1), np.corrcoef(pre_djf_cru_samz, pre_djf_gcm_samz)[0,1], 'DJF', 's', 'g'],
                          [pre_mam_cru_samz.std(ddof=1), np.corrcoef(pre_mam_cru_samz, pre_mam_gcm_samz)[0,1], 'MAM', 's', 'b'],
                          [pre_jja_cru_samz.std(ddof=1), np.corrcoef(pre_jja_cru_samz, pre_jja_gcm_samz)[0,1], 'JJA', 's', 'y'],
                          [pre_son_cru_samz.std(ddof=1), np.corrcoef(pre_son_cru_samz, pre_son_gcm_samz)[0,1], 'SON', 's', 'c'],
                          [pre_annual_cru_samz.std(ddof=1), np.corrcoef(pre_annual_cru_samz, pre_annual_gcm_samz)[0,1], 'ANN', 's', 'r']],
                   SAMZ2=[[tas_djf_cru_samz.std(ddof=1), np.corrcoef(tas_djf_cru_samz, np.nanmean(tas_djf_rcm_samz, axis=1))[0,1], 'DJF', 'o', 'g'],
                          [tas_mam_cru_samz.std(ddof=1), np.corrcoef(tas_mam_cru_samz, np.nanmean(tas_mam_rcm_samz, axis=1))[0,1], 'MAM', 'o', 'b'],
                          [tas_jja_cru_samz.std(ddof=1), np.corrcoef(tas_jja_cru_samz, np.nanmean(tas_jja_rcm_samz, axis=1))[0,1], 'JJA', 'o', 'y'],
                          [tas_son_cru_samz.std(ddof=1), np.corrcoef(tas_son_cru_samz, np.nanmean(tas_son_rcm_samz, axis=1))[0,1], 'SON', 'o', 'c'],
                          [tas_annual_cru_samz.std(ddof=1), np.corrcoef(tas_annual_cru_samz, np.nanmean(tas_annual_rcm_samz, axis=1))[0,1], 'ANN', 'o', 'r'],
                          [tas_djf_cru_samz.std(ddof=1), np.corrcoef(tas_djf_cru_samz, tas_djf_gcm_samz)[0,1], 'DJF', 's', 'g'],
                          [tas_mam_cru_samz.std(ddof=1), np.corrcoef(tas_mam_cru_samz, tas_mam_gcm_samz)[0,1], 'MAM', 's', 'b'],
                          [tas_jja_cru_samz.std(ddof=1), np.corrcoef(tas_jja_cru_samz, tas_jja_gcm_samz)[0,1], 'JJA', 's', 'y'],
                          [tas_son_cru_samz.std(ddof=1), np.corrcoef(tas_son_cru_samz, tas_son_gcm_samz)[0,1], 'SON', 's', 'c'],
                          [tas_annual_cru_samz.std(ddof=1), np.corrcoef(tas_annual_cru_samz, tas_annual_gcm_samz)[0,1], 'ANN', 's', 'r']],
                   ENEB1=[[pre_djf_cru_eneb.std(ddof=1), np.corrcoef(pre_djf_cru_eneb, pre_djf_rcm_eneb)[0,1], 'DJF', 'o', 'g'],
                          [pre_mam_cru_eneb.std(ddof=1), np.corrcoef(pre_mam_cru_eneb, pre_mam_rcm_eneb)[0,1], 'MAM', 'o', 'b'],
                          [pre_jja_cru_eneb.std(ddof=1), np.corrcoef(pre_jja_cru_eneb, pre_jja_rcm_eneb)[0,1], 'JJA', 'o', 'y'],
                          [pre_son_cru_eneb.std(ddof=1), np.corrcoef(pre_son_cru_eneb, pre_son_rcm_eneb)[0,1], 'SON', 'o', 'c'],
                          [pre_annual_cru_eneb.std(ddof=1), np.corrcoef(pre_annual_cru_eneb, pre_annual_rcm_eneb)[0,1], 'ANN', 'o', 'r'],
                          [pre_djf_cru_eneb.std(ddof=1), np.corrcoef(pre_djf_cru_eneb, pre_djf_gcm_eneb)[0,1], 'DJF', 's', 'g'],
                          [pre_mam_cru_eneb.std(ddof=1), np.corrcoef(pre_mam_cru_eneb, pre_mam_gcm_eneb)[0,1], 'MAM', 's', 'b'],
                          [pre_jja_cru_eneb.std(ddof=1), np.corrcoef(pre_jja_cru_eneb, pre_jja_gcm_eneb)[0,1], 'JJA', 's', 'y'],
                          [pre_son_cru_eneb.std(ddof=1), np.corrcoef(pre_son_cru_eneb, pre_son_gcm_eneb)[0,1], 'SON', 's', 'c'],
                          [pre_annual_cru_eneb.std(ddof=1), np.corrcoef(pre_annual_cru_eneb, pre_annual_gcm_eneb)[0,1], 'ANN', 's', 'r']],                  
                   ENEB2=[[tas_djf_cru_eneb.std(ddof=1), np.corrcoef(tas_djf_cru_eneb, np.nanmean(tas_djf_rcm_eneb, axis=1))[0,1], 'DJF', 'o', 'g'],
                          [tas_mam_cru_eneb.std(ddof=1), np.corrcoef(tas_mam_cru_eneb, np.nanmean(tas_mam_rcm_eneb, axis=1))[0,1], 'MAM', 'o', 'b'],
                          [tas_jja_cru_eneb.std(ddof=1), np.corrcoef(tas_jja_cru_eneb, np.nanmean(tas_jja_rcm_eneb, axis=1))[0,1], 'JJA', 'o', 'y'],
                          [tas_son_cru_eneb.std(ddof=1), np.corrcoef(tas_son_cru_eneb, np.nanmean(tas_son_rcm_eneb, axis=1))[0,1], 'SON', 'o', 'c'],
                          [tas_annual_cru_eneb.std(ddof=1), np.corrcoef(tas_annual_cru_eneb, np.nanmean(tas_annual_rcm_eneb, axis=1))[0,1], 'ANN', 'o', 'r'],
                          [tas_djf_cru_eneb.std(ddof=1), np.corrcoef(tas_djf_cru_eneb, tas_djf_gcm_eneb)[0,1], 'DJF', 's', 'g'],
                          [tas_mam_cru_eneb.std(ddof=1), np.corrcoef(tas_mam_cru_eneb, tas_mam_gcm_eneb)[0,1], 'MAM', 's', 'b'],
                          [tas_jja_cru_eneb.std(ddof=1), np.corrcoef(tas_jja_cru_eneb, tas_jja_gcm_eneb)[0,1], 'JJA', 's', 'y'],
                          [tas_son_cru_eneb.std(ddof=1), np.corrcoef(tas_son_cru_eneb, tas_son_gcm_eneb)[0,1], 'SON', 's', 'c'],
                          [tas_annual_cru_eneb.std(ddof=1), np.corrcoef(tas_annual_cru_eneb, tas_annual_gcm_eneb)[0,1], 'ANN', 's', 'r']],
               MATOPIBA1=[[pre_djf_cru_matopiba.std(ddof=1), np.corrcoef(pre_djf_cru_matopiba, pre_djf_rcm_matopiba)[0,1], 'DJF', 'o', 'g'],
                          [pre_mam_cru_matopiba.std(ddof=1), np.corrcoef(pre_mam_cru_matopiba, pre_mam_rcm_matopiba)[0,1], 'MAM', 'o', 'b'],
                          [pre_jja_cru_matopiba.std(ddof=1), np.corrcoef(pre_jja_cru_matopiba, pre_jja_rcm_matopiba)[0,1], 'JJA', 'o', 'y'],
                          [pre_son_cru_matopiba.std(ddof=1), np.corrcoef(pre_son_cru_matopiba, pre_son_rcm_matopiba)[0,1], 'SON', 'o', 'c'],
                          [pre_annual_cru_matopiba.std(ddof=1), np.corrcoef(pre_annual_cru_matopiba, pre_annual_rcm_matopiba)[0,1], 'ANN', 'o', 'r'],
                          [pre_djf_cru_matopiba.std(ddof=1), np.corrcoef(pre_djf_cru_matopiba, pre_djf_gcm_matopiba)[0,1], 'DJF', 's', 'g'],
                          [pre_mam_cru_matopiba.std(ddof=1), np.corrcoef(pre_mam_cru_matopiba, pre_mam_gcm_matopiba)[0,1], 'MAM', 's', 'b'],
                          [pre_jja_cru_matopiba.std(ddof=1), np.corrcoef(pre_jja_cru_matopiba, pre_jja_gcm_matopiba)[0,1], 'JJA', 's', 'y'],
                          [pre_son_cru_matopiba.std(ddof=1), np.corrcoef(pre_son_cru_matopiba, pre_son_gcm_matopiba)[0,1], 'SON', 's', 'c'],
                          [pre_annual_cru_matopiba.std(ddof=1), np.corrcoef(pre_annual_cru_matopiba, pre_annual_gcm_matopiba)[0,1], 'ANN', 's', 'r']],                  
               MATOPIBA2=[[tas_djf_cru_matopiba.std(ddof=1), np.corrcoef(tas_djf_cru_matopiba, np.nanmean(tas_djf_rcm_matopiba, axis=1))[0,1], 'DJF', 'o', 'g'],
                          [tas_mam_cru_matopiba.std(ddof=1), np.corrcoef(tas_mam_cru_matopiba, np.nanmean(tas_mam_rcm_matopiba, axis=1))[0,1], 'MAM', 'o', 'b'],
                          [tas_jja_cru_matopiba.std(ddof=1), np.corrcoef(tas_jja_cru_matopiba, np.nanmean(tas_jja_rcm_matopiba, axis=1))[0,1], 'JJA', 'o', 'y'],
                          [tas_son_cru_matopiba.std(ddof=1), np.corrcoef(tas_son_cru_matopiba, np.nanmean(tas_son_rcm_matopiba, axis=1))[0,1], 'SON', 'o', 'c'],
                          [tas_annual_cru_matopiba.std(ddof=1), np.corrcoef(tas_annual_cru_matopiba, np.nanmean(tas_annual_rcm_matopiba, axis=1))[0,1], 'ANN', 'o', 'r'],
                          [tas_djf_cru_matopiba.std(ddof=1), np.corrcoef(tas_djf_cru_matopiba, tas_djf_gcm_matopiba)[0,1], 'DJF', 's', 'g'],
                          [tas_mam_cru_matopiba.std(ddof=1), np.corrcoef(tas_mam_cru_matopiba, tas_mam_gcm_matopiba)[0,1], 'MAM', 's', 'b'],
                          [tas_jja_cru_matopiba.std(ddof=1), np.corrcoef(tas_jja_cru_matopiba, tas_jja_gcm_matopiba)[0,1], 'JJA', 's', 'y'],
                          [tas_son_cru_matopiba.std(ddof=1), np.corrcoef(tas_son_cru_matopiba, tas_son_gcm_matopiba)[0,1], 'SON', 's', 'c'],
                          [tas_annual_cru_matopiba.std(ddof=1), np.corrcoef(tas_annual_cru_matopiba, tas_annual_gcm_matopiba)[0,1], 'ANN', 's', 'r']])

	x95 = [0.01, 8.5] # For Tair, this is for 95th level (r = 0.195)
	y95 = [0.0, 3]
	x99 = [0.01, 0.9] # For Tair, this is for 99th level (r = 0.254)
	y99 = [0.0, 3]

	rects = dict(SAMZ1=321,
				 SAMZ2=322,
				 ENEB1=323,
				 ENEB2=324,
				 MATOPIBA1=325,
				 MATOPIBA2=326)

	# Plot models and obs database taylor diagram
	fig = plt.figure(figsize=(6, 7))
	
	for var in ['SAMZ1', 'SAMZ2', 'ENEB1', 'ENEB2', 'MATOPIBA1', 'MATOPIBA2']:

		dia = TaylorDiagram(stdrefs[var], fig=fig, rect=rects[var], label=u'Referência', srange=(0., 3.), extend=True)
		dia.samplePoints[0].set_color('r')
		dia.ax.plot(x95,y95,color='black')
		dia.ax.plot(x99,y99,color='black')
				
		# Add samples to Taylor diagram
		for i, (stddev,corrcoef,name,mark,cor) in enumerate(samples[var]):
			dia.add_sample(stddev, corrcoef,
						   label=name, marker=mark, mfc=cor, color='black', ms=8, ls='')			   
			plt.text(-3.3, 3.7, text1[var], fontweight='bold')
			
		# Add RMS contours, and label them
		contours = dia.add_contours(levels=5, colors='0.5')
		plt.clabel(contours, inline=1, fontsize=8, fmt='%.1f')

	# Add a figure legend
	fig.legend(dia.samplePoints, 
			   [ p.get_label() for p in dia.samplePoints ], 
			   prop=dict(size=8), ncol=6, numpoints=1, loc='lower center')
		
	# Path out to save figure
	path_out = '/home/nice/Downloads'
	name_out = 'pyplt_taylor_diagram_reg_had_obs_1986-2005.png'
	if not os.path.exists(path_out):
		create_path(path_out)
	plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
	plt.show()
	exit()
