# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "03/25/2019"
__description__ = "This script plot correlation Portrait Diagram from CMIP5 models"

import os
import netCDF4
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from comp_statist_indices import compute_bias


def import_cmip5(model):
	
	param = 'pr' # pr or tas
	area  = 'neb' # amz or neb
	exp   = 'historical_r1i1p1'
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/cmip5_hist'
	arq   = '{0}/{1}_{2}_Amon_{3}_{4}_{5}.nc'.format(path, param, area,	model, exp, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_mdl = np.nanmean(np.nanmean(value, axis=1), axis=1) 
	sea_mdl = np.nanmean(np.nanmean(value[0:360:3,:,:], axis=1), axis=1)
	djf_mdl = sea_mdl[0:120:4]
	mam_mdl = sea_mdl[1:120:4]
	jja_mdl = sea_mdl[2:120:4]
	son_mdl = sea_mdl[3:120:4]

	return annual_mdl, sea_mdl, djf_mdl, mam_mdl, jja_mdl, son_mdl


def import_obs(database):
	
	param = 'pre' # pre or tmp
	area  = 'neb' # amz or neb
	date  = '197512-200511'

	path  = '/home/nice/Documentos/ufrn/PhD_project/datas/obs_data'
	arq   = '{0}/{1}_{2}_{3}_obs_mon_{4}.nc'.format(path, param, area, database, date)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	annual_obs = np.nanmean(np.nanmean(value, axis=1), axis=1) 
	sea_obs = np.nanmean(np.nanmean(value[0:360:3,:,:], axis=1), axis=1)
	djf_obs = sea_obs[0:120:4]
	mam_obs = sea_obs[1:120:4]
	jja_obs = sea_obs[2:120:4]
	son_obs = sea_obs[3:120:4]

	return annual_obs, sea_obs, djf_obs, mam_obs, jja_obs, son_obs
	
	
def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels
                     for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, fontsize=20, va="center")

    # We want to show all ticks.
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # And label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor", fontsize=20)
             
    plt.setp(ax.get_yticklabels(), ha="right", rotation_mode="anchor", 
			 fontsize=20)
             
    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="white", linestyle='-', linewidth=5)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["gray", "black"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Arguments:
        im         : The AxesImage to be labeled.
    Optional arguments:
        data       : Data used to annotate. If None, the image's data is used.
        valfmt     : The format of the annotations inside the heatmap.
                     This should either use the string format method, e.g.
                     "$ {x:.2f}", or be a :class:`matplotlib.ticker.Formatter`.
        textcolors : A list or array of two color specifications. The first is
                     used for values below a threshold, the second for those
                     above.
        threshold  : Value in data units according to which the colors from
                     textcolors are applied. If None (the default) uses the
                     middle of the colormap as separation.

    Further arguments are passed on to the created text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[im.norm(data[i, j]) > threshold])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)
	
    return texts


seasons = ['DJF','MAM','JJA','SON','All seasons', 'Annual']
models = ['BCC-CSM1.1','BCC-CSM1.1M','BNU-ESM','CanESM2','CNRM-CM5','CSIRO-ACCESS-1','CSIRO-ACCESS-3','CSIRO-MK36',
'FIO-ESM','GISS-E2-H-CC','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','INMCM4','IPSL-CM5A-LR',
'IPSL-CM5A-MR','IPSL-CM5B-LR','LASG-FGOALS-G2','LASG-FGOALS-S2','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR',
'MPI-ESM-MR','MRI-CGCM3','NCAR-CCSM4','NCAR-CESM1-BGC','NCAR-CESM1-CAM5','NorESM1-ME','NorESM1-M','ensmean_cmip5']

djf = []
mam = []
jja = []
son = []
season = []
annual = []

for mdl in models:
	print 'CMIP5 Model:', mdl
	
	# Import cmip5 model end obs database monthly	
	annual_sim, sea_sim, djf_sim, mam_sim, jja_sim, son_sim= import_cmip5(mdl)
		
	obs  = u'cru_ts4.02'
	annual_cru, sea_cru, djf_cru, mam_cru, jja_cru, son_cru = import_obs(obs)
	
	bias_djf = compute_bias(djf_sim, djf_cru)
	bias_mam = compute_bias(mam_sim, mam_cru)
	bias_jja = compute_bias(jja_sim, jja_cru)
	bias_son = compute_bias(son_sim, son_cru)
	bias_sea = compute_bias(sea_sim, sea_cru)
	bias_anual = compute_bias(annual_sim, annual_cru)
	
	djf.append(bias_djf)
	mam.append(bias_mam)
	jja.append(bias_jja)
	son.append(bias_son)
	season.append(bias_sea)
	annual.append(bias_anual)
	
harvest = np.array([djf, mam, jja, son, season, annual])
print harvest
                    
fig, ax = plt.subplots(figsize=(28,10))

# Choice variable: Rainfall (AMZ and AMZ) or Temperature (AMZ and AMZ) 
out_var    = u'pre' # pre or tmp
out_area   = u'neb' # amz or neb
area_name  = u'NEB (Lat:15S 2N, Lon:46W 34W)' # AMZ (Lat:16S 4N, Lon:74W 48W) or NEB (Lat:15S 2N, Lon:46W 34W)

if out_var == 'pre':
	var_name   = u'Rainfall'
else:
	var_name   = u'Temperature' 

fig.suptitle(u'{0} Bias - {1}  \n CMIP5-hist x CRU-ts4.02 - 1975-2005 (Reference period: 1850-2005)'.format(var_name, area_name), fontsize=30, x=0.45, y=0.80)
im, cbar = heatmap(harvest, seasons, models, vmin=-6, vmax=6, ax=ax, cmap='bwr', cbarlabel='Bias')
texts = annotate_heatmap(im, valfmt="{x:.1f}", fontsize=20)

# Save figure
path_out = '/home/nice/Documentos/ufrn/PhD_project/results/cmip5'
name_out = 'pyplt_portrait_diagram_{0}_{1}_cmip5_cru_1975-2005.png'.format(out_var, out_area)

if not os.path.exists(path_out):
	create_path(path_out)
	
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

fig.tight_layout()
plt.show()

exit()


