import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from pathlib import Path
import colorcet as cc
import xarray as xr
import netCDF4 as nc4
import filebrowser as fb
import rioxarray
import h5py
import pygmt
import matplotlib.ticker as tkr

import AirGravQC.utility.utility as util
import AirGravQC.gridFiles.read_ers as ers
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.config as config

from AirGravQC.gridFiles.graphicsShaded import graphicsShaded

groupName = config.groupName
projectName = config.projectName


def xdImage(data_array, mytitle, colormap=cc.m_CET_L9, cmap_norm='nonorm', 
        minClip=np.nan, maxClip=np.nan, gridlines=True, cb_ticks='stats', nSigma=2,
        hs=True, azdeg=45, ax=None, clipTo3Std = True):
    """
    Uses `graphicsShaded()` to display the gridded data in data_array. All
    parameters after the name of the whizzFile are just passed through
    to `graphicsShaded()`.

    Assumes only one DataArray in the xarray dataset

    Parameters
    ----------
    data_array : 2D xarray
        The data to be imaged.
    mytitle : String
        The figure title.
    colormap : Colormap, optional
        A colour map, eg cc.m_CET_L9. The default is cc.m_CET_L9.
    cmap_norm : String, optional
        Must be one of 'nonorm' (no normalisation, ie linear stretch); 'equalize'
        (equalization stretch); 'auto'. The default is 'nonorm'.
    minClip : Float, optional
        z -> z < minClip : minClip: z. The default is np.nan - no clipping.
    maxClip : Float, optional
        z -> z > maxClip : maxClip: z. The default is np.nan - no clipping.
    cb_ticks : TYPE, optional
        DESCRIPTION. The default is 'stats'.
    nSigma : TYPE, optional
        Not currently used. The default is 2.
    hs : Bool, optional
        hill-shading. The default is True.

    Returns
    -------
    None.

    """
    # if dataname == '':
    #     dataArray = 'data'

    vmin = np.nan
    vmax = np.nan
    
    fd_mean = data_array.mean()
    fd_std = data_array.std()
    if clipTo3Std:
        a = fd_mean - 3.0 * fd_std
        vmin = a.data * 1
        a = fd_mean + 3.0 * fd_std
        vmax = 1 * a.data
    elif ~np.isnan(minClip + maxClip):
        vmin = minClip
        vmax = maxClip
    
    graphicsShaded(data_array.x, data_array.y, data_array, mytitle, colormap, cmap_norm, minClip=vmin, maxClip=vmax, gridlines=gridlines, 
                   cb_ticks=cb_ticks, nSigma=nSigma, hs=hs, azdeg=azdeg, ax=ax, origin='lower')
