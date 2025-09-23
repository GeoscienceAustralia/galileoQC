#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make a shaded image of grid data.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as ml
import h5py
from pathlib import Path
import matplotlib.ticker as tkr

import pegasusQC.gridFiles.graphics as graphics
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.utility.utility as util
import pegasusQC.config as config

groupName = config.groupName
projectName = config.projectName


def graphicsShaded(e, n, z, mytitle, colormap=config.qc_colormap, cmap_norm='nonorm', 
                   minClip=np.nan, maxClip=np.nan, gridlines=True, cb_ticks='stats', nSigma=2,
                   hs=True, azdeg=45, ax=None, origin='upper', cb_title='',
                   whizzfile=None, e_chan = '', n_chan=''):
    """
    Creates a colour image of a data array, with colour bar and grid-lines. The
    shape of z must be $shape(e) \times shape(n)$. The (e, n, z) typically are the 
    output of erm.read_ers_image(). If `whizzfile` is provided, then a flight-line
    map is drawn over the image.

    Parameters
    ----------
    e : np.array(Float, 1D
        The easting vector.
    n : np.array(Float, 1D)
        The northing vector.
    z : np.array(Float, 2D)
        The data to be imaged (referenced to the easting and northing positions).
    mytitle : String
        The figure title.
    colormap : Colormap, optional
        A colour map, eg cc.m_CET_L9. The default is config.qc_colormap.
    cmap_norm : String, optional
        Must be one of 'nonorm' (no normalisation, ie linear stretch); 'equalize'
        (equlaization stretch); 'auto'. The default is 'nonorm'.
    minClip : Float, optional
        z -> z < minClip : minClip: z. The default is np.nan - no clipping.
    maxClip : Float, optional
        z -> z > maxClip : maxClip: z. The default is np.nan - no clipping.
    gridlines : Bool, optional
        If True (the default), then grid lines are drawn on the image, else not.
    cb_ticks : TYPE, optional
        DESCRIPTION. The default is 'stats'.
    nSigma : Float, optional
        Not currently used. The default is 2.
    hs : Bool, optional
        hill-shading. The default is True.
    azdeg: Float, optional
        The shading azimuth in degrees from north, defaults to 45 deg.
    ax : Axis, optional
        The Matplotlib figure axis to be plotted to. Default None, in which case a new
        figure is made.
    origin : String, optional
        {'upper', 'lower'} Place the [0, 0] index of the array in the upper left or lower left corner
        of the Axes. The convention (the default) 'upper' is typically used for
        matrices and images.
    whizzfile : pathlib Path, optional
        If provided, the path to the whizz survey file. Default None.
    e_chan : String, optional
        The name of the field containing eastings. The default is the name
        stored in the Coordinates attribute XChannel.
    n_chan : String, optional
        The name of the field containing northings. The default is the name
        stored in the Coordinates attribute YChannel.

    Returns
    -------
    None.

    """
    if not np.isnan(minClip) and not np.isnan(maxClip):
        z = np.clip(z, minClip, maxClip)
    elif np.isnan(minClip) and (not np.isnan(maxClip)):
        z = np.clip(z, np.min(z), maxClip)
    elif (not np.isnan(minClip)) and np.isnan(maxClip):
        z = np.clip(z, minClip, np.max(z))
    
    # register the supplied colormap for access via its name
    # Somewhat dodgy if-elif code to cope with matplotlib deprecation

    # if not 'myCmap' in plt.colormaps():
    #     if "register_cmap" in dir(plt):
    #         plt.register_cmap('myCmap', colormap)
    #     elif "colormaps" in dir(ml) and "register" in dir(ml.colormaps):
    #         ml.colormaps.register(colormap, name='myCmap')
    
    if ax == None:
        fig, ax = plt.subplots()#figsize=(12,6))
    thou_format = tkr.FuncFormatter(util._space_thou)
    fig.suptitle(mytitle)#, fontsize=10)
    fig.subplots_adjust(top=0.85)
    
    ax.set_xlabel('Eastings [m]')#, fontsize=8)
    ax.set_ylabel('Northings [m]')#, fontsize=8)
    ax.grid(gridlines)
    ax.axes.set_aspect('equal')
    plt.tight_layout()
    ax.xaxis.set_major_formatter(thou_format)
    ax.yaxis.set_major_formatter(thou_format)
    # for label in ax.get_xticklabels(): label.set_fontsize(6)
    # for label in ax.get_yticklabels(): label.set_fontsize(6)
    graphics.imshow_hs(z, ax, cmap=colormap,  cmap_norm=cmap_norm, hs=hs,
                   azdeg=azdeg, altdeg=45, blend_mode='alpha', alpha=0.7, cb_title=cb_title,
                   extent=(e[0], e[-1], n[0], n[-1]), origin=origin)

    if not whizzfile == None:
        with h5py.File(whizzfile, 'r') as f:
            if e_chan == '':
                e_chan = f[groupName]['CoordinateFrame'].attrs['XChannel']
            if n_chan == '':
                n_chan = f[groupName]['CoordinateFrame'].attrs['YChannel']
            g = f[groupName]['Lines']
            for line in list(g.keys()):
                lX = rd.getLineData(g[line], e_chan)[0:]
                lY = rd.getLineData(g[line], n_chan)[0:]
                flownline, = ax.plot(lX, lY, color='blue', lw=0.6, alpha=0.7)

    return fig
    
