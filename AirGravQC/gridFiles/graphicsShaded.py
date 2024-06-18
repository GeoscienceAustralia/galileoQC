import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from pathlib import Path
import colorcet as cc
import matplotlib.ticker as tkr

import AirGravQC.gridFiles.graphics as graphics
import AirGravQC.utility.utility as util
import AirGravQC.config as config

groupName = config.groupName
projectName = config.projectName


def graphicsShaded(e, n, z, mytitle, colormap=cc.m_CET_L9, cmap_norm='nonorm', 
                   minClip=np.nan, maxClip=np.nan, gridlines=True, cb_ticks='stats', nSigma=2,
                   hs=True, azdeg=45, ax=None, origin='upper'):
    """
    Creates a colour image of a data array, with colour bar and grid-lines. The
    shape of z must be $shape(e) \times shape(n)$. The (e, n, z) typically are the 
    output of erm.read_ers_image().

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
        A colour map, eg cc.m_CET_L9. The default is cc.m_CET_L9.
    cmap_norm : String, optional
        Must be one of 'nonorm' (no normalisation, ie linear stretch); 'equalize'
        (equlaization stretch); 'auto'. The default is 'nonorm'.
    minClip : Float, optional
        z -> z < minClip : minClip: z. The default is np.nan - no clipping.
    maxClip : Float, optional
        z -> z > maxClip : maxClip: z. The default is np.nan - no clipping.
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
    if not 'myCmap' in plt.colormaps():
        plt.register_cmap('myCmap', colormap)
    
    if ax == None:
        fig, ax = plt.subplots(figsize=(12,6))
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
    graphics.imshow_hs(z, ax, cmap='myCmap',  cmap_norm=cmap_norm, hs=hs, colorbar=True,
                   azdeg=45, altdeg=45, blend_mode='alpha', alpha=0.7,
                   extent=(e[0], e[-1], n[0], n[-1]),origin=origin)
    return fig
    
