import numpy as np
import colorcet as cc

import AirGravQC.utility.utility as util
import AirGravQC.gridFiles.read_ers as ers
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.config as config
import AirGravQC.gridFiles.gridutility as gut

from AirGravQC.gridFiles.graphicsShaded import graphicsShaded

groupName = config.groupName
projectName = config.projectName


def xdImage(data_array, mytitle, colormap=cc.m_CET_L9, cmap_norm='nonorm', 
        minClip=np.nan, maxClip=np.nan, gridlines=True, cb_ticks='stats', nSigma=2,
        hs=True, azdeg=45, ax=None, clipTo3Std=True, mask_polygon=[]):
    """
    Uses `graphicsShaded()` to display the gridded data in data_array. All
    parameters after the name of the whizzFile are just passed through
    to `graphicsShaded()`.

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
    gridlines : Bool, optional
        If True (the default), then grid lines are drawn on the image, else not.
    cb_ticks : TYPE, optional
        DESCRIPTION. The default is 'stats'.
    nSigma : TYPE, optional
        Not currently used. The default is 2.
    hs : Bool, optional
        hill-shading. The default is True.
    azdeg: Float, optional
        The shading azimuth in degrees from north, defaults to 45 deg.
    ax : Axis, optional
        The Matplotlib figure axis to be plotted to. Default None, in which case a new
        figure is made.
    clipTo3Std : Boolean, optional
        If True (the default), the data are clipped to +/- 3 standard deviations from
        the mean before imaging. This over-rides minClip and maxClip.
    mask_polygon : numpy 2D array, optional
        If the size of mask_polygon > 0, then data_array will be masked to the area
        within the polygon defined by it.

    Returns
    -------
    None.

    """

    if np.array(mask_polygon).size > 0:
        data_array = gut.maskGridByPolygon(data_array, mask_polygon, x_chan='x', y_chan='y')
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

    if 'E' in list(data_array.coords):
        data_array = data_array.rename({'E': 'x','N': 'y'})
    
    graphicsShaded(data_array.x, data_array.y, data_array, mytitle, colormap, cmap_norm, minClip=vmin, maxClip=vmax, gridlines=gridlines, 
                   cb_ticks=cb_ticks, nSigma=nSigma, hs=hs, azdeg=azdeg, ax=ax, origin='lower')


def xdsImage(data_set, mytitle, colormap=cc.m_CET_L9, cmap_norm='nonorm', 
        minClip=np.nan, maxClip=np.nan, gridlines=True, cb_ticks='stats', nSigma=2,
        hs=True, azdeg=45, ax=None, clipTo3Std = True, mask_polygon=[]):
    """
    Uses `graphicsShaded()` to display the gridded data in data_array. All
    parameters after the name of the whizzFile are just passed through
    to `graphicsShaded()`.

    Assumes only one DataArray in the xarray dataset

    Parameters
    ----------
    data_set : 2D xarray dataset
        The first dataArray in dataSet will be imaged.
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
    nSigma : TYPE, optional
        Not currently used. The default is 2.
    hs : Bool, optional
        hill-shading. The default is True.
    azdeg: Float, optional
        The shading azimuth in degrees from north, defaults to 45 deg.
    ax : Axis, optional
        The Matplotlib figure axis to be plotted to. Default None, in which case a new
        figure is made.
    clipTo3Std : Boolean, optional
        If True (the default), the data are clipped to +/- 3 standard deviations from
        the mean before imaging. This over-rides minClip and maxClip.
    mask_polygon : numpy 2D array, optional
        If the size of mask_polygon > 0, then data_array will be masked to the area
        within the polygon defined by it.

    Returns
    -------
    None.

    """

    first_data = data_set.to_array()
    first_data = first_data.squeeze('variable')
    xdImage(first_data, mytitle, colormap=colormap, cmap_norm=cmap_norm, 
        minClip=minClip, maxClip=maxClip, gridlines=gridlines, cb_ticks=cb_ticks, nSigma=nSigma,
        hs=hs, azdeg=azdeg, ax=ax, clipTo3Std=clipTo3Std, mask_polygon=mask_polygon)



