#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Grid and image many channels of data from a `geoWhizz` file.
"""
import numpy as np

import AirGravQC.utility.utility as util
import AirGravQC.gridFiles.read_ers as ers
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.config as config

from AirGravQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from AirGravQC.gridFiles.xarray_to_grid import xarray_to_grid
from AirGravQC.gridFiles.xdImage import xdImage
import AirGravQC.gridFiles.gridutility as gut

groupName = config.groupName
projectName = config.projectName


def grid_n_image(whizz_file, z_chans, grid_space, *, lines=[], e_chan='', n_chan='', mr_chans=[], d1_chans=[], sh_chans=[], minClip=np.nan, maxClip=np.nan, 
    gridlines=True, method='neighbours', mask_polygon=[], mask_pixels=1, numneighbours=1):
    """
    Every channel in `z_chans` from `whizz_file` is interpolated onto a grid and imaged.
    Channels listed in `mr_chans` have the mean value of each survey line subtracted first.
    Channels listed in `d1_chans` are first differenced along each survey line first.

    Parameters
    ----------
    whizzFile : Path or String
        The Path to, or String name of, the whizz file in HDF5 format.
    z_chans : [String]
        An array of names of channels in `whizz_file` to be interpolated to a regular grid and imaged.
    grid_space : Float
        The distance between grid cell centres in grid distance units.
    lines : String Array, optional
        List of lines to be gridded. Default all lines.
    e_chan : String, optional
        The name of the x (easting) channel in `whizz_file`. Default is to use
        the `XChannel` attribute.
    n_chan : String, optional
        The name of the y (northing) channel in `whizz_file`. Default is to use
        the `YChannel` attribute.
    mr_chans : [String]
        An array of names of channels from `z_chans` whose mean along each survey line should be subtracted before gridding and imaging.
    d1_chans : [String]
        An array of names of channels from `z_chans` whose first difference along each survey line should be gridded and imaged.
    sh_chans : [String]
        An array of names of channels from `z_chans` whose imaged grid will be shaded.
    gridlines : Bool, optional
        If True (the default), then grid lines are drawn on the image, else not.
    method : string, optional
        The gridding algorithm to use in interpolating the data. Available is the Verde nearest
        neighbour method - "neighbours" and the SciPy GridData "linear" method. "neighbours" is
        much faster if `pykdtree` is installed. Default `neighbours` method.
    mask_polygon : numpy 2D array, optional
        If the size of mask_polygon > 0, then data_array will be masked to the area
        within the polygon defined by it.
    mask_pixels : Integer, optional
        If mask_pixels > 0, then all pixels further than `mask_pixels * grid_space` from a data
        location will be masked out. Default 1.
    numneighbours : Integer, optional
        If method='neighbours', then this is the number of neighbours to average. Default 1.

    Returns
    -------
    Nothing.

    """
    # I changed the order of the inputs, so now I need to give feedback to users.
    if not isinstance(z_chans, list):
        print('ERROR - the 2nd argument to grid_n_image should be a list (z_chans). Check docstring for details.')
        return
    if not isinstance(grid_space, float):
        print('ERROR - the 3rd argument to grid_n_image should be a float (grid_space). Check docstring for details.')
        return
    if not isinstance(lines, list):
        print('ERROR - the lines argument to grid_n_image should be a list. Check docstring for details.')
        return
    if not isinstance(mr_chans, list):
        print('ERROR - the mr_chans argument to grid_n_image should be a list. Check docstring for details.')
        return
    if not isinstance(d1_chans, list):
        print('ERROR - the d1_chans argument to grid_n_image should be a list. Check docstring for details.')
        return
    if not isinstance(sh_chans, list):
        print('ERROR - the sh_chans argument to grid_n_image should be a list. Check docstring for details.')
        return
    for z_chan in z_chans:
        remove_mean = False
        diff_one = False
        shaded = False
        if z_chan in mr_chans:
            remove_mean = True
        if z_chan in d1_chans:
            diff_one = True
        if z_chan in sh_chans:
            shaded = True

        print(f'Gridding and imaging {z_chan}')
        my_data = whizz_to_xarray(whizz_file, z_chan, n_chan=n_chan, e_chan=e_chan, lines=lines, remove_mean=remove_mean, diff_one=diff_one)
        if len(my_data.attrs) == 0:
            continue
        my_grid, my_region = xarray_to_grid(my_data, grid_space, region=[], method=method, mask_polygon=mask_polygon, 
            mask_pixels=mask_pixels, numneighbours=numneighbours)
        xdImage(my_grid, f'{my_grid.attrs["title"]}', minClip=minClip, maxClip=maxClip, gridlines=gridlines, hs=shaded)
        gut.report_gridStats(my_grid, mask_polygon=mask_polygon)
        
