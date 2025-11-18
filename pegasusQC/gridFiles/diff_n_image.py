#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Image the grid of the difference between two channels.

Author: Mark Helm Dransfield

Created: Sat Jul 18 16:43:31 2020

License: CC BY-SA
"""

import numpy as np

from pegasusQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from pegasusQC.gridFiles.xarray_to_grid import xarray_to_grid
from pegasusQC.gridFiles.xdImage import xdImage
import pegasusQC.gridFiles.gridutility as gut


def diff_n_image(whizz_file, channel1, channel2, grid_space, *, method='neighbours', mask_polygon=[], mask_pixels=1, numneighbours=1):
    """
    Subtracts the data in `channel2` from those in `channel1`, then grids and images that difference.

    Parameters
    ----------
    whizzFile : Path or String

        The Path to, or String name of, the whizz file in HDF5 format.

    channel1 : String

        A name of a channel in `whizz_file`.

    channel2 : String

        A name of a channel in `whizz_file`.

    grid_space : Float

        The distance between grid cell centres in grid distance units.

    mask_polygon : numpy 2D array, optional

        If the size of mask_polygon > 0, then data_array will be masked to the area
        within the polygon defined by it.

    mask_pixels : Integer, optional

        If mask_pixels > 0, then all pixels further than `mask_pixels * grid_space` from a data
        location will be masked out. Default 1.

    numneighbours : Integer, optional

        If method='neighbours', then this is the number of neighbours to average. Default 5.

    Returns
    -------
    None.

    """
    print(f'Gridding and imaging {channel1} - {channel2}.')
    my_data1 = whizz_to_xarray(whizz_file, channel1, remove_mean=False, diff_one=False)
    my_data2 = whizz_to_xarray(whizz_file, channel2, remove_mean=False, diff_one=False)
    my_data1[channel1] = my_data1[channel1] - my_data2[channel2]
    my_data1.attrs['title'] = f'{channel1} - {channel2}'
    if 'units' in my_data2[channel2].attrs:
        my_data1[channel1].attrs['units'] =  my_data2[channel2].attrs['units']
    my_grid, my_region = xarray_to_grid(my_data1, grid_space, method=method, mask_polygon=mask_polygon, 
        mask_pixels=mask_pixels, numneighbours=numneighbours)

    xdImage(my_grid, my_grid.attrs['title'], cmap_norm='nonorm', 
        minClip=np.nan, maxClip=np.nan, gridlines=True, cb_ticks='stats', nSigma=2,
        hs=True, azdeg=45, ax=None, clipTo3Std = True, mask_polygon=mask_polygon)

    gut.report_gridStats(my_grid, mask_polygon=mask_polygon)

