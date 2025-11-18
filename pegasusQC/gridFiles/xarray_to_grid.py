#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Interpolate point-located data from a 1-D to a 2-D xArray DataArray.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import xarray as xr
import verde as vd
from scipy.interpolate import griddata

import pegasusQC.utility.utility as util
import pegasusQC.gridFiles.read_ers as ers
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.gridFiles.gridutility as gut
import pegasusQC.config as config
from pegasusQC.gridFiles.minc import minc

groupName = config.groupName
projectName = config.projectName


def xarray_to_grid(my_data, grid_space, region=None, method='neighbours', mask_polygon=[], mask_pixels=1, numneighbours=5, bdist=None, maxiters=100):
    """
    Interpolates `my_data` onto a regular grid.

    Parameters
    ----------
    my_data : xArray Dataset

        Contains x, y, and z data dimensioned by fiducial.

    grid_space : Float

        The distance between grid cell centres in grid distance units.

    region : Array of float, optional

        Coordinates of the corners of the bounding rectangle [West, East, South, North].
        Default uses the minimum and maximum coordinates.

    method : string, optional

        The gridding algorithm to use in interpolating the data. Available are the Verde methods:
        "neighbours", "bicubic", and "biharmonic", the pygmi method "minc"," and the SciPy
        GridData "linear" method. The "neighbours" method is much faster if `pykdtree` is installed.
        Default "neighbours" method.

    mask_polygon : numpy 2D array, optional

        If the size of mask_polygon > 0, then data_array will be masked to the area
        within the polygon defined by it.

    mask_pixels : Integer, optional

        If mask_pixels > 0, then all pixels further than `mask_pixels * grid_space` from a data
        location will be masked out. Default 1.

    numneighbours : Integer, optional

        If method='neighbours', then this is the number of neighbours to average. Default 5.

    bdist : float, optional

        If method is "minc", then this is the blanking distance in units of cell. Default None.

    maxiters : int, optional

        Maximum number of iterations for the minimum curvature method. Default 100.

    Returns
    -------
    grid : xArray DataArray

        Contains the gridded data.

    region : tuple

        Four floating point values for xmin, xmax, ymin, ymax.

    """
    x_chan = my_data.attrs['x_channel']
    y_chan = my_data.attrs['y_channel']
    z_chan = my_data.attrs['z_channel']
    myunits = "unknown units"
    if 'units' in my_data[z_chan].attrs:
        if len(my_data[z_chan].attrs['units']) > 0 and not my_data[z_chan].attrs['units'].isspace():
            myunits = my_data[z_chan].attrs['units']
    print(f'Processing (x, y, z) = ({x_chan}, {y_chan}, {z_chan}).\n {z_chan} in {myunits}.')

    if method == 'scipy':

        if region is None:
            region = [
                        min(my_data[x_chan].data), 
                        max(my_data[x_chan].data),
                        min(my_data[y_chan].data),
                        max(my_data[y_chan].data)
                    ]
        X = np.arange(region[0], region[1], grid_space)
        Y = np.arange(region[2], region[3], grid_space)

        tmpgrid = xr.DataArray(np.zeros((len(X), len(Y))), coords=[X, Y], dims=['x', 'y'])
        grid = tmpgrid.transpose()

        east, north = np.meshgrid(X, Y)  # 2D grid for interpolation

        grid.values = griddata(list(zip(my_data[x_chan].data, my_data[y_chan].data)), my_data[z_chan].data, (east, north), method='linear')

    # method is minimum curvature
    elif method == 'minc':
        # mask_pixels = 0
        # mask_polygon = []

        myvalues, myx, myy = minc(my_data[x_chan].data, my_data[y_chan].data, my_data[z_chan].data, grid_space, extent=region, bdist=bdist,
                                maxiters=maxiters)
        X = np.array(myx)
        Y = np.array(myy)
        east, north = np.meshgrid(X, Y)  # 2D grid for masking

        tmpgrid = xr.DataArray(np.zeros((len(X), len(Y))), coords=[X, Y], dims=['x', 'y'])
        grid = tmpgrid.transpose()
        grid.values = myvalues[::-1]

    # method is one of the Verde methods
    else:
        # 1. Decimate the data along line
        spacing = grid_space / 2.0
        reducer = vd.BlockReduce(reduction=np.median, spacing=spacing)
        coordinates, z_data = reducer.filter(
            (my_data[x_chan].data, my_data[y_chan].data), my_data[z_chan].data
        )
        # 2. Define the grid region and coordinates
        if region is None:
            region = vd.get_region([my_data[x_chan], my_data[y_chan]])
            myreg = float(region[0].values), float(region[1].values), float(region[2].values), float(region[3].values)
        else:
            myreg = float(region[0]), float(region[1]), float(region[2]), float(region[3])

        east, north = vd.grid_coordinates(region=myreg, spacing=grid_space)
        
        # 3. Grid by the chosen method
        if method == 'neighbours':
            grd = vd.KNeighbors(k=numneighbours)
            grd.fit((my_data[x_chan], my_data[y_chan]), my_data[z_chan])
        elif method == 'bicubic':
            grd = vd.Cubic().fit((my_data[x_chan], my_data[y_chan]), my_data[z_chan])
        elif method == 'biharmonic':
            grd = vd.Spline(damping=1e-10)
            grd.fit(coordinates, z_data)
        else:
            print('Error - method must be "scipy", "neighbours", "bicubic", "biharmonic".')
            return
        
        # 4. Tidy up ...
        my_grid_ds = grd.grid(coordinates=(east, north) )
        my_grid = my_grid_ds.to_array()
        my_grid = my_grid.squeeze('variable')
        grid = my_grid.rename({'easting': 'x','northing': 'y'})

    if np.array(mask_polygon).size > 0:
        grid = gut.maskGridByPolygon(grid, mask_polygon, x_chan='x', y_chan='y')

    if mask_pixels > 0:
        mymask = vd.distance_mask(
            (my_data[x_chan], my_data[y_chan]),
            maxdist = grid_space * mask_pixels,
            coordinates = (east, north),
        )
        grid = grid.where(mymask)

    grid.attrs['units'] = myunits
    grid.attrs['long_name'] = z_chan
    grid.attrs['title'] = my_data.attrs['title']
    grid.attrs['x_channel'] = x_chan
    grid.attrs['y_channel'] = y_chan

    return grid, region
