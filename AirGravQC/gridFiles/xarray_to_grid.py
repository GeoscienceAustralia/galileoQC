#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Interpolate point-located data from a 1-D to a 2-D xArray DataArray.
"""
import numpy as np
import xarray as xr
import verde as vd
from scipy.interpolate import griddata

import AirGravQC.utility.utility as util
import AirGravQC.gridFiles.read_ers as ers
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.gridFiles.gridutility as gut
import AirGravQC.config as config

groupName = config.groupName
projectName = config.projectName


def xarray_to_grid(my_data, grid_space, region=[], method='scipy', mask_polygon=[], numneighbours=5):
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
        The gridding algorithm to use in interpolating the data. Available is the Verde nearest
        neighbour method - "neighbours" and the SciPy GridData "linear" method. "neighbours" is
        much faster if `pykdtree` is installed. Default SciPy `linear` method.
    mask_polygon : numpy 2D array, optional
        If the size of mask_polygon > 0, then data_array will be masked to the area
        within the polygon defined by it.
    numneighbours : Integer, optional
        If method='neighbours', then this is the number of neighbours to average. Default 5.

    Returns
    -------
    grid : xArray DataArray
        Contains the gridded data.
    region : tuple
        Four floating point values for xmin, xmax, ymin, ymax.

    """

    # grid_space = 500.0

    x_chan = my_data.attrs['x_channel']
    y_chan = my_data.attrs['y_channel']
    z_chan = my_data.attrs['z_channel']
    if 'units' in my_data[z_chan].attrs:
        myunits = my_data[z_chan].attrs['units']
    else:
        myunits = "no units"
    print(f'Processing (x, y, z) = ({x_chan}, {y_chan}, {z_chan}). {z_chan} in {myunits}.')

    if method == 'pygmt':
        print('"pygmt" method has been deprecated. Please use "scipy".')
        # grid spacing and search radius
        # inspc = f'{grid_space:.0f}+e'
        # maxradius = 5.0 * grid_space

        # if region == []:
        #     region = pygmt.info(data=my_data[[x_chan, y_chan]], spacing=1)  # West, East, South, North
        # print(f"Data points cover region: {region}")

        # # Preprocess z_chan data using blockmedian
        # data_trm = pygmt.blockmedian(
        #     data=my_data[[x_chan, y_chan, z_chan]],
        #     spacing=inspc,
        #     region=region,
        # )

        # # make grid
        # grid = pygmt.surface(
        #     x=data_trm[0],
        #     y=data_trm[1],
        #     z=data_trm[2],
        #     spacing=inspc,
        #     M=maxradius, # maxradius parameter (use long form name for later version of pygmt)
        #     region=region,  # xmin, xmax, ymin, ymax
        #     T=0.35,  # tension factor
        # )
    elif method == 'scipy':

        if region == []:
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

        X, Y = np.meshgrid(X, Y)  # 2D grid for interpolation

        grid.values = griddata(list(zip(my_data[x_chan].data, my_data[y_chan].data)), my_data[z_chan].data, (X, Y), method='linear')

    elif method == 'neighbours':
        if region == []:
            region = vd.get_region([my_data[x_chan], my_data[y_chan]])
            myreg = float(region[0].values), float(region[1].values), float(region[2].values), float(region[3].values)
        else:
            myreg = float(region[0]), float(region[1]), float(region[2]), float(region[3])

        east, north = vd.grid_coordinates(region=myreg, spacing=grid_space)
        grd = vd.KNeighbors(k=numneighbours)
        grd.fit((my_data[x_chan], my_data[y_chan]), my_data[z_chan])
        my_grid_ds = grd.grid(coordinates=(east, north) )
        my_grid = my_grid_ds.to_array()
        my_grid = my_grid.squeeze('variable')
        grid = my_grid.rename({'easting': 'x','northing': 'y'})
    else:
        print('Error - method must be "scipy".')
        return

    if np.array(mask_polygon).size > 0:
        grid = gut.maskGridByPolygon(grid, mask_polygon, x_chan='x', y_chan='y')

    grid.attrs['units'] = myunits
    grid.attrs['long_name'] = z_chan
    grid.attrs['title'] = my_data.attrs['title']
    grid['x'].attrs['orig_name'] = x_chan
    grid['y'].attrs['orig_name'] = y_chan

    return grid, region
