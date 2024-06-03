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

groupName = config.groupName
projectName = config.projectName



def xarray_to_grid(my_data, grid_space, region=[]):
    """
    Uses `PyGMT` to interpolate `my_data` onto a regular grid. Method
    uses data to 5 x the grid spacing to focus on QC issues.

    Parameters
    ----------
    my_data : xArray Dataset
        Contains x, y, and z data dimensioned by fiducial.

    Returns
    -------
    grid : xArray DataArray
        Contains the gridded data.
    region : tuple
        Four floating point values for xmin, xmax, ymin, ymax.
    grid_space : Float
        The distance between grid cell centres in grid distance units.

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

    # grid spacing and search radius
    inspc = f'{grid_space:.0f}+e'
    maxradius = 5.0 * grid_space

    if region == []:
        region = pygmt.info(data=my_data[[x_chan, y_chan]], spacing=1)  # West, East, South, North
    print(f"Data points cover region: {region}")

    # Preprocess z_chan data using blockmedian
    data_trm = pygmt.blockmedian(
        data=my_data[[x_chan, y_chan, z_chan]],
        spacing=inspc,
        region=region,
    )

    # make grid
    grid = pygmt.surface(
        x=data_trm[0],
        y=data_trm[1],
        z=data_trm[2],
        spacing=inspc,
        M=maxradius, # maxradius parameter (use long form name for later version of pygmt)
        region=region,  # xmin, xmax, ymin, ymax
        T=0.35,  # tension factor
    )

    grid.attrs['units'] = myunits
    grid.attrs['long_name'] = z_chan
    grid.attrs['title'] = my_data.attrs['title']
    grid['x'].attrs['orig_name'] = x_chan
    grid['y'].attrs['orig_name'] = y_chan

    return grid, region
