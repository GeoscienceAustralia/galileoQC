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

from AirGravQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from AirGravQC.gridFiles.xarray_to_grid import xarray_to_grid
from AirGravQC.gridFiles.xdImage import xdImage

groupName = config.groupName
projectName = config.projectName


def grid_n_image(whizz_file, z_chans, mr_chans, d1_chans, grid_space, lines=[], n_chan='', e_chan='', sh_chans=[]):
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
    mr_chans : [String]
        An array of names of channels from `z_chans` whose mean along each survey line should be subtracted before gridding and imaging.
    d1_chans : [String]
        An array of names of channels from `z_chans` whose first difference along each survey line should be gridded and imaged.
    sh_chans : [String]
        An array of names of channels from `z_chans` whose imaged grid will be shaded.
    grid_space : Float
        The distance between grid cell centres in grid distance units.

    Returns
    -------
    Nothing.

    """
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
        my_grid, my_region = xarray_to_grid(my_data, grid_space)
        xdImage(my_grid, f'{my_grid.attrs["title"]}', hs=shaded)
