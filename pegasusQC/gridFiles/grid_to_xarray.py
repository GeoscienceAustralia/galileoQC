#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads a grid file into xarray DataSet.

Author: Mark Helm Dransfield

Created: Sat Jul 18 16:43:31 2020

License: CC BY-SA
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from pathlib import Path
import colorcet as cc
import xarray as xr
import netCDF4 as nc4
import filebrowser as fb
# import rioxarray
import h5py
import matplotlib.ticker as tkr

import pegasusQC.gridFiles.graphics as graphics
import pegasusQC.utility.utility as util
import pegasusQC.gridFiles.read_ers as ers
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.config as config
from pegasusQC.gridFiles.graphicsShaded import graphicsShaded
from pegasusQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from pegasusQC.gridFiles.xarray_to_grid import xarray_to_grid
from pegasusQC.gridFiles.xdImage import xdImage
import pegasusQC.gridFiles.gridutility as gut

groupName = config.groupName
projectName = config.projectName

def gridfile_to_xr(whizzFile='', bandout=0):
    """
    Returns an xarray Dataset containing the geographically located data from a
    gridfile in either ERS or NC format.

    Parameters
    ----------
    whizzFile : Path

        The Path to the grid file, must have extension `ers` or `nc`.

    bandout : Int, optional

        The band to be read if the grid file is `ERS`. The default is 0.

    Returns
    -------
    xd : xarray Dataset

        The data from `whizzFile`.

    """
    if whizzFile == '':
        whizzFile = Path(fb.get_grid_filename(filetypes=(('NetCdf4 grid', '*.nc'),
                                                         ('ERMapper grid', '*.ers'))))
    if whizzFile.suffix.upper() == '.ERS':
        e, n, z, datum, projection = ers.read_ers_image(whizzFile, bandout=bandout)
        xa = xr.DataArray(data=z,#np.flip(z, 0), # DANGER!!!
                          dims=["N", "E"],
                          coords={"N": n, "E": e})
        xa.dropna(dim='N',how='all')
        xa.dropna(dim='E',how='all')
        fname = whizzFile.with_suffix('').name
        xd = xr.Dataset(data_vars={fname: xa})
        xd.attrs["long_name"] = fname
        xd.attrs["datum"] = datum
        xd.attrs["projection"] = projection
        if datum == 'WGS84' and projection == 'SUTM55':
            xd.rio.write_crs("epsg:32755", inplace=True)
    elif whizzFile.suffix.upper() == '.NC':
        nc = nc4.Dataset(str(whizzFile), mode='r')
        xd = xr.open_dataset(xr.backends.NetCDF4DataStore(nc))
    else:
        print(f'ERROR: whizzFile has suffix {Path(whizzFile).suffix.upper()} but must be `.nc` or `.ers`.')
    return xd, whizzFile


def gridfile_to_xa(whizzFile='', bandout=0):
    """
    Returns an xarray DataArray containing the geographically located data from a
    gridfile in either ERS or NC format.

    Parameters
    ----------
    whizzFile : Path

        The Path to the grid file, must have extension `ers` or `nc`.

    bandout : Int, optional

        The band to be read if the grid file is `ERS`. The default is 0.

    Returns
    -------
    xd : xarray Dataset

        The data from `whizzFile`.

    """
    if whizzFile == '':
        whizzFile = Path(fb.get_grid_filename(filetypes=(('NetCdf4 grid', '*.nc'),
                                                         ('ERMapper grid', '*.ers'))))
    if whizzFile.suffix.upper() == '.ERS':
        e, n, z, datum, projection = ers.read_ers_image(whizzFile, bandout=bandout)
        xa = xr.DataArray(data=np.flip(z, 0), # DANGER!!!
                          dims=["y", "x"],
                          coords={"y": n, "x": e})
        xa.dropna(dim='y',how='all')
        xa.dropna(dim='x',how='all')
        fname = whizzFile.with_suffix('').name

        xa.attrs['x_channel'] = 'E'
        xa.attrs['y_channel'] = 'N'
        xa.attrs['long_name'] = fname
        xa.attrs['datum'] = datum
        xa.attrs['projection'] = projection
        if datum == 'WGS84' and projection == 'SUTM55':
            xa.rio.write_crs("epsg:32755", inplace=True)
    elif whizzFile.suffix.upper() == '.NC':
        # nc = nc4.Dataset(str(whizzFile), mode='r')
        xa = xr.load_dataarray(str(whizzFile))#xr.backends.NetCDF4DataStore(nc))
    else:
        print(f'ERROR: whizzFile has suffix {Path(whizzFile).suffix.upper()} but must be `.nc` or `.ers`.')
    return xa, whizzFile


