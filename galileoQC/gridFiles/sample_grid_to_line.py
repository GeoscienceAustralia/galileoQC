#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sample gridded data onto a located line whizzFile.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

# import necessary modules

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy import interpolate
import xarray
import pathlib
from scipy.interpolate import RegularGridInterpolator

import galileoQC.gridFiles.read_ers as grd
import galileoQC.whizzFiles.retrieveData as rd
import galileoQC.utility.utility as util
import galileoQC.config as config

groupName = config.groupName
projectName = config.projectName


def sample_grid_to_line(grid, hdfPath, lines=[]):
    """
    Interpolates the data in a grid file onto a new channel, with the same name
    as the grid stem, in the geoWhizz HDF5 file. Fills empty channel samples
    with nans.

    Parameters
    ----------
    grid : String or pathlib Path

        Name of a ERMapper file, including path and extension.

    hdfPath : String or pathlib Path

        Name of a HDF5 Whizz file, including path and extension.

    lines : Array{String}, optional

        Array of line numbers as strings. Default = [], meaning all lines are checked.

    Returns
    -------
    None.

    """
    hdfFile = str(hdfPath)
    # If we were passed a filename, get data from file
    if isinstance(grid, (pathlib.PurePath, str)):
        eg, ng, zg, datum, projection = grd.read_ers_image(grid)
        zg = zg[::-1, :]
        newChannelName = grid.stem

    # if we were passed a DataArray, extract its data
    elif type(grid) is xarray.core.dataarray.DataArray:
        # First, fill any holes in the rectangular grid array
        # rioxarray does this really well but it wants a CRS first
        grid.rio.write_crs(4283, inplace=True)
        grid.rio.write_nodata(np.nan, encoded=True, inplace=True)
        grid_fill = grid.rio.interpolate_na()

        # Then expand it to ensure flight-line data is inside grid
        pad_cells = 1
        grid_mean = np.nanmean(grid_fill.values)
        parr = grid_fill.pad(x=(pad_cells, pad_cells),
                     y=(pad_cells, pad_cells), 
                     keep_attrs=True,
                     mode="mean",
                     fill_value=grid_mean,
                     stat_length=None)

        eg = grid.x.data
        ng = grid.y.data
        zg = np.transpose(grid_fill.data)

        if 'name' in grid.attrs:
            newChannelName = grid.attrs['name']
        elif 'Name' in grid.attrs:
            newChannelName = grid.attrs['Name']
        elif 'long_name' in grid.attrs:
            newChannelName = grid.attrs['long_name']
        else:
            newChannelName = 'from_grid'
    else:
        print('\nERROR - type of grid not recognised.')
        print('    No grid sampled to database.')
        return

    ngmin = np.min(ng)
    ngmax = np.max(ng)
    ngd = np.abs(ng[1] - ng[0])
    egmin = np.min(eg)
    egmax = np.max(eg)
    egd = np.abs(eg[1] - eg[0])
    
    print('\nGrid file read for channel ', newChannelName)
    print("  corners (", egmin, ngmin, "), (", egmax,  ngmax, ")")
    print("  spacings ", ngd, egd)
    print("  shape ", eg.shape, ng.shape, zg.shape)

    # create interpolator for gridded data

    itp = RegularGridInterpolator((eg, ng), zg)#, method='cubic', bounds_error=False, fill_value=np.nan)

    with h5py.File(hdfFile, 'r+') as f:
        g = f[groupName]['Lines']
        lineGroups = list(g.values())
        channelNames = list(lineGroups[0].keys())
        x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        y = f[groupName]['CoordinateFrame'].attrs['YChannel']
        t = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
        height = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        if lines == []:
            lines = list(g.keys())

        for line in lines:
            lineText = 'Line ' + str(line)
            if newChannelName in g[line]:
                print(f'{lineText}: {newChannelName} already in whizz_file, not written.')
                continue

            # retrieve line positions
            em = rd.getLineData(g[line], x)
            nm = rd.getLineData(g[line], y)

            # The interpolator cannot extrapolate, so only use it inside the grid area
            condition_nmin = nm > ngmin 
            condition_nmax = nm < ngmax
            condition_nlim = condition_nmax * condition_nmin
            condition_emin = em > egmin 
            condition_emax = em < egmax
            condition_enlim = (condition_emax * condition_emin * condition_nlim).flatten()

            zm_interped = itp((em[condition_enlim], nm[condition_enlim]))

            # Then move the answers to a full size array, filling gaps with NaNs
            zm = np.zeros(em.shape)
            zm[np.where(~condition_enlim)] = np.nan
            zm[np.where(condition_enlim)] = zm_interped

            # and store in line sub-group
            dd = g[line].create_dataset(newChannelName, data = zm, compression="gzip", compression_opts=4)
            dd.attrs['Name'] = newChannelName

    return


