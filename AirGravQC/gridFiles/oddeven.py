#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Perform odd-even analysis of a channel of data in a `whizz_file`.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from pathlib import Path
import xarray as xr
import netCDF4 as nc4
import filebrowser as fb
# import rioxarray
import h5py
import matplotlib.ticker as tkr
import collections

import AirGravQC.gridFiles.graphics as graphics
import AirGravQC.utility.utility as util
import AirGravQC.gridFiles.read_ers as ers
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.config as config
from AirGravQC.qualitycontrol.psdChannelDiff import _time_frequency
from AirGravQC.gridFiles.graphicsShaded import graphicsShaded
from AirGravQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from AirGravQC.gridFiles.xarray_to_grid import xarray_to_grid
from AirGravQC.gridFiles.xdImage import xdImage
import AirGravQC.gridFiles.gridutility as gut


groupName = config.groupName
projectName = config.projectName


def updateLineTracks(whizzFile, lines=[], x='', y='', trackError=5.0, verbose=False):
    """
    Writes the mean line track direction as an attribute for each flight-line,
    writes it as a line attribute and sets the `LineVariety` attribute to 'Control'
    or 'Traverse'.

    Parameters
    ----------
    whizzFile : Path or String
        The Path to, or String name of, the whizz file in HDF5 format.
    lines : String Array, optional
        List of lines to have their track set. Default all lines.
    x : String, optional
        The name of the x (easting) channel in `whizz_file`. Default is to use
        the `XChannel` attribute.
    y : String, optional
        The name of the y (northing) channel in `whizz_file`. Default is to use
        the `YChannel` attribute.
    trackError : Float, optional IGNORED FOR NOW!
        Any track with +/- trackError degrees of another track is assumed to
        have the same nominal direction. Default is 5.0 degrees.
    verbose : Bool, optional
        The verbosity of output. Default False.

    Returns
    -------
    Nothing.

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r+') as f:
        g = f[groupName]['Lines']
        if lines == []:
            lines = list(g.keys())
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']

        numlines = len(lines)
        tracks = np.zeros((numlines,))

        # assign mean track to line `MeanTrack` attribute
        for idx, line in enumerate(lines):
            mean_track = calcMeanTrack(g[line], x, y)
            g[line].attrs['MeanTrack'] = mean_track
            tracks[idx] = mean_track

        # Traverses are more common than controls
        tracksDict = collections.Counter(np.round(tracks, -1))
        traverseTrack = tracksDict.most_common()[0][0] % 180.0
        print(f'Traverse track = {traverseTrack:.1f} degrees (relative East).')
        controlsTrack = (traverseTrack - 90) % 180.0
        print(f'Controls track = {controlsTrack:.1f} degrees (relative East).')

        # set Control / Traverse attribute
        # WARNING - assigning all non-traverses to 'Control'
        for line in lines:
            gg = g[line]
            if np.round(gg.attrs['MeanTrack'] - traverseTrack, -1) == 0.0: # refine this test using trackError
                gg.attrs['LineVariety'] = 'Traverse'
            elif np.round(gg.attrs['MeanTrack'] - traverseTrack, -1) == 180.0: # refine this test using trackError
                gg.attrs['LineVariety'] = 'Traverse'
            elif np.round(gg.attrs['MeanTrack'] - controlsTrack, -1) == 0.0:
                gg.attrs['LineVariety'] = 'Control'
            elif np.round(gg.attrs['MeanTrack'] - controlsTrack, -1) == 180.0:
                gg.attrs['LineVariety'] = 'Control'
            else:
                gg.attrs['LineVariety'] = 'Unknown'
                mean_track = gg.attrs['MeanTrack']
                print(f'Line {line} has unknown line purpose (track = {mean_track:.1f}).')
            if verbose:
                line_purpose = gg.attrs['LineVariety']
                print(f'Line {line} is {line_purpose} on track {gg.attrs["MeanTrack"]:.1f}')


def calcMeanTrack(lineGroup, easting, northing):
    """
    Return the aircraft track for the flight-line. The track is the angle
    (in degrees) from North in a clockwise direction.

    Parameters
    ----------
    linegroup : HDF5 group
        A whizzFile line group.
    easting : String
        The name of the x (easting) channel in `linegroup`.
    northing : String
        The name of the y (northing) channel in `linegroup`.

    Returns
    -------
    track : float
        The mean aircraft track angle east of north in [0, 180] degrees.

    """

    # # calculate the by-sample track direction
    # dx = np.diff(rd.getLineData(lineGroup, easting))
    # dy = np.diff(rd.getLineData(lineGroup, northing))
    # # arctan returns angle north from east and [-pi, pi] range
    # theta = np.arctan2(dy, dx) * 180.0 / np.pi
    # # use east from north and [0, 180] range.
    # track = (90.0 - theta) % 180.0

    # # trim start and end of vector by 5% to remove any residual
    # # aircraft manoeuvres at start or end of flight-line.
    # idx1 = int(len(track) / 20.0)
    # idx2 = len(track) - idx1
    # return np.mean(track[idx1:idx2])


    # calculate the by-sample track direction
    dx = np.mean(np.diff(rd.getLineData(lineGroup, easting)))
    dy = np.mean(np.diff(rd.getLineData(lineGroup, northing)))
    # arctan returns angle north from east and [-pi, pi] range
    theta = np.arctan2(dy, dx) * 180.0 / np.pi
    # use east from north and [0, 180] range.
    track = (90.0 - theta) % 180.0
    return track


def oddevenlines(whizz_file, channel, grid_space, oddlines=[], evenlines=[], method='neighbours', mask_polygon=[], mask_pixels=1, numneighbours=1, bdist=None, maxiters=100, hs=True):
    """
    Performs odd-even analysis of the `channel` data in `whizz_file`. The data are
    sorted into two sets of odd and even lines. Each set is gridded and the difference
    of the grids is imaged, and its RMS value reported as an estimate of the error in
    the data.

    Parameters
    ----------
    whizz_file : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    channel : String
        The channel or field name to analyse. Must exist in `whizz_file`.
    grid_space : Float
        The width of the grid cell to be used in gridding. Recommend: 1/5 - 1/4 line spacing.
    oddlines : Array of String, optional
        An array of line numbers that will constitute the odd lines. The default is NOT WORKING! ...to take the first,
        and then every second traverse thereafter.
    evenlines : Array of String, optional
        An array of line numbers that will constitute the even lines. The default is NOT WORKING! ...to take every
        second traverse (alternates to the oddlines).
    method : string, optional
        The gridding algorithm to use in interpolating the data. Available are the Verde methods:
        "neighbours", "bicubic", and "biharmonic" as well as "minc" and the SciPy GridData "linear" method.
        The "neighbours" method is much faster if `pykdtree` is installed. Default `neighbours` method.
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
        Maximum number of iterations for minc method. The default is 100.
    hs : Bool, optional
        If True (the default), then the image is hill-shaded.
    Returns
    -------
    None.

    """

    if np.array(oddlines).size == 0 or np.array(evenlines).size == 0:
        oddlines, evenlines = _getOddEvenLines(str(whizz_file))

    if len(oddlines) == 0 or len(evenlines) == 0:
        print("ERROR - cannot proceed without odd and even lines!")
        return
    # Read data to xarrays.
    even_data = whizz_to_xarray(whizz_file, channel, n_chan='', e_chan='', lines=evenlines, remove_mean=False, diff_one=False)
    odd_data = whizz_to_xarray(whizz_file, channel, n_chan='', e_chan='', lines=oddlines, remove_mean=False, diff_one=False)

    # Find the coverage for each gridded data set
    e_chan = even_data.attrs['x_channel']
    n_chan = even_data.attrs['y_channel']
    even_region = [
                    np.min(even_data[e_chan].values),
                    np.max(even_data[e_chan].values),
                    np.min(even_data[n_chan].values),
                    np.max(even_data[n_chan].values)
                ]
    odd_region = [
                    np.min(odd_data[e_chan].values),
                    np.max(odd_data[e_chan].values),
                    np.min(odd_data[n_chan].values),
                    np.max(odd_data[n_chan].values)
                ]

    # We are only interested in the statistics over the intersection of the regions.
    intersectregion = [
                        max(even_region[0], odd_region[0]),
                        min(even_region[1], odd_region[1]),
                        max(even_region[2], odd_region[2]),
                        min(even_region[3], odd_region[3]),
                        ]
                        
    # Grid and difference the data sets
    even_grid, even_region = xarray_to_grid(even_data, grid_space, region=intersectregion, method=method, 
        mask_polygon=mask_polygon, mask_pixels=mask_pixels, numneighbours=numneighbours, bdist=bdist, maxiters=maxiters)
    odd_grid, odd_region = xarray_to_grid(odd_data, grid_space, region=intersectregion, method=method, 
        mask_polygon=mask_polygon, mask_pixels=mask_pixels, numneighbours=numneighbours, bdist=bdist, maxiters=maxiters)
    d_grid = even_grid - odd_grid

    # Subtraction does not preserve attributes
    d_grid.attrs['units'] = even_grid.attrs['units']
    d_grid.attrs['long_name'] = even_grid.attrs['long_name']
    d_grid.attrs['title'] = f"even minus odd : {even_grid.attrs['title']}"
    d_grid['x'].attrs['orig_name'] = even_grid['x'].attrs['orig_name']
    d_grid['y'].attrs['orig_name'] = even_grid['y'].attrs['orig_name']

    # Image and report statistics

    xdImage(d_grid, d_grid.attrs['title'], colormap=config.qc_colormap, cmap_norm='nonorm', 
        minClip=np.nan, maxClip=np.nan, gridlines=True, cb_ticks='stats', nSigma=2,
        hs=hs, azdeg=45, ax=None, clipTo3Std = True)

    gut.report_gridStats(d_grid, mask_polygon=mask_polygon)


def _getOddEvenLines(whizz_file):
    """
    Returns the flight-lines in `whizz_file` sorted into odd lines and even lines,
    according to the line attribute `Parity`.

    Parameters
    ----------
    whizz_file : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.

    Returns
    -------
    oddlines : list of string
        The Odd survey lines.
    evenlines : list of string
        The Even survey lines.

    """

    filename = str(whizz_file)
    numevens = 0
    evenlines = []
    numodds = 0
    oddlines = []
    
    with h5py.File(filename, 'r') as f:
        lines_group = f[groupName]['Lines']
        lines = lines_group.keys()
        
        for line in lines:
            # Only want Traverse lines
            lineIsTraverse = False
            try:
                if lines_group[line].attrs['LineVariety'] == 'Traverse':
                    lineIsTraverse = True
                else:
                    continue
            except:
                continue

            # 
            if "Parity" in lines_group[line].attrs.keys():
                if lines_group[line].attrs["Parity"]:
                    oddlines.append(line)
                    numodds += 1
                else:
                    evenlines.append(line)
                    numevens += 1

    print(f'{numevens} even lines, {numodds} odd lines.')
    return oddlines, evenlines


def _getPlannedLines(whizz_file):
    """
    Returns the flight-lines in `whizz_file` sorted into those with,
    and those without, a PlannedLine attribute set.

    Parameters
    ----------
    whizz_file : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.

    Returns
    -------
    oddlines : list of string
        The Odd survey lines.
    evenlines : list of string
        The Even survey lines.

    """

    filename = str(whizz_file)
    numplanned = 0
    plannedlines = []
    numunplanned = 0
    unplannedlines = []
    
    with h5py.File(filename, 'r') as f:
        lines_group = f[groupName]['Lines']
        lines = lines_group.keys()
        
        for line in lines:
            try:
                if 'PlannedLine' in lines_group[line].attrs:#.keys():
                    plannedlines.append(line)
                    numplanned += 1
                else:
                    unplannedlines.append(line)
                    numunplanned += 1
            except:
                continue

    print(f'{numplanned} planned lines, {numunplanned} unplanned lines, total {numplanned + numunplanned}.')
    return plannedlines, unplannedlines


def _getTravCtrlLines(whizz_file):
    """
    Returns the flight-lines in `whizz_file` sorted into traverse lines and control lines,
    according to the line attribute `LineVariety`.

    Parameters
    ----------
    whizz_file : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.

    Returns
    -------
    travlines : list of string
        The Traverse survey lines.
    ctrllines : list of string
        The Control survey lines.

    """

    filename = str(whizz_file)
    numtravs = 0
    travlines = []
    numctrls = 0
    ctrllines = []
    
    with h5py.File(filename, 'r') as f:
        lines_group = f[groupName]['Lines']
        lines = lines_group.keys()
        numlines = len(lines)
        
        for line in lines:
            if "LineVariety" in lines_group[line].attrs.keys():
                if lines_group[line].attrs["LineVariety"] == "Traverse":
                    travlines.append(line)
                    numtravs += 1
                elif lines_group[line].attrs["LineVariety"] == "Control":
                    ctrllines.append(line)
                    numctrls += 1

    print(f'{numtravs} traverse lines, {numctrls} control lines, {numlines - numctrls - numtravs} not classified.')
    return travlines, ctrllines


def altsample_grid(whizz_file, channel, filter_length, grid_space, method='neighbours', mask_polygon=[], mask_pixels=1, numneighbours=1, hs=True):
    """
    Performs odd-even analysis of the `channel` data in `whizz_file`. The data are
    sorted into two sets of odd and even lines. Each set is gridded and the difference
    of the grids is imaged, and its RMS value reported as an estimate of the error in
    the data.

    Parameters
    ----------
    whizz_file : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    channel : String
        The channel or field name to analyse. Must exist in `whizz_file`.
    filter_length : Float
        The length in seconds of the filter applied to channel.
    grid_space : Float
        The width of the grid cell to be used in gridding. Recommend: 1/5 - 1/4 line spacing.
    method : string, optional
        The gridding algorithm to use in interpolating the data. Available are the Verde methods:
        "neighbours", "bicubic", and "biharmonic" and the SciPy GridData "linear" method. "neighbours"
        is much faster if `pykdtree` is installed. Default `neighbours` method.
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

    filename = str(whizz_file)
    with h5py.File(filename, 'r') as f:
        fs = _time_frequency(f[groupName])
    half_fids = int(0.5 * fs * filter_length)
    odd_start = int(0.5 * half_fids)

    # Read data to xarrays.
    data = whizz_to_xarray(whizz_file, channel, n_chan='', e_chan='', lines=[], remove_mean=False, diff_one=False)
    odd_data = data.isel(fiducials=slice(odd_start,-1,half_fids))
    even_data = data.isel(fiducials=slice(0,-1,half_fids))

    # Find the coverage for each gridded data set
    e_chan = even_data.attrs['x_channel']
    n_chan = even_data.attrs['y_channel']
    even_region = [
                    np.min(even_data[e_chan].values),
                    np.max(even_data[e_chan].values),
                    np.min(even_data[n_chan].values),
                    np.max(even_data[n_chan].values)
                ]
    odd_region = [
                    np.min(odd_data[e_chan].values),
                    np.max(odd_data[e_chan].values),
                    np.min(odd_data[n_chan].values),
                    np.max(odd_data[n_chan].values)
                ]

    # We are only interested in the statistics over the intersection of the regions.
    intersectregion = [
                        max(even_region[0], odd_region[0]),
                        min(even_region[1], odd_region[1]),
                        max(even_region[2], odd_region[2]),
                        min(even_region[3], odd_region[3]),
                        ]
                        
    # Grid and difference the data sets
    even_grid, even_region = xarray_to_grid(even_data, grid_space, region=intersectregion, method=method, 
        mask_polygon=mask_polygon, mask_pixels=mask_pixels, numneighbours=numneighbours)
    odd_grid, odd_region = xarray_to_grid(odd_data, grid_space, region=intersectregion, method=method, 
        mask_polygon=mask_polygon, mask_pixels=mask_pixels, numneighbours=numneighbours)
    d_grid = even_grid - odd_grid

    # Subtraction does not preserve attributes
    d_grid.attrs['units'] = even_grid.attrs['units']
    d_grid.attrs['long_name'] = even_grid.attrs['long_name']
    d_grid.attrs['title'] = f"Grid difference of alternates : {even_grid.attrs['title']}"
    d_grid['x'].attrs['orig_name'] = even_grid['x'].attrs['orig_name']
    d_grid['y'].attrs['orig_name'] = even_grid['y'].attrs['orig_name']

    # Image and report statistics

    xdImage(d_grid, d_grid.attrs['title'], colormap=config.qc_colormap, cmap_norm='nonorm', 
        minClip=np.nan, maxClip=np.nan, gridlines=True, cb_ticks='stats', nSigma=2,
        hs=hs, azdeg=45, ax=None, clipTo3Std = True)

    gut.report_gridStats(d_grid, mask_polygon=mask_polygon)




