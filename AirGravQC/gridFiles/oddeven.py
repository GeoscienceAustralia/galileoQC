import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from pathlib import Path
import xarray as xr
import netCDF4 as nc4
import filebrowser as fb
import rioxarray
import h5py
import pygmt
import matplotlib.ticker as tkr
import collections
import colorcet as cc

import AirGravQC.gridFiles.graphics as graphics
import AirGravQC.utility.utility as util
import AirGravQC.gridFiles.read_ers as ers
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.config as config
from AirGravQC.gridFiles.graphicsShaded import graphicsShaded
from AirGravQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from AirGravQC.gridFiles.xarray_to_grid import xarray_to_grid
from AirGravQC.gridFiles.grids_gmt import image_pygmt
from AirGravQC.gridFiles.xdImage import xdImage

groupName = config.groupName
projectName = config.projectName


"""
ODDEVENLINES_v2()
- spacing = calcLineSpacing()
- for trav in traverses
    meany = mean(y)
    spaces = int(meany / spacing)
    if spaces % 2 == 0:
        - assign trav to evens
    - else:
        - assign trav to odds
- analyse_odds_evens()

"""


def updateOddOrEven(whizzFile, lines=[], x='', y='', verbose=False):
    """
    Sorts traverse lines as 'odd' or 'even' and writes result as an attribute for each flight-line.

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
    verbose : Bool, optional
        The verbosity of output. Default False.

    Returns
    -------
    Nothing.

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r+') as f:
        lines_group = f[groupName]['Lines']
        if lines == []:
            lines = list(g.keys())
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']

        trav_oddness = True
        first = True

        for line in lines:
            linegroup = lines_group[line]
            if not (linegroup.attrs['LinePurpose'] == 'Traverse'):
                continue

            if first:
                first = False
                linegroup.attrs['OddNotEven'] = trav_oddness
                lastgroup = linegroup

            trav_spacing = linegroup.attrs['TravSpacing']
            if approx_spacing(linegroup, lastgroup, trav_spacing):
                trav_oddness = not trav_oddness
            elif not approx_spacing(linegroup, lastgroup, 0.0):
                if verbose:
                    print(f'Skipped line {line}: spacing to previous inconsistent with {trav_spacing}.')
                continue

            linegroup.attrs['OddNotEven'] = trav_oddness
            if verbose:
                print(f'Line {line}: OddNotEven = {trav_oddness}.')


def updateLineTracks(whizzFile, lines=[], x='', y='', trackError=5.0, verbose=False):
    """
    Writes the mean line track direction as an attribute for each flight-line,
    writes it as a line attribute and sets the `LinePurpose` attribute to 'Control'
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
                gg.attrs['LinePurpose'] = 'Traverse'
            elif np.round(gg.attrs['MeanTrack'] - traverseTrack, -1) == 180.0: # refine this test using trackError
                gg.attrs['LinePurpose'] = 'Traverse'
            elif np.round(gg.attrs['MeanTrack'] - controlsTrack, -1) == 0.0:
                gg.attrs['LinePurpose'] = 'Control'
            elif np.round(gg.attrs['MeanTrack'] - controlsTrack, -1) == 180.0:
                gg.attrs['LinePurpose'] = 'Control'
            else:
                gg.attrs['LinePurpose'] = 'Unknown'
                mean_track = gg.attrs['MeanTrack']
                print(f'Line {line} has unknown line purpose (track = {mean_track:.1f}).')
            if verbose:
                line_purpose = gg.attrs['LinePurpose']
                print(f'Line {line}) is {line_purpose} on track {mean_track:.1f}')


def calcMeanTrack(lineGroup, easting, northing):
    """
    Return the aircraft track for the flight-line. The track is the angle
    (in degrees) from North in a clockwise direction.

    Parameters
    ----------
    linegroup : HDF5 group
        A whizzFile line group.heir track set. Default all lines.
    easting : String
        The name of the x (easting) channel in `linegroup`.
    northing : String
        The name of the y (northing) channel in `linegroup`.

    Returns
    -------
    track : float
        The mean aircraft track angle east of north in [0, 180] degrees.

    """
    # calculate the by-sample track direction
    dx = np.diff(rd.getLineData(lineGroup, easting))
    dy = np.diff(rd.getLineData(lineGroup, northing))
    # arctan returns angle north from east and [-pi, pi] range
    theta = np.arctan2(dy, dx) * 180.0 / np.pi
    # use east from north and [0, 180] range.
    track = (90.0 - theta) % 180.0

    # trim start and end of vector by 5% to remove any residual
    # aircraft manoeuvres at start or end of flight-line.
    idx1 = int(len(track) / 20.0)
    idx2 = len(track) - idx1
    return np.mean(track[idx1:idx2])


def oddevenlines(whizz_file, channel, grid_space, oddlines=[], evenlines=[]):
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

    Returns
    -------
    None.

    """

    if oddlines.size == 0 or evenlines.size == 0:
        print('ERROR. Please specify the odd and even lines. Automatic sorting of the lines is not yet functional.')
        # oddlines, evenlines = _getOddEvenLines(whizz_file)
        return

    # Read data to xarrays.
    even_data = whizz_to_xarray(whizz_file, channel, n_chan='', e_chan='', lines=evenlines, remove_mean=False, diff_one=False)
    odd_data = whizz_to_xarray(whizz_file, channel, n_chan='', e_chan='', lines=oddlines, remove_mean=False, diff_one=False)

    # Find the coverage for each gridded data set
    _, even_region = xarray_to_grid(even_data, grid_space)
    _, odd_region = xarray_to_grid(odd_data, grid_space)

    # We are only interested in the statistics over the intersection of the regions.
    intersectregion = [
                        max(even_region[0], odd_region[0]),
                        min(even_region[1], odd_region[1]),
                        max(even_region[2], odd_region[2]),
                        min(even_region[3], odd_region[3]),
                        ]
    # Grid and difference the data sets
    even_grid, even_region = xarray_to_grid(even_data, grid_space, region=intersectregion)
    odd_grid, odd_region = xarray_to_grid(odd_data, grid_space, region=intersectregion)
    d_grid = even_grid - odd_grid

    # Subtraction does not preserve attributes
    d_grid.attrs['units'] = even_grid.attrs['units']
    d_grid.attrs['long_name'] = even_grid.attrs['long_name']
    d_grid.attrs['title'] = f"even minus odd : {even_grid.attrs['title']}"
    d_grid['x'].attrs['orig_name'] = even_grid['x'].attrs['orig_name']
    d_grid['y'].attrs['orig_name'] = even_grid['y'].attrs['orig_name']

    # Image and report statistics
    # image_pygmt(d_grid, intersectregion)

    xdImage(d_grid, d_grid.attrs['title'], colormap=cc.m_CET_L9, cmap_norm='nonorm', 
        minClip=np.nan, maxClip=np.nan, gridlines=True, cb_ticks='stats', nSigma=2,
        hs=True, azdeg=45, ax=None, clipTo3Std = True)

    print(f'RMS of odd-even difference = {d_grid.std().data.item():.2f} {d_grid.attrs["units"]}')
    print(f'Estimated noise in the data = {d_grid.std().data.item()/2.0:.2f} {d_grid.attrs["units"]}')


def _getOddEvenLines(whizz_file):


    filename = str(whizz_file)
    numevens = 0
    numodds = 0
    
    with h5py.File(filename, 'r') as f:
        lines_group = f[groupName]['Lines']
        lines = lines_group.keys()
        
        for line in lines:
            # Only want Traverse lines
            lineIsTraverse = False
            try:
                if lines_group[line].attrs['LinePurpose'] == 'Traverse':
                    lineIsTraverse = True
                else:
                    continue
            except:
                continue

            # 
    print(f'{numevens} even lines, {numodds} odd lines.')
    return odds, evens