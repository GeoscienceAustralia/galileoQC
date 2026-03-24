#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Display airborne and ground gravity data along a flight-line.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import galileoQC.config as config
import galileoQC.whizzFiles.retrieveData as rd
import galileoQC.gridFiles.read_ers as grd
import galileoQC.whizzPlots.whizzPlot as wpl
import galileoQC.utility.utility as util

groupName = config.groupName
    
                   
def diffGroundGrid(whizzFile, whizzChannel, whizzLine, gridPath, plot_title='Ground & Airborne Comparison'):
    """
    Samples data from the grid file onto the line from whizzFile
    and displays those data and the whizzChannel on a plot.

    Parameters
    ----------
    whizzFile : String or pathlib Path

        Name of a HDF5 Whizz file, including path and extension.

    whizzChannel : String

        Name of the channel of line data to be compared.

    whizzLine : String

        The identifier (eg '10030.0') for the line from which the
        whizzChannel data are taken for comparison.

    gridPath : Path

        The path (PathLib) to the ERS gridfile from which data is
        to be extracted for comparison.

    plot_title : String, optional

        A title for the plot. Default is 'Ground & Airborne Comparison'.

    Returns
    -------
    None.

    """

    # retrieve the grid information
    eg, ng, zg, datum, projection = grd.read_ers_image(gridPath)
    ngmin = np.min(ng)
    ngmax = np.max(ng)
    ngd = np.abs(ng[1] - ng[0])
    egmin = np.min(eg)
    egmax = np.max(eg)
    egd = np.abs(eg[1] - eg[0])
    zg = zg[::-1, :]
    newChannelName = gridPath.stem
    
    print('\nGrid file read for channel ', newChannelName)
    print('  SW Corner = ',ngmin, egmin, '  NE Corner = ',ngmax, egmax, '. Spacings = ', ngd, egd)

    # retrieve the line information
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        y = f[groupName]['CoordinateFrame'].attrs['YChannel']
        em = np.array(g[whizzLine][x])
        nm = np.array(g[whizzLine][y])
        zm = np.array(g[whizzLine][whizzChannel])

    # Some useful summary info, particularly if the line is not inside
    # the grid area.
    lineNo = whizzLine
    lineText = 'Line ' + lineNo
    print(lineText)
    print('  North min, max ', np.min(nm), np.max(nm), '  East min, max ', np.min(em), np.max(em))

    # Check that the endpoints of the line are both within the grid area.
    i_float, i_disp = divmod(np.nanmin(nm) - ngmin, ngd)
    j_float, j_disp = divmod(np.nanmin(em) - egmin, egd)
    if i_float < 0 or j_float < 0 or i_float > len(ng) or j_float > len(eg):
        print('ERROR - the flight line is not wholly within the area of the grid')
        return
    i_float, i_disp = divmod(np.nanmax(nm) - ngmin, ngd)
    j_float, j_disp = divmod(np.nanmax(em) - egmin, egd)
    if i_float < 0 or j_float < 0 or i_float > len(ng) or j_float > len(eg):
        print('ERROR - the flight line is not wholly within the area of the grid')
        return

    # A point on the line is close to the projection of a grid cell centre if
    # it is closer than eps (half the sample spacing along the line).
    sample_spacing = util._distance(em[1] - em[0], nm[1] - nm[0])
    eps = sample_spacing / 2.0

    # initialise the arrays to store the sampled data
    num_grid_samples = int(float(len(em)) * np.sqrt(2) * sample_spacing / max(ngd, egd))
    ns = np.zeros((num_grid_samples,))
    es = np.zeros((num_grid_samples,))
    zs = np.zeros((num_grid_samples,))
    diffs = np.zeros((num_grid_samples,))
    m = 0

    # traverse along the line extracting grid data as it occurs
    for k in range(0, len(em)):
        # find the closest grid cell centres (i,j), (i+1,j), (i,j+1), (i+1,j+1)
        i_float, i_disp = divmod(nm[k] - ngmin, ngd)
        j_float, j_disp = divmod(em[k] - egmin, egd)

        if i_float < 0 or j_float < 0 or i_float > len(ng) or j_float > len(eg):
            print('ERROR - the flight line is not wholly within the area of the grid')
            return

        # if the nearest grid cell centre is on the line, use the grid data
        if i_disp < eps and j_disp < eps:
            i0 = int(i_float)
            j0 = int(j_float)
            ns[m] = nm[k]
            es[m] = em[k]
            zs[m] = zg[i0,j0]
            diffs[m] = zs[m] - zm[k]
            m += 1
        # if the nearest grid cell centre is on the same northing, interpolate the grid data over easting
        elif i_disp < eps:
            i0 = int(i_float)
            j0 = int(j_float)
            j1 = j0 + 1
            ns[m] = nm[k]
            es[m] = em[k]
            zs[m] = zg[i0,j0] + (zg[i0,j1] - zg[i0,j0]) * (em[k] - eg[j0]) / (eg[j1] - eg[j0])
            diffs[m] = zs[m] - zm[k]
            m += 1
        # similarly ... for same easting, interpolate over northing
        elif j_disp < eps:
            i0 = int(i_float)
            i1 = i0 + 1
            j0 = int(j_float)
            ns[m] = nm[k]
            es[m] = em[k]
            zs[m] = zg[i0,j0] + (zg[i1,j0] - zg[i0,j0]) * (nm[k] - ng[i0]) / (ng[i1] - ng[i0])
            diffs[m] = zs[m] - zm[k]
            m += 1

    # clean up a couple of things before plotting
    if (np.ptp(em) > np.ptp(nm)):
        plotx_s = es[:m-1]
        plotx_m = em
    else:
        plotx_s = ns[:m-1]
        plotx_m = nm
    plotz_s = zs[:m-1]
    plotz_m = zm * 10.0
    # results ...
    print(f'Stdev(diff) = {np.nanstd(diffs):.1f}')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(plotx_s, plotz_s, 'b', label=newChannelName)
    ax.plot(plotx_m, plotz_m, 'g', label=lineText + ' ' + whizzChannel)
    ax.set_title(plot_title)
    ax.set_xlabel('Easting [m]')
    ax.set_ylabel('Bouguer Gravity [um/s/s]')
    ax.legend()
    ax.grid()
