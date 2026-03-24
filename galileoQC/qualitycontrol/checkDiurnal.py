#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check magnetic diurnal data against specification.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt

import galileoQC.config as config
import galileoQC.whizzFiles.retrieveData as rd
import galileoQC.utility.utility as util

groupName = config.groupName


def checkDiurnal(whizzFile, basemag, lines=[], rangeLimit = 5.0, nSamples = 3000, xChannel='', yChannel='', tChannel='', 
    maxDuration=0.0, maxDistance=0.0, plot_flag=False, verbose=False):
    """
    Checks the `basemag` data for diurnal exceedances. These occur at any data
    value whose difference from a chord `nSamples` long is greater than `rangeLimit`.

    If `maxDuration` is greater than 1.0, then that number of seconds
    will be used instead of `nSamples`. If `maxDuration` is less than 1.0
    and `maxDistance` is greater than 1.0, then that distance in metres
    will be used instead of `nSamples`.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path

        The pathlib Path to the Whizz HDF5 file containing the survey line data.

    basemag : String

        The name of the channel in whizzFile containing the mag data to be checked.

    lines : Array{String}, optional

        Array of line numbers. Default = [], meaning all lines are checked.

    rangeLimit : Float, optional

        The maximum allowed deviation from the straight line chord. Default = 5.0 nT

    nSamples : Integer, optional

        The number of samples (moving window) over which the test is applied.
        Default = 3000

    xChannel : String, optional

        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.

    yChannel : String, optional

        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.

    tChannel : String, optional

        The name of the geoWhizz field or channel containing the measured times. The
        default is to read the tChannel field name from the Coordinate Frame.

    maxDuration : Float, optional

        The time in seconds (moving window) over which the test is applied.
        The default is to ignore this parameter.

    maxDistance : Float, optional

        The distance in metres (moving window) over which the test is applied.
        The default is to ignore this parameter.

    plot_flag : Bool, optional

        If True, all plots are generated.

    verbose : Bool, optional

        If True, a more verbose output is provided.

    Returns
    -------
    None

    """
    
    filename = str(whizzFile)
    report = ''
    num_failed_lines = 0

    if maxDistance < 1.0 and maxDuration < 1.0:
        measure = "samples"
    elif maxDuration > 0.0:
        measure = "duration"
    elif maxDistance > 0.0:
        measure = "distance"
    else:
        print('ERROR - problem with nSamples.')
        return

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if lines == []:
            lines = g.keys()
        numLines = len(lines)
        if measure == "distance":
            if xChannel == '':
                xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
            if yChannel == '':
                yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if measure == "duration":
            if tChannel == '':
                tChannel = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
            t = rd.getLineData(g[lines[0]], tChannel)
            nSamples = int(maxDuration / (t[1] - t[0]))

        for line in lines:
            # Wind speed means different lines are at different speeds so
            # calculate nSamples from duration for each line.
            if measure == "distance":
                x_deltas = np.diff(rd.getLineData(gLines[line], xChannel))
                y_deltas = np.diff(rd.getLineData(gLines[line], yChannel))
                dd = util._distance(x_deltas, y_deltas)
                nSamples = int(maxDistance / np.mean(dd))

            diurnalExceeded = False
            failedSample = 0
            bigExtremum = 0.0
            data = rd.getLineData(g[line], basemag)
            data = data[np.logical_not(np.isnan(data))]
                
            if nSamples > len(data):
                print(f'\n  Short line: {len(data)} < {nSamples}.')
                nSamples = len(data)
            if (nSamples % 2) == 0:
                nSamples = nSamples - 1
            nSam = (nSamples - 1) // 2

            for ii in range(nSam, len(data)-nSam):
                localData = data[ii-nSam:ii+nSam]
                localSlope = (localData[-1] - localData[0]) / nSamples
                deviation = localData - localSlope * range(0, len(localData)) - localData[0]
                extremum = np.max(deviation) if np.max(deviation) > -np.min(deviation) else -np.min(deviation)
                if extremum > rangeLimit:
                    diurnalExceeded = True
                    if extremum > bigExtremum:
                        bigExtremum = extremum
                        failedSample = ii
                    
            if diurnalExceeded:
                report += f'\n  L {line}: Diurnal for {basemag} at sample number {failedSample} diverges from chord by {bigExtremum:.2f},'
                report += f'\n  exceeding {rangeLimit:.1f} - FAIL'
                num_failed_lines += 1
                if plot_flag:
                    fig = plt.figure()
                    ax = fig.add_subplot(1,1,1)
                    ax.plot(data)
                    plotTitle = f'Line {line} Channel {basemag}: reaches {bigExtremum:.2f} at {failedSample}, exceeding {rangeLimit} - FAIL'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax.get_xticklabels(): label.set_fontsize(6)
                    for label in ax.get_yticklabels(): label.set_fontsize(6)
                    fig.tight_layout()

    print(f'  Checked {numLines} lines, {num_failed_lines} failed.\n')
    print(report)
    if plot_flag and num_failed_lines > 0:
        plt.show()


