#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check that the fourth difference of the magnetics is within specification.
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
# from scipy.signal import butter, lfilter

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.utility.utility as util

groupName = config.groupName
                            

def checkTCDiff4(whizzFile, TCDiff4='', rawMag='', lines=[], limit=0.02, nSamples=3000, xChannel='', yChannel='', tChannel='', 
    maxDuration=0.0, maxDistance=0.0, plot_flag=False, verbose=False):
    """
    Checks the total magnetic field fourth difference channel in a whizzFile
    against the specification that the peak to peak variation over a set
    number of samples must not exceed some peak value. If `TCDiff4` is not
    available, then the rawMag channel may be used.

    If `maxDuration` is greater than 1.0, then that number of seconds
    will be used instead of `nSamples`. If `maxDuration` is less than 1.0
    and `maxDistance` is greater than 1.0, then that distance in metres
    will be used instead of `nSamples`.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    TCDiff4 : String, optional
        The name of the channel in whizzFile containing the 4th difference mag data.
        Default is '', in which case rawMag is used.
    rawMag : String, optional
        The name of the channel in whizzFile containing the raw magnetic data data.
        Default is ''; if both rawMag and TCDiff4 are '', then an error is reported.
    lines : Array{String}, optional
        Array of line numbers as strings. Default = [], meaning all lines are checked.
    limit : Float, optional
        The maximum allowed peak to peak variation. Default = 0.02
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
            lines = list(g.keys())
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
                x_deltas = np.diff(rd.getLineData(g[line], xChannel))
                y_deltas = np.diff(rd.getLineData(g[line], yChannel))
                dd = util._distance(x_deltas, y_deltas)
                nSamples = int(maxDistance / np.mean(dd))

            if TCDiff4 == '':
                if rawMag == '':
                    print('ERROR - no rawmag or 4th difference channel name supplied.')
                else:
                    mag = rd.getLineData(g[line], rawMag)
                    md4 = np.diff(mag, n=4)
                    data = np.append(np.append(md4[0:2],md4),md4[-3:-1])
                    plotTitle = line + ' 4th difference of ' + rawMag + ' Range'
            else:
                data = rd.getLineData(g[line], TCDiff4)
                plotTitle = line + ' ' + TCDiff4 + ' Range'
            rangeTooHigh = False
            
            rangeTooHigh, numExceedances = util._failsDeviation(data, limit, nSamples)
            
            if rangeTooHigh:
                num_failed_lines += 1
                report += f'\n  4th difference on line {line} fails.'
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                ax.plot(data)
                plt.title(plotTitle)
                plt.grid(True)
                plt.show()
                report += f'\nLine {line}: exceedances = {numExceedances} > {nSamples} - FAIL'
              
            elif numExceedances > 0:
                report += f'\nLine {line}: exceedances = {numExceedances} < {nSamples} - PASS'
                if plot_flag:
                    fig = plt.figure()
                    ax = fig.add_subplot(1,1,1)
                    ax.plot(data)
                    plt.title(plotTitle)
                    plt.grid(True)
                    plt.show()
                
            else:
                aaa = 1

    print(f'  Checked {numLines} lines, {num_failed_lines} failed.\n')
    if verbose:
        print(report)


