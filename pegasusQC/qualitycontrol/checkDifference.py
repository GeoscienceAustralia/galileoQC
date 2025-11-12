#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check that the fourth difference of the magnetics is within specification.
Author: Mark Helm Dransfield
Created: ca 2023
License: CC BY-SA
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
# from scipy.signal import butter, lfilter

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.utility.utility as util

groupName = config.groupName
                            

def checkDiff(whizzFile, diff_chan='', raw_chan='', 
    lines=[], xChannel='', yChannel='', tChannel='', 
    limit=0.02, nSamples=3000, maxDuration=0.0, maxDistance=0.0, num_diff=4,
    plot_flag=False, verbose=False
):
    """
    Checks the (usually: total magnetic field fourth) difference channel in
    a whizzFile against the specification that its peak to peak variation over
    a set number of samples must not exceed some peak value. If `diff_chan` is not
    available, then the raw_chan channel may be used.

    If `maxDuration` is greater than 1.0, then that number of seconds
    will be used instead of `nSamples`. If `maxDuration` is less than 1.0
    and `maxDistance` is greater than 1.0, then that distance in metres
    will be used instead of `nSamples`.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    diff_chan : String, optional
        The name of the channel in whizzFile containing the differenced mag data.
        Default is '', in which case raw_chan is used.
    raw_chan : String, optional
        The name of the channel in whizzFile containing the raw magnetic data data.
        Default is ''; if both raw_chan and diff_chan are '', then an error is reported.
    lines : Array{String}, optional
        Array of line numbers as strings. Default = [], meaning all lines are checked.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    tChannel : String, optional
        The name of the geoWhizz field or channel containing the measured times. The
        default is to read the tChannel field name from the Coordinate Frame.
    limit : Float, optional
        The maximum allowed peak to peak variation. Default = 0.02
    nSamples : Integer, optional
        The number of samples (moving window) over which the test is applied.
        Default = 3000
    maxDuration : Float, optional
        The time in seconds (moving window) over which the test is applied.
        The default is to ignore this parameter.
    maxDistance : Float, optional
        The distance in metres (moving window) over which the test is applied.
        The default is to ignore this parameter.
    num_diff : Int, optional
        The number of differences to apply to `raw_chan` if used. The most common
        values are 4 or 8. Default 4.
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
    num_exceed_lines = 0
    num_failures = 0

    if maxDistance < 1.0 and maxDuration < 1.0:
        measure = "samples"
    elif maxDuration > 0.0:
        measure = "duration"
    elif maxDistance > 0.0:
        measure = "distance"
    else:
        print('ERROR - problem with nSamples.')
        return

    if diff_chan == '' and raw_chan == '':
        print('ERROR - no raw_chan, or difference channel, name supplied.')

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
            if verbose:
                report += f'{maxDuration} seconds evaluated as {nSamples} fids.'
        for line in lines:
            # Wind speed means different lines are at different speeds so
            # calculate nSamples from duration for each line.
            if measure == "distance":
                x_deltas = np.diff(rd.getLineData(g[line], xChannel))
                y_deltas = np.diff(rd.getLineData(g[line], yChannel))
                dd = util._distance(x_deltas, y_deltas)
                nSamples = int(maxDistance / np.nanmean(dd))
                if verbose:
                    report += f''

            if diff_chan == '':
                mag = rd.getLineData(g[line], raw_chan)
                data = _expanded_diff(mag, num_diff)
                # data = np.append(np.append(md4[0:2],md4),md4[-3:-1])
                plotTitle = line + f' {num_diff}-th difference of ' + raw_chan + ' Range'
            else:
                data = rd.getLineData(g[line], diff_chan)
                plotTitle = line + ' ' + diff_chan + ' Range'
            
            num_exceedances_in_line, num_failures_on_line, report = _failsDiffTest(
                f'L {line}', num_diff, data, limit, nSamples, num_exceed_lines, report)

            if num_failures_on_line > 0:
                num_failed_lines += 1
                num_failures += num_failures_on_line
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                ax.plot(data)
                plt.title(plotTitle)
                plt.grid(True)
              
            elif num_exceedances_in_line > 0:
                if verbose:
                    report += f'\nLine {line}: {num_exceedances_in_line} exceedances,'
                    report += f' less than {nSamples} ({maxDistance} m -> {nSamples} fids) consecutively - PASS'
                
            else:
                aaa = 1

    report = f'  Checked {numLines} line(s), {num_failures} failure(s) over {num_failed_lines} failed line(s).\n' + report
    print(report)
    if plot_flag:
        plt.show()


def _expanded_diff(data, num):
    """
    Return the `num`-th difference of `data`. In a function to allow padding to `len(data)` if needed.
    """
    return np.diff(data, n=num)


def _odd(number):
    """
    True if `number` is odd, else False.
    """
    if number % 2 == 0:
        return False
    else:
        return True


def _failsDiffTest(lineName, num_diff, data, limit, nSamples, num_exceed_lines, report):
    num_fids_in_exceedance = 0
    num_failures_on_line = 0
    num_exceedances_in_line = 0
    exceedance_in_line = False

    for fid in range(0, len(data)):
        # There is an exceedance ...
        if np.abs(data[fid]) > limit:
            # If a new exceedance, then initialise variables;
            if num_fids_in_exceedance == 0:
                num_exceed_lines += 1
                num_fids_in_exceedance = 1
                diff_extreme = data[fid]
            # Else increment and update on the current exceedance.
            else:
                num_fids_in_exceedance += 1
                if np.abs(diff_extreme) < abs(data[fid]):
                    diff_extreme = data[fid]
        else:
            if num_fids_in_exceedance > 0: # the current exceedance has ended
                if num_fids_in_exceedance > nSamples:
                    report += f'\n{lineName}: {num_diff}-th difference > {limit}'
                    report += f' for {num_fids_in_exceedance} > {nSamples} fids;'
                    report += f' peak exceedance = {diff_extreme:.3f} - FAIL'
                    exceedance_in_line = True
                    num_failures_on_line += 1
                num_exceedances_in_line += num_fids_in_exceedance
                num_fids_in_exceedance = 0
    return num_exceedances_in_line, num_failures_on_line, report

