#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check for spikes in the data (not recommended).
"""
import numpy as np
import h5py

import AirGravQC.config as config

groupName = config.groupName

def checkSpikes(whizzFile, channels=[], lines=[], numStd = 8.0, window=0, verbose=False):
    """
    Checks for spikes in all the given channels of data.
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    channels : String List
        List of field names from the database to be checked.
    lines : Array{String}, optional
        Array of line numbers as strings. Default = [], meaning all lines are checked.
    numStd : Float, optional
        maximum allowed number of standard deviations allowed. The default is 8.0.

    Returns
    -------
    None

    """
    
    filename = str(whizzFile)
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if lines == []:
            lines = list(g.keys())
        numLines = len(lines)
        if channels == []:
            lineGroups = list(g.values())
            channels = list(lineGroups[0].keys())
        numChannels = len(channels)
        report = ''
        for line in lines:
            for channel in channels:
                report += spike_test2(g[line][channel], numStd, f'\n  {line}; {channel}', window=window, verbose=verbose)

        if report == '':
            report = f'Checked {numChannels} channels in {numLines} lines and no spikes found.'
    print(report)
    return


def spike_test2(data, threshold, report_start, window=0, verbose=False):
    """
    If, in a window of width `window`, a data value is greater than threshold x standard deviation, then we have a spike.
    """
    if window == 0:
        window = len(data)

    report = ''
    if len(data) < window:
        if verbose:
            return f'{report_start} data length {len(data)+1} is too short for {window} window.'
        else:
            return ''

    data_range = np.max(data) - np.min(data)
    if data_range == 0.0:
        if verbose:
            return f'{report_start} data range {data_range} is too small for analysis.'
        else:
            return ''

    too_small = data_range / 1000.0
    myStd = np.std(data)
    if myStd < too_small:
        if verbose:
            return f'{report_start} data standard deviation {myStd} is too small for analysis.'
        else:
            return ''

    step = window // 2
    for ii in range(0, len(data)-window, step):

        deriv = data[ii:ii+window] - np.mean(data[ii:ii+window])
        extremum = np.max(deriv) if np.max(deriv) > -np.min(deriv) else -np.min(deriv)
        if extremum > threshold * myStd:
            report += f'{report_start} Extremum: {extremum:.3g} > {(threshold * myStd):.3g} = {threshold:.3g} x STD of {myStd:.3g}'

    return report


def spike_test1(data, threshold, report_start, verbose=False):
    """
    If the extremum of the first difference of any pair of samples > threshold x the standard deviation of the first differences,
    then we have a spike (or step).
    """
    report = ''
    deriv = np.diff(data, n = 1)
    deriv = deriv - np.mean(deriv)
    if len(deriv) > 10:
        myStd = np.std(deriv)
        extremum = np.max(deriv) if np.max(deriv) > -np.min(deriv) else -np.min(deriv)
        if extremum > threshold * myStd:
            report += f'{report_start} Extremum: {extremum:.3g} > {(threshold * myStd):.3g} = {threshold:.3g} x STD of {myStd:.3g}'
    else:
        report = f'{report_start} data length {len(deriv)+1} is too short for analysis.'
    return report


