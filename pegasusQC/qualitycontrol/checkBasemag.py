#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check base station magnetometer values for excess variation.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.utility.utility as util

groupName = config.groupName


def checkBasemag(whizzFile, basemag, peak = 0.5, nSamples = 3000, verbose=False):
    """
    Checks the basemag channel in a whizzFile against the specification that
    the peak to peak variation over a set number of samples must not exceed
    some peak value.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path

        The pathlib Path to the Whizz HDF5 file containing the survey line data.

    basemag : String

        The name of the channel in whizzFile containing the basemag data.

    peak : Float

        The maximum allowed peak to peak variation.

    nSamples : Integer

        The number of samples (moving window) over which the test is applied.

    Returns
    -------
    None

    """
    filename = str(whizzFile)
    report = ''

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        lines = list(g.keys())
        numLines = len(lines)
        numfail = 0

        for line in lines:
            dataFail = False
            # plotTitle = line + ' ' + basemag + ' Peak-to-peak'
            data = rd.getLineData(g[line], basemag)
            data = data[np.logical_not(np.isnan(data))]
            if len(data) < 2:
                print(line, ' insufficient data')
                continue

            if _peakToPeak(data) < 2.0 * peak:
                if verbose:
                    report += f'{line} passed easily.\n'
                continue
            
            if len(data) < nSamples:
                    if _peakToPeak(data) > 2.0 * peak:
                        dataFail = True
                        report += f'{line}  FAIL, peak to peak range = {_peakToPeak(data)}\n'
                        continue
                    else:
                        print(line, ' passed on one segment.')
                        continue
            else:
                for ii in range(0, len(data) - nSamples):
                    if _peakToPeak(data[ii:ii+nSamples]) > 2.0 * peak:
                        dataFail = True
                        print(line, ' FAIL, peak to peak range = ', _peakToPeak(data[ii:ii+nSamples]))
                        break
                if dataFail == False:
                    print(line, ' passed on many segments.')
                    continue

    print(report)
    return
            

def _peakToPeak(data):
    return np.max(data) - np.min(data)
    
