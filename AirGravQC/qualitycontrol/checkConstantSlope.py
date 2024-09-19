#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the required channels have a constant slope with fiducial.
"""
import numpy as np
import h5py

import AirGravQC.config as config
import AirGravQC.whizzFiles.retrieveData as rd

groupName = config.groupName


def checkConstantSlope(whizzFile, *, lines=[], channels=[]):
    """
    Checks for constant slope (`np.diff`) in all the given channels of data.
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    lines : Array{String}, optional
        Array of line numbers as strings. Default = [], meaning all lines are checked.
    channels : String List
        List of channel names from the database to be checked.

    Returns
    -------
    None

    """
    filename = str(whizzFile)
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if channels == []:
            lineGroups = list(g.values())
            channels = list(lineGroups[0].keys())
        data_is_good = True
        report = ''

        if lines == []:
            lines = list(g.keys())

        for line in lines:
            for channel in channels:
                data = rd.getLineData(g[line], channel)
                deriv = np.diff(data, n = 1)
                mean_deriv = np.mean(deriv)
                deriv = deriv - mean_deriv
                if len(deriv) > 10:
                    extremum = np.max(deriv) if np.max(deriv) > -np.min(deriv) else -np.min(deriv)
                    if extremum > mean_deriv / 1000.0:
                        report += f'\n  {line}; {channel} Largest difference (= {extremum:.3g}) > 0.1% of mean difference (= {(mean_deriv / 100.0):.3g})'
                        data_is_good = False
                else:
                    report += f'\n  {line}; {channel} data length {len(deriv)+1} is too short for analysis.'
                    data_is_good = False
    if data_is_good:
        report += 'All channels tested were either constant or of constant slope for all lines tested.'
    print(report)
    return


