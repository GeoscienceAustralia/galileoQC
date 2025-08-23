#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Report the Falcon difference noise between A and B complement data.
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.utility.utility as util

groupName = config.groupName

def diffGravVturb(whizzFile, turbulence, aD, bD, error_spec=5.0, low_cut=0.001, measX='', measY=''):
    """
    For a Falcon AGG. For each line, reports the gD difference noise,
    stdev(aD-bD)/2, and plots this as a scatter plot against the mean turbulence.
    All lines with difference noise greater than `error_spec` are reported.
    
    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    turbulence : String
        The name of the channel containing the turbulence field.
    aD : String
        The name of the channel containing the A complement gravity data.
    bD : String
        The name of the channel containing the B complement gravity data.
    error_spec : Float, optional
        The value above which the difference noise is excessive and is
        reported. Default 5.0.
    low_cut : Float, optional
        The low frequency cut-off frequency (in Hz) for the band-pass filtering
        applied before differencing. Default 0.001 (ie 1 mHz).

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    num_lines_failed = 0
    report = 'diffGravVturb() estimates the noise in (A+B)/2 as the stdev(A-B)/2.\n'
    period = 1.0 / low_cut
    wavelength = 60.0 * period / 1000.0 # km
    report += f'The input data are band-pass filtered at [{period}, 1.0] sec or [{wavelength}, 0.06] km at 60m/s.\n'
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        numLines = len(g.items())
        turbMean = np.zeros((numLines,))
        errmean = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        line_strs = list(g.keys())
        flightLine = line_strs[0]
        turb_units = g[flightLine][turbulence].attrs['Units']
        err_units = g[flightLine][aD].attrs['Units']

        if measX == '':
            measX = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if measY == '':
            measY = f[groupName]['CoordinateFrame'].attrs['YChannel']


        failed_lines = 0
        for line in g.keys():
            xM = rd.getLineData(g[line], measX)
            yM = rd.getLineData(g[line], measY)
            line_length = util._displacement2(xM[0], xM[-1], yM[0], yM[-1]) / 1000.0

            turb = rd.getLineData(g[line], turbulence)
            A_d = rd.getLineData(g[line], aD)
            B_d = rd.getLineData(g[line], bD)
            idx = np.where(~np.isnan(A_d + B_d))
            Ad = _butter_bandpass_filter(A_d[idx], low_cut, 1.0, 8.0, order=3)
            Bd = _butter_bandpass_filter(B_d[idx], low_cut, 1.0, 8.0, order=3)
            turb_clean = turb[idx]
            err_data = (Ad - Bd)/2.0
            err_data = err_data - np.mean(err_data)

            turbMean[count] = np.mean(turb_clean)
            errmean[count] = np.std(err_data)
            
            if errmean[count] > error_spec:
                report += f'{line} (length {line_length:6.1f} km) fails with noise {errmean[count]:.1f} > {error_spec}.\n'
                num_lines_failed += 1

            count += 1

        fig = plt.figure()
        fig.suptitle(f'Noise vs Turbulence - {projName}', fontsize=12)
        fig.subplots_adjust(top=0.85)
        ax = fig.add_subplot(1,1,1)
        neplot, = ax.plot(turbMean, errmean, 'bo', label='$E_{D}$')
        for ii in range(0, len(turbMean)):
            ax.text(turbMean[ii], errmean[ii], f'{line_strs[ii]}', fontsize=8)
        plt.ylabel(f'Difference Noise [{err_units}]', fontsize = 8)
        plt.xlabel(f'Turbulence [{turb_units}]', fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(8)
        for label in ax.get_yticklabels(): label.set_fontsize(8)
        fig.tight_layout()
        plt.show()
        report = f'{num_lines_failed} failed of {count} lines total.' + report
    print(report)
