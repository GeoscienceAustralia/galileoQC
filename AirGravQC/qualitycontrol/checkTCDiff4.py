#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check that the fourth difference of the magnetics is within specification.
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
# from scipy.signal import butter, lfilter

import AirGravQC.config as config
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.utility.utility as util

groupName = config.groupName
                            

def checkTCDiff4(whizzFile, TCDiff4='', rawMag='', lines=[], limit = 0.02, nSamples = 3000, plot_flag = False, verbose=False):
    """
    Checks the total magnetic field fourth difference channel in a whizzFile
    against the specification that the peak to peak variation over a set
    number of samples must not exceed some peak value. If `TCDiff4` is not
    available, then the rawMag channel may be used.

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
    limit : Float, optional
        The maximum allowed peak to peak variation. Default = 0.02
    nSamples : Integer, optional
        The number of samples (moving window) over which the test is applied.
        Default = 3000
    plot_flag : Bool, optional
        If True, all plots are generated.

    Returns
    -------
    None

    """
    
    filename = str(whizzFile)
    report = ''
    num_failed_lines = 0

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if lines == []:
            lines = g.keys()
        numLines = len(lines)
        for line in lines:
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


