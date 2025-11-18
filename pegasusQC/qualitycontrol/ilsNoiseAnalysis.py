#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Report the FTG in-line sum noise.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.utility.utility as util
from pegasusQC.qualitycontrol.psdChannelDiff import _time_frequency
from pegasusQC.whizzPlots.plotBoxWhisker import plotBoxWhisker


groupName = config.groupName
    

def ilsNoiseAnalysis(whizzFile, diagComponent1, diagComponent2, diagComponent3, noiseSpec=17.0, vertaccel='', vertvelocity='', vertdispl='', lines=[], labelLines=False, lowcut=0.03, highcut=0.1, dontfilter=True, verbose=False):
    """
    For an FTG. For each line, reports the standard deviation of the in-line sums,
    and plots these as a scatter plot against the standard deviation of the vertical
    acceleration (if the acceleration is not supplied, it is estimated as the time difference
    of the vertical velocity, or if that is not supplied, the second time difference of the
    vertical displacement). Also plots the min, max, mean and stdev (units of eotvos) for each survey line
    as a function of line.

    Parameters
    ----------
    whizzFile : String or pathlib Path

        Name of a HDF5 Whizz file, including path and extension.

    diagComponent1 : String

        The name of the channel containing the first tensor diagonal (in-line) component.

    diagComponent2 : String

        The name of the channel containing the second tensor diagonal (in-line) component.

    diagComponent3 : String

        The name of the channel containing the third tensor diagonal (in-line) component.

    noiseSpec : Float, optional

        The noise specification (largest allowed in-line sum for any flight line). Default 17.0 E.

    vertaccel : String, optional

        The name of the channel containing the vertical velocity field. Default ''.

    vertvelocity : String, optional

        The name of the channel containing the vertical velocity field. Default ''.

    vertdispl : String, optional

        The name of the channel containing the vertical velocity field. Default ''.

    lines : String list, optional

        The line numbers to be checked. Default is all lines in the whizzFile.

    labelLines : Bool, optional

        if True, label (with the line number) all points on the plot where the
        line failed the specification. Default False

    lowcut : Float, optional

        The low-pass frequency in Hz of the Butterworth filter applied if
        dontfilter == False. Default 0.03 Hz.

    highcut : Float, optional

        The high-pass frequency in Hz of the Butterworth filter applied if
        dontfilter == False. Default 0.1 Hz.

    dontfilter : Bool, optional

        If True, do not filter the data before calculating the in-line sum. Default True.

    verbose : Bool, optional

        If True, provide verbose printed reporting. Default False

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    report = ''

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        f_sample = _time_frequency(f[groupName])

        if lines == []:
            lines = list(g.keys())

        numLines = len(lines)
        chMin = np.zeros((numLines,))
        chMax = np.zeros((numLines,))
        chMean = np.zeros((numLines,))
        chStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        accStd = np.zeros((numLines,))
        ilsStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        labelx = []
        labely = []
        labelt = []
        failed_lines = 0
        count = 0

        for line in lines:
            linegroup = g[line]

            plot_ILSTurb = True
            if vertaccel != '':
                accel = rd.getLineData(linegroup, vertaccel)
            elif vertvelocity != '':
                data = rd.getLineData(linegroup, vertvelocity)
                accel = np.diff(data, n = 1)
            elif vertdispl != '':
                data = rd.getLineData(linegroup, vertdispl)
                accel = np.diff(data, n = 2)
            else:
                print("WARNING - no vertical acceleration, velocity or displacement (height/altitude) channel name supplied.")
                print("    Plot of in-line sum versus turbulence will not be produced.")
                plot_ILSTurb = False
                return

            accStd[count] = np.std(accel)

            data1 = rd.getLineData(linegroup, diagComponent1)
            data2 = rd.getLineData(linegroup, diagComponent2)
            data3 = rd.getLineData(linegroup, diagComponent3)
            ils = util._inLineSum(data1, data2, data3, fs=f_sample, lowcut=lowcut, highcut=highcut, dontfilter=dontfilter)

            lineNo[count] = line
            chMin[count] = np.min(ils)
            chMax[count] = np.max(ils)
            chMean[count] = np.mean(ils)
            chStd[count] = np.std(ils)
            if verbose:
                report += f'  Line {line}, standard deviation of band-pass filtered in-line sum = {chStd[count]:.2g}.\n'
            if plot_ILSTurb and chStd[count] > noiseSpec:
                if labelLines:
                    labelx.append(accStd[count])
                    labely.append(chStd[count])
                    labelt.append(line)
                failed_lines += 1
                report += f'Line {line}: in-line sum = {chStd[count]:.1f} exceeds specification of {noiseSpec}.\n'
            count += 1
        
        fig = plt.figure()
        fig.suptitle(f'In-line Sum Noise Analysis - {projName}', fontsize=10)
        fig.subplots_adjust(top=0.85)
        ax = fig.add_subplot(2,1,1)
        ax.plot(accStd, chStd, 'go')
        if labelLines:
            for ii in range(failed_lines):
                plt.text(labelx[ii], labely[ii], labelt[ii], va='top', ha='right', size=6.0)
        plt.ylabel(f'Inline Sum [E]', fontsize = 8)
        plt.xlabel(f'Turbulence [m/s/s]', fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(8)
        for label in ax.get_yticklabels(): label.set_fontsize(8)

        plotBoxWhisker(chMin, chMax, chMean, chStd, lineNo, '', 'In-line Sum Statistics', xlabelstr='Line Number', ylabelstr='Inline Sum [E]', xaxis='linenumber')

        fig.tight_layout()
        plt.show()
    print(report)

