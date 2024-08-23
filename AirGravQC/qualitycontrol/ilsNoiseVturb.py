#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Report the FTG in-line sum noise.
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import AirGravQC.config as config
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.utility.utility as util

groupName = config.groupName
    

def ilsNoiseVturb(whizzFile, diagComponent1, diagComponent2, diagComponent3, noiseSpec=17.0, vertaccel='', vertvelocity='', vertdispl='', labelLines=False):
    """
    For a Bell Air-FTG. For each line, reports the standard deviation of the in-line sums,
    and plots these as a scatter plot against the standard deviation of the vertical
    acceleration (if the acceleration is not supplied, it is estimated as the time difference
    of the vertical velocity, or if that is not supplied, the second time difference of the
    vertical displacement).

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
    labelLines : Bool, optional
        if True, label (with the line number) all points on the plot where the
        line failed the specification. Defaults to False

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    report = ''

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        numLines = len(g.items())
        accStd = np.zeros((numLines,))
        ilsStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        labelx = []
        labely = []
        labelt = []
        failed_lines = 0
        count = 0

        for line in g.keys():
            linegroup = g[line]
            if vertaccel != '':
                accel = rd.getLineData(linegroup, vertaccel)
            elif vertvelocity != '':
                data = rd.getLineData(linegroup, vertvelocity)
                accel = np.diff(data, n = 1)
            elif vertdispl != '':
                data = rd.getLineData(linegroup, vertdispl)
                accel = np.diff(data, n = 2)
            else:
                print("ERROR - need one of vertical acceleration, velocity or displacement (height/altitude).")
                return

            accStd[count] = np.std(accel)

            data1 = rd.getLineData(linegroup, diagComponent1)
            data2 = rd.getLineData(linegroup, diagComponent2)
            data3 = rd.getLineData(linegroup, diagComponent3)
            ilsStd[count] = np.std(util._inLineSum(data1, data2, data3))
            if ilsStd[count] > noiseSpec:
                if labelLines:
                    labelx.append(accStd[count])
                    labely.append(ilsStd[count])
                    labelt.append(line)
                failed_lines += 1
                report += f'Line {line}: in-line sum = {ilsStd[count]:.1f} exceeds specification of {noiseSpec}.\n'
                # print(f'Line {line}: in-line sum = {ilsStd[count]:.1f} exceeds specification of {noiseSpec}.')
            lineNo[count] = line
            count += 1
        
        fig = plt.figure()
        fig.suptitle(f'In-line Sum Noise vs Turbulence - {projName}', fontsize=12)
        fig.subplots_adjust(top=0.85)
        ax = fig.add_subplot(1,1,1)
        ax.plot(accStd, ilsStd, 'go')
        if labelLines:
            for ii in range(failed_lines):
                plt.text(labelx[ii], labely[ii], labelt[ii], va='top', ha='right', size=6.0)
        plt.ylabel(f'Inline Sum [E]', fontsize = 8)
        plt.xlabel(f'Turbulence [m/s/s]', fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(8)
        for label in ax.get_yticklabels(): label.set_fontsize(8)
        fig.tight_layout()
        plt.show()
    print(report)

