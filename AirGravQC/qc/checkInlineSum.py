import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import AirGravQC.config as config
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.utility.utility as util

groupName = config.groupName


def checkInlineSum(whizzFile, inline1='', inline2='', inline3='', dontfilter=False, verbose=False):
    """
    Estimates the inline sum for each sample on each line in an FTG whizzFile.
    Plots the min, max, mean and stdev (units of eotvos) for each survey line
    as a function of line.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.

    Returns
    -------
    None

    """
    if inline1 == '':
        inline1 = 'Inline1_raw'
    if inline2 == '':
        inline2 = 'Inline2_raw'
    if inline3 == '':
        inline3 = 'Inline3_raw'

    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']

        numLines = len(g.items())
        chMin = np.zeros((numLines,))
        chMax = np.zeros((numLines,))
        chMean = np.zeros((numLines,))
        chStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            data1 = rd.getLineData(g[line], inline1)
            data2 = rd.getLineData(g[line], inline2)
            data3 = rd.getLineData(g[line], inline3)
            ils_BP = util._inLineSum(data1, data2, data3, dontfilter=dontfilter)
            if verbose:
                print(f'Line {line}, standard deviation of band-pass filtered in-line sum = {np.std(ils_BP):.2g}')
            lineNo[count] = line
            chMin[count] = np.min(ils_BP)
            chMax[count] = np.max(ils_BP)
            chMean[count] = np.mean(ils_BP)
            chStd[count] = np.std(ils_BP)
            count += 1

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(lineNo, chMin, 'bo', mfc='w')
        ax.plot(lineNo, chMax, 'bo', mfc='w')
        ax.errorbar(lineNo, chMean, chStd, capsize=3, marker='s', c='blue')
        plt.ylabel('Inline Sum', fontsize = 6)
        plotTitle = groupName + ': ' + ' Inline Sum Stats'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        plt.show()

