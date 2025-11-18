#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate statistical plots for the difference between two channels.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
# from scipy.signal import butter, lfilter
from matplotlib.ticker import StrMethodFormatter
import matplotlib.ticker as tkr

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.gridFiles.read_ers as ers
import pegasusQC.gridFiles.gridfiles as grd
import pegasusQC.utility.utility as util
from pegasusQC.whizzPlots.plotBoxWhisker import plotBoxWhisker
import pegasusQC.whizzPlots.whizzPlot as wpl

groupName = config.groupName


def diffChanStats(whizzFile, channel1, channel2, lines=[], verbose=True):
    """
    Generate statistical plots for the difference between two channels across all lines. The plots show
    the min, mean, max and stdev for each channel as a function of line number.

    Parameters
    ----------
    whizzFile : String or pathlib Path

        Name of a HDF5 Whizz file, including path and extension.

    channel1 : String

        A channel to be subtracted from.

    channel2 : String

        A channel to subtract.

    lines : Array{String}, optional

        Array of line numbers as strings. Default = [], meaning all lines are checked.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        gProject = f[groupName]
        g = gProject['Lines']
        if lines == []:
            lines = list(g.keys())

        numLines = len(lines)        
        allmean = 0.0
        numsamp = 0
        chMin = np.zeros((numLines,))
        chMax = np.zeros((numLines,))
        chMean = np.zeros((numLines,))
        chStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0
        # build a y label
        dd = g[lines[0]][channel1]
        xlabelstr = 'Line Number'
        ylabelstr = channel1 + ' - ' + channel2
        if 'Units' in dd.attrs.keys():
            ylabelstr += ' ' + dd.attrs['Units']
    
        report = ''

        for line in lines:
            # if line != 'CoordinateFrame':
            lineNo[count] = line
            mydata = g[line][channel1][:] - g[line][channel2][:]
            if np.sum(~np.isnan(mydata)) > 3:
                chMin[count] = np.nanmin(mydata)
                chMax[count] = np.nanmax(mydata)
                chMean[count] = np.nanmean(mydata)
                chStd[count] = np.nanstd(mydata)
                allmean += np.nansum(mydata)
                numsamp += len(mydata)
            else:
                chMin[count] = 0.0
                chMax[count] = 0.0
                chMean[count] = 0.0
                chStd[count] = 0.0
                report += f'Less than three real values in {lineNo[count]:.2f} for data, no statistics.\n'
            count += 1

        print(f'Overall mean difference is {allmean / numsamp :.2f}.')
        if verbose:
            print(report)
        figtitle = wpl.make_plot_title(gProject)
        titlestr = channel1 + ' - ' + channel2 + ' Stats'
        plotBoxWhisker(chMin, chMax, chMean, chStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)
    return
    
 