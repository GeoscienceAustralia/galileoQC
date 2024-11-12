#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot a statistical analysis.
"""
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr

import AirGravQC.utility.utility as util


def plotBoxWhisker(chMin, chMax, chMean, chStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr='Line Number', xaxis='linenumber'):
    """
	Creates a box whisker plot showing the range, standard deviation and mean for each flight-line.

    Parameters
    ----------
    chMin : Numpy array
        Contains the minimum value of some channel for all lines.
    chMax : Numpy array
        Contains the maximum value of some channel for all lines.
    chMean : Numpy array
        Contains the mean value of some channel for all lines.
    chStd : Numpy array
        Contains the standard deviation value of some channel for all lines.
    lineNo : Numpy array
        Contains line numbers.
    figtitle : String
        Contains figure title.
    titlestr : String
        Contains plot title.
    xlabelstr : String
        Contains x label.
    ylabelstr : String, optional
        Contains y label. Defaults to 'Line Number'.
    xaxis : String, optional
        If "linenumber" then the x axis variable is line number;
        if "ordinal" then it is just ordinal.

    Returns
    -------
    None.


    """
    if not xaxis in ['linenumber', 'ordinal']:
        print("Error. xaxis must be either 'linenumber' or 'ordinal' but was {xaxis}.")
        return
    fig = plt.figure()
    fig.suptitle(figtitle, fontsize=10)
    fig.subplots_adjust(top=0.85)
    thou_format = tkr.FuncFormatter(util._space_thou)

    ax = fig.add_subplot(1,1,1)
    if xaxis == 'linenumber':
        ax.plot(lineNo, chMin, 'bo', mfc='w')
        ax.plot(lineNo, chMax, 'bo', mfc='w')
        ax.errorbar(lineNo, chMean, chStd, capsize=3, marker='s', c='blue', linestyle='')
    else:
        ordinal = range(0,len(lineNo))
        ax.plot(ordinal, chMin, 'bo', mfc='w')
        ax.plot(ordinal, chMax, 'bo', mfc='w')
        ax.errorbar(ordinal, chMean, chStd, capsize=3, marker='s', c='blue', linestyle='')
        xlabelstr = 'line number ordinal'

    ax.xaxis.set_major_formatter(thou_format)
    plt.title(titlestr)
    plt.xlabel(xlabelstr, fontsize = 6)
    plt.ylabel(ylabelstr, fontsize = 6)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    plt.show()
    return        
    
 