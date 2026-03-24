#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the values of a channel for each survey flight-line.

Author: Mark Helm Dransfield

Created: 2023

License: CC BY-SA
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import xarray as xr
import verde as vd
import pooch

import galileoQC.config as config
# import galileoQC.qualitycontrol.qualityAnalysis as qc
import galileoQC.whizzFiles.retrieveData as rd
import galileoQC.whizzPlots.whizzPlot as wpl
import galileoQC.utility.utility as util
import matplotlib.ticker as tkr

groupName = config.groupName


def plotChannelsLines(whizzFile, channels, flightLines, x='', mean_remove=False, xOffset=True):
    """
    For the given channel in the geoWhizz file filename, plot the values for all
    specified flightlines versus x.

    Parameters
    ----------
    whizzFile : String or pathlib Path

        Name of a HDF5 Whizz file, including path and extension.

    channel : String

        The name of the channel or field to plot.

    flightLines : [String]

        A list of flightlines, e.g. ['1000110.0', '1000210.0', '1000310.0'] .

    x : String

        The name of the independent variable for the plot.

    xOffset : Bool, optional

        If True, map x to x - x[0] before plotting. The default is True.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
        plotPreTitle = wpl.make_plot_title(f[groupName])
        xDel = 0.0
        
        for line in flightLines:
            xData = rd.getLineData(g[line], x)
            xUnits = rd.getChannelAttrs(g[line], x)
            if xUnits == '':
                xlabelstr = f'{x}'
            else:
                xlabelstr = f'{x} [{xUnits}]'
            for channel in channels:
                fig = plt.figure(figsize=(6,9))
                thou_format = tkr.FuncFormatter(util._space_thou)
                ax = fig.add_subplot(1,1,1)
                yData = rd.getLineData(g[line], channel)
                yUnits = rd.getChannelAttrs(g[line], channel)
                if yUnits == '':
                    ylabelstr = f'{channel}'
                else:
                    ylabelstr = f'{channel} [{yUnits}]'
                if xOffset and line == flightLines[0]:
                    xDel = xData[0]
                xData = xData - xDel

                myPlot, = ax.plot(xData, yData, color='blue', lw=0.3)
                plotTitle = f'{plotPreTitle}, L{line}'
                ax.xaxis.set_major_formatter(thou_format)
                # ax.yaxis.set_major_formatter(thou_format)
                plt.xlabel(xlabelstr, fontsize = 6)
                plt.ylabel(ylabelstr, fontsize = 6)
                plt.title(plotTitle, fontsize = 8)
                plt.grid(True)
                for label in ax.get_xticklabels(): label.set_fontsize(6)
                for label in ax.get_yticklabels(): label.set_fontsize(6)
                plt.show()
