#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot a survey line map.
Author: Mark Helm Dransfield
Created: 2023
License: CC BY-SA
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
from matplotlib.patches import Polygon

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.whizzPlots.whizzPlot as wpl
import pegasusQC.utility.utility as util
import matplotlib.ticker as tkr

groupName = config.groupName


def linesMap(whizzFiles=[], easting='', northing='', flight_lines=[], whizzPlanFile='', planLines=[], planEast='', planNorth='', blockpolygon=[], colourchan='', colourvalue=None):
    """
    Plots a line map of the flown survey over the planned survey lines.

    Parameters
    ----------
    whizzFiles : Array of String or pathlib Path
        Each element is the name of a HDF5 Whizz file, including path and extension.
    easting : String, optional
        The name of the field containing eastings. The default is the name
        stored in the Coordinates attribute XChannel.
    northing : String, optional
        The name of the field containing eastings. The default is the name
        stored in the Coordinates attribute YChannel.
    whizzPlanFile : String or pathlib Path
        The name of a HDF5 Whizz file, including path and extension.
    planLines : Array of String
        Each element is the name of a survey line in whizzPlanFile, flown lines will
        only be plotted if their planned line number is in this list. Default [] ignores
        this limit.
    planEast : String, optional
        The name of the field containing eastings. The default is the name
        stored in the Coordinates attribute XChannel.
    planNorth : String, optional
        The name of the field containing eastings. The default is the name
        stored in the Coordinates attribute YChannel.
    blockpolygon : Array of lists of easting, northing, optional
        If provided, then the polygon is drawn.
    colourchan : String, optional
        The name of a channel in the whizzFiles.
    colourvalue : Float, optional
        If the colourchan and colourvalue are both provided, then the lines in the map
        will be coloured when the value of colourchan = colourvalue.

    Returns
    -------
    None.

    """
    observed_files = True
    if whizzFiles == []:
        print("No files of observed data provided.")
        observed_files = False
    planned_file = True
    if whizzPlanFile == '':
        print("No file of planned data provided.")
        planned_file = False
    if not (observed_files or planned_file):
        print("Need at least one file input.")
        return

    plot_subtitle = ''
    if observed_files and planned_file:
        plot_subtitle = '[planned (red); flown (blue)]'

    fig = plt.figure(figsize=(8, 6))
    thou_format = tkr.FuncFormatter(util._space_thou)
    ax = fig.add_subplot(1,1,1)
    plotTitle = 'Line Map: '

    if planned_file:
        planname = str(whizzPlanFile)

        with h5py.File(planname, 'r') as f:
            if planEast == '':
                planEast = f[groupName]['CoordinateFrame'].attrs['XChannel']
            if planNorth == '':
                planNorth = f[groupName]['CoordinateFrame'].attrs['YChannel']
            g = f[groupName]['Lines']
            plotTitle += wpl.make_plot_title(f[groupName])
                    
            if planLines == []:
                planLines = list(g.keys())
            for line in planLines:
                lX = g[line][planEast][0:]
                lY = g[line][planNorth][0:]
                planline, = ax.plot(lX, lY, color='red', lw=0.6, alpha=0.7)

    if observed_files:
        for file in whizzFiles:
            filename = str(file)

            with h5py.File(filename, 'r') as f:
                if easting == '':
                    easting = f[groupName]['CoordinateFrame'].attrs['XChannel']
                if northing == '':
                    northing = f[groupName]['CoordinateFrame'].attrs['YChannel']
                g = f[groupName]['Lines']
                if flight_lines == []:
                    lines_to_plot = list(g.keys())
                else:
                    lines_to_plot = flight_lines
                if not planned_file:
                    plotTitle += wpl.make_plot_title(f[groupName])
                for line in lines_to_plot:
                    if planned_file:
                        if 'PlannedLine' in g[line].attrs.keys(): #whizzAttrExists(g[line], 'PlannedLine'):
                            planned_line = f"{g[line].attrs['PlannedLine']:.3f}" ## AAARGH HACK
                            # When comparing with a plan, only show the planned lines ...
                            if planned_line in planLines:
                                lX = rd.getLineData(g[line], easting)[0:]
                                lY = rd.getLineData(g[line], northing)[0:]
                                flownline, = ax.plot(lX, lY, color='blue', lw=0.6, alpha=0.7)
                                if not (colourchan == '' or colourvalue is None):
                                    coldata = rd.getLineData(g[line], colourchan)[0:]
                                    coly = np.ma.masked_where(coldata != colourvalue, lY)
                                    colline, = ax.plot(lX, coly, color='orange', lw=0.8, alpha=0.9)
                    else:
                        # ... otherwise, show all observed lines.
                        lX = rd.getLineData(g[line], easting)[0:]
                        lY = rd.getLineData(g[line], northing)[0:]
                        flownline, = ax.plot(lX, lY, color='blue', lw=0.6, alpha=0.7)
                        if not (colourchan == '' or colourvalue is None):
                            coldata = rd.getLineData(g[line], colourchan)[0:]
                            coly = np.ma.masked_where(coldata != colourvalue, lY)
                            colline, = ax.plot(lX, coly, color='orange', lw=0.8, alpha=0.9)
            
    ax.set_aspect('equal')
    ax.xaxis.set_major_formatter(thou_format)
    ax.yaxis.set_major_formatter(thou_format)
    plt.xlabel(f'{easting} [m]', fontsize = 10)
    plt.ylabel(f'{northing} [m]', fontsize = 10)
    plt.suptitle(plotTitle, fontsize = 12)
    plt.title(plot_subtitle, fontsize = 10)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    if len(blockpolygon) > 0:
        plotsurveyboundary(ax, blockpolygon)
    plt.show()


def plotsurveyboundary(ax, blockpolygon):
    polygon = Polygon(blockpolygon, edgecolor='green', fill=False)
    ax.add_patch(polygon)  
    return