#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A collection of possibly useful plot functions, not in common use.

Created on Sat Aug 14 18:21:15 2021

@author: markdransfield
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import xarray as xr
import verde as vd
import pooch

import AirGravQC.config as config
import AirGravQC.qualitycontrol.qualityAnalysis as qc
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.utility.utility as util
import matplotlib.ticker as tkr

groupName = config.groupName


def plotWsLineChannel(whizzFile1, flightLine1, channel1,  
                      whizzFile2, flightLine2, channel2, 
                      x1='', x2='', y1='', y2='', h1='', h2='',
                      plotTitle = '', xOffset=False):
    """
    This is a one-off for the Vic/SA project(5371). Given flightLines from the Otway
    project and this project that are close to each spatially, compare them in plots.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    flightLine : String
        A flightline, e.g. '1000110.0'.
    channel : String
        The name of the channel or field to plot.
    x : String, optional
        The name of the independent variable for the plot. The default is ''.
    plotTitle : String, optional
        A title for the plot. The default is '' in which case the title will be Project Line Channel.
    xOffset : Bool, optional
        If True, map x to x - x[0] before plotting. The default is True.

    Returns
    -------
    None.

    """
    
    x1Data, x1Units, y1Data, y1Units, h1Data, h1Units, z1Data, z1Units, proj1Name = _get_data(
        whizzFile1, flightLine1, channel1, x=x1, y=y1, h=h1)
    x2Data, x2Units, y2Data, y2Units, h2Data, h2Units, z2Data, z2Units, proj2Name = _get_data(
        whizzFile2, flightLine2, channel2, x=x2, y=y2, h=h2)
    
    clipx0 = 142.7
    clipx1 = 143.35
    z2Scale = 1.0
    if z2Units == 'mGal':
        z2Scale = 10.0
        z2Units = 'gu'
    z1Scale = 1.0
    if z1Units == 'mGal':
        z1Scale = 10.0
        z1Units = 'gu'
    z1D = z1Scale * z1Data[(clipx0 < x1Data) & (x1Data  < clipx1)]
    y1D = y1Data[(clipx0 < x1Data) & (x1Data  < clipx1)]
    h1D = h1Data[(clipx0 < x1Data) & (x1Data  < clipx1)]
    x1D = x1Data[(clipx0 < x1Data) & (x1Data  < clipx1)]
    z2D = z2Scale * z2Data[(clipx0 < x2Data) & (x2Data  < clipx1)]
    y2D = y2Data[(clipx0 < x2Data) & (x2Data  < clipx1)]
    h2D = h2Data[(clipx0 < x2Data) & (x2Data  < clipx1)]
    x2D = x2Data[(clipx0 < x2Data) & (x2Data  < clipx1)]
    
    xLabel = f'[{x1Units}] / [{x1Units}]'
    yLabel = f'{y1} [{y1Units}] / {y2} [{y2Units}]'
    hLabel = f'{h1} [{h1Units}] / {h2} [{h2Units}]'
    zLabel = f'Bouguer Gravity [gu]'
    
    fig = plt.figure()
    fig.suptitle(f'Lines {flightLine1} & {flightLine2}', fontsize=10)
    # fig.text(0.8, 0.95, f'acq: {acqDate:.0f}', fontsize=7)
    fig.subplots_adjust(top=0.90)
    
    if plotTitle == '':
        plotTitle = f'{channel1} / {channel2}'
    if xOffset:
        xLabel += f' - {x1Data[0]}'
        x1Data = x1Data - x1Data[0]
        x2Data = x2Data - x1Data[0]

    ax1 = fig.add_subplot(3,1,1)
    _subplotCompare(ax1, x1D, z1D, x2D, z2D, channel1, channel2, '', 'Gravity [gu]', 'Gravity')
    ax2 = fig.add_subplot(3,1,2)
    _subplotCompare(ax2, x1D, h1D, x2D, h2D, flightLine1, flightLine2, '', 'Height [m]', 'Height')
    ax3 = fig.add_subplot(3,1,3)
    _subplotCompare(ax3, x1D, y1D, x2D, y2D, flightLine1, flightLine2, 'Longitude [degree]', 'Latitude [degree]', 'Plan View')

    plt.show()
    
    
def _subplotCompare(ax, x1, y1, x2, y2, c1, c2, xLabel, yLabel, plotTitle):
    line1, = ax.plot(x1, y1, color='blue', lw=0.6, label=c1)
    line2, = ax.plot(x2, y2, color='green', lw=0.6, label=c2)
    ax.set_xlabel(xLabel, fontsize = 8)
    ax.set_ylabel(yLabel, fontsize = 8)
    plt.title(plotTitle, fontsize = 9)
    ax.xaxis.grid(visible=True, color='black', lw=0.3, ls='--')
    ax.yaxis.grid(visible=True, color='black', lw=0.3, ls='--')
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    ax.legend(fontsize=7)


def _get_data(whizzFile, flightLine, channel, x='', y='', h=''):
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        
        xName = x
        if x == '':
            xName = f[groupName]['CoordinateFrame'].attrs['XChannel']
        xUnits = getChannelAttrs(g[flightLine], xName)#g[flightLine][xName].attrs['Units']
        xData = np.array(g[flightLine][xName])

        yName = y
        if y == '':
            yName = f[groupName]['CoordinateFrame'].attrs['YChannel']
        yUnits = g[flightLine][yName].attrs['Units']
        yData = np.array(g[flightLine][yName])

        hName = h
        if h == '':
            hName = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        hUnits = g[flightLine][hName].attrs['Units']
        hData = np.array(g[flightLine][hName])

        zData = np.array(g[flightLine][channel])
        zUnits = g[flightLine][channel].attrs['Units']

    return xData, xUnits, yData, yUnits, hData, hUnits, zData, zUnits, projName


def plot_xcohere(whizzFile, flightLine, xchannel, ychannel):
    """
    Plot coherence between `xchannel` and `ychannel`.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    flightLine : String
        A flightline, e.g. '1000110.0'.
    xchannel : String
        The name of the first channel or field to analyse.
    ychannel : String
        The name of the second channel or field to analyse.

    Returns
    -------
    None.

    """
    import scipy.signal as sig

    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        xdata = np.array(g[flightLine][xchannel])
        ydata = np.array(g[flightLine][ychannel])


    (f, c) = sig.coherence(xdata, ydata, fs = 2.0, nfft=1024) 
    period = 1.0 / f[1:]

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.plot(period, c[1:]) 
    plt.xlabel('Period [s]', fontsize = 6)
    plt.ylabel('Coherence', fontsize = 6)
    plotTitle = f'{projName} : {xchannel} v {ychannel} L{flightLine}'
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    plt.show() 
    
    
def make_plot_title(group):
    """
    Make a title string to be used in plots.

    Parameters
    ----------
    group : HDF5 Group
        Name of a HDF5 group that contains `ProjectName` and `BlockID` attributes.

    Returns
    -------
    plotTitle : String
        The plot title.

    """
    plotTitle = ''
    if 'ProjectName' in group.attrs:
        plotTitle += group.attrs['ProjectName']
    if 'BlockID' in group.attrs:
        plotTitle += ' (' + group.attrs['BlockID'] + ')'
    return plotTitle


def specificLinesMap(whizzFile, lines, easting='', northing=''):
    """
    Plots a line map of the survey contained in the HDF5 Whizz file.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        The name of a HDF5 Whizz file, including path and extension.
    lines : Array of String
        The lines to plot to the map
    easting : String, optional
        The name of the field containing eastings. The default is the name
        stored in the Coordinates attribute XChannel.
    northing : String, optional
        The name of the field containing eastings. The default is the name
        stored in the Coordinates attribute YChannel.

    Returns
    -------
    None.

    """
    from matplotlib.ticker import StrMethodFormatter

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        if easting == '':
            easting = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if northing == '':
            northing = f[groupName]['CoordinateFrame'].attrs['YChannel']
        g = f[groupName]['Lines']
        plotTitle = f[groupName].attrs['ProjectName'] + ': Line Map'
        
        mycolours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',
        'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
        myc = 0
        for line in lines:
            lX = rd.getLineData(g[line], easting)[0:]
            lY = rd.getLineData(g[line], northing)[0:]
            flownline, = ax.plot(lX, lY, color=mycolours[myc], label=line, lw=3, alpha=0.7)
            myc += 1
            
    ax.set_aspect('equal')
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    plt.xlabel('X [m]', fontsize = 10)
    plt.ylabel('Y [m]', fontsize = 10)
    plt.suptitle(plotTitle, fontsize = 12)
    # plt.title('[planned (red); flown (blue)]', fontsize = 10)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(10)
    for label in ax.get_yticklabels(): label.set_fontsize(10)
    plt.legend(fontsize=8)
    plt.show()


def statusMap(planFile='', planEast='', planNorth='', plotTitle=''):
    """
    Plots a line map of the accepted segments of the survey as updated in
    the HDF5 Whizz plan file.

    Parameters
    ----------
    planFile : Path
        The name of a HDF5 Whizz file, including path and extension, containg the plan.
    easting : String, optional
        The name of the field containing eastings. The default is the name
        stored in the Coordinates attribute XChannel.
    northing : String, optional
        The name of the field containing eastings. The default is the name
        stored in the Coordinates attribute YChannel.

    Returns
    -------
    None.

    """
    from matplotlib.ticker import StrMethodFormatter

    # set up the figure for plotting
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    if plotTitle == '':
        plotTitle = 'Line QC Status Map'
    
    planname = str(planFile)

    with h5py.File(planname, 'r') as f:
        if planEast == '':
            planEast = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if planNorth == '':
            planNorth = f[groupName]['CoordinateFrame'].attrs['YChannel']
        statusChannel = 'AcceptedReflight'
        g = f[groupName]['Lines']
        
        # for each line that has been flown, plot the parts where the status >= 0
        for line in list(g.keys()):
            line_path = g[line]
            if line_path.attrs['HasBeenFlown']:
                status = rd.getLineData(g[line], statusChannel).astype(np.float32)
                lX = rd.getLineData(g[line], planEast)[status >= 0]
                lY = rd.getLineData(g[line], planNorth)[status >= 0]
                lS = status[status >= 0]
                ax.plot(lX, lY, lw=0.2, color='blue')

    ax.set_aspect('equal')
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    plt.xlabel('X [m]', fontsize = 10)
    plt.ylabel('Y [m]', fontsize = 10)
    plt.title(plotTitle, fontsize = 12)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(10)
    for label in ax.get_yticklabels(): label.set_fontsize(10)
    plt.show()


def _plot_speed(t, t_label, speed, min_speed=54, max_speed=66, plot_title=''):
    """
    Plots the speed as a function of t and adds limit lines showing the allowed
    tolerance

    Parameters
    ----------
    t : TYPE
        Time (sec) data vector.
    speed : TYPE
        The speed data vector.
    min_speed : TYPE, optional
        DESCRIPTION. The default is 54.
    max_speed : TYPE, optional
        DESCRIPTION. The default is 66.
    plot_title : TYPE, optional
        DESCRIPTION. The default is ''.

    Returns
    -------
    None.

    """
    fig = plt.figure()
    fig.suptitle(plot_title, fontsize=10)
    fig.subplots_adjust(top=0.85)
    
    ax = fig.add_subplot(1,1,1)
    thou_format = tkr.FuncFormatter(util._space_thou)
    ax.plot(t, speed, 'b', t, np.ones(t.size) * min_speed, 'r', t, np.ones(t.size) * max_speed, 'r', mfc='w')
    ax.xaxis.set_major_formatter(thou_format)
    plt.xlabel(t_label, fontsize = 6)
    plt.ylabel('Speed [m/s]', fontsize = 6)
    plotTitle = 'Speed' + ' Stats'
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)

   
def _plotcheckSafeClearance(projName, line, distance, clearance_chan='', altitude_chan='', terrain_chan='', alt=[], dtm=[], clearance=[]):
    fig = plt.figure()

    ax = fig.add_subplot(2,1,1)
    thou_format = tkr.FuncFormatter(util._space_thou)
    if clearance_chan == '':
        ax.plot(distance, alt)
        ax.plot(distance, dtm)
        plt.legend([altitude_chan, terrain_chan], fontsize=8)
    else:
        ax.plot(distance, clearance)
        plt.legend([clearance_chan], fontsize=8)
    ax.xaxis.set_major_formatter(thou_format)
    plt.xlabel('distance along line [m]', fontsize = 6)
    plt.ylabel('height', fontsize = 6)
    plotTitle = projName + ': Minimum Clearance Check ' + line
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)


