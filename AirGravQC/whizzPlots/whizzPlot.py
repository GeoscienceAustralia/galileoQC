#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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
import AirGravQC.qc.qualityAnalysis as qc
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.utility.utility as util
import matplotlib.ticker as tkr

groupName = config.groupName


def plot_grid(whizzFile, channel, cellsize, x='', y=''):
    """
    WARNING - FAILS ON LARGE DATASETS
    Grids a channel in the whizzFile, and displays it as an image.
    The grid interpolation is performed by Verde's spline function.
    (Another option would be to use pygmt's xyz2grd method.)
    The image is displayed as a simple pcolormesh.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    channel : String
        The name of the first channel or field to analyse.
    cellsize : Float
        The width of the grid cells in the same units as `x` and `y`.
    x : String, optional
        The name of the easting channel. Default is to use 'XChannel' 
        from the whizzFile's attributes.
    y : String, optional
        The name of the northing channel. Default is to use 'XChannel' 
        from the whizzFile's attributes.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    xpos = np.zeros((0,))
    ypos = np.zeros((0,))
    z = np.zeros((0,))
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']
        for line in list(g.keys()):
            xpos = np.append(xpos, rd.getLineData(g[line], x))
            ypos = np.append(ypos, rd.getLineData(g[line], y))
            z = np.append(z, rd.getLineData(g[line], channel))
            
        x_min = np.nanmin(xpos)
        x_max = np.nanmax(xpos)

        y_min = np.nanmin(ypos)
        y_max = np.nanmax(ypos)

        spline = vd.Spline()
        spline.fit((xpos, ypos), z)
        grid_coords = vd.grid_coordinates(region=(x_min, x_max, y_min, y_max), spacing=cellsize)
        z_reg = spline.predict(grid_coords)
        
        ax.pcolormesh(grid_coords[0], grid_coords[1], z_reg)
        ax.set_title(f'{projName} gridded {channel}')
        plt.colorbar()
        plt.tight_layout()
        plt.show()
        

def plotxy(y, x='', plotTitle = '', xOffset=True, plot_symbol=''):
    """
    

    Parameters
    ----------
    y : DICT - y{'label': axis_label, 'data': numpy 1D array}
        The y axis label and data.
    x : DICT - x{'label': axis_label, 'data': numpy 1D array}, optional
        The x axis label and data. The default is ''.
    plotTitle : String, optional
        Usually includes project name. The default is ''.
    xOffset : Bool, optional
        If True, map x to x - x[0] before plotting. The default is True.
    plot_symbol : String, optional
        The symbol to use in plotting. The default is '' giving a continuous line.

    Returns
    -------
    None.

    """
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    
    yData = y['data']
    yLabel = y['label']
    if x != '':
        xData = x['data']
        xLabel = x['label']
    else:
        xData = range(0, len(yData))
        xLabel = 'sample'
    
    fig.suptitle(plotTitle, fontsize=10)
    fig.subplots_adjust(top=0.85)
    
    if xOffset and x != '':
        xLabel += f' - {xData[0]}'
        xData = xData - xData[0]

    line, = ax.plot(xData, yData, color='blue', lw=0.3, )
    if plot_symbol != '':
        line.set_linestyle('None')
        line.set_marker(plot_symbol)
    plt.xlabel(xLabel, fontsize = 7)
    plt.ylabel(yLabel, fontsize = 7)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    plt.show()


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


def plotInlineSum(whizzFile, line):
    """
    Plots the raw, and filtered, inline sum for a line in the whizzFile.
    The routine relies on too many hard-coded values, and should be
    improved or removed.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    line : String
        The flight line identifier.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        acqDate = g[line].attrs['Date_Local']

        time = g[line]['time']
        data1 = rd.getLineData(g[line], 'Inline1_raw')
        data2 = rd.getLineData(g[line], 'Inline2_raw')
        data3 = rd.getLineData(g[line], 'Inline3_raw')
        ils = gw.inLineSum(data1, data2, data3)
        
        fig = plt.figure()
        fig.suptitle(f'{projName}: Line {line}', fontsize=10)
        fig.text(0.8, 0.95, f'acq: {acqDate:.0f}', fontsize=7)
        fig.subplots_adjust(top=0.85)
        
        ax = fig.add_subplot(2,1,1)
        ax.plot(time, ils, ms=2)
        ax.set_ylim(-100.0, 100.0)
        plt.ylabel('Inline Sum [E]', fontsize = 6)
        plotTitle = 'Inline Sum'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        
        ax2 = fig.add_subplot(2,1,2)
        ax2.plot(time, data1+data2+data3, ms=2)
        #ax2.set_ylim(-100.0, 100.0)
        plt.ylabel('Raw Inline Sum [E]', fontsize = 6)
        plotTitle = 'Inline Sum'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        plt.show()
        
    
def plotAllRepeatLines(filename, flightLines, x='', channels=[], xOffset=True):
    """
    For all lines in flightLines (assumed to be repeats), plot (x, channel) and
    report stats of differences to mean.
    This will require trimming to [minX, maxX] and interpolating to common x.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    flightLines : [String]
        An array of flightline, e.g. ['1000110.0'].
    x : String, optional
        The name of the independent variable for the plot. Default is the `whizzFile`
        `XChannel`.
    channels : [String]
        The names of the channels to analyse and plot.
    xOffset : Bool, optional
        If True, map x to x - x[0] before plotting. The default is True.

    Returns
    -------
    None.

    """

    # nSamples = 0
    # nLines = len(flightLines)
    # minBigX = 1.0E12
    # maxSmallX = -1.0E12
    baseLine = float(flightLines[0])
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if channels == []:
            lineGroups = list(g.values())
            channels = list(lineGroups[0].keys())

        # nSamples is the array width for data storage
        # for line in flightLines:
        #     xs = np.array(g[line][x])
        #     nSamples = max(nSamples, xs.size)
        #     minBigX = min(max(xs), minBigX)
        #     maxSmallX = max(min(xs), maxSmallX)
            
        # xData = np.empty((nLines, nSamples))
        # xData[:] = np.nan
        # yData = np.empty((nLines, nSamples))
        # yData[:] = np.nan
        
        # read the data into the arrays
        for channel in channels:
            fig = plt.figure(figsize=(6,9))
            ax = fig.add_subplot(1,1,1)
            for line in flightLines:
                xd = rd.getLineData(g[line], x)
                yd = rd.getLineData(g[line], channel)
                #x1 = x1[np.logical_not(np.isnan(x1))]
                #y1 = y1[np.logical_not(np.isnan(y1))]
                myPlot, = ax.plot(xd, yd, lw=0.5)
                ax.text(xd[0], yd[0], f'Line {line}', fontsize=6)
            plt.xlabel('x', fontsize = 8)
            plt.ylabel(channel, fontsize = 8)
            plotTitle = f'{baseLine:.0f} Repeat Lines {channel}'
            plt.title(plotTitle, fontsize = 10)
            plt.grid(True)
            for label in ax.get_xticklabels(): label.set_fontsize(6)
            for label in ax.get_yticklabels(): label.set_fontsize(6)
            
            plt.show()
            
    return


def plotLineXd4Channels(whizzFile, flightLine, x, channel1, channel2, xOffset=True):
    """
    For the given flightLine in the whizzFile, plot the 4th difference of both
    channel1 and channel2 against x.
    

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    flightLine : String
        A flightline, e.g. '1000110.0'.
    x : String
        The name of the independent variable for the plot.
    channel : String
        The name of the first channel or field to analyse.
    channel : String
        The name of the second channel or field to analyse.
    xOffset : Bool, optional
        If True, map x to x - x[0] before plotting. The default is True.

    Returns
    -------
    None.

    """
    
    filename = str(whizzFile)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        
        xData = np.array(g[flightLine][x])
        y1Data = np.diff(np.array(g[flightLine][channel1]), n = 4)
        y2Data = np.diff(np.array(g[flightLine][channel2]), n = 4)
    
    if xOffset:
        xData = xData - xData[0]
    ax.plot(xData[2:-2], y1Data, xData[2:-2], y2Data, lw=0.5)
    plt.xlabel(x, fontsize = 6)
    plt.ylabel(channel1, fontsize = 6)
    plotTitle = f'{projName}: L{flightLine}, {x} vs D4({channel1}, {channel2})'
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    plt.show()


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
    
    
def psdLineChannel(whizzFile, flightLine, channel, time='', plotTitle = ''):
    """
    Plot the PSD (log-log Sqrt(Power) from welch method) of channel in flightLine. 

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    flightLine : String
        A flightline, e.g. '1000110.0'.
    channel : String
        The name of the channel or field to plot.
    plotTitle : String, optional
        A title for the plot. The default is '' in which case the title will be Project Line Channel.

    Returns
    -------
    None.

    """
    import scipy.signal as sig
    
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if time == '':
            time = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
        data = np.array(g[flightLine][channel])
        t = np.array(g[flightLine][time])
        f_sample = 1.0 / (t[1] - t[0])

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    freq, Pxx = sig.welch(data, nfft=2048, fs = f_sample)
    period = 1.0 / freq[1:]
    plt.plot(period, np.sqrt(Pxx[1:]), color='blue', lw=0.3)
    # plt.semilogx(freq, np.sqrt(Pxx), color='blue', lw=0.3)

    plt.xlim([0, 200])
    plt.xlabel('Period [s]', fontsize = 6)
    # plt.xlabel('Frequency [Hz]', fontsize = 6)
    plt.ylabel(channel, fontsize = 6)
    if plotTitle == '':
        plotTitle = f'{projName} : {channel} L{flightLine}'
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    plt.show()


def psdLineChannels(whizzFile, flightLine, channel1, channel2, time='', plotTitle = ''):
    """
    Plot the PSD (log-log Sqrt(Power) from welch method) of channel in flightLine. 

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    flightLine : String
        A flightline, e.g. '1000110.0'.
    channel : String
        The name of the channel or field to plot.
    plotTitle : String, optional
        A title for the plot. The default is '' in which case the title will be Project Line Channel.

    Returns
    -------
    None.

    """
    import scipy.signal as sig
    
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if time == '':
            time = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
        data1 = np.array(g[flightLine][channel1])
        data2 = np.array(g[flightLine][channel2])
        t = np.array(g[flightLine][time])
        f_sample = 1.0 / (t[1] - t[0])

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    freq, Pxx1 = sig.welch(data1, nfft=2048, fs = f_sample)
    freq, Pxx2 = sig.welch(data2, nfft=2048, fs = f_sample)
    period = 1.0 / freq[1:]
    plt.plot(period, np.sqrt(Pxx2[1:]) - np.sqrt(Pxx1[1:]), 'g', lw=0.6)
    plt.xlim([0, 200])
    
    plt.xlabel('Period [s]', fontsize = 6)
    plt.ylabel(f'{channel2} / {channel1}', fontsize = 6)
    if plotTitle == '':
        plotTitle = f'{projName} L{flightLine}: sqrt(Pwr({channel2})) - sqrt(Pwr({channel2}))'
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    plt.show()


def make_plot_title(group):
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
    
    
def plotBoxWhisker(chMin, chMax, chMean, chStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr='Line Number'):
    """

    """
    fig = plt.figure()
    fig.suptitle(figtitle, fontsize=10)
    fig.subplots_adjust(top=0.85)
    thou_format = tkr.FuncFormatter(util._space_thou)

    ax = fig.add_subplot(1,1,1)
    ax.plot(lineNo, chMin, 'bo', mfc='w')
    ax.plot(lineNo, chMax, 'bo', mfc='w')
    ax.errorbar(lineNo, chMean, chStd, capsize=3, marker='s', c='blue', linestyle='')
    ax.xaxis.set_major_formatter(thou_format)
    plt.title(titlestr)
    plt.xlabel(xlabelstr, fontsize = 6)
    plt.ylabel(ylabelstr, fontsize = 6)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    plt.show()
    return        
    
    
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


def _plotRepeatAnalysis(xBase, xOffset, nLines, xData, yData, zData, channel, flightLines, baseLine, z, chan_z_units, chan_y_label, chan_y_units):
    xPlot = xBase
    if xOffset:
            xPlot = xPlot - xPlot[0]
        
    fig = plt.figure(figsize=(12,9))
    thou_format = tkr.FuncFormatter(util._space_thou)
    
    #plot the y data
    ax = fig.add_subplot(2,2,1)
    ax.xaxis.set_major_formatter(thou_format)
    for line in range(0, nLines):
        y1 = yData[line,:]
        y1 = y1[np.logical_not(np.isnan(y1))]
        x1 = xData[line,:]
        x1 = x1[np.logical_not(np.isnan(x1))]
        myPlot, = ax.plot(x1, y1, lw=0.5, label=f'Line {flightLines[line]}')
            
    ax.legend(fontsize=8)
    plt.xlabel('x', fontsize = 10)
    plt.ylabel(chan_y_label, fontsize = 10)
    plotTitle = f'{baseLine:.0f} Repeat Lines {channel}'
    plt.title(plotTitle, fontsize = 12)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    
    #plot the z data
    ax = fig.add_subplot(2,2,2)
    ax.xaxis.set_major_formatter(thou_format)
    for line in range(0, nLines):
        z1 = zData[line,:]
        z1 = z1[np.logical_not(np.isnan(z1))]
        x1 = xData[line,:]
        x1 = x1[np.logical_not(np.isnan(x1))]
        myPlot, = ax.plot(x1, z1, lw=0.5, label=f'Line {flightLines[line]}')
        ax.legend(fontsize=8)
            
    plt.xlabel('x', fontsize = 10)
    plt.ylabel(f'z {chan_z_units}', fontsize = 10)
    plotTitle = f'{baseLine:.0f} Repeat Lines {z}'
    plt.title(plotTitle, fontsize = 12)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    
    # plot the y differences and RMS
    ax = fig.add_subplot(2,2,3)
    ax.xaxis.set_major_formatter(thou_format)
    yMean = np.mean(yData, axis=0)
    ySum = np.zeros(yMean.shape)

    if yData.shape[0] > 2:
        yOverallStd = np.nanmean(np.nanstd(yData, axis=0))
        print(f'All lines: stdev({channel}) = {yOverallStd:.2f} {chan_y_units}')
    for line in range(0, nLines):
        ySum = ySum + yData[line,:] - yMean
        yStd = np.nanstd(yData[line,:]-yMean)
        y1 = yData[line,:] - yMean
        y1 = y1[np.logical_not(np.isnan(y1))]
        x1 = xData[line,:]
        x1 = x1[np.logical_not(np.isnan(x1))]

        myPlot, = ax.plot(x1, y1, lw=0.5, label=f'RMS = {yStd:.1f}')
        ax.legend(fontsize=8)
        print(f'Line {flightLines[line]}: stdev({channel}) = {yStd:.2f} {chan_y_units}')
            
    plt.xlabel('x', fontsize = 10)
    plt.ylabel(chan_y_label, fontsize = 10)
    plotTitle = f'Differences to mean and RMS {channel}'
    plt.title(plotTitle, fontsize = 12)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    
    # plot the z differences and RMS
    ax = fig.add_subplot(2,2,4)
    ax.xaxis.set_major_formatter(thou_format)
    zMean = np.mean(zData, axis=0)
    zSum = np.zeros(zMean.shape)

    if yData.shape[0] > 2:
        zOverallStd = np.nanmean(np.nanstd(zData, axis=0))
        print(f'All lines: stdev({z}) = {zOverallStd:.1f} {chan_z_units}')
    for line in range(0, nLines):
        zSum = zSum + zData[line,:] - zMean
        zStd = np.nanstd(zData[line,:]-zMean)
        z1 = zData[line,:]-zMean
        z1 = z1[np.logical_not(np.isnan(z1))]
        x1 = xData[line,:]
        x1 = x1[np.logical_not(np.isnan(x1))]

        myPlot, = ax.plot(x1, z1, lw=0.5, label=f'RMS = {zStd:.2f}')
        ax.legend(fontsize=8)
        # ax.text(x1[0], z1[0], f'RMS = {zStd:.2f}', fontsize=8)
        print(f'Line {flightLines[line]}: stdev({z}) = {zStd:.1f} {chan_z_units}')
            
    plt.xlabel('x', fontsize = 10)
    # plt.ylabel(f'z {chan_z_units}', fontsize = 10)
    plotTitle = f'Differences to mean and RMS {z}'
    plt.title(plotTitle, fontsize = 12)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    
    plt.subplots_adjust(0.1, 0.05, 0.9, 0.95)
    plt.show()
            
    return


