#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check statistics of repeat flight-lines.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter
import matplotlib.ticker as tkr
from scipy import interpolate

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.gridFiles.read_ers as grd
import pegasusQC.utility.utility as util

groupName = config.groupName


def checkRepeatLines(whizzFiles, channel, repeatLines, x='', z='', xOffset=True, verbose=False):
    """
    For all repeatLines, plot (x, channel) and report stats of differences to mean.
    Each line is trimmed to [minX, maxX] and interpolated to common x.
    The analysis is repeated for the `z` channel (height).

    If there are more than 2 lines, then the standard deviation is calculated for
    each x over all the lines, and its mean is reported as the "All lines stdev".

    Parameters
    ----------
    whizzFiles : array of HDF5 Whizz file pathlib Paths

        The pathlib Paths to the Whizz HDF5 files containing the survey repeat line data.

    channel : String

        The name of the channel or field to analyse and plot. Usually a gravity channel.

    repeatLines : [String], optional

        A list of flightlines, e.g. ['1000110.0', '1000210.0', '1000310.0']. 

    x : String, optional

        The name of the independent variable for the analysis. Defaults to whichever of `XChannel`
        or 'YChannel' has the largest range of values.

    y : String, optional

        The name of the independent variable for the analysis. Defaults to the `YChannel`.

    z : String, optional

        The name of the height variable for the analysis and plot. Defaults to the `AltitudeChannel`.

    xOffset : Bool, optional

        If True, map x to x - x[0] before plotting. The default is True.

    Returns
    -------
    None

    """

    # build the arrays to store the data
    if not hasattr(whizzFiles, "__len__"):
        print('ERROR - whizzFiles not an array.')
        return
    if x == '':
        x = _selectXchannel(whizzFiles[0], repeatLines)
    temp_repeats = repeatLines.copy()
    try:
        xBase, xData, gData, zData, minBigX, maxSmallX, deltaX  = _xBaseInterpolant(whizzFiles, channel, temp_repeats, x, z, verbose=verbose)
    except:
        return
    temp_repeats = repeatLines.copy()
    vec_len = len(xBase) - 1 # interpolateLine has lost a datapoint in outputs

    # Interpolate the data to common x and store in arrays
    lineCount = 0
    baseLine = -1.0
    for whizzFile in whizzFiles:

        filename = str(whizzFile)
        with h5py.File(filename, 'r') as f:
            g = f[groupName]['Lines']
            if z == '':
                z = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
            all_flightLines = list(g.keys())

            # if the channel has an attribute 'Units' 
            dd = g[all_flightLines[0]][channel]
            chan_g_label = channel
            chan_g_units = rd.getChannelAttrs(g[all_flightLines[0]], channel)
            chan_g_label += ' ' + chan_g_units

            chan_z_units = rd.getChannelAttrs(g[all_flightLines[0]], z)
            chan_x_units = rd.getChannelAttrs(g[all_flightLines[0]], x)

            # read the data into the arrays
            for line in all_flightLines:
                if line in temp_repeats:
                    baseLine = -1.0
                    if 'PlannedLine' in g[line].attrs.keys():
                        baseLine = g[line].attrs['PlannedLine']
                    xd = rd.getLineData(g[line], x)
                    gd = rd.getLineData(g[line], channel)
                    zd = rd.getLineData(g[line], z)

                    # Get the heading TODO: use this to check RMS(mean difference vs heading direction)
                    try:
                        reportTrackDirection(f[groupName], line)
                    except:
                        print(f'    Cannot report heading for line {line}.')

                    # clean out 'nan's'
                    good = ~np.isnan(xd + gd + zd)
                    xd = xd[good]
                    gd = gd[good]
                    zd = zd[good]
                    if xd.size < 10:
                        print(f'ERROR - after trimming NaNs, the data vectors in {line} are too short for analysis. Stopping.')
                        return

                    # ensure ordered in increasing x
                    if xd[1] < xd[0]:
                        xd = xd[::-1]
                        gd = gd[::-1]
                        zd = zd[::-1]
                    
                    # trim data and store
                    keepsml = xd < xBase[-1]
                    keepbig = xd > xBase[0]
                    keep = keepsml & keepbig

                    # old debug code
                    # print(len(keep), len(xd), len(gd), xBase[0], xd[0], xd[-1])
                    # fig = plt.figure(figsize=(12,9))
                    # ax = fig.add_subplot(1,1,1)
                    # ax.plot(xd[keep]-xBase[0], gd[keep], '.', lw=0.5)
                    # ax.plot(xBase-xBase[0], range(0, len(xBase)), '.', lw=0.5)
                    # plt.show()

                    # interpolate data (last xBase argument is just a dummy)
                    spl = interpolate.splrep(xd[keep]-xBase[0], gd[keep])#, k=3, s=0)
                    gOut = interpolate.splev(xBase-xBase[0], spl)
                    spl = interpolate.splrep(xd[keep]-xBase[0], zd[keep])#, k=3, s=0)
                    zOut = interpolate.splev(xBase-xBase[0], spl)
                    vec_len = xBase.shape[0]

                    xData[lineCount, 0:vec_len] = xBase
                    gData[lineCount, 0:vec_len] = gOut
                    zData[lineCount, 0:vec_len] = zOut

                    lineCount += 1
                    # In case the line is in more than one geoWhizz file
                    temp_repeats.remove(line)
        
    # analyse statistics and report with plots
    _plotRepeatAnalysis(xBase, xOffset, lineCount, xData, gData, zData, channel, repeatLines, z, chan_z_units, chan_g_label, chan_g_units, x, chan_x_units)#, baseLine=baseLine)
            
    return

def reportTrackDirection(surveygroup, line, east='', north=''):
    if east == '':
        east = surveygroup['CoordinateFrame'].attrs['XChannel']
    if north == '':
        north = surveygroup['CoordinateFrame'].attrs['YChannel']

    g = surveygroup['Lines']
    linegroup = g[line]
    dx = np.diff(rd.getLineData(linegroup, east))
    dy = np.diff(rd.getLineData(linegroup, north))
    heading = np.arctan2(dx, dy) * 180.0 / np.pi
    mean_heading = np.mean(heading)
    print(f'    Line {line} heading = {mean_heading:.1f} deg.')


def _xBaseInterpolant(whizzFiles, channel, repeatLines, x='', z='', verbose=False):
    """
    For all `repeatLines` found in the `whizzFiles`, set up the data vectors
    and parameters for interpolations of both `channel` and `z` over `x` to a 
    single pseudo- flight-line.
    
    Parameters
    ----------
    whizzFiles : Array of String or pathlib Path

        Names of HDF5 Whizz files, including path and extension.

    channel : String

        The name of the channel which forms the first dependent variable.
        Its start, end and length after removal of NaNs and dummies is required.

    repeatLines : Array of String

        The flight-lines to be analysed.

    x : String, optional

        The name of the channel containing the `x` data to form the independent
        variable. Defaults to the `XChannel` attribute of each whizzFile.

    z : String, optional

        The name of the channel containing the altitude data which form the
        second dependent variable. Defaults to the `AltitudeChannel` attribute
        of each whizzFile.

    Returns
    -------
    xBase : numpy 1D array of float

        The uniformly sampled `x` locations onto which all `z` and `channel` 
        data will be interpolated

    xData : numpy 2D array of float

        Returned as NaNs but of the correct size to write the data to. 

    gData : numpy 2D array of float

        Returned as NaNs but of the correct size to write the data to. 

    zData : numpy 2D array of float

        Returned as NaNs but of the correct size to write the data to. 

    minBigX : Float

        For all flight-lines, the minimum of the largest value of `x`.

    maxSmallX : Float

        For all flight-lines, the maximum of the smallest value of `x`.

    deltaX : Float

        The sampling spacing of the interpolant, `x`.

    """

    nSamples = 1000000000
    minBigX = 1.0E12
    maxSmallX = -1.0E12
    nLines = len(repeatLines)
    linecount = 0
    
    print('Analysing line data to work out size of interpolated output arrays.')
    for whizzFile in whizzFiles:
        filename = str(whizzFile)
        with h5py.File(filename, 'r') as f:
            g = f[groupName]['Lines']
            if x == '':
                x = f[groupName]['CoordinateFrame'].attrs['XChannel']
            if z == '':
                z = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
            all_flightLines = list(g.keys())

            # nSamples is the array width for data storage
            for line in all_flightLines:
                if line in repeatLines:
                    linecount += 1
                    xs = rd.getLineData(g[line], x)
                    zs = rd.getLineData(g[line], z)
                    cs = rd.getLineData(g[line], channel)

                    xs = xs[~np.isnan(xs + zs + cs)]
                    if xs.size < 10:
                        print(f'ERROR - after trimming NaNs, the data vectors in {line} are too short for analysis. Stopping.')
                        return
                    nSamples = min(nSamples, xs.size)
                    if verbose:
                        print(f'    Line {line}, Shape of x array = {xs.shape}, number of samples for interpolant = {nSamples}.')
                    minBigX = min(max(xs), minBigX)
                    maxSmallX = max(min(xs), maxSmallX)
                    repeatLines.remove(line)
                
    if minBigX < maxSmallX:
        print(f'ERROR: stopping because minBigX ({minBigX}) < maxSmallX ({maxSmallX})')
        return 0.0
    deltaX = (minBigX - maxSmallX) / (nSamples - 1)
    xBase = np.linspace(maxSmallX, minBigX, num=nSamples, endpoint=True)
    print(f'{linecount} of {nLines} lines analysed, interpolant length set to {nSamples} samples.')
    xData = np.empty((linecount, nSamples))
    xData[:] = np.nan
    gData = np.empty((linecount, nSamples))
    gData[:] = np.nan
    zData = np.empty((linecount, nSamples))
    zData[:] = np.nan
    if verbose:
        print(f'  x range = ({maxSmallX:.2f}, {minBigX:.2f}) , deltaX = {deltaX:.2f}.')
        print(f'  output array shapes: x {xBase.shape}, y {gData.shape}, z {zData.shape}.')
    print('')
    
    return xBase, xData, gData, zData, minBigX, maxSmallX, deltaX


def _plotRepeatAnalysis(xBase, xOffset, nLines, xData, gData, zData, channel, flightLines, z, chan_z_units, chan_g_label, chan_g_units, chan_x_label, chan_x_units, baseLine=-1.0):
    xPlot = xBase
    if xOffset:
            xPlot = xPlot - xPlot[0]
    if baseLine < 0.0:
        baseLineStr = ''
    else:
        baseLineStr = f'{baseLine:.0f}'
    fig = plt.figure(figsize=(12,9))
    thou_format = tkr.FuncFormatter(util._space_thou)
    
    #plot the y data
    ax = fig.add_subplot(2,2,1)
    ax.xaxis.set_major_formatter(thou_format)
    for line in range(0, nLines):
        g1 = gData[line,:]
        g1 = g1[np.logical_not(np.isnan(g1))]
        x1 = xData[line,:]
        x1 = x1[np.logical_not(np.isnan(x1))]
        myPlot, = ax.plot(x1, g1, lw=0.5, label=f'Line {flightLines[line]}')
            
    ax.legend(fontsize=8)
    plt.xlabel(f'{chan_x_label} {chan_x_units}', fontsize = 10)
    plt.ylabel(chan_g_label, fontsize = 10)        
    plotTitle = f'{baseLineStr} Repeat Lines {channel}'
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
            
    plt.xlabel(f'{chan_x_label} {chan_x_units}', fontsize = 10)
    plt.ylabel(f'{z} {chan_z_units}', fontsize = 10)
    plotTitle = f'{baseLineStr} Repeat Lines {z}'
    plt.title(plotTitle, fontsize = 12)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    
    # plot the y differences and RMS
    ax = fig.add_subplot(2,2,3)
    ax.xaxis.set_major_formatter(thou_format)
    gMean = np.mean(gData, axis=0)
    gSum = np.zeros(gMean.shape)

    if gData.shape[0] > 2:
        yOverallStd = np.nanmean(np.nanstd(gData, axis=0))
        print('\nSummary Statistics')
        print(f'All lines: stdev({channel}) = {yOverallStd:.2f} {chan_g_units}')
    for line in range(0, nLines):
        gSum = gSum + gData[line,:] - gMean
        gStd = np.nanstd(gData[line,:]-gMean)
        g1 = gData[line,:] - gMean
        g1 = g1[np.logical_not(np.isnan(g1))]
        x1 = xData[line,:]
        x1 = x1[np.logical_not(np.isnan(x1))]

        myPlot, = ax.plot(x1, g1, lw=0.5, label=f'RMS = {gStd:.1f}')
        ax.legend(fontsize=8)
        print(f'    Line {flightLines[line]}: stdev({channel}) = {gStd:.2f} {chan_g_units}')
            
    plt.xlabel(f'{chan_x_label} {chan_x_units}', fontsize = 10)
    plt.ylabel(f'$\Delta$ {chan_g_label}', fontsize = 10)
    plotTitle = f'{channel}: Differences to mean'
    plt.title(plotTitle, fontsize = 12)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    
    # plot the z differences and RMS
    ax = fig.add_subplot(2,2,4)
    ax.xaxis.set_major_formatter(thou_format)
    zMean = np.mean(zData, axis=0)
    zSum = np.zeros(zMean.shape)

    if gData.shape[0] > 2:
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
        print(f'    Line {flightLines[line]}: stdev({z}) = {zStd:.1f} {chan_z_units}')
            
    plt.xlabel(f'{chan_x_label} {chan_x_units}', fontsize = 10)
    plt.ylabel(f'$\Delta$ {z} {chan_z_units}', fontsize = 10)
    plotTitle = f'{z}: Differences to mean'
    plt.title(plotTitle, fontsize = 12)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    
    plt.subplots_adjust(0.1, 0.05, 0.9, 0.95)
    plt.show()
            
    return


def _selectXchannel(whizzFile, repeatLines):
    """
    For all `repeatLines` found in the `whizzFiles`, set up the data vectors
    and parameters for interpolations of both `channel` and `z` over `x` to a 
    single pseudo- flight-line.
    
    Parameters
    ----------
    whizzFile : String or pathlib Path

        Name of HDF5 Whizz file, including path and extension.

    repeatLines : Array of String

        The flight-lines to be analysed.

    Returns
    -------
    x : String

        The name of the channel to be used as the independent variable.

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        y = f[groupName]['CoordinateFrame'].attrs['YChannel']

        all_flightLines = list(g.keys())
        # read the data from the first repeat line into the arrays
        for line in all_flightLines:
            if line in repeatLines:
                xd = rd.getLineData(g[line], x)
                yd = rd.getLineData(g[line], y)
                break
        if np.abs(xd[-1] - xd[0]) > np.abs(yd[-1] - yd[0]):
            return x
        else:
            return y















