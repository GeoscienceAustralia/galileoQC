import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import AirGravQC.config as config
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.whizzFiles.pointfiles as gw
import AirGravQC.gridFiles.read_ers as grd
import AirGravQC.utility.utility as util

groupName = config.groupName


def checkRepeatLines(whizzFiles, channel, repeatLines, x='', z='', xOffset=True):
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
        The name of the independent variable for the plot. Defaults to the `XChannel`.
    z : String, optional
        The name of the height variable for the analysis and plot. Defaults to the `XChannel`.
    xOffset : Bool, optional
        If True, map x to x - x[0] before plotting. The default is True.

    Returns
    -------
    None

    """

    # build the arrays to store the data
    temp_repeats = repeatLines.copy()
    xBase, xData, yData, zData, minBigX, maxSmallX, deltaX  = _xBaseInterpolant(whizzFiles, channel, temp_repeats, x, z)
    temp_repeats = repeatLines.copy()

    # Interpolate the data to common x and store in arrays
    lineCount = 0
    for whizzFile in whizzFiles:

        filename = str(whizzFile)
        with h5py.File(filename, 'r') as f:
            g = f[groupName]['Lines']
            if x == '':
                x = f[groupName]['CoordinateFrame'].attrs['XChannel']
                north = f[groupName]['CoordinateFrame'].attrs['YChannel']
            if z == '':
                z = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
            all_flightLines = list(g.keys())

            # if the channel has an attribute 'Units'
            dd = g[all_flightLines[0]][channel]
            chan_y_label = channel
            if 'Units' in dd.attrs.keys():
                chan_y_units = dd.attrs['Units']
                chan_y_label += ' ' + chan_y_units
            ddz = g[all_flightLines[0]][z]
            if 'Units' in ddz.attrs.keys():
                chan_z_units = ddz.attrs['Units']

            # read the data into the arrays
            for line in all_flightLines:
                if line in temp_repeats:
                    if 'PlannedLine' in g[line].attrs.keys():
                        baseLine = g[line].attrs['PlannedLine']
                    else:
                        baseline = ''
                    xd = rd.getLineData(g[line], x)
                    yd = rd.getLineData(g[line], channel)
                    zd = rd.getLineData(g[line], z)

                    # Get the heading TODO: use this to check RMS(mean difference vs heading direction)
                    dx = np.diff(xd)
                    dy = np.diff(rd.getLineData(g[line], north))
                    heading = np.arctan2(dx, dy) * 180.0 / np.pi
                    mean_heading = np.mean(heading)
                    print(f'Line {line} heading = {mean_heading:.1f} deg.')

                    # ensure ordered in increasing x
                    if xd[1] < xd[0]:
                        xd = xd[::-1]
                        yd = yd[::-1]
                        zd = zd[::-1]

                    xStart = 0
                    xEnd = xd.size - 1
                    
                    # trim data and store
                    for xSample in range(0, xd.size):
                        if xd[xSample] < (maxSmallX - deltaX / 2.0):
                            xStart = max(xSample, xStart)
                        else:
                            break
                    for xSample in range(xd.size-1, 0, -1):
                        if xd[xSample] > (minBigX + deltaX / 2.0):
                            xEnd = min(xSample, xEnd)
                        else:
                            break
                            
                    # interpolate data
                    (yOut, _) = gw.interpolateLine(xd-xBase[0], yd, xBase-xBase[0])
                    (zOut, _) = gw.interpolateLine(xd-xBase[0], zd, xBase-xBase[0])

                    vec_len = len(xBase)-1 # interpolateLine has lost a datapoint in outputs
                    # print(f'line {line}, shapes: xBase {xBase.shape}, xData {xData.shape}')
                    xData[lineCount, 0:vec_len] = xBase[1:]
                    yData[lineCount, 0:vec_len] = yOut
                    zData[lineCount, 0:vec_len] = zOut
                    lineCount += 1
                    # In case the line is in more than one geoWhizz file
                    temp_repeats.remove(line)
        
    # analyse statistics and report with plots
    _plotRepeatAnalysis(xBase, xOffset, lineCount, xData, yData, zData, channel, repeatLines, baseLine, z, chan_z_units, chan_y_label, chan_y_units)
            
    return


def _xBaseInterpolant(whizzFiles, channel, repeatLines, x='', z=''):

    nSamples = 0
    minBigX = 1.0E12
    maxSmallX = -1.0E12
    nLines = len(repeatLines)
    linecount = 0
    
    for whizzFile in whizzFiles:
        filename = str(whizzFile)
        with h5py.File(filename, 'r') as f:
            g = f[groupName]['Lines']
            if x == '':
                x = f[groupName]['CoordinateFrame'].attrs['XChannel']
                north = f[groupName]['CoordinateFrame'].attrs['YChannel']
            if z == '':
                z = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
            all_flightLines = list(g.keys())

            # nSamples is the array width for data storage
            for line in all_flightLines:
                if line in repeatLines:
                    linecount += 1
                    xs = rd.getLineData(g[line], x)
                    nSamples = max(nSamples, xs.size)
                    minBigX = min(max(xs), minBigX)
                    maxSmallX = max(min(xs), maxSmallX)
                    deltaX = np.abs(xs[1] - xs[0])
                    repeatLines.remove(line)
                
    if minBigX < maxSmallX:
        return 0.0
    xBase = np.linspace(maxSmallX, minBigX, num=nSamples, endpoint=True)
    print(f'{linecount} of {nLines} lines analysed, each with {nSamples} samples.')
    xData = np.empty((nLines, nSamples))
    xData[:] = np.nan
    yData = np.empty((nLines, nSamples))
    yData[:] = np.nan
    zData = np.empty((nLines, nSamples))
    zData[:] = np.nan

    return xBase, xData, yData, zData, minBigX, maxSmallX, deltaX


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



