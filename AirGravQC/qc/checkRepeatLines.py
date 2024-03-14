import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import AirGravQC.config as config
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.whizzFiles.pointfiles as gw
import AirGravQC.gridFiles.read_ers as grd
import AirGravQC.whizzPlots.whizzPlot as wpl
import AirGravQC.utility.utility as util

groupName = config.groupName


def checkRepeatLines(whizzFiles, channel, repeatLines, x='', z='', xOffset=True):
    """
    For all repeatLines, plot (x, channel) and report stats of differences to mean.
    This will require trimming to [minX, maxX] and interpolating to common x.
    Repeat the analysis for the `z` channel (height).

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
    wpl._plotRepeatAnalysis(xBase, xOffset, lineCount, xData, yData, zData, channel, repeatLines, baseLine, z, chan_z_units, chan_y_label, chan_y_units)
            
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

