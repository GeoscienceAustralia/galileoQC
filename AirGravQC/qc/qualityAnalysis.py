#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Aug 14 20:28:20 2021

@author: markdransfield


"""


# import necessary modules
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter
from matplotlib.ticker import StrMethodFormatter
import matplotlib.ticker as tkr

import AirGravQC.config as config
import AirGravQC.whizzFiles.pointfiles as gw
import AirGravQC.gridFiles.read_ers as ers
import AirGravQC.gridFiles.gridfiles as grd
import AirGravQC.utility.utility as util
import AirGravQC.whizzPlots.whizzPlot as wpl

groupName = config.groupName

def checkPhase(filename, channel1, channel2):
    """
    For every survey line in a geoWhizz HDF5 file, given two channels, calculate
    the phase shift (in number of samples) required to maximise the correlation
    between the two channels

    Parameters
    ----------
    filename : String
        The name of a geoWhizz HDF5 file.
    channel1 : String
        The name of the first channel.
    channel2 : String
        The name of the second channel.
    
    Returns
    -------
    None.

    """
    from scipy.signal import correlate
    with h5py.File(filename, 'r') as f:
        g = f[groupName]
        numLines = len(g.items())
        offsets = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            linegroup = g[line]
            A = gw.getLineData(linegroup, channel1)
            B = gw.getLineData(linegroup, channel2)
            nsamples = A.size

            # regularize datasets by subtracting mean and dividing by s.d.
            A -= A.mean(); A /= A.std()
            B -= B.mean(); B /= B.std()

            # Put in an artificial time shift between the two datasets
            #time_shift = 20
            #A = numpy.roll(A, time_shift)

            # Find cross-correlation
            xcorr = correlate(A, B)
            
            # delta time array to match xcorr
            dt = np.arange(1-nsamples, nsamples)
            recovered_time_shift = dt[xcorr.argmax()]
            offsets[count] = recovered_time_shift
            count += 1

            time = np.arange(dt[0], dt[-1], 0.1)
            # Now interpolate through gaps by cubic spline
            (xcorrInt, _) = gw.interpolateLine(dt, xcorr, time)
            recovered_time_shift2 = time[xcorrInt.argmax()]
            print(f'Line {line}: Recovered time shift = {recovered_time_shift2:.1f}')
            
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(time, xcorrInt)
        plt.ylabel('Correlation [arbitrary units]', fontsize = 6)
        plt.xlabel('FID difference [s]', fontsize = 6)
        plotTitle = f'Line {line}: Correlation of {channel1} v {channel2}'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        plt.show()
        print(f'Offset MEAN = {np.mean(offsets):.2f}; STD = {np.std(offsets):.3f}')
        
        
def calcDrift(whizzFile, time, gradient):
    """
    NOT USED
    """
    filename = str(whizzFile)
    from sklearn.linear_model import LinearRegression
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        
        for line in g.keys():
            t = gw.getLineData(g[line], time).reshape((-1,1)) 
            gamma = gw.getLineData(g[line], gradient)#g[line][gradient]
            
            model = LinearRegression().fit(t, gamma)
            r_sq = model.score(t, gamma)
            print('coefficient of determination:', r_sq, 'intercept:', model.intercept_, 'slope:', model.coef_)


def checkVertAcc(whizzFile, vertvelocity):
    """
    Uses numpy.diff to estimate the vertical acceleration from the vertical velocity
    data (units of m/s) in the channel called `vertvelocity` in whizzFile. Plots the
    result (units of milli-g) for each survey line.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.

    Returns
    -------
    None

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']

        numLines = len(g.items())
        chStd = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            data = gw.getLineData(g[line], vertvelocity)#g[line]['vertvelocity']
            accel = np.diff(data, n = 1)
            chStd[count] = int(np.std(accel) * 100.0)
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            ax.plot(accel, '.', ms=2)
            plt.ylabel('Acceleration', fontsize = 6)
            plotTitle = groupName + ': Vert Accel [m/s/s] ' + str(chStd[count]) + 'mg ' + line
            plt.title(plotTitle, fontsize = 8)
            plt.grid(True)
            for label in ax.get_xticklabels(): label.set_fontsize(6)
            for label in ax.get_yticklabels(): label.set_fontsize(6)
            plt.show()
            count += 1


def checkVertAccStats(whizzFile):
    """
    Uses numpy.diff to estimate the vertical acceleration from the vertical velocity
    data (units of m/s) in the channel called `vertvelocity` in whizzFile. Plots the
    min, max, mean and stdev (units of milli-g) for each survey line as a function of
    line.
    Assumes samples are 1 sec apart.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.

    Returns
    -------
    None

    """
    outlierRatio = 6.0
    specAcc = 1.0
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]

        numLines = len(g.items())
        chMin = np.zeros((numLines,))
        chMax = np.zeros((numLines,))
        chMean = np.zeros((numLines,))
        chStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            data = g[line]['vertvelocity']
            accel = np.diff(data, n = 1)
            lineNo[count] = line
            chMin[count] = np.min(accel)
            chMax[count] = np.max(accel)
            chMean[count] = np.mean(accel)
            chStd[count] = np.std(accel)
            maxDevRatio = (chMax[count] - chMean[count]) / chStd[count]
            minDevRatio = (chMean[count] - chMin[count]) / chStd[count]
            if chStd[count] > specAcc:
                print(lineNo[count], ': 100 mill-g exceedance ', chStd[count])
            if maxDevRatio > outlierRatio:
                print(lineNo[count], ': max outlier ratio ', maxDevRatio)
            if minDevRatio > outlierRatio:
                print(lineNo[count], ': min outlier ratio', minDevRatio)
            count += 1

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(lineNo, chMin, 'bo', mfc='w')
    ax.plot(lineNo, chMax, 'bo', mfc='w')
    ax.errorbar(lineNo, chMean, chStd, capsize=3, marker='s', c='blue')
    plt.ylabel('Acceleration [m/s/s]', fontsize = 6)
    plotTitle = groupName + ': Acceleration Stats'
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    plt.show()


def checkStatcor(whizzFile, statcor, flight=''):
    """
    Plots `statcor` vs `flight` as a scatter plot. Used to compare with static gravimeter readings.
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    statcor : String
        The name of the channel containing the static corrections.
    flight : String, optional
        The name of the channel containing the flight number. The default is to use the line attribute.

    Returns
    -------
    None

    """
    
    filename = str(whizzFile)
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']

        num_lines = len(g.keys())
        flight_num = np.zeros((num_lines,))
        static_num = np.zeros((num_lines,))
        count = 0
        for line in g.keys():
            linegroup = g[line]
            if flight == '':
                flight_num[count] = linegroup.attrs['Flight']
            else:
                flight_num[count] = linegroup[flight][0]
            static_num[count] = linegroup[statcor][0]
            count += 1

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(flight_num, static_num, 'bo', mfc='w')
        plt.xlabel(flight, fontsize = 6)
        plt.ylabel(statcor, fontsize = 6)
        plotTitle = groupName + ': ' + ' Static Correction Analysis'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        plt.show()
    return


def checkConstantSlope(whizzFile, channels=[]):
    """
    Checks for constant slope (`np.diff`) in all the given channels of data.
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    channels : String List
        List of channel names from the database to be checked.

    Returns
    -------
    None

    """
    filename = str(whizzFile)
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if channels == []:
            lineGroups = list(g.values())
            channels = list(lineGroups[0].keys())
        data_is_good = True
        report = ''
        for line in g.keys():
            for channel in channels:
                data = gw.getLineData(g[line], channel)
                deriv = np.diff(data, n = 1)
                mean_deriv = np.mean(deriv)
                deriv = deriv - mean_deriv
                if len(deriv) > 10:
                    extremum = np.max(deriv) if np.max(deriv) > -np.min(deriv) else -np.min(deriv)
                    if extremum > mean_deriv / 1000.0:
                        report += f'\n  {line}; {channel} Largest difference (= {extremum:.3g}) > 0.1% of mean difference (= {(mean_deriv / 100.0):.3g})'
                        data_is_good = False
                else:
                    report += f'\n  {line}; {channel} data length {len(deriv)+1} is too short for analysis.'
                    data_is_good = False
    if data_is_good:
        report += 'All channels tested were either constant or of constant slope for all lines tested.'
    print(report)
    return


def diffChanStats(whizzFile, channel1, channel2):
    """
    Generate statistical plots for the difference between two channels across all lines. The plots show
    the min, mean, max and stdev for each channel as a function of line number.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    channel1 : String
        A channel to be subtracted from.
    channel2 : String
        A channel to subtract.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        gProject = f[groupName]
        g = gProject['Lines']
        lines = list(g.keys())
        numLines = len(g.items())
        
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
    
        for line in g.keys():
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
                print(f'Less than three real values in {lineNo[count]:.2f} for data, no statistics.')
            count += 1

        print(f'Overall mean difference is {allmean / numsamp :.2f}.')
        figtitle = wpl.make_plot_title(gProject)
        titlestr = channel1 + ' - ' + channel2 + ' Stats'
        wpl.plotBoxWhisker(chMin, chMax, chMean, chStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)
    return
    
    
def lineStats(whizzFile, channel):
    """
    Generate a statistical plot for the channel across all lines. The plot shows
    the min, mean, max and stdev for the channel as a function of line number.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    channel : String
        The name of the channel to plot.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    outlierRatio = 5.0

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        
        numLines = len(g.items())
        chMin = np.zeros((numLines,))
        chMax = np.zeros((numLines,))
        chMean = np.zeros((numLines,))
        chStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0
                
        for line in g.keys():
            #acqDate = g[line].attrs['Date_Local']
            lineNo[count] = line
            chMin[count] = np.min(g[line][channel])
            chMax[count] = np.max(g[line][channel])
            chMean[count] = np.mean(g[line][channel])
            chStd[count] = np.std(g[line][channel])
            maxDevRatio = (chMax[count] - chMean[count]) / chStd[count]
            minDevRatio = (chMean[count] - chMin[count]) / chStd[count]
            if maxDevRatio > outlierRatio:
                print(lineNo[count], ': max outlier ratio ', maxDevRatio)
            if minDevRatio > outlierRatio:
                print(lineNo[count], ': min outlier ratio', minDevRatio)
            count += 1
    figtitle = wpl.make_plot_title(f[groupName])
    titlestr = channel + ' Stats'
    wpl.plotBoxWhisker(chMin, chMax, chMean, chStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)
    

def statsChannelDiff(whizzFile, channel1, channel2, flightLines=[]):
    """
    Plot the statistics of channel1 - channel2 in each flightLine. 

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    flightLines : String List, optional
        A list of flightline, e.g. ['1000110.0']. Default is all lines in whizzFile.
    channel1 : String
        The name of a channel.
    channel2 : String
        The name of a channel.

    Returns
    -------
    None.

    """
    
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if flightLines == []:
            flightLines = list(g.keys())
        corr_units = g[flightLines[0]][channel1].attrs['Units']
        if not (g[flightLines[0]][channel1].attrs['Units'] == corr_units):
            print('Error: {channel1} and {channel2} do not have the same units.')
            return

        y_label = f'{channel1} - {channel2} [{corr_units}]'
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
            
        # initialise variables
        chMin = np.zeros((numLines,))
        chMax = np.zeros((numLines,))
        chMean = np.zeros((numLines,))
        chStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in flightLines:
            if line != 'CoordinateFrame':
                lineNo[count] = line
                dd = gw.getLineData(g[line], channel1) - gw.getLineData(g[line], channel2)

                if np.sum(~np.isnan(dd)) > 3:
                    chMin[count] = np.nanmin(dd)
                    chMax[count] = np.nanmax(dd)
                    chMean[count] = np.nanmean(dd)
                    chStd[count] = np.nanstd(dd)
                else:
                    chMin[count] = 0.0
                    chMax[count] = 0.0
                    chMean[count] = 0.0
                    chStd[count] = 0.0
                    print(f'Less than three real values in {lineNo[count]:.2f} for {channel1}, {channel2}, no statistics.')
                count += 1

        figtitle = wpl.make_plot_title(gProject)
        titlestr = channel1 + ' - ' + channel2 + ' Stats'
        wpl.plotBoxWhisker(chMin, chMax, chMean, chStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)
    return
