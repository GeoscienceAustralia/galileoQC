#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Aug 14 20:28:20 2021

@author: markdransfield

# GENERAL
checkSpikes()
allChanStats()
checkConstantSlope()
checkVertAcc()
checkVertAccStats()()
failsDeviation()
checkGaps()
checkPhase()
calcDrift()
checkErsHeaders()
lineStats()
psdChannelDiff()

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
import AirGravQC.gridFiles.gridfiles as grd
import AirGravQC.utility.utility as util
# import AirGravQC.gridfiles as grd
import AirGravQC.whizzPlots.whizzPlot as wpl

groupName = config.groupName

def checkPhase(filename, channel1, channel2):
    """
    For every survey line in a geoWhizz HDF5 file, given two fields, calculate
    the phase shift (in number of samples) required to maximise the correlation
    between the two fields

    Parameters
    ----------
    filename : String
        The name of a geoWhizz HDF5 file.
    channel1 : String
        The name of the first field.
    channel2 : String
        The name of the second field.
    
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


def failsDeviation(data, limit, nSamples):
    """
    Checks data, testing for more than `nSamples` consecutive instances where
    the data values are outside the range [-limit, limit].
    
    Parameters
    ----------
    data : Numpy 1D Float array
        The data to be checked.
    limit : Float
        The limit, assumed > 0.
    nSamples : Float
        Maximum allowed number of deviations allowed.

    Returns
    -------
    Bool
        True if the data failed.
    Float
        The number of exceedances (regardless of whether consecutive).
    
    """
    
    if np.max(data) < limit and np.min(data) > -limit:
        print('all ok')
        return False, 0
    
    indxSpikes = np.argwhere(np.logical_or(data > limit, data < -limit))
    numExceedances = indxSpikes.shape[0]
    if numExceedances < nSamples:
        print(f'only {numExceedances} exceedances')
        return False, numExceedances
    
    else:
        count = 1
        for jj in range(1, len(indxSpikes)):
            if indxSpikes[jj][0] == indxSpikes[jj-1][0] + 1:
                count += 1
            else:
                if count >= nSamples:
                    return True, count
                else:
                    count = 1
        if count >= nSamples:
            return True, count
        print(f'{numExceedances} but not consecutive, most {count}')
        return False, numExceedances
        

def checkSpikes(whizzFile, fields = [], numStd = 8.0):
    """
    Checks for spikes in all the given fields of data.
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    fields : String List
        List of field names from the database to be checked.
    numStd : Float, optional
        maximum allowed number of standard deviations allowed. The default is 8.0.

    Returns
    -------
    None

    """
    
    filename = str(whizzFile)
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if fields == []:
            lineGroups = list(g.values())
            fields = list(lineGroups[0].keys())
        noSpikes = True
        report = ''
        for line in g.keys():
            for field in fields:
                deriv = np.diff(g[line][field], n = 1)
                deriv = deriv - np.mean(deriv)
                if len(deriv) > 10:
                    myStd = np.std(deriv)
                    extremum = np.max(deriv) if np.max(deriv) > -np.min(deriv) else -np.min(deriv)
                    if extremum > numStd * myStd:
                        report += f'\n  {line}; {field} Extremum: {extremum:.3g} > {(numStd * myStd):.3g} = {numStd:.3g} x STD of {myStd:.3g}'
                        noSpikes = False
                else:
                    report += f'\n  {line}; {field} data length {len(deriv)+1} is too short for analysis.'
                    noSpikes = False
    print(report)
    return


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


def checkConstantSlope(whizzFile, fields=[]):
    """
    Checks for constant slope (`np.diff`) in all the given fields of data.
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    fields : String List
        List of field names from the database to be checked.

    Returns
    -------
    None

    """
    filename = str(whizzFile)
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if fields == []:
            lineGroups = list(g.values())
            fields = list(lineGroups[0].keys())
        data_is_good = True
        report = ''
        for line in g.keys():
            for field in fields:
                data = gw.getLineData(g[line], field)
                deriv = np.diff(data, n = 1)
                mean_deriv = np.mean(deriv)
                deriv = deriv - mean_deriv
                if len(deriv) > 10:
                    extremum = np.max(deriv) if np.max(deriv) > -np.min(deriv) else -np.min(deriv)
                    if extremum > mean_deriv / 1000.0:
                        report += f'\n  {line}; {field} Largest difference (= {extremum:.3g}) > 0.1% of mean difference (= {(mean_deriv / 100.0):.3g})'
                        data_is_good = False
                else:
                    report += f'\n  {line}; {field} data length {len(deriv)+1} is too short for analysis.'
                    data_is_good = False
    if data_is_good:
        report += 'All fields tested were either constant or of constant slope for all lines tested.'
    print(report)
    return


def checkGaps(whizzFile, maxGapSec=0.0, maxNumGaps=0):
    """
    Checks every dataset for each channel and each survey line in filePath for
    gaps, and reports all gaps found.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    maxGapSec :  Float, optional
        The largest allowed gap measured in seconds. Default 0.0
    maxNumGaps : Integer, optional
        The maximum number of gaps allowed on any survey line. Default 0

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        _reportGaps(g, maxGapSec, maxNumGaps)
        
        
def _reportGaps(group, maxGapSec=0.0, maxNumGaps=0):
    """
    Checks every dataset for each channel and each survey line in the HDF5 group
    for gaps, and reports all gaps found.

    Parameters
    ----------
    group : HDF5 Whizz file 'Lines' group
        The group containing the survey line data.
    maxGapSec :  Float, optional
        The largest allowed gap measured in seconds. Default 0.0
    maxNumGaps : Integer, optional
        The maximum number of gaps allowed on any survey line. Default 0

    Returns
    -------
    None.

    """
    lineGroups = list(group.values())
    channelNames = list(lineGroups[0].keys())
    num_channels = len(channelNames)
    num_lines_failed = 0
    total_num_lines = 0
    message = ''

    for line in group.keys():
        total_num_lines += 1
        gaps_on_line = 0
        lineNo = line
        lineText = f'Line {lineNo}'
        for channel in channelNames:
            numberMissing = np.count_nonzero(np.isnan(group[line][channel]))
            if numberMissing > 0:
                lineText += f'\n    {channel}, nans: {numberMissing}'
                gaps_on_line += 1
        if gaps_on_line > 0:
            num_lines_failed += 1
            message += lineText + '\n'
    print(f'Checking for all gaps in all {num_channels} channels on all {total_num_lines} lines.')
    print(message)
    print(f'{num_lines_failed} lines failed.')


def allChanStats(whizzFile, allChannels=[], lines=[], d1_chans=[], mr_chans=[], sin_chans=[]):
    """
    Generate statistical plots for the channels across all lines. The plots show
    the min, mean, max and stdev for each channel as a function of line number.

    The statistics can be optionally calculated on the data after first differencing,
    or mean removal, or both. In the latter case, first differencing is done first.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath.
        Name of a HDF5 Whizz file, including path and extension.
    allChannels : [String], optional.
        A list of the channels or fields to plot. Default is all in whizzFile.
    lines : String list, optional.
        The line numbers to be checked. Default is all lines in the whizzFile.
    d1_chans : [String], optional.
        An array of names of channels from `allChannels` whose first difference
        along each survey line should be calculated before the statistics.
    mr_chans : [String], optional.
        An array of names of channels from `allChannels` whose mean along each
        survey line should be subtracted before calculating statistics.
    sin_chans : [String], optional.
        An array of names of channels from `allChannels` whose sine
        along each survey line should be calculated before the statistics.

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
        if allChannels == []:
            lineGroups = list(g.values())
            allChannels = list(lineGroups[0].keys())
        numLines = len(lines)
        
        for channel in allChannels:
            remove_mean = False
            diff_one = False
            remove_sine = False
            if channel in mr_chans:
                remove_mean = True
            if channel in d1_chans:
                diff_one = True
            if channel in sin_chans:
                remove_sine = True

            # initialise variables
            chMin = np.zeros((numLines,))
            chMax = np.zeros((numLines,))
            chMean = np.zeros((numLines,))
            chStd = np.zeros((numLines,))
            lineNo = np.zeros((numLines,))
            count = 0

            # get the units for the y axis label
            dd = g[lines[0]][channel]
            xlabelstr = 'Line number'
            #if the channel has an attribute 'Units'
            ylabelstr = channel
            if 'Units' in dd.attrs.keys():
                ylabelstr += ' ' + dd.attrs['Units']
        
            for line in lines:
                if line != 'CoordinateFrame':
                    lineNo[count] = line
                    dd = gw.getLineData(g[line], channel)
                    if remove_sine:
                        dd = np.sin(dd * (np.pi / 180.0))
                    if diff_one:
                        dd = np.append(np.diff(dd), dd[-1]-dd[-2])
                    if remove_mean:
                        dd = dd - np.mean(dd)

                    # if np.sum(~np.isnan(g[line][channel])) > 3:
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
                        print(f'Less than three real values in {lineNo[count]:.2f} for {channel}, no statistics.')
                    count += 1

            figtitle = wpl.make_plot_title(gProject)
            titlestr = ''
            if remove_mean:
                titlestr += 'mr('
            if diff_one:
                titlestr += ' d1('
            if remove_sine:
                titlestr += 'sin('

            titlestr += channel

            if remove_mean:
                titlestr += ')'
            if diff_one:
                titlestr += ')'
            if remove_sine:
                titlestr += ')'
            titlestr += ' Stats'
            wpl.plotBoxWhisker(chMin, chMax, chMean, chStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)
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
        The name of the channel or field to plot.

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


def psdChannelDiff(whizzFile, channel1, channel2, flightLines=[]):
    """
    Plot the PSD (log-log Sqrt(Power) from welch method) of
    channel1 - channel2 in each flightLine. 

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
    import scipy.signal as sig
    global mean_speed
    
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
            
        for line in flightLines:
            mean_speed = _mean_line_speed(f[groupName], line)
            f_sample = _time_frequency(f[groupName])
            data = gw.getLineData(g[line], channel1) - gw.getLineData(g[line], channel2)#

            freq, Pxx = sig.welch(data, nfft=4*4096, fs = f_sample)
            period = 1.0 / freq[1:]
            rootPwr = np.sqrt(Pxx[1:]) / freq[1:]
            maxPwr = np.max(rootPwr)
            #print(f'{line} - low-f limit: {rootPwr[0]:.2f}, max: {maxPwr:.2f}')
            plt.loglog(period, rootPwr, color='blue', lw=0.3)
            if rootPwr[0] > 3.0:
                ax.text(period[0], rootPwr[0], f'{line}', fontsize=6)
    
        plt.ylim([1,1E5])
        plt.xlabel('Period [s]', fontsize = 6)
        plt.ylabel(y_label, fontsize = 6)
        plotTitle = f'{projName} : {y_label}'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        secax = ax.secondary_xaxis('top',
            functions=(_period_to_dist, _dist_to_period))
        secax.set_xlabel('wavelength [m]', fontsize=6)
        for label in secax.get_xticklabels(): label.set_fontsize(6)
        plt.show()


def _period_to_dist(p):
    """
    Converts period (sec) to distance by multiplying by the GLOBAL
    mean speed (m/s).
    For use in creating a secondary wavelength axis for `psdChannelDiff`.

    Parameters
    ----------
    p : Float
        The period in seconds.

    Returns
    -------
    Float. The distance in metres.

    """
    global mean_speed
    return mean_speed * p


def _dist_to_period(x):
    """
    Converts distance (metres) to period(sec) by dividing by the GLOBAL
    mean speed (m/s).
    For use in creating a secondary wavelength axis for `psdChannelDiff`.

    Parameters
    ----------
    x : Float
        The distance in metres.

    Returns
    -------
    Float. The period in seconds.

    """
    global mean_speed
    return x / mean_speed


def _time_frequency(group):
    """
    Returns the sample frequency of the data in the project group.
    Simply calculated as the inverse of the difference of the first
    two sample times in the first line. Relies on the TimeChannel
    attribute being set. 

    Parameters
    ----------
    group : HDF5 group
        Must be the project group (top level of hierarchy in whizz File).

    Returns
    -------
    Float. The sample frequency in Hz.

    """
    flightLines = list(group['Lines'].keys())
    time = group['CoordinateFrame'].attrs['TimeChannel']
    t = np.array(group['Lines'][flightLines[0]][time])
    return 1.0 / np.abs(t[1] - t[0])


def _mean_line_speed(group, line):
    """
    Returns the mean line speed of the data in the line (which is
    in the given project group). Instantaneous speeds are calculated
    for each sample along the line and their average value returned.
    Relies on the XChaneel, YChannel and TimeChannel attributes being set. 

    Parameters
    ----------
    group : HDF5 group
        Must be the project group (top level of hierarchy in whizz File).
    line : String
        The line identifier.

    Returns
    -------
    Float. The mean speed along the line.

    """
    xChannel = group['CoordinateFrame'].attrs['XChannel']
    yChannel = group['CoordinateFrame'].attrs['YChannel']
    tChannel = group['CoordinateFrame'].attrs['TimeChannel']
    xPos = np.array(group['Lines'][line][xChannel])
    yPos = np.array(group['Lines'][line][yChannel])
    sampleTime = np.array(group['Lines'][line][tChannel])
    dt = np.gradient(sampleTime)
    xVel = np.gradient(xPos) / dt
    yVel = np.gradient(yPos) / dt
    sample_speed = np.sqrt(xVel * xVel + yVel * yVel)
    return np.mean(sample_speed)


def checkErsHeaders(folderPath='\.'):
    """
    Compares all .ers files in the folder for certain key parameters and reports
    any that are different from the first file found in any parameter value.

    Parameters
    ----------
    folderPath : Path, optional
        The folder or directory containing the .ers files. The default is '\.'.

    Returns
    -------
    None.

    """
    reportString = ''
    allOk = True
    fileOK = True
    
    # get a list of the ers file paths
    file_count = 0
    ersFiles = []
    folderFiles = folderPath.iterdir()
    for aFile in folderFiles:
        if aFile.is_file() and (aFile.suffix == '.ERS' or aFile.suffix == '.ers') and (aFile.name[0] != '.'):
            ersFiles.append(aFile)
            file_count += 1
    print(f'Found {file_count} .ers files ...')
    print(f'in: {str(folderPath)}')

    # get the header dict for each ers file
    ersDicts = []
    for aFile in ersFiles:
            [ncells, nrows, nbands, eastings, northings, nullcell,
                    precision, headerbytes, originalnullcell, byteorder, datum,
                    projection] = grd.read_ers_header(str(aFile))
            commonOK, reportStr = _commonErsHdrErrors(ncells, nrows, nbands, nullcell,
                    precision, headerbytes, originalnullcell, byteorder, datum,
                    projection)
            ersDicts.append([ncells, nrows, nbands, nullcell,
                    precision, headerbytes, originalnullcell, byteorder, datum,
                    projection])
            
    # compare the contents line-by-line
    print(f'Comparing ERS files against {ersFiles[0].name}.')
    firstDict = ersDicts[0]
    print(firstDict)
    for jj in range(1, len(ersDicts)):
        fileOK = True
        print(f'Checking file {ersFiles[jj].name}')
        for ii in range(0, len(firstDict)):
            if firstDict[ii] != ersDicts[jj][ii]:
                allOk = False
                fileOK = False
                print(f'  Different element {ii}. {firstDict[ii]} != {ersDicts[jj][ii]}')
        if fileOK:
            print('  Checked OK.')
    
    return


def _commonErsHdrErrors(ncells, nrows, nbands, nullcell, precision, headerbytes, originalnullcell, byteorder, datum, projection):
    """
    Checks a set of ERS header fields for various common errors.

    Parameters
    ----------
    ncells : Integer
        DESCRIPTION.
    nrows : Integer
        DESCRIPTION.
    nbands : Integer
        DESCRIPTION.
    nullcell : Integer
        DESCRIPTION.
    precision : Integer
        DESCRIPTION.
    headerbytes : Integer
        DESCRIPTION.
    originalnullcell : Integer
        DESCRIPTION.
    byteorder : Integer
        DESCRIPTION.
    datum : Integer
        DESCRIPTION.
    projection : Integer
        DESCRIPTION.

    Returns
    -------
    allOk : Bool
        True if no errors found, else False.
    reportStr : String
        A report listing all errors found.

    """
    allOk = True
    reportStr = ''
    
    if projection == 'WGS84':
        reportStr += 'ERROR: projection cannot be WGS84'
        allOk = False
    
    return allOk, reportStr




#=========================
#
# NAV
#
#=========================


def checkGNSS(whizzFile, num_sats, pdop, vdop, hdop, nsats_min=4, max_pdop=6, max_vdop=4, max_hdop=4, lines=[]):
    """
    Checks that the data in a whizzFile meets the requirements for the minimum
    number of satellites, and maximum PDOP, VDOP and HDOP.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    num_sats : String
        Name of the channel containing the number of satellites visible for
        each measurement.
    pdop : String
        Name of the channel containing the PDOP for each measurement.
    vdop : String
        Name of the channel containing the VDOP for each measurement.
    hdop : String
        Name of the channel containing the HDOP for each measurement.
    nsats_min : Integer, optional
        The minimum number of satellites required, default 4.
    max_pdop : Integer, optional
        The maximum PDOP allowed, default 6
    max_vdop : Integer, optional
        The maximum VDOP allowed, default 6
    max_hdop : Integer, optional
        The maximum HDOP allowed, default 6

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    report = ''
    num_failed_lines = 0

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        xchan = f[groupName]['CoordinateFrame'].attrs['XChannel']
        ychan = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if lines == []:
            lines = g.keys()

        error_count = 0
        for line in lines:
            x = gw.getLineData(g[line], xchan)
            y = gw.getLineData(g[line], ychan)

            nsats_data = gw.getLineData(g[line], num_sats)
            min_nsats_data = np.nanmin(nsats_data)
            if min_nsats_data < nsats_min:
                xmin_fail = np.nanmin(x[nsats_data < nsats_min])
                ymin_fail = np.nanmin(y[nsats_data < nsats_min])
                xmax_fail = np.nanmax(x[nsats_data < nsats_min])
                ymax_fail = np.nanmax(y[nsats_data < nsats_min])
                numfail = np.count_nonzero(nsats_data < nsats_min)
                report += f'Line {line} failed for {numfail} fids: min sats = {min_nsats_data:.0f} < {nsats_min}.'
                report += f' Easting [{xmin_fail:.0f}, {xmax_fail:.0f}],'
                report += f' Northing [{ymin_fail:.0f}, {ymax_fail:.0f}].\n'
                error_count += 1

            pdop_data = gw.getLineData(g[line], pdop)
            max_pdop_data = np.nanmax(pdop_data)
            if max_pdop_data > max_pdop:
                xmin_fail = np.nanmin(x[pdop_data > max_pdop])
                ymin_fail = np.nanmin(y[pdop_data > max_pdop])
                xmax_fail = np.nanmax(x[pdop_data > max_pdop])
                ymax_fail = np.nanmax(y[pdop_data > max_pdop])
                numfail = np.count_nonzero(pdop_data > max_pdop)
                report += f'Line {line} failed for {numfail} fids: max pdop = {max_pdop_data:.1f} > {max_pdop}.'
                report += f' Easting [{xmin_fail:.0f}, {xmax_fail:.0f}],'
                report += f' Northing [{ymin_fail:.0f}, {ymax_fail:.0f}].\n'
                error_count += 1

            vdop_data = gw.getLineData(g[line], vdop)
            max_vdop_data = np.nanmax(vdop_data)
            if max_vdop_data > max_vdop:
                xmin_fail = np.nanmin(x[vdop_data > max_vdop])
                ymin_fail = np.nanmin(y[vdop_data > max_vdop])
                xmax_fail = np.nanmax(x[vdop_data > max_vdop])
                ymax_fail = np.nanmax(y[vdop_data > max_vdop])
                numfail = np.count_nonzero(vdop_data > max_vdop)
                report += f'Line {line} failed for {numfail} fids: max vdop = {max_vdop_data:.1f} > {max_vdop}.'
                report += f' Easting [{xmin_fail:.0f}, {xmax_fail:.0f}],'
                report += f' Northing [{ymin_fail:.0f}, {ymax_fail:.0f}].\n'
                error_count += 1

            hdop_data = gw.getLineData(g[line], hdop)
            max_hdop_data = np.nanmax(hdop_data)
            if max_hdop_data > max_hdop:
                xmin_fail = np.nanmin(x[hdop_data > max_hdop])
                ymin_fail = np.nanmin(y[hdop_data > max_hdop])
                xmax_fail = np.nanmax(x[hdop_data > max_hdop])
                ymax_fail = np.nanmax(y[hdop_data > max_hdop])
                numfail = np.count_nonzero(hdop_data > max_hdop)
                report += f'Line {line} failed for {numfail} fids: max hdop = {max_hdop_data:.1f} > {max_hdop}.'
                report += f' Easting [{xmin_fail:.0f}, {xmax_fail:.0f}],'
                report += f' Northing [{ymin_fail:.0f}, {ymax_fail:.0f}].\n'
                error_count += 1

        if error_count == 1:
            errstr = f'Found 1 error.'
        elif error_count > 1:
            errstr = f'Found {error_count} errors.'
        else:
            errstr = f'Found 0 errors.'
        report = f'In {projName}, checked num sats, PDOP, VDOP and HDOP. {errstr}\n' + report
        print(report)
        

def checkHeading(whizzFile, nominalHeadings, lines = [], x='', y='', tolerance=10.0, plot_flag=False):
    """
    Checks heading in degrees is within +/- tolerance (in degrees) of nominal (in degrees). Actually
    checks against `sin(nominalHeading +/- tolerance)`.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    nominalHeadings : [Float]
        The desired headings in degrees from north.
    x : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    y : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    tolerance : Float, optional
        Headings within +/- tolerance degrees of nominalHeading are ok.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    report = ''
    num_failed_lines = 0

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if lines == []:
            lines = g.keys()
        numLines = len(lines)

        for line in lines:
            dx = np.diff(gw.getLineData(g[line], x))
            dy = np.diff(gw.getLineData(g[line], y))
            allok = True
            for nomhead in nominalHeadings:
                tol_1 = np.cos(np.pi * (nomhead + tolerance) / 180.0)
                tol_2 = np.cos(np.pi * (nomhead - tolerance) / 180.0)
                upper_limit = max(tol_1, tol_2)
                lower_limit = min(tol_1, tol_2)
                heading = np.arctan2(dx, dy) * 180.0 / np.pi
                min_heading = np.nanmin(heading)
                max_heading = np.nanmax(heading)
                mean_heading = np.mean(heading)
                cosheading = np.cos(np.pi * heading / 180.0)
                allok = all(h <= upper_limit for h in cosheading) and all(h >= lower_limit for h in cosheading)
                if allok:
                    break
            
            if not allok:
                num_failed_lines += 1
                report += f'Line {line}: heading range exceeded. Mean {mean_heading:.2f}, '
                report += f'Min {min_heading:.2f}, Max {max_heading:.2f} deg.\n'
                if plot_flag:
                    fig = plt.figure()
                    fig.suptitle(f'Heading Check Line {line}', fontsize=10)
                    fig.subplots_adjust(top=0.85)
                    
                    ax = fig.add_subplot(1,1,1)
                    thou_format = tkr.FuncFormatter(util._space_thou)
                    ax.plot(gw.getLineData(g[line], x)[1:], heading, 'b', mfc='w')
                    ax.xaxis.set_major_formatter(thou_format)
                    plt.ylabel('Estimated heading [deg]', fontsize = 6)
                    plt.xlabel(f'{x} [m]', fontsize = 6)
                    plt.grid(True)
                    for label in ax.get_xticklabels(): label.set_fontsize(6)
                    for label in ax.get_yticklabels(): label.set_fontsize(6)

    # print(f'Heading limits: [{nominalHeading}, +/-{tolerance}] deg or equivalent.')
    print(f'  Checked {numLines} lines, {num_failed_lines} failed.\n')
    print(report)
    if plot_flag:
        plt.show()
    return
    

def checkLineLengths(whizzFile, min_len=50.0, measX='', measY=''):
    """
    Checks that all lines in whizzFile are at least min_len km long.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    min_len : TYPE, optional
        The minimum allowed line length in km. The default is 50.0.
    measX : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    measY : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.

    Returns
    -------
    None.

    """
    measFile = str(whizzFile)
    
    with h5py.File(measFile, 'r') as f:
        gMeas = f[groupName]['Lines']
        if measX == '':
            measX = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if measY == '':
            measY = f[groupName]['CoordinateFrame'].attrs['YChannel']
        
        num_failed_lines = 0
        for line in gMeas.keys():
            xM = np.array(gMeas[line][measX])
            yM = np.array(gMeas[line][measY])
            line_length = _displacement2(xM[0], xM[-1], yM[0], yM[-1])
            if line_length < min_len * 1000.0:
                num_failed_lines += 1
                print(f'Line {line} length = {line_length:.1f} less than allowed min {min_len*1000.0:.1f}')
        print(f'Number failed lines = {num_failed_lines}')
            

def checkOverlaps(whizzFile, min_overlap = 7.6, lines = [], verbose=False, plot_flag=False):
    """
    For every line in the file whizzFile, calculate the overlap with each other
    line that has the same prefix. Plot a map. Report overlaps.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    min_overlap : Float, optional
        The minimum overlap distance in km, default 7.6 km.
    lines : String list, optional.
        The line numbers to be checked. Default is all lines in the whizzFile.
    verbose : Bool, optional
        If True, report status of all overlaps, else only report errors. Default False.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line.

    Returns
    -------
    None.

    """
    from matplotlib.ticker import StrMethodFormatter
    measFile = str(whizzFile)
    if plot_flag:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    
    with h5py.File(measFile, 'r') as f:
        gMeas = f[groupName]['Lines']
        east = f[groupName]['CoordinateFrame'].attrs['XChannel']
        nrth = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if lines == []:
            lines = list(gMeas.keys())
        num_coinc_lines = 0
        report = ''

        for idx1 in range(0, len(lines)): #line1 in lines: #gMeas.keys():
            line1 = lines[idx1]
            line1_plan = gMeas[line1].attrs['PlannedLine']

            for idx2 in range (idx1 + 1, len(lines)):#line2 in gMeas.keys():
                line2 = lines[idx2]
                line2_plan = gMeas[line2].attrs['PlannedLine']
                # if the second line isn't the first line but has the same planned line no.
                if line1 != line2 and line1_plan == line2_plan:
                    num_coinc_lines += 1
                    # extract positions
                    n1 = np.array(gMeas[line1][nrth])
                    e1 = np.array(gMeas[line1][east])
                    n2 = np.array(gMeas[line2][nrth])
                    e2 = np.array(gMeas[line2][east])
                    # get line direction in radians
                    dirn = np.arctan2((e1[-1] - e1[0]), (n1[-1] - n1[0]))
                    
                    # make life easier by transforming to a 2D problem
                    n0 = n1[0]
                    e0 = e1[0]
                    (y1, x1) = _rotateCoords(e1 - e0, n1 - n0, -dirn)
                    (y2, x2) = _rotateCoords(e2 - e0, n2 - n0, -dirn)

                    # find ends of lines
                    min1 = np.nanmin(x1)
                    max1 = np.nanmax(x1)
                    min2 = np.nanmin(x2)
                    max2 = np.nanmax(x2)
                    overlap = 0.0
                    whole_line = 0
                    
                    # find overlap length if any
                    if (min2 < min1 and min1 < max2) and not (min2 < max1 and max1 < max2):
                        overlap = abs(max2 - min1)
                    elif (min2 < max1 and max1 < max2) and not (min2 < min1 and min1 < max2):
                        overlap = abs(min2 - max1)
                    elif (min2 < max1 and max1 < max2) and (min2 < min1 and min1 < max2):
                        overlap = abs(max1 - min1)
                        whole_line = 1
                    elif (min1 < max2 and max2 < max1) and (min1 < min2 and min2 < max1):
                        overlap = abs(max2 - min2)
                        whole_line = 2

                    if plot_flag:
                        plotline1, = ax.plot(e1, n1, color='blue', lw=0.2)
                    
                    if whole_line == 1:
                        if verbose:
                            report += f'  OK: All of line {line1} in {line2}: length {overlap:.0f} m.\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='green', lw=0.2)
                    elif whole_line == 2:
                        if verbose:
                            report += f'  OK: All of line {line2} in {line1}: length {overlap:.0f} m.\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='green', lw=0.2)
                    elif overlap < 1:
                        if verbose:
                            report += f'  OK: Non-overlapping lines {line1}, {line2}: overlap {overlap:.0f} m < {min_overlap * 1000}.\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='green', lw=0.2)
                    elif overlap < min_overlap * 1000.0:
                        report += f'\n  ERROR: Repeat line {line1}, {line2}: overlap {overlap:.0f} m < {min_overlap * 1000}.\n\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='red', lw=0.2)
                    else:
                        if verbose:
                            report += f'  OK: Repeat line {line1}, {line2}: overlap by {overlap:.0f} m.\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='green', lw=0.2)

    print(f'{num_coinc_lines} coincident lines found.')
    if report == '':
        report = f'All overlaps meet requirement (>{min_overlap:.1f} km).'

    print(report)
    if num_coinc_lines > 0 and plot_flag:
        ax.set_aspect('equal')
        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        plt.xlabel('X [m]', fontsize = 10)
        plt.ylabel('Y [m]', fontsize = 10)
        plt.suptitle('Overlap Map', fontsize = 12)
        plt.title('[1st line blue, accepted line green, error line red]', fontsize = 10)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(10)
        for label in ax.get_yticklabels(): label.set_fontsize(10)
        plt.show()


def checkXYPlan(planPath, measPath, lines=[], planX='', planY='', measX='', measY='', allowance=200.0, maxCounter=14, maxDistance=0, known='', plot_flag=False, verbose=False):
    """
    Reports exceedances of actual horizontal position from planned horizontal
    positions for an airborne survey Whizz database.
    The positions (planX, planY) of the start and end of each planned survey line
    are read from planPath. The measured positions (measX, measY) are read from
    measPath and the perpendicular distance of each from the planned line is
    calculated. If this distance exceeds allowance for maxCounter consecutive
    fiducial, or maxDistance consecutive metres, then an out-of-specification
    exceedance is reported for that line.

    Parameters
    ----------
    planPath : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, with the survey
        positions plan.
    planX : String, optional
        The name of the geoWhizz field or channel containing the planned x position. The
        default is to read the xChannel field name from the Coordinate Frame.
    planY : String, optional
        The name of the geoWhizz field or channel containing the planned y position. The
        default is to read the yChannel field name from the Coordinate Frame.
    measPath : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, with the survey
        measured data.
    measX : String, optional
        The name of the geoWhizz field or channel containing the measured x position. The
        default is to read the xChannel field name from the Coordinate Frame.
    measY : String, optional
        The name of the geoWhizz field or channel containing the measured y position. The
        default is to read the yChannel field name from the Coordinate Frame.
    allowance : Float, optional
        The allowed horizontal distance for the measured line from the planned
        line. If any portion of a measured line is further than this from the
        planned position for more than the maximum allowed number of fids or the
        maximum allowed distance, then the line fails. The default is 200.0.
    maxCounter : Int, optional
        The maximum number of consecutive fids for which an exceedance
        greater than allowance is permitted. If 0, then the constraint is
        ignored. The default is 14.
    maxDistance : Float, optional
        The maximum number of consecutive metres for which an exceedance
        greater than allowance is permitted. If 0, then the constraint is
        ignored. The default is 0.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line.

    Returns
    -------
    None.

    """
    start_x = 0.0
    end_x = 0.0
    start_y = 0.0
    end_y = 0.0

    exceedances_known = False
    this_exc_known = False
    number_known = 0

    planfile = str(planPath)
    measFile = str(measPath)
    
    with h5py.File(planfile, 'r') as fp:
        gPlan = fp[groupName]['Lines']
        if planX == '':
            planX = fp[groupName]['CoordinateFrame'].attrs['XChannel']
        if planY == '':
            planY = fp[groupName]['CoordinateFrame'].attrs['YChannel']

        with h5py.File(measFile, 'r+') as fm:
            gMeas = fm[groupName]['Lines']
            if measX == '':
                measX = fm[groupName]['CoordinateFrame'].attrs['XChannel']
            if measY == '':
                measY = fm[groupName]['CoordinateFrame'].attrs['YChannel']
            numLines = len(gMeas.items())

            message = ''
            num_lines_exceeded = 0
            total_num_excs = 0
            num_lines_unplanned = 0


            if lines == []:
                lines = gMeas.keys()

            for line in lines:
                planLine = f"{gMeas[line].attrs['PlannedLine']:.1f}"
                lineName = _get_lineName(gMeas[line])
                exceedance_in_line = False
                if planLine in gPlan:
                    xP = np.array(gPlan[planLine][planX])
                    yP = np.array(gPlan[planLine][planY])
                    xM = np.array(gMeas[line][measX])
                    yM = np.array(gMeas[line][measY])
                    max_deviation = 0.0

                    if known != '':
                        exceedances_known = True
                        exc_known = np.array(gMeas[line][known])
                        report_known = -1

                    # rotate to line of x ~ 0 using line direction in radians
                    # if abs(yP[-1] - yP[0]) > epsilon:
                    dirn = np.arctan2((yP[-1] - yP[0]), (xP[-1] - xP[0]))
                    # else:
                    #     dirn = np.pi / 2.0    

                    [x, y] = _rotateCoords(xM - xP[0], yM - yP[0], dirn)

                    num_fids_in_exceedance = 0
                    exceedance_in_line = False

                    for fid in range(0, len(x)):
                        # There is an exceedance ...
                        if np.abs(y[fid]) > allowance: #x
                            # If a new exceedance, then initialise variables;
                            if num_fids_in_exceedance == 0:
                                start_x = xM[fid]
                                start_y = yM[fid]
                                num_fids_in_exceedance = 1
                                max_deviation = abs(y[fid]) - allowance #x
                                this_exc_known = False

                            # Else increment and update on the current exceedance.
                            else:
                                num_fids_in_exceedance += 1
                                max_deviation = max(np.abs(y[fid]) - allowance, max_deviation) #x
                                if exceedances_known:
                                    if exc_known[fid] > 0:
                                        report_known = exc_known[fid]
                                        this_exc_known = True

                        else:
                            if num_fids_in_exceedance > 0: # the current exceedance has ended
                                end_x = xM[fid]
                                end_y = yM[fid]
                                len_exceedance = util._length([start_x, end_x], [start_y, end_y])[1]
                                if _exceedance_fail(num_fids_in_exceedance, len_exceedance, maxCounter, maxDistance):
                                    message += f'\nL {lineName} deviates more than {allowance} m for '
                                    message += f'{num_fids_in_exceedance} fids ({len_exceedance:.0f} m), max exceedance = {max_deviation:.0f} m.'
                                    message += f'\n  From ({start_x:.0f} E {start_y:.0f} N) to ({end_x:.0f} E {end_y:.0f} N).'
                                    total_num_excs += 1
                                    if this_exc_known:
                                        number_known += 1
                                        message += f' Known exceedance {report_known}.'
                                        this_exc_known = False
                                    exceedance_in_line = True
                                num_fids_in_exceedance = 0
                else:
                    print(f'Line {lineName} / {planLine} not in plan.')
                    num_lines_unplanned += 1
                if exceedance_in_line:
                    num_lines_exceeded += 1
                    if plot_flag:
                        if abs(np.cos(dirn)) > 0.5:
                            _plot_exceeding_line(x, y, xP, yP, xM, yM, measX, measY, allowance, line, planLine, dirn)
                        else:
                            _plot_exceeding_line(x, y, yP, xP, yM, xM, measY, measX, allowance, line, planLine, dirn)

            message = f'\n{num_lines_exceeded} lines with horizontal exceedances.\n' + message # 5 DEC
            message = f'\n{total_num_excs} horizontal exceedances.\n' + message # 5 DEC
            message = f'\n{num_lines_unplanned} lines not in plan and not checked.\n' + message # 5 DEC
            message = f'\n{number_known} exceedances known in the database.\n' + message # 5 DEC
            print(message)
            if plot_flag:
                plt.show()


def _get_lineName(linegroup):
    lineNo = linegroup.attrs['LineNumber']
    lineName = f'{lineNo:.2f}'
    if 'Flight' in linegroup.attrs:
        flightNo = linegroup.attrs['Flight']
        lineName += ":" + f'{flightNo:.0f}'
    return lineName


def _exceedance_fail(num_fids_in_exceedance, len_exceedance, maxCounter, maxDistance):
    """
    Given a specification, maxCounter, on number of fids and, maxDistance,
    on distance, checks to see if either  num_fids_in_exceedance > maxCounter
    or len_exceedance > maxDistance. If either is True, returns True.

    Parameters
    ----------
    num_fids_in_exceedance : Integer
        DESCRIPTION.
    len_exceedance : Float
        DESCRIPTION.
    maxCounter : Integer
        DESCRIPTION.
    maxDistance : Float
        DESCRIPTION.

    Returns
    -------
    Bool.

    """
    if maxCounter < 1 and maxDistance < 1:
        return False
    if maxCounter > 0:
        if num_fids_in_exceedance > maxCounter:
            return True
    if maxDistance > 0:
        if len_exceedance > maxDistance:
            return True
    return False


def _plot_exceeding_line(x, y, xP, yP, xM, yM, measX, measY, allowance, line, planLine, dirn):
    """
    Plots a standard exceedance figure for checkXYPlan().

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    allowance : TYPE
        DESCRIPTION.
    line : TYPE
        DESCRIPTION.
    planLine : TYPE
        DESCRIPTION.
    dirn : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    fig = plt.figure()
    plot_title = f'Line {line} (plan: {planLine}); Bearing {dirn * 180 / np.pi:.1f} deg E'
    fig.suptitle(plot_title, fontsize=10)

    ax = fig.add_subplot(2,1,1)
    ax.plot(x, y, 'b')
    ax.plot(x, -allowance * np.ones(y.shape), 'r')
    ax.plot(x, allowance * np.ones(y.shape), 'r')
    ax.set_xlim(x[0], x[-1])
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax.set_xlabel('deviation from planned line [m]', fontsize = 8)
    ax.set_ylabel('distance along line', fontsize = 8)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    ax.grid()

    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(xP, yP, color='darkorange', label='Plan', lw=1.5, alpha=0.7)
    ax2.plot(xM, yM, color='blue', label='Measured', lw=1.5, alpha=0.7)
    ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax2.set_xlim(xM[0], xM[-1])
    ax2.set_ylabel(f'{measY} [m]', fontsize=8)
    ax2.set_xlabel(f'{measX} [m]', fontsize=8)
    for label in ax2.get_xticklabels(): label.set_fontsize(6)
    for label in ax2.get_yticklabels(): label.set_fontsize(6)
    ax2.legend(fontsize=8)
    ax2.grid()
    return

   
def checkVertPlan(planPath, measPath, lines=[], planX='', planY='', planZ='', measX='', measY='', measZ='', allowance=30.0, maxCounter=13, maxDistance=0.0, known='', plot_flag=False):
    """
    Reports exceedances of actual vertical position from planned vertical positions
    (stored in a plan file) for an airborne survey Whizz database.

    The positions (`planX`, `planY`, `planY`) of each planned survey line
    are read from `planPath`. The measured positions (`measX`, `measY`, `measZ`) are read
    from `measPath` and the vertical distance of each from the planned line is
    calculated. If this distance exceeds `allowance` for a distance greater than
    `maxDistance`, then an out-of-specification exceedance is reported for that line.
    If `maxDistance` is less than 1.0, then the test is instead against `maxCounter`
    consecutive positions. The default for `maxCounter` is 13.

    See also `checkClearance`, `checkSafeClearance`, `checkDrape`.

    Parameters
    ----------
    planPath : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, of survey plan.
    measPath : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, of measured data.
    lines : Array{String}, optional
        Array of line numbers as strings. Default = [], meaning all lines are checked.
    planX : String, optional
        The name of the geoWhizz field or channel containing the planned x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    planY : String, optional
        The name of the geoWhizz field or channel containing the planned y positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    measX : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    measY : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    allowance : Float, optional
        The maximum allowed deviation in height. The default is 30.0.
    maxCounter : Int, optional
        The maximum number of consecutive fids for which an exceedance
        greater than allowance is permitted. The default is 13.
    maxDistance : Float, optional
        The maximum number of consecutive metres for which an exceedance
        greater than allowance is permitted. If 0, then `maxCounter` is
        used instead of `maxDistance`. The default is 0.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line.

    Returns
    -------
    None.

    TODO: If maxCounter unspecified, and maxDistance is specified, test against maxDistance.
          Use _exceedance_fail() to do this.

    """
    planfile = str(planPath)
    measFile = str(measPath)
    if maxDistance > 1.0:
        counting = False
    elif maxCounter > 0:
        counting = True
    else:
        print(f'ERROR. Require maxDistance > 1.0 m or maxCounter > 0')
        print(f'  maxDistance = {maxDistance}, maxCounter = {maxCounter}.')
    
    exceedances_known = False
    this_exc_known = False
    number_exc_known = 0

    with h5py.File(planfile, 'r') as fp:
        gPlan = fp[groupName]['Lines']
        if planX == '':
            planX = fp[groupName]['CoordinateFrame'].attrs['XChannel']
        if planY == '':
            planY = fp[groupName]['CoordinateFrame'].attrs['YChannel']
        if planZ == '':
            planZ = fp[groupName]['CoordinateFrame'].attrs['AltitudeChannel']

        with h5py.File(measFile, 'r') as fm:
            gMeas = fm[groupName]['Lines']
            if measX == '':
                measX = fm[groupName]['CoordinateFrame'].attrs['XChannel']
            if measY == '':
                measY = fm[groupName]['CoordinateFrame'].attrs['YChannel']
            if measZ == '':
                measZ = fm[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
            numErrors = int(0)
            numErrLines = 0
            num_lines_unplanned = 0
            report = ''

            if lines == []:
                lines = gMeas.keys()

            for line in lines:
                line_flagged = False
                lineName = _get_lineName(gMeas[line])
                planLine = f"{gMeas[line].attrs['PlannedLine']:.1f}"
                if planLine in gPlan: # 5 DEC
                    xP = np.array(gPlan[planLine][planX])
                    yP = np.array(gPlan[planLine][planY])
                    zP = np.array(gPlan[planLine][planZ])
                    xM = np.array(gMeas[line][measX])
                    yM = np.array(gMeas[line][measY])
                    zM = np.array(gMeas[line][measZ])

                    if known != '':
                        exceedances_known = True
                        exc_known = np.array(gMeas[line][known])
                        report_known = -1
                    
                    # make life easier by transforming to a 2D problem
                    dirn = np.arctan2((yM[-1] - yM[0]), (xM[-1] - xM[0]))
                    (xm, ym) = _rotateCoords(xM - xP[0], yM - yP[0], dirn)
                    (xp0, yp0) = _rotateCoords(xP - xP[0], yP - yP[0], dirn)
                    xp = xp0
                    yp = yp0
                    zp = zP
                    
                    # interpolate (xm, zM) onto (xp, zmp)
                    if abs(xm[-1] - xm[0]) < abs(ym[-1] - ym[0]):
                        print('ERROR - expect xms > yms but this is not so.')
                    (zmp, zM_trim) = gw.interpolateLine(xp, zp, xm, zM, plot_flag=False)
                    # calculate the deviation vector
                    z_dev = zM_trim - zmp
                
                    # check vertical deviations
                    # initialise fiducial counter and error report for line, assume no exceedances
                    fid = int(0)
                    exceeding = False
                
                    for one_x in z_dev:
                        fid += 1
                        in_spec = abs(one_x) < allowance
                        
                        if not (exceeding or in_spec):
                            start_fid = fid
                            exceeding = True
                            exc_fids = 0
                            this_exc_known = False
                        elif exceeding and not in_spec:
                            exc_fids += 1
                            if exceedances_known:
                                if exc_known[fid] > 0:
                                    report_known = exc_known[fid]
                                    this_exc_known = True
                        elif exceeding and in_spec:
                            num_fids = exc_fids
                            exceeding = False
                            ex0 = xM[start_fid]
                            ex1 = xM[start_fid + num_fids]
                            ey0 = yM[start_fid]
                            ey1 = yM[start_fid + num_fids]
                            exc_dist = _displacement2(ex0, ex1, ey0, ey1)
                            if (counting and num_fids > maxCounter) or (not counting and exc_dist > maxDistance):
                                if not line_flagged:
                                    numErrLines += 1
                                    line_flagged = True
                                numErrors += 1
                                max_dev = np.nanmax(abs(z_dev[start_fid:start_fid + num_fids]))
                                report += f'\nL {lineName} deviates more than {allowance:.1f} m for'
                                report += f' {num_fids} fids ({exc_dist:.0f} m),'
                                report += f' max exceedance = {max_dev - allowance:.1f} m.'
                                report += f'\n  From ({ex0:.0f} E, {ey0:.0f} N) to ({ex1:.0f} E, {ey1:.0f} N).'
                                if this_exc_known:
                                    number_exc_known += 1
                                    report += f' Known exceedance {report_known}.'
                                    this_exc_known = False
                    if plot_flag and line_flagged:
                        if abs(np.cos(dirn)) > 0.5:
                            _plot_vert_exceedance(xm, z_dev, xP, zP, xM, zM, measX, measZ, allowance, line, planLine, dirn)
                        else:
                            _plot_vert_exceedance(xm, z_dev, yP, zP, yM, zM, measY, measZ, allowance, line, planLine, dirn)
                else:
                    print(f'Line {lineName} not in plan.')
                    num_lines_unplanned += 1

        print(f'\n{num_lines_unplanned} lines not in plan and not checked.\n')
        print(f'Total number of exceedances = {numErrors} over {numErrLines} erroneous lines.')
        print(f'\n{number_exc_known} exceedances known in the database.\n')
        print(report)
        if plot_flag:
            plt.show()
                                       
   
def _rotateCoords(x, y, angle):
    """
    Rotates Cartesian vectors x and y, by `angle` to xr and yr.

    Parameters
    ----------
    x : numpy Float 1D array
        DESCRIPTION.
    y : numpy Float 1D array
        DESCRIPTION.
    angle : Float
        Angle in radians by which coordinates are to be rotated.

    Returns
    -------
    xr : numpy Float 1D array
        DESCRIPTION.
    yr : numpy Float 1D array
        DESCRIPTION.

    """
    xr = np.cos(angle) * x + np.sin(angle) * y
    yr = np.sin(angle) * x - np.cos(angle) * y
    return xr, yr
 

def _plot_vert_exceedance(xm, z_dev, xP, zP, xM, zM, measX, measZ, allowance, line, planLine, dirn):
    fig = plt.figure()
    plot_title = f'Line {line} (plan: {planLine}); Bearing {dirn * 180 / np.pi:.1f} deg E'
    fig.suptitle(plot_title, fontsize=10)
    
    ax = fig.add_subplot(2,1,1)
    ax.plot(xm[1:], z_dev, 'b', lw=0.6)
    ax.plot(xm[1:], -allowance * np.ones(z_dev.shape), 'r')
    ax.plot(xm[1:], allowance * np.ones(z_dev.shape), 'r')
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax.set_xlim(xm[0], xm[-1])
    ax.set_ylabel('deviation from planned drape [m]', fontsize=8)
    ax.set_xlabel('distance along line [m]', fontsize=8)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    ax.grid()

    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(xP, zP, color='darkorange', label='Plan', lw=1.5, alpha=0.7)
    ax2.plot(xM, zM, color='blue', label='Measured', lw=1.5, alpha=0.7)
    ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax2.set_xlim(xM[0], xM[-1])
    ax2.set_ylabel(f'{measZ} [m]', fontsize=8)
    ax2.set_xlabel(f'{measX} [m]', fontsize=8)
    for label in ax2.get_xticklabels(): label.set_fontsize(6)
    for label in ax2.get_yticklabels(): label.set_fontsize(6)
    ax2.legend(fontsize=8)
    ax2.grid()
    return


def checkSafeClearance(whizzFile, minimumAllowedClearance, clearance_chan='', altitude_chan='', terrain_chan='', xChannel='', yChannel='', plot_flag=False):
    """
    Checks the data from Whizz HDF5 file for low clearance above the terrain (as a safety check).
    
    The clearance (calculated as altitude - terrain if clearance_chan=''), must be
    greater than minimumAllowedClearance or a warning is reported.
    If any samples along the line fail the test, then profiles are plotted
    for visual analysis.
    The positions (xChannel and yChannel) are only used for plotting.

    See also `checkClearance`, `checkDrape`, `checkVertPlan`.    

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    minimumAllowedClearance : Float
        The planned clearance between above the terrain in metres.
    clearance_chan : String, optional
        The name of the terrain clearance field in the Lines group
        of the Whizz HDF5 file.
    altitude_chan : String, optional
        The name of the absolute altitude or height field in the Lines group
        of the Whizz HDF5 file.
    terrain_chan : String, optional
        The name of the absolute terrain or DTM height field in the Lines group
        of the Whizz HDF5 file.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    plot_flag : Bool, optional
        If True, plot a map of each pair of intersecting lines where the `zChannel`
        values differ by more than `max_allowed_deltaZ`. Default False.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    if clearance_chan == '':
        if altitude_chan == '' or terrain_chan == '':
            print('ERROR - either clearance_chan, or both altitude_chan and terrain_chan must be specified.')

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        num_failed_lines = 0
        report = ''

        for line in g.keys():
            lineName = _get_lineName(g[line])
            if clearance_chan == '':
                alt = gw.getLineData(g[line], altitude_chan)
                dtm = gw.getLineData(g[line], terrain_chan)
                clearance = alt - dtm
            else:
                alt = []
                dtm = []
                clearance = gw.getLineData(g[line], clearance_chan)
            minActualClearance = np.min(clearance)
            # print(f'Line {line}: Max deviation from {nominalClearance:.0f} m clearance = {maxDeviation:.0f} m.')

            if minActualClearance < minimumAllowedClearance:
                num_failed_lines += 1
                report += f'\nClearance too low at {minActualClearance:.0f} m on line {lineName}'
                x = gw.getLineData(g[line], xChannel)
                y = gw.getLineData(g[line], yChannel)
                distance = util._length(x, y)
                if plot_flag:
                    wpl._plotcheckSafeClearance(projName, line, distance, clearance_chan, altitude_chan, terrain_chan, alt, dtm, clearance)
        print(f'Number of failed lines = {num_failed_lines}.')
        print(report)
        if num_failed_lines > 0 and plot_flag:
            plt.show()



def checkClearance(whizzFile, nominalClearance, clearance_chan='', altitude_chan='', terrain_chan='', allowance=20.0, maxDistance=1000.0, xChannel='', yChannel=''):
    """
    Reports exceedances from contract specifications of the height above terrain 
    for an airborne survey Whizz database.

    The specification requires heights to be within a relative range (allowance) about a nominal
    value (nominalClearance) over a particular distance (maxDistance).
    
    The clearance (calculated as altitude - terrain if clearance_chan=''), must be
    within +/- allowance of nominalClearance.
    For each line, the maximum absolute deviation of clearance from nominalClearance is
    reported. 
    If any samples along the line exceed the allowance, then profiles are plotted
    for visual analysis. No use is made of maxDistance (yet).
    The positions (xChannel and yChannel) are only used for plotting.

    See also `checkSafeClearance`, `checkDrape`, `checkVertPlan`.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    clearance_chan : String, optional
        The name of the terrain clearance field in the Lines group
        of the Whizz HDF5 file.
    altitude_chan : String, optional
        The name of the absolute altitude or height field in the Lines group
        of the Whizz HDF5 file.
    terrain_chan : String, optional
        The name of the absolute terrain or DTM height field in the Lines group
        of the Whizz HDF5 file.
    nominalClearance : Float
        The planned clearance between above the terrain in metres.
    allowance : Float, optional
        The absolute maximum deviation allowed from the planned clearance in
        metres. The default is 20.0.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.

    Returns
    -------
    None.

    """
    
    filename = str(whizzFile)
    if clearance_chan == '':
        if altitude_chan == '' or terrain_chan == '':
            print('ERROR - either clearance_chan, or both altitude_chan and terrain_chan must be specified.')

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        num_failed_lines = 0
        report = ''

        for line in g.keys():
            lineName = _get_lineName(g[line])
            if clearance_chan == '':
                alt = gw.getLineData(g[line], altitude_chan)
                dtm = gw.getLineData(g[line], terrain_chan)
                clearance = alt - dtm
            else:
                clearance = gw.getLineData(g[line], clearance_chan)
            deviation = nominalClearance - clearance
            maxDeviation = np.max(abs(deviation))
            # print(f'Line {line}: Max deviation from {nominalClearance:.0f} m clearance = {maxDeviation:.0f} m.')
            if maxDeviation > allowance:
                num_failed_lines += 1
                report += f'\nClearance deviation of {maxDeviation:.0f} m on line {lineName}'
                x = gw.getLineData(g[line], xChannel)
                y = gw.getLineData(g[line], yChannel)
                distance = util._length(x, y)
                fig = plt.figure()

                ax = fig.add_subplot(2,1,1)
                thou_format = tkr.FuncFormatter(util._space_thou)
                if clearance_chan == '':
                    ax.plot(distance, alt)
                    ax.plot(distance, dtm)
                    ax.xaxis.set_major_formatter(thou_format)
                    plt.legend([altitude_chan, terrain_chan], fontsize=8)
                else:
                    ax.plot(distance, clearance)
                    ax.xaxis.set_major_formatter(thou_format)
                    plt.legend([clearance_chan], fontsize=8)
                plt.xlabel('distance along line [m]', fontsize = 6)
                plt.ylabel('height', fontsize = 6)
                plotTitle = projName + ': Absolute Clearance Check ' + line
                plt.title(plotTitle, fontsize = 8)
                plt.grid(True)
                for label in ax.get_xticklabels(): label.set_fontsize(6)
                for label in ax.get_yticklabels(): label.set_fontsize(6)
                ax2 = fig.add_subplot(2,1,2)
                ax2.plot(distance, deviation)
                ax2.plot(distance, allowance * np.ones(y.shape), 'r')
                ax2.plot(distance, -allowance * np.ones(y.shape), 'r')
                ax2.xaxis.set_major_formatter(thou_format)
                plt.ylabel('deviation', fontsize = 6)
                plt.xlabel('distance along line [m]', fontsize = 6)
                plt.grid(True)
                for label in ax2.get_xticklabels(): label.set_fontsize(6)
                for label in ax2.get_yticklabels(): label.set_fontsize(6)
        print(f'Number of failed lines = {num_failed_lines}.')
        print(report)
        if num_failed_lines > 0:
            plt.show()

        
def checkDrape(whizzFile, altitude, drapeChannel, warningClearance = 20.0, xChannel = '', yChannel = '', plot_flag=True):
    """
    Checks actual altitude flown against planned drape in drape channel.

    See also `checkClearance`, `checkSafeClearance`, `checkVertPlan`.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    altitude : String
        The name of the geoWhizz field or channel containing the measured altitudes.
    drapeChannel : String
        The name of the geoWhizz field or channel containing the measured drape heights.
    warningClearance : Float, optional
        DESCRIPTION. The default is 20.0.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    plot_flag : Bool, optional
        If True, plot a map of each pair of intersecting lines where the `zChannel`
        values differ by more than `max_allowed_deltaZ`. Default False.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        
        for line in g.keys():
            lineName = _get_lineName(g[line])
            alt = g[line][altitude]
            drape = g[line][drapeChannel]
            clearance = np.abs(alt[()] - drape[()])
            maxDeviation = np.max(clearance)
            if maxDeviation > warningClearance:
                print(f'Exceedance: maximum drape deviation of {maxDeviation:.1f} on line {lineName}')
                if plot_flag:
                    x = g[line][xChannel]
                    y = g[line][yChannel]
                    distance = util._length(x, y)
                    fig = plt.figure()
                    ax = fig.add_subplot(2,1,1)
                    thou_format = tkr.FuncFormatter(util._space_thou)
                    ax.plot(distance, alt[()])
                    ax.plot(distance, drape[()])
                    ax.xaxis.set_major_formatter(thou_format)
                    plt.xlabel('distance along line [m]', fontsize = 6)
                    plt.ylabel('height', fontsize = 6)
                    plotTitle = projName + ': Clearance Check ' + lineName
                    plt.title(plotTitle, fontsize = 8)
                    plt.legend([altitude, drapeChannel])
                    plt.grid(True)
                    for label in ax.get_xticklabels(): label.set_fontsize(6)
                    for label in ax.get_yticklabels(): label.set_fontsize(6)
                    ax2 = fig.add_subplot(2,1,2)
                    ax2.plot(distance, clearance)
                    ax2.xaxis.set_major_formatter(thou_format)
                    plt.ylabel('deviation', fontsize = 6)
                    plt.xlabel('distance along line [m]', fontsize = 6)
                    plt.grid(True)
                    for label in ax2.get_xticklabels(): label.set_fontsize(6)
                    for label in ax2.get_yticklabels(): label.set_fontsize(6)
                    plt.show()
        

def checkIntersectionZ(whizzFile, controls=[], travs=[], xChannel='', yChannel='', zChannel='', max_allowed_deltaZ=10.0, plot_flag=False):
    """
    Checks that the values of the data in `zChannel` at the intersection of traverse and control
    lines are different by no more than the maximum allowed delta.
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    controls : [String]
        A list of control flightlines, e.g. ['1000110.0', '1000210.0', '1000310.0'].
    travs : [String], optional
        A list of traverse flightlines, e.g. ['1000110.0', '1000210.0', '1000310.0']. Defaults
        to all flight lines in the `whizzFile`.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    zChannel : String, optional
        The name of the geoWhizz field or channel containing the data to be tested. The
        default is to read the AltitudeChannel field name from the Coordinate Frame.
    max_allowed_deltaZ : Float, optional
        The maximum allowed difference in `zChannel` between the traverse and control lines
        at each intersection point. Defaults to 10.0.
    plot_flag : Bool, optional
        If True, plot a map of each pair of intersecting lines where the `zChannel`
        values differ by more than `max_allowed_deltaZ`. Default False.

    Returns
    -------
    None

    """
    data_is_good = True
    report = ''
    filename = str(whizzFile)
    num_intersections_checked = 0
    num_failed_intersections = 0
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if zChannel == '':
            zChannel = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        all_lines = list(g.keys())
        alltravs, allcontrols = controls_lessthan_1000(all_lines)
        if controls == []:
            controls = allcontrols
        if travs == []:
            lines = alltravs
        else:
            lines = travs
        #if the channel has an attribute 'Units'
        dd = g[lines[0]][zChannel]
        z_units = ''
        if 'Units' in dd.attrs.keys():
            z_units = dd.attrs['Units']

        for linec in controls:
            x_ctrl = np.array(g[linec][xChannel])
            y_ctrl = np.array(g[linec][yChannel])
            z_ctrl = np.array(g[linec][zChannel])

            bear_ctrl = _calc_bearing(x_ctrl, y_ctrl)
            (y_ctrl1, x_ctrl1) = _rotateCoords(x_ctrl-x_ctrl[0], y_ctrl-y_ctrl[0], -bear_ctrl)
            for linet in lines:
                # if linet == linec, then it is a control line, not a traverse.
                # TODO: compare the PlannedLine for linet and linec rather than the lines themselves,
                #       since we don't want to compare different segments of the same line!
                if linet == linec:
                    continue
                x_trav = np.array(g[linet][xChannel])
                y_trav = np.array(g[linet][yChannel])
                z_trav = np.array(g[linet][zChannel])
                (y_trav1, x_trav1) = _rotateCoords(x_trav-x_ctrl[0], y_trav-y_ctrl[0], -bear_ctrl)
                # if _lines_cross(x_ctrl, y_ctrl, x_trav, y_trav):
                if _intersect(x_ctrl[0], y_ctrl[0], x_ctrl[-1], y_ctrl[-1], x_trav[0], y_trav[0], x_trav[-1], y_trav[-1]):
                    # print(f'bearings: {bearingt}, {bearingt} -- {np.abs(np.cos(bearingc - bearingt))}')
                    dh = _intersection_height(x_trav1, y_trav1, z_trav, x_ctrl1, y_ctrl1, z_ctrl, bear_ctrl)
                    num_intersections_checked += 1
                    if dh > max_allowed_deltaZ:
                        num_failed_intersections += 1
                        # print(dh)
                        bc_deg = bear_ctrl * 180.0 / np.pi
                        report += f'\n  {linet} : {linec} [bearing={bc_deg:.1f}] intersection {zChannel} difference = {dh:.1f} > {max_allowed_deltaZ:.1f}'
                        if z_units != '':
                            report += ' ' + z_units
                        report += '.'
                        data_is_good = False
                        if plot_flag:
                            fig = plt.figure()
                            fig.suptitle(f'Title {linet}', fontsize=10)
                            fig.subplots_adjust(top=0.85)
                            
                            ax = fig.add_subplot(1,2,1)
                            ax.plot(x_trav, y_trav, x_ctrl, y_ctrl)
                            plt.ylabel('y_trav', fontsize = 6)
                            plt.grid(True)
                            for label in ax.get_xticklabels(): label.set_fontsize(6)
                            for label in ax.get_yticklabels(): label.set_fontsize(6)

                            ax = fig.add_subplot(1,2,2)
                            ax.plot(x_trav1, y_trav1, x_ctrl1, y_ctrl1)
                            plt.ylabel('y_ctrl', fontsize = 6)
                            plt.grid(True)
                            for label in ax.get_xticklabels(): label.set_fontsize(6)
                            for label in ax.get_yticklabels(): label.set_fontsize(6)
                            plt.show()
                # else:
                #     report += f'\n  {linet} : {linec} un-tested since not perpendicular.'
    if data_is_good:
        report += f'All {num_intersections_checked} intersection {zChannel} differences were less than {max_allowed_deltaZ:.1f}'
        if z_units != '':
            report += ' ' + z_units
        report += '.'
    else:
        tmpstr = f'Of {num_intersections_checked} intersections checked'
        tmpstr += f', {num_failed_intersections} exceeded the '
        tmpstr += f'{max_allowed_deltaZ} allowed {zChannel} difference.\n'
        report = tmpstr + report
    print(report)
    return


def controls_lessthan_1000(all_lines):
    ctrl_strs = []
    trav_strs = []
    for line in all_lines:
        if float(line) < 999.99:
            ctrl_strs.append(line)
        else:
            trav_strs.append(line)
    return trav_strs, ctrl_strs


def _ccw(x1, y1, x2, y2, x3, y3):
    return (y3-y1)*(x2-x1) > (y2-y1)*(x3-x1)


def _intersect(cx1, cy1, cx2, cy2, tx1, ty1, tx2, ty2):
    return _ccw(cx1, cy1, tx1, ty1, tx2, ty2) != _ccw(cx2, cy2, tx1, ty1, tx2, ty2) and _ccw(cx1, cy1, cx2, cy2, tx1, ty1) != _ccw(cx1, cy1,cx2, cy2, tx2, ty2)


def _intersection_height(x_trav, y_trav, z_trav, x_ctrl, y_ctrl, z_ctrl, bearingc):
    """
    """

    y = np.abs(y_trav - _mean_1std(y_ctrl))
    it_arr = np.where(y == y.min())
    it = it_arr[0]
    if len(it) > 1:
        print('\nBug in _intersection_height - this intersection height cant be calculated')
        print(it, type(it), len(it))
        return 0.0

    x = np.abs(x_ctrl - x_trav[it])
    ic_arr = np.where(x == x.min())
    ic = ic_arr[0]

    # print(it, z_trav[it], ic, z_ctrl[ic])
    return np.abs(z_trav[it] - z_ctrl[ic])[0]


def _mean_1std(x):
    """
    Calculate the mean of the values in x that fall in the range of +/- stdev(x).
    This is a simplistic "mean excluding outliers".
    """
    mean1 = np.mean(x)
    std1 = np.std(x)
    idx = np.argwhere(np.logical_and(x > (mean1 - std1), x < (mean1 + std1)))
    return np.mean(x[idx])


def _lines_cross(x_ctrl, y_ctrl, x_trav, y_trav):
    """
    Assumes x,y coordinates rotated so that y_ctrl is approximately constant (y-axis
    approximately parallel to the control line).
    """

    # The traverse is 'north' or 'south' of the control line
    if y_trav.min() > y_ctrl.max() or y_trav.max() < y_ctrl.min():
        return False

    # We don't allow shallow angle crossovers (and won't count lines that were supposed to be parallel!).
    min_cosine = 0.1
    bear_ctrl = _calc_bearing(x_ctrl, y_ctrl)
    bear_trav = _calc_bearing(x_trav, y_trav)
    if np.abs(np.cos(bear_ctrl - bear_trav)) < min_cosine:
        return False

    # Find the index to the traverse sample whose y coordinate is closest to the mean control y coordinate
    y = np.abs(y_trav - _mean_1std(y_ctrl))
    it_arr = np.where(y == y.min())
    if np.size(it_arr):
        itc = it_arr[0]
    else:
        return False

    # Now can check if the traverse is 'east' or 'west' of the control line.
    if x_ctrl.min() > x_trav[itc] or x_ctrl.max() < x_trav[itc]:
        return False

    return True


def _calc_bearing(x, y):
    """
    arctan(mean(diff(x) / mean(diff(y))))
    """
    return np.arctan2(_mean_1std(np.diff(x)), _mean_1std(np.diff(y)))


def checkSpeeds(whizzFile, xChannel='', yChannel='', tChannel='', vel_north='', vel_east='', nominalSpeed=60.0, 
    maxDuration=0.0, maxDistance=0.0, allowance=0.1, allowed_range=[], minSafeSpeed=42.0, known='', plot_flag=False):
    """
    Checks the data from Whizz HDF5 file for speed exceedances against a specification
    requiring ground speeds to be within a relative range (allowance) about a nominal
    value (nominalSpeed) over a particular distance (maxDistance).
    
    The positions (xChannel and yChannel) are assumed to be sampled uniformly in time
    and the first two time (tChannel) values for the first line in the file are 
    differenced to obtain the sampling interval. From this, the number of samples, N, 
    for the aircraft to travel maxDistance is calculated. The algorithm compares the
    actual distance flown in N samples with maxDistance. If the actual distance is
    greater than

    *maxDistance + allowance * maxDistance*
    
    or less than 

    *maxDistance - allowance * maxDistance*
    
    an error is printed.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    tChannel : String, optional
        The name of the geoWhizz field or channel containing the measured times. The
        default is to read the tChannel field name from the Coordinate Frame.
    nominalSpeed : Float, optional
        The specified ground speed in m/s. The default is 60.0.
    maxDuration : Float, optional
        The time in seconds over which the speed estimate is determined.
        The default is 13.3.
    maxDistance : Float, optional
        The distance in metres over which the speed estimate is determined.
        The default is 1000.0.
    allowance : Float, optional
        The magnitude, relative to the nominalSpeed, of the range of allowed
        speeds. The default is to use `allowed_range` or, failing that, 0.1 (i.e. +/- 10% of nominal).
    allowed_range : [Float], optional
        The minimum and maximum allowed speeds as an array, `[min_allowed, max_allowed]`. The
        default is to use `allowance`.
    minSafeSpeed : Float, optional
        The minimum allowed instantaneous safe speed in m/s. The default is 42.0.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line. Default False.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    if maxDistance > 1.0:
        check_by_time = False
    elif maxDuration > 1.0:
        check_by_time = True
    else:
        print(f'ERROR. Require maxDistance > 1.0 m or maxDuration > 1.0')
        print(f'  maxDistance = {maxDistance}, maxDuration = {maxDuration}.')
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if tChannel == '':
            tChannel = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
        title_str = ''
        project = f[groupName].attrs['ProjectName']
        block = f[groupName].attrs['BlockID']
        if project != '':
            title_str += f'{project}'
        if block != '':
            title_str += f' {block}'
        _reportSpeeds(g, maxDuration=maxDuration, maxDistance=maxDistance, xChannel=xChannel, yChannel=yChannel,
                tChannel=tChannel, vel_north=vel_north, vel_east=vel_east, nominalSpeed=nominalSpeed, 
                     allowance=allowance, allowed_range=allowed_range, minSafeSpeed=minSafeSpeed, title_str=title_str, known=known,
                     plot_flag=plot_flag)
        

def _reportSpeeds(group, maxDuration=0.0, maxDistance=0.0, xChannel='X', yChannel='Y',
                tChannel='time', vel_north='', vel_east='', nominalSpeed=60.0,
                allowance=0.1, allowed_range=[], minSafeSpeed=42.0, title_str='', known='', plot_flag=False):
    """
    Checks the data from Whizz Line group for speed exceedances against a specification
    requiring ground speeds to be within a relative range (`allowance`) about a nominal
    value (`nominalSpeed`) over a specified distance (`maxDistance`). If `maxDistance` is
    not provided, and `maxDuration` is provided, then the check is instead over the
    specified duration.
    
    If `vel_north` and `vel_east` are both provided, then the data in those channels
    is used to calculate speed along the line. If not, then speeds are calculated
    by differencing the positions (`xChannel` and `yChannel`) and dividing by the
    difference of the first two time (`tChannel`) values for the first line in the
    group. This assumes that data are sampled uniformly in time.
    
    From this, the number of samples, N, 
    for the aircraft to travel maxDistance is calculated. The algorithm compares the
    actual distance flown in N samples with maxDistance. If the actual distance is
    greater than
        maxDistance + allowance * maxDistance
    or less than 
        maxDistance - allowance * maxDistance
    an error is printed.

    Parameters
    ----------
    group : HDF5 Group
        The Whizz line group containing the survey line data.
    maxDuration : Float, optional
        The time in seconds over which the speed estimate is determined.
        The default is 0.0.
    maxDistance : Float, optional
        The distance in metres over which the speed estimate is determined.
        The default is 0.0.
    xChannel : String, optional
        The field in the line containing X positions. The default is 'X'.
    yChannel : String, optional
        The field in the line containing Y positions. The default is 'Y'.
    tChannel : String, optional
        The field in the line containing sample times. The default is 'time'.
    vel_north : String, optional
        The field in the line containing velocity north. The default is ''.
    vel_east : String, optional
        The field in the line containing velocity east. The default is ''.
    nominalSpeed : Float, optional
        The specified ground speed in m/s. The default is 60.0.
    allowance : Float, optional
        The magnitude, relative to the nominalSpeed, of the range of allowed
        speeds. The default is 0.1 (i.e. +/- 10% of nominal).
    allowed_range : [Float], optional
        The minimum and maximum allowed speeds as an array, `[min_allowed, max_allowed]`. The
        default is to use `allowance`.
    minSafeSpeed : Float, optional
        The minimum allowed instantaneous safe speed in m/s. The default is 42.0.
    title_str : String, optional
        A title string for the plots. Default ''.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line. Default False.

    Returns
    -------
    None.

    """
    calc_from_pos = False
    if vel_north == '' or vel_east == '':
        calc_from_pos = True
        print('Velocities not known - will calculate from positions')

    check_against_dist = False
    max_allowed_str = f'{maxDuration:.1f} s'
    if maxDuration < 1.0:
        check_against_dist = True
        max_allowed_str = f'{maxDistance:.1f} m'

    num_failed_lines = 0
    num_failures = 0
    num_exceed_lines = 0
    total_num_lines = 0
    report = ''

    lines = list(group.keys())
    total_num_lines = len(lines)

    exceedances_known = False
    this_exc_known = False
    number_exc_known = 0

    if len(allowed_range) == 2:
        min_allowance = allowed_range[0]
        max_allowance = allowed_range[1]
    else:
        min_allowance = nominalSpeed * (1.0 - allowance)
        max_allowance = nominalSpeed * (1.0 + allowance)

    settings = f'Nominal ground speed {nominalSpeed:.1f} m/s; '
    settings += f'allowed {min_allowance:.1f} : '
    settings += f'{max_allowance:.1f} for < {max_allowed_str}.\n'
    print(settings)
    
    for line in lines:
        lineName = _get_lineName(group[line])
        if title_str == '':
            plot_title = lineName
        else:
            plot_title = f'{title_str}: {lineName}'

        x, y, dist, t, speed = _get_data(group[line], xChannel, yChannel, tChannel, vel_north, vel_east)

        if speed[speed < minSafeSpeed].size > 0:
            lineUnsafeSlow = True
            print(f' For at least one reading in L{lineName}, the ground speed was < {minSafeSpeed} (might be unsafe).')
        
        speed_extreme = 0.0
        num_fids_in_exceedance = 0
        exceedance_in_line = False
        too_slow = False
        line_fails = False

        if known != '':
            exceedances_known = True
            exc_known = np.array(group[line][known])
            report_known = -1
        
        for fid in range(0, len(x)):
            # There is a speed exceedance ...
            if speed[fid] > max_allowance or speed[fid] < min_allowance:
                # If a new exceedance, then initialise variables;
                if num_fids_in_exceedance == 0:
                    num_exceed_lines += 1
                    if speed[fid] < min_allowance:
                        too_slow = True
                    else:
                        too_slow = False
                    start_x = x[fid]
                    start_y = y[fid]
                    start_t = t[fid]
                    num_fids_in_exceedance = 1
                    speed_extreme = speed[fid]
                    this_exc_known = False
                # Else increment and update on the current exceedance.
                else:
                    # check we haven't swapped from too slow to too fast or vice versa
                    if (too_slow and speed[fid] > max_allowance) or ((not too_slow) and speed[fid] < min_allowance):
                        print(f'WARNING: Exceedance reversed speed in one fid. NOT BELIEVABLE. At time={t[fid]:.3f}')
                    num_fids_in_exceedance += 1
                    if exceedances_known:
                        if exc_known[fid] > 0:
                            report_known = exc_known[fid]
                            this_exc_known = True
                    if too_slow:
                        speed_extreme = min(speed[fid], speed_extreme)
                    else:
                        speed_extreme = max(speed[fid], speed_extreme)
            else:
                if num_fids_in_exceedance > 0: # the current exceedance has ended
                    end_x = x[fid]
                    end_y = y[fid]
                    end_t = t[fid]
                    dist_exceedance = util._length([start_x, end_x], [start_y, end_y])[1]
                    durn_exceedance = end_t - start_t
                    if too_slow:
                        speed_msg = "too slow"
                    else:
                        speed_msg = "too fast"
                    if _exceedance_fail(durn_exceedance, dist_exceedance, maxDuration, maxDistance):
                        report += f'\nL {lineName} {speed_msg} for {durn_exceedance} sec '
                        report += f'({dist_exceedance:.0f} m), peak exceedance = {speed_extreme:.0f} m/s.'
                        report += f'\n  From ({start_x:.0f} E {start_y:.0f} N) to ({end_x:.0f} E {end_y:.0f} N).'
                        exceedance_in_line = True
                        num_failures += 1
                        if this_exc_known:
                            number_exc_known += 1
                            report += f' Known exceedance {report_known}.'
                            this_exc_known = False
                    num_fids_in_exceedance = 0
                    too_slow = False
        if exceedance_in_line:
            num_failed_lines += 1
            if plot_flag:
                if check_against_dist:
                    wpl._plot_speed(dist, 'Distance from start of line [m]', speed, min_allowance, max_allowance, plot_title=plot_title)
                else:
                    wpl._plot_speed(t-t[0], 'Time from start of line [sec]', speed, min_allowance, max_allowance, plot_title=plot_title)

    print(f' Checked {total_num_lines} lines and {num_exceed_lines} had some short exceedance(s).')
    print(f' {num_failed_lines} lines failed for exceedance > allowed.')
    print(f' Total number of full exceedances = {num_failures}.')
    print(f'\n{number_exc_known} exceedances known in the database.\n')
    print(report)
    if plot_flag:
        plt.show()


def _get_data(line_group, xChannel, yChannel, tChannel, vel_north, vel_east):
    """
    """
    x = np.array(line_group[xChannel])
    y = np.array(line_group[yChannel])
    t = np.array(line_group[tChannel])
    distance = util._length(x, y)
    if vel_north == '' or vel_east == '':
        sampleTime = t[1] - t[0]
        xVel = np.diff(x) / sampleTime
        yVel = np.diff(y) / sampleTime
        temp = np.sqrt(xVel * xVel + yVel * yVel)
        speed = np.append(temp, np.mean(temp))
    else:
        xVel = np.array(line_group[vel_east])
        yVel = np.array(line_group[vel_north])
        speed = np.sqrt(xVel * xVel + yVel * yVel)

    return x, y, distance, t, speed


def _displacement2(x0, x1, y0, y1):
    """
    The displacement distance from (xo, y0) to (x1, y1).

    Parameters
    ----------
    x0 : TYPE
        DESCRIPTION.
    x1 : TYPE
        DESCRIPTION.
    y0 : TYPE
        DESCRIPTION.
    y1 : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return np.sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1))
            

