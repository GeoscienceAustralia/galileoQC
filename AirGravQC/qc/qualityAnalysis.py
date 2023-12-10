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
                    projection] = ers.read_ers_header(str(aFile))
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



