#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot power spectral density of the difference between two channels.
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.fft import rfft, rfftfreq
from scipy.signal.windows import hann

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.utility.utility as util
import pegasusQC.whizzFiles.reportData as rp

groupName = config.groupName


def psdChannelDiff(whizzFile, channel1, channel2, flightLines=[]):
    """
    Plot the PSD (log-log Sqrt(Power) from welch method) of
    channel1 - channel2 in each flightLine. 

    Parameters
    ----------
    whizzFile : String or pathlib Path
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
            data = rd.getLineData(g[line], channel1) - rd.getLineData(g[line], channel2)#

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


def psdChannel(whizzFile, channel, flightLines=[], shortestPeriod=0.0, minlinelenkm=None, verbose=False):
    """
    Plot the PSD of channel averaged over the flightLines. 
        
    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    channel : String
        The name of the channel to be analysed.
    flightLines : String List, optional
        A list of flightline, e.g. ['1000110.0']. Default is all lines in whizzFile.
    shortestPeriod : Float, optional
        The left hand limit of the x (period) axis of the plot in seconds. Default is 0.0.
    minlinelenkm : Float, optional
        Flightlines shorter than this number of kilometres will be ignored in the calculation. Default is None.
    verbose : Bool, optional
        If True, more information is printed. Default is False.

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

        myscaling = 'density'
        corr_units = rd.getChannelAttrs(g[flightLines[0]], channel)
        if myscaling == 'density':
            y_label = f'{channel} [{corr_units}^2/Hz]'
        else:
            y_label = f'{channel} [{corr_units}^2]'
        shortestN = pow(2, 30)

        # It is clearer code if we first count the number of lines and number of samples
        numlines = 0
        numsamples = pow(2, 30)
        linelist = []
        for line in flightLines:
            # ignore lines that are too short
            if not minlinelenkm is None:
                x = f[groupName]['CoordinateFrame'].attrs['XChannel']
                y = f[groupName]['CoordinateFrame'].attrs['YChannel']
                xPos = rd.getLineData(g[line], x)
                yPos = rd.getLineData(g[line], y)
                linelen = util._displacement2(xPos[0], xPos[-1], yPos[0], yPos[-1]) / 1000.0
                if linelen < minlinelenkm:
                    if verbose:
                        print(f'line {line} skipped - too short at {linelen} km.')
                    continue

            # ignore lines where raw and filtered data vectors are not the same length
            rawdata = rd.getLineData(g[line], channel)

            # count the good lines and find the shortest length
            numlines += 1
            linelist.append(line)
            if rawdata.size < numsamples:
                numsamples = rawdata.size

        print(f'Of {len(flightLines)} lines in database, {numlines} will be processed.')
        print(f'The first {numsamples} samples from each line will be used.')
        if numlines == 1:
            print('Only one flight line in analysis, so results are not reliable.')
        elif numlines < 1:
            print('WARNING - no lines found. Perhaps try with a shorter min line length.')
            return

        f_sample = _time_frequency(f[groupName])
        q = 3  # over-sampling factor
        N = numsamples
        w = hann(N)

        # initialise the vectors by analysing the first line
        firstline = linelist[0]

        rawdata = rd.getLineData(g[firstline], channel)[0:N]
        freq, Pxx = sig.welch((rawdata - np.mean(rawdata)) * w, nfft=N*q, fs = f_sample)
        period = 1.0 / freq[1:]
        sumWelch = Pxx
        if verbose:    
            print(f'Line {line}; Num samples = {N}; f_sample = {f_sample} Hz; est line length = {63 * N / f_sample} m. period shape {period.shape}')

        for line in linelist[1:]:

            rawdata = rd.getLineData(g[line], channel)[0:N]
            freq, Pxx = sig.welch((rawdata - np.mean(rawdata)) * w, nfft=N*q, fs = f_sample, scaling=myscaling)
            sumWelch = sumWelch + Pxx
            
            if verbose:    
                print(f'Line {line}; Num samples = {N}; f_sample = {f_sample} Hz; est line length = {63 * N / f_sample} m. period shape {period.shape}')

    avgWelch = sumWelch / numlines

    fig, ax = plt.subplots()

    ax.loglog(period, avgWelch[1:], color='blue', lw=0.3)
    
    if shortestPeriod > 0.0:
        ax.set_xlim(left=shortestPeriod)

    ax.yaxis.grid(True, which='major')
    ax.set_xlabel(f'Period [s] at sample rate {f_sample:.2f} Hz', fontsize = 6)
    ax.set_ylabel(y_label, fontsize = 6)
    plotTitle = f'{projName} : {y_label}'
    ax.set_title(plotTitle, fontsize = 8)
    ax.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)


def psdChannelGain(whizzFile, rawchan, filchan, flightLines=[], nominalPeriod=0.0, shortestPeriod=0.0, minlinelenkm=None, verbose=False):
    """
    Plot the FFT of filchan / rawchan averaged over the flightLines. 
        
    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    rawchan : String
        The name of the channel to be the denominator in the gain.
    filchan : String
        The name of the channel to be the numerator in the gain.
    flightLines : String List, optional
        A list of flightline, e.g. ['1000110.0']. Default is all lines in whizzFile.
    nominalPeriod : Float, optional
        At this period in seconds, a vertical red line is drawn. Default (0.0) is to not draw the line.
    shortestPeriod : Float, optional
        The left hand limit of the x (period) axis of the plot in seconds. Default is 0.0.
    minlinelenkm : Float, optional
        Flightlines shorter than this number of kilometres will be ignored in the calculation. Default is None.
    verbose : Bool, optional
        If True, more information is printed. Default is False.

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
        corr_units = g[flightLines[0]][rawchan].attrs['Units']
        if not (g[flightLines[0]][filchan].attrs['Units'] == corr_units):
            print('Error: {rawchan} and {filchan} do not have the same units.')
            return

        y_label = f'{filchan} / {rawchan}' #' [{corr_units}]'
        shortestN = pow(2, 30)

        # It is clearer code if we first count the number of lines and number of samples
        numlines = 0
        numsamples = pow(2, 30)
        linelist = []
        for line in flightLines:
            # ignore lines that are too short
            if not minlinelenkm is None:
                x = f[groupName]['CoordinateFrame'].attrs['XChannel']
                y = f[groupName]['CoordinateFrame'].attrs['YChannel']
                xPos = rd.getLineData(g[line], x)
                yPos = rd.getLineData(g[line], y)
                linelen = util._displacement2(xPos[0], xPos[-1], yPos[0], yPos[-1]) / 1000.0
                if linelen < minlinelenkm:
                    if verbose:
                        print(f'line {line} skipped - too short at {linelen} km.')
                    continue

            # ignore lines where raw and filtered data vectors are not the same length
            rawdata = rd.getLineData(g[line], rawchan)
            fildata = rd.getLineData(g[line], filchan)
            if rawdata.size != fildata.size:
                print('ERROR line {line} - vector lengths unmatched.')
                continue

            # count the good lines and find the shortest length
            numlines += 1
            linelist.append(line)
            if rawdata.size < numsamples:
                numsamples = rawdata.size

        print(f'Of {len(flightLines)} lines in database, {numlines} will be processed.')
        print(f'The first {numsamples} samples from each line will be used.')
        if numlines == 1:
            print('Only one flight line in analysis, so results are not reliable.')

        f_sample = _time_frequency(f[groupName])
        q = 3  # over-sampling factor
        N = numsamples
        w = hann(N)

        # initialise the vectors by analysing the first line
        firstline = linelist[0]
        rawdata = rd.getLineData(g[firstline], rawchan)[0:N]
        fildata = rd.getLineData(g[firstline], filchan)[0:N]

        Rxx = rfft((rawdata - np.mean(rawdata)) * w, n=N*q)
        freq = rfftfreq(N * q, 1.0 / f_sample)
        Fxx = rfft((fildata - np.mean(fildata)) * w, n=N*q)
        freq = rfftfreq(N * q, 1.0 / f_sample)

        if nominalPeriod > 0.0:
            gain = _costaper(freq, nominalPeriod)
            Rxx_filtered = Rxx * gain
            ratioFftcalc = np.abs((Rxx_filtered[1:] / Rxx[1:]))

        period = 1.0 / freq[1:]
        ratioFftmeas = np.abs((Fxx[1:] / Rxx[1:]))

        for line in linelist[1:]:

            rawdata = rd.getLineData(g[line], rawchan)[0:N]
            fildata = rd.getLineData(g[line], filchan)[0:N]

            Rxx = rfft((rawdata - np.mean(rawdata)) * w, n=N*q)
            freq = rfftfreq(N * q, 1.0 / f_sample)
            Fxx = rfft((fildata - np.mean(fildata)) * w, n=N*q)
            freq = rfftfreq(N * q, 1.0 / f_sample)

            if nominalPeriod > 0.0:
                gain = _costaper(freq, nominalPeriod)
                Rxx_filtered = Rxx * gain
                ratioFftcalc = ratioFftcalc + np.abs((Rxx_filtered[1:] / Rxx[1:]))

            ratioFftmeas = ratioFftmeas + np.abs((Fxx[1:] / Rxx[1:]))
            
            if verbose:    
                print(f'Line {line}; Num samples = {N}; f_sample = {f_sample} Hz; est line length = {63 * N / f_sample} m. period shape {period.shape}')

    ratioFftmeas = ratioFftmeas / numlines
    if nominalPeriod > 0.0:
        ratioFftcalc = ratioFftcalc / numlines
        transfer = gain[1:]

    fig, ax = plt.subplots()

    ax.semilogx(period, ratioFftmeas, color='blue', lw=0.3)
    if nominalPeriod > 0.0:
        ax.semilogx(period, ratioFftcalc, color='green', lw=0.5)
        ax.semilogx(period, transfer, color='black', lw=0.5)
        ax.semilogx(nominalPeriod, 0.5, '+', color="red")
    
    ax.set_ylim([1E-2,1])
    if shortestPeriod > 0.0:
        ax.set_xlim(left=shortestPeriod)

    ax.legend(["data gain", "check gain", "filter"], loc="lower left")
    ax.set_yticks([0.5], minor=True)
    ax.yaxis.grid(True, which='major')
    ax.yaxis.grid(True, which='minor', color='r')
    ax.set_xlabel(f'Period [s] at sample rate {f_sample} Hz', fontsize = 6)
    ax.set_ylabel(y_label, fontsize = 6)
    plotTitle = f'{projName} : {y_label}'
    ax.set_title(plotTitle, fontsize = 8)
    ax.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)


def _costaper(freq, nominalPeriod):
    """
    Return the gain vector corresponding to the freq vector
    for a cosine taper low-pass filter with centre frequency
    corresponding to `nominalperiod` and pass frequency 2/3rd
    the 50% pass frequency; cut frequency 4/3rds.
    """

    # Apply a cosine taper filter
    passfreq = 2.0 / (3.0 * nominalPeriod)
    cutfreq = 2.0 * passfreq
    numfreqs = len(freq)

    gain = np.ones(numfreqs, dtype=float)
    for i in range(numfreqs):
        gain[i] = 0.5 - 0.5 * np.cos(np.pi * (freq[i] - cutfreq) / (passfreq - cutfreq))
    gain[freq > cutfreq] = 0.0
    gain[freq < passfreq] = 1.0

    return gain


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
