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

import AirGravQC.config as config
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.utility.utility as util

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


def psdChannelGain(whizzFile, rawchan, filchan, flightLines=[], nominalPeriod=0.0, verbose=False):
    """
    Plot the FFT of filchan / rawchan in each flightLine. 
    
    ToDo: average the spectrum over many flight lines to achieve a less noisy result.
    
    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    flightLines : String List, optional
        A list of flightline, e.g. ['1000110.0']. Default is all lines in whizzFile.
    rawchan : String
        The name of a channel.
    filchan : String
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
        corr_units = g[flightLines[0]][rawchan].attrs['Units']
        if not (g[flightLines[0]][filchan].attrs['Units'] == corr_units):
            print('Error: {rawchan} and {filchan} do not have the same units.')
            return

        y_label = f'{filchan} / {rawchan}' #' [{corr_units}]'
        fig, ax = plt.subplots()
        shortestN = pow(2, 30)
            
        for line in flightLines:
            f_sample = _time_frequency(f[groupName])
            rawdata = rd.getLineData(g[line], rawchan)
            fildata = rd.getLineData(g[line], filchan)

            q = 3  # over-sampling factor
            N = rawdata.size
            if N != fildata.size:
                print('ERROR - vector lengths unmatched.')
                continue
            w = hann(N)

            Rxx = rfft((rawdata - np.mean(rawdata)) * w, n=N*q)
            freq = rfftfreq(N * q, 1.0 / f_sample)
            Fxx = rfft((fildata - np.mean(fildata)) * w, n=N*q)
            freq = rfftfreq(N * q, 1.0 / f_sample)

            periodone = 1.0 / freq[1:]
            ratioFftone = np.abs((Fxx[1:] / Rxx[1:]))

            if shortestN < pow(2, 30):
                shortestN = min(periodone.size, shortestN)
                period = np.column_stack([period[:shortestN], periodone[:shortestN]])
                ratioFft = np.column_stack([ratioFft[:shortestN], ratioFftone[:shortestN]])
            else:
                period = periodone
                ratioFft = ratioFftone
                shortestN = min(periodone.size, shortestN)
            
            if verbose:    
                print(f'Line {line}; Num samples = {N}; f_sample = {f_sample} Hz; est line length = {63 * N / f_sample} m. period shape {period.shape}')

        perioda = np.mean(period, axis=1)
        ratioFfta = np.mean(ratioFft , axis=1)
        plt.semilogx(perioda, ratioFfta, color='blue', lw=0.3)
    
        plt.ylim([1E-2,1])

        ax.set_yticks([0.5], minor=True)
        ax.yaxis.grid(True, which='major')
        ax.yaxis.grid(True, which='minor', color='r')
        if nominalPeriod > 0.0:
            ax.set_xticks([nominalPeriod], minor=True)
            ax.xaxis.grid(True, which='minor', color='r')

        plt.xlabel(f'Period [s] at sample rate {f_sample} Hz', fontsize = 6)
        plt.ylabel(y_label, fontsize = 6)
        plotTitle = f'{projName} : {y_label}'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
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
