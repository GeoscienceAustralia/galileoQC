#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the calculation of the Eotvos correction.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import galileoQC.config as config
import galileoQC.whizzFiles.retrieveData as rd
import galileoQC.gridFiles.read_ers as grd
from galileoQC.whizzPlots.plotBoxWhisker import plotBoxWhisker
import galileoQC.whizzPlots.whizzPlot as wpl
import galileoQC.utility.utility as util

groupName = config.groupName
    

def checkEotvosCorr(whizzFile, eotCorr, latitude='', x='', y='', GRS80_height='', time='', east_vel='', north_vel='', lines=[], changesign=False, plot_flag=False):
    """
    Subtracts the eotvos correction in the data file from one calculated using
    Jekeli (2016, slide 51) and the latitude and position data in the data file.
    The min, max, mean and standard deviation of the difference is calculated
    for each line and presented in an extended whisker plot.

    Parameters
    ----------
    whizzFile : String or pathlib Path

        Name of a HDF5 Whizz file, including path and extension.

    eotCorr : String

        The name of the geoWhizz field or channel containing the eotvos correction.

    latitude : String, optional

        The name of the geoWhizz field or channel containing the latitude. The
        default is to read the latitude field name from the Coordinate Frame.

    x : String, optional

        The name of the geoWhizz field or channel containing the x position. The
        default is to read the xChannel field name from the Coordinate Frame.

    y : String, optional

        The name of the geoWhizz field or channel containing the y position. The
        default is to read the yChannel field name from the Coordinate Frame.

    GRS80_height : String, optional

        The name of the geoWhizz field or channel containing the GRS80_height. The
        default is to read the altitudeChannel field name from the Coordinate Frame.

    time : String, optional

        The name of the geoWhizz field or channel containing the time. The
        default is to read the timeChannel field name from the Coordinate Frame.

    east_vel : String, optional

        The name of the geoWhizz field or channel containing the velocity in the
        east direction. The default ('') is to calculate this from the x and time
        channels.

    north_vel : String, optional

        The name of the geoWhizz field or channel containing the velocity in the
        north direction. The default ('') is to calculate this from the y and time
        channels.

    lines : Array{String}, optional

        Array of line numbers as strings. Default = [], meaning all lines are checked.

    changesign : Bool, optional

        If True, adds the calculated result to the data channel, instead of subtracting.
        Default False.

    plot_flag : Bool, optional

        If true, plot details where the difference is large. Summary plot is always done. Default False.

    Returns
    -------
    None.

    """
    if (east_vel == '')  | (north_vel == ''):
        mywarning = 'WARNING - no velocity channel names provided. Calculation of the \n'
        mywarning += 'Eotvos effect will proceed using velocities estimated from \n'
        mywarning += 'positions. This ignores aircraft heading changes due to cross \n'
        mywarning += 'winds and thus can be in error by up to about 10,000 um/s/s).\n'
        print(mywarning)

    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if latitude == '':
            latitude = f[groupName]['CoordinateFrame'].attrs['LatitudeChannel']
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if GRS80_height == '':
            GRS80_height = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        if time == '':
            time = f[groupName]['CoordinateFrame'].attrs['TimeChannel']            

        flightLine = list(g.keys())[0]
        corr_units = rd.getChannelAttrs(g[flightLine], eotCorr, myattribute='Units')
        if corr_units == 'mGal':
            unit_scale = 10.0
        elif corr_units in ['gu', 'µm/s/s', 'um/s/s', 'um/s2']:
            unit_scale = 1.0
        else:
            print(f'ERROR - correction units {corr_units} not recognised, expected mGal or µm/s/s or um/s/s or gu')
            return

        if lines == []:
            lines = list(g.keys())

        numLines = len(lines)
        diffMin = np.zeros((numLines,))
        diffMax = np.zeros((numLines,))
        diffMean = np.zeros((numLines,))
        diffStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in lines:
            lat_data = rd.getLineData(g[line], latitude)
            x_data = rd.getLineData(g[line], x)
            y_data = rd.getLineData(g[line], y)
            ht_data = rd.getLineData(g[line], GRS80_height)
            time_data = rd.getLineData(g[line], time)
            cor_data = rd.getLineData(g[line], eotCorr)
            if (east_vel == '')  | (north_vel == ''):
                (e_speed, n_speed) = _calc_speed(x_data, y_data, time_data)
            else:
                n_speed = rd.getLineData(g[line], north_vel)
                e_speed = rd.getLineData(g[line], east_vel)
            cal_data = _eotvosCorrection(e_speed, n_speed, lat_data, ht_data)
            if changesign:
                err_data = cor_data * unit_scale + cal_data
            else:
                err_data = cor_data * unit_scale - cal_data 
            diffMin[count] = np.min(err_data)
            diffMax[count] = np.max(err_data)
            diffMean[count] = np.mean(err_data)
            diffStd[count] = np.std(err_data)
            lineNo[count] = line
            
            if plot_flag > 10.0 and np.abs(diffMin[count] - diffMax[count]):
                fig = plt.figure()
                fig.suptitle(f'{projName} L{lineNo[count]}', fontsize=10)
                fig.subplots_adjust(top=0.85)
                ax = fig.add_subplot(2,1,1)
                ax.plot(x_data, cor_data*unit_scale, x_data, -cal_data)
                ax = fig.add_subplot(2,1,2)
                ax.plot(x_data, err_data)
                
            count += 1
        
        figtitle = wpl.make_plot_title(f[groupName]) + ' Eotvos Correction Check'
        titlestr = eotCorr + ' Stats'
        xlabelstr = 'Line number'
        ylabelstr = 'Eotvos Correction Error [um/s/s]'
        plotBoxWhisker(diffMin, diffMax, diffMean, diffStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)
    return


def _calc_speed(e, n, t):
    """
    Returns the (e_speed, n_speed) velocity given the (e, n) positions as
    a function of t (time). Uses a simple numpy.diff approach.
 
    Parameters
    ----------
    e : Numpy 1D array

        The easting data.

    n : Numpy 1D array

        The northing data.

    t : Numpy 1D array

        The time data.

    Returns
    -------
    e_speed : Numpy 1D array

        The velocity in the east direction.

    n_speed : Numpy 1D array

        The velocity in the north direction.

    """
    n_speed = np.diff(n) / np.diff(t)
    last = n_speed[-1]
    n_speed = np.append(n_speed, last)
    e_speed = np.diff(e) / np.diff(t)
    last = e_speed[-1]
    e_speed = np.append(e_speed, last)
    return (e_speed, n_speed)
    

def _eotvosCorrection(eSpeed, nSpeed, latitude, height=0):
    """
    Calculates the Eotvos correction for moving-base gravimetry. Uses the exact
    equation (see for example Jekeli (2016), slide 51)

    Parameters
    ----------
    eSpeed : float

        the aircraft speed in the East direction in m/s/s.
    
    nSpeed : float

        the aircraft speed in the North direction in m/s/s.
    
    latitude : float

        latitude in degrees (N pos, S neg).
    
    height : float

        the height in metres above the GRS80 ellipsoid; optional, default 0.0.

    Returns
    -------
    eotvos : float

        the eotvos correction for gravity in um/s/s.

    """
    
    radius = 6378137.0
    ellipticity = 0.0818191908426
    ellSquared = ellipticity * ellipticity
    angularVelocity = 7.2921158553E-5
    cosLat = np.cos(latitude * np.pi / 180)
    sinLat = np.sin(latitude * np.pi / 180)
    ellSinLatSquared = ellSquared * sinLat * sinLat
    
    angTerm = -2.0 * angularVelocity * eSpeed * cosLat
    eastTerm = - (eSpeed * eSpeed) / (height + radius / (np.sqrt(1 - ellSinLatSquared)))
    northTerm = - (nSpeed * nSpeed) / (height + (radius * (1 - ellSquared)) / ((1 - ellSinLatSquared) * np.sqrt(1 - ellSinLatSquared)))
    
    eotvos = 1.0E6 * (angTerm + eastTerm + northTerm)
    return eotvos
