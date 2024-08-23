#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the calculation of the Eotvos correction.
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import AirGravQC.config as config
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.gridFiles.read_ers as grd
from AirGravQC.whizzPlots.plotBoxWhisker import plotBoxWhisker
import AirGravQC.whizzPlots.whizzPlot as wpl
import AirGravQC.utility.utility as util

groupName = config.groupName
    

def checkEotvosCorr(whizzFile, eotCorr, latitude='', x='', y='', GRS80_height='', time='', east_vel='', north_vel='', plot_flag=False):
    """
    Subtracts the eotvos correction in the data file from one calculated using
    Hinze et al (2005) and the latitude and position data in the data file.
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

    Returns
    -------
    None.

    """
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
        corr_units = g[flightLine][eotCorr].attrs['Units']
        if corr_units == 'mGal':
            unit_scale = 10.0
        elif corr_units == 'gu' or corr_units == 'µm/s/s':
            unit_scale = 1.0
        else:
            print('ERROR - correction units not recognised, expected mGal or µm/s/s (gu)')
            return

        numLines = len(g.items())
        diffMin = np.zeros((numLines,))
        diffMax = np.zeros((numLines,))
        diffMean = np.zeros((numLines,))
        diffStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            lat_data = rd.getLineData(g[line], latitude)
            x_data = rd.getLineData(g[line], x)
            y_data = rd.getLineData(g[line], y)
            ht_data = rd.getLineData(g[line], GRS80_height)
            time_data = rd.getLineData(g[line], time)
            cor_data = rd.getLineData(g[line], eotCorr)
            if (east_vel == '')  | (north_vel == ''):
                (n_speed, e_speed) = _calc_speed(x_data, y_data, time_data)
            else:
                n_speed = rd.getLineData(g[line], north_vel)
                e_speed = rd.getLineData(g[line], east_vel)
            cal_data = _eotvosCorrection(e_speed, n_speed, lat_data, ht_data)
            err_data = cor_data * unit_scale + cal_data
            diffMin[count] = np.min(err_data)
            diffMax[count] = np.max(err_data)
            diffMean[count] = np.mean(err_data)
            diffStd[count] = np.std(err_data)
            lineNo[count] = line
            
            if np.abs(diffMin[count] - diffMax[count]) > 10.0 and plot_flag:
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
    eSpeed (Float) : the aircraft speed in the East direction in m/s/s.
    
    nSpeed (Float) : the aircraft speed in the North direction in m/s/s.
    
    latitude (Float) : latitude in degrees (N pos, S neg).
    
    height (Float) : the height in metres above the GRS80 ellipsoid; optional, default 0.0.

    Returns
    -------
    eotvos (Float) : the eotvos correction for gravity in um/s/s.

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
