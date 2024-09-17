#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the calculation of the latitude (normal) correction.
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

    
def checkLatCorr(whizzFile, latCorr, latitude='', lines=[]):
    """
    Subtracts the latitude correction in the data file from one calculated using
    Hinze et al (2005) and the latitude data in the data file. The min, max, mean
    and standard deviation of the difference is calculated for each line and
    presented in an extended whisker plot.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    latCorr : String
        The name of the geoWhizz field or channel containing the latitude correction
        (sometimes called "normal gravity").
    latitude : String, optional
        The name of the geoWhizz field or channel containing the latitude. The
        default is to read the latitude field name from the Coordinate Frame.
    lines : Array{String}, optional
        Array of line numbers as strings. Default = [], meaning all lines are checked.

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

        flightLine = list(g.keys())[0]
        corr_units = g[flightLine][latCorr].attrs['Units']
        if corr_units == 'mGal':
            unit_scale = 10.0
        elif corr_units == 'gu' or corr_units == 'µm/s/s' or corr_units == 'um/s/s':
            unit_scale = 1.0
        else:
            print('ERROR - correction units not recognised, expected mGal or µm/s/s (gu)')
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
            cor_data = rd.getLineData(g[line], latCorr)
            if line == '8474.0':
                fig = plt.figure()
                ax = fig.add_subplot(3,1,1)
                ax.plot(cor_data * 10)
                plt.title('lat corr in data')
                ax = fig.add_subplot(3,1,2)
                ax.plot(lat_data)
                plt.title('Latitude')
                ax = fig.add_subplot(3,1,3)
                ax.plot(_normalGravity(lat_data))
                plt.title('My Estimate')
                
            err_data = cor_data * unit_scale + _normalGravity(lat_data)
            diffMin[count] = np.min(err_data)
            diffMax[count] = np.max(err_data)
            diffMean[count] = np.mean(err_data)
            diffStd[count] = np.std(err_data)
            lineNo[count] = line
            count += 1
        
        figtitle = wpl.make_plot_title(f[groupName]) + ' Latitude Correction Check'
        titlestr = latCorr + ' Stats'
        xlabelstr = 'Line number'
        ylabelstr = 'Latitude Correction Error [um/s/s]'
        plotBoxWhisker(diffMin, diffMax, diffMean, diffStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)
    return
    

def _normalGravity(latitude):
    """
    Returns the ellipsoid theoretical gravity in µm/s/s of Hinze et al (2005), eqn 2,
    based on the GRS80 ellipsoid.    

    Parameters
    ----------
    latitude (Float) : latitude in degrees (N pos, S neg).
    

    Returns
    -------
    latitude correction (Float) : ellipsoid theoretical gravity in um/s/s.

    """
    
    gNormal = 9780326.7715
    k = 0.001931851353
    eSquared = 0.0066943800229
    sinLatSquared = np.sin(latitude * np.pi / 180.0)
    sinLatSquared = sinLatSquared * sinLatSquared
    
    
    ratio = (1 + k * sinLatSquared) / np.sqrt(1 - eSquared * sinLatSquared)
    return gNormal * (ratio)# - 1.0)
