#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the free-air correction calculation.
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
                    

def checkFreeAirCorr(whizzFile, faCorr, latitude='', GRS80_height='', lines=[], changesign=False):
    """
    Subtracts the free-air correction in the data file from one calculated using
    Hinze et al (2005) and the latitude and height data in the data file.
    The min, max, mean and standard deviation of the difference is calculated
    for each line and presented in an extended whisker plot.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    faCorr : String
        The name of the geoWhizz field or channel containing the free-air correction.
    latitude : String, optional
        The name of the geoWhizz field or channel containing the latitude. The
        default is to read the latitude field name from the Coordinate Frame.
    GRS80_height : String, optional
        The name of the geoWhizz field or channel containing the GRS80_height. The
        default is to read the altitudeChannel field name from the Coordinate Frame.
    lines : Array{String}, optional
        Array of line numbers as strings. Default = [], meaning all lines are checked.
    changesign : Bool, optional
        If True, adds the calculated result to the data channel, instead of subtracting.
        Default False.

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
        if GRS80_height == '':
            GRS80_height = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']

        flightLine = list(g.keys())[0]
        corr_units = g[flightLine][faCorr].attrs['Units']
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
            ht_data = rd.getLineData(g[line], GRS80_height)
            cor_data = rd.getLineData(g[line], faCorr)
           
            cal_data = _freeAirCorrection(ht_data, lat_data)
            if changesign:
                err_data = cor_data * unit_scale + cal_data
            else:
                err_data = cor_data * unit_scale - cal_data 
            diffMin[count] = np.min(err_data)
            diffMax[count] = np.max(err_data)
            diffMean[count] = np.mean(err_data)
            diffStd[count] = np.std(err_data)
            lineNo[count] = line
            count += 1
            
            if line == '8334.0':
                fig = plt.figure()
                fig.suptitle(f'{projName}', fontsize=10)
                fig.subplots_adjust(top=0.85)
                ax = fig.add_subplot(3,1,1)
                ax.plot(cor_data*1e6)
                ax = fig.add_subplot(3,1,2)
                ax.plot(cal_data)
                ax = fig.add_subplot(3,1,3)
                ax.plot(err_data)
        
        figtitle = wpl.make_plot_title(f[groupName]) + ' Free Air Correction Check'
        titlestr = faCorr + ' Stats'
        xlabelstr = 'Line number'
        ylabelstr = 'Free Air Correction Error [um/s/s]'
        plotBoxWhisker(diffMin, diffMax, diffMean, diffStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)
    return


def _freeAirCorrection(height, latitude):
    """
    Returns the free-air gravity correction in um/s/s of Hinze et al (2005), eqn 3.
    
    Parameters
    ----------
    height (Float) : the height in metres above the GRS80 ellipsoid; optional, default 0.0.

    latitude (Float) : latitude in degrees (N pos, S neg).

    Returns
    -------
    free-air gravity correction in um/s/s.

    """
    
    sinLatSquared = np.sin(latitude * np.pi / 180.0)
    sinLatSquared = sinLatSquared * sinLatSquared
    freeAir = -(0.3087691 - 0.0004398 * sinLatSquared) * height
    freeAir += 7.2125E-8 * height * height
    # Hinze et al work in mGal but we want um/s/s so x 10
    return 10.0 * freeAir    
