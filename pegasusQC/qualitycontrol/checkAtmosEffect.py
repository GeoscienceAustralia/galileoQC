#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check atmospheric gravity effect calculation.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.gridFiles.read_ers as grd
from pegasusQC.whizzPlots.plotBoxWhisker import plotBoxWhisker
import pegasusQC.whizzPlots.whizzPlot as wpl
import pegasusQC.utility.utility as util

groupName = config.groupName

def checkAtmosEffect(whizzFile, atmosCorr, lines=[], GRS80_height=''):
    """
    Subtracts the atomspheric correction in the data file from one calculated using
    Hinze et al (2005) and the height data in the data file.
    The min, max, mean and standard deviation of the difference is calculated
    for each line and presented in an extended whisker plot.

    Parameters
    ----------
    whizzFile : String or pathlib Path

        Name of a HDF5 Whizz file, including path and extension.

    atmosCorr : String

        The name of the geoWhizz field or channel containing the atmospheric correction. These
        data must be in units of 'gu', 'µm/s/s', 'um/s/s' or 'um/s2'.

    lines : Array{String}, optional

        Array of line numbers as strings. Default = [], meaning all lines are checked.

    GRS80_height : String, optional

        The name of the geoWhizz field or channel containing the GRS80_height. The
        default is to read the altitudeChannel field name from the Coordinate Frame.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if GRS80_height == '':
            GRS80_height = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']

        flightLine = list(g.keys())[0]
        corr_units = g[flightLine][atmosCorr].attrs['Units']
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
            ht_data = rd.getLineData(g[line], GRS80_height)
            cor_data = rd.getLineData(g[line], atmosCorr)
           
            cal_data = _atmosEffect(ht_data)
            err_data = cor_data * unit_scale - cal_data 
            diffMin[count] = np.min(err_data)
            diffMax[count] = np.max(err_data)
            diffMean[count] = np.mean(err_data)
            diffStd[count] = np.std(err_data)
            lineNo[count] = line
            count += 1
            
        figtitle = wpl.make_plot_title(f[groupName]) + ' Atmospheric Correction Check'
        titlestr = atmosCorr + ' Stats'
        xlabelstr = 'Line number'
        ylabelstr = 'Atmospheric Correction Error [um/s/s]'
        plotBoxWhisker(diffMin, diffMax, diffMean, diffStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)

    return


def _atmosEffect(height):
    """
    Returns the atmospheric gravity correction in um/s/s using Hinze et al (2005), eqn 3.

    Parameters
    ----------
    height : float

        the height in metres above the GRS80 ellipsoid; optional, default 0.0.

    Returns
    -------
    atmospheric correction (Float).

    """
    return 10.0 * (0.874 - 9.9E-5 * height + 3.56E-9 * height * height)

