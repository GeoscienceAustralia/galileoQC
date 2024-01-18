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


