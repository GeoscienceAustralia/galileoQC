#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check aircraft ground clearances.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.utility.utility as util
import pegasusQC.whizzPlots.whizzPlot as wpl

groupName = config.groupName
                                       
   
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
            lineName = util._get_lineName(g[line])
            if clearance_chan == '':
                alt = rd.getLineData(g[line], altitude_chan)
                dtm = rd.getLineData(g[line], terrain_chan)
                clearance = alt - dtm
            else:
                alt = []
                dtm = []
                clearance = rd.getLineData(g[line], clearance_chan)
            minActualClearance = np.min(clearance)
            # print(f'Line {line}: Max deviation from {nominalClearance:.0f} m clearance = {maxDeviation:.0f} m.')

            if minActualClearance < minimumAllowedClearance:
                num_failed_lines += 1
                report += f'\nClearance too low at {minActualClearance:.0f} m on line {lineName}'
                x = rd.getLineData(g[line], xChannel)
                y = rd.getLineData(g[line], yChannel)
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
            lineName = util._get_lineName(g[line])
            if clearance_chan == '':
                alt = rd.getLineData(g[line], altitude_chan)
                dtm = rd.getLineData(g[line], terrain_chan)
                clearance = alt - dtm
            else:
                clearance = rd.getLineData(g[line], clearance_chan)
            deviation = nominalClearance - clearance
            maxDeviation = np.max(abs(deviation))
            # print(f'Line {line}: Max deviation from {nominalClearance:.0f} m clearance = {maxDeviation:.0f} m.')
            if maxDeviation > allowance:
                num_failed_lines += 1
                report += f'\nClearance deviation of {maxDeviation:.0f} m on line {lineName}'
                x = rd.getLineData(g[line], xChannel)
                y = rd.getLineData(g[line], yChannel)
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
            lineName = util._get_lineName(g[line])
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



