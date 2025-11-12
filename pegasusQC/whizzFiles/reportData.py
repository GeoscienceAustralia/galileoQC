#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Summarise contents of a `geoWhizz` file.
Author: Mark Helm Dransfield
Created: 2023
License: CC BY-SA
"""

import numpy as np
import h5py
from pathlib import Path
import pathlib

import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.utility.utility as util
import pegasusQC.config as config
from pegasusQC.gridFiles.oddeven import (_getOddEvenLines, _getTravCtrlLines, _getPlannedLines)

groupName = config.groupName
projectName = config.projectName

def reportWhizz(whizzFile, line='', channel=''):
    """
    Prints a short summary of the data in a HDF5 Whizz file.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    line : String, optional
        The line number, formatted as a string, to report in detail. The default is '' and no line details.
    channel : String, optional
        The channel or field name to report in detail. The default is '' and no channel details.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
        
    with h5py.File(filename, 'r') as f:
        
        whizzHeader = list(f.keys())[0]
        g = f[whizzHeader]
        
        gAttributeNames = list(g.attrs)
        
        gCoord = g['CoordinateFrame']
        coordAttrs = list(gCoord.attrs)
        
        gLines = g['Lines']
        lineNames = list(gLines.keys())
        numLines = len(lineNames)
        lineGroups = list(gLines.values())
        channelNames = list(lineGroups[0].keys())
        numChannels = len(channelNames)
        
        print(whizzHeader)
        for attribute in gAttributeNames:
            print(f'    {attribute:.20}: {g.attrs[attribute]}')
        print('\nCoordinates')
        for attribute in coordAttrs:
            print(f'    {attribute}: {gCoord.attrs[attribute]}')
        if line != '':
            myLineGroup = gLines[line]
            lineAttrs = list(myLineGroup.attrs)
            print(f'\nLine {myLineGroup}')
            _ = _distanceFlown(whizzFile, lines=[line])
            for attribute in lineAttrs:
                print(f'    {attribute}: {myLineGroup.attrs[attribute]}')
        else:
            _ = _distanceFlown(whizzFile)
            print(f'\n{numLines} lines:')
        util._print_wrappedlist(lineNames)

        print(f'\n{numChannels} channels:')
        util._print_wrappedlist(channelNames)
        if channel != '':
            if line == '':
                line = lineNames[0]
                
            myChanGroup = gLines[line][channel]
            chanAttrs = list(myChanGroup.attrs)
            print(f'\nChannel {myChanGroup}')
            for attribute in chanAttrs:
                print(f'    {attribute}: {myChanGroup.attrs[attribute]}')


def reportLines(whizzFile):
    """
    Prints a short summary of the flight-lines in a HDF5 Whizz file.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    oddlines, evenlines = _getOddEvenLines(whizzFile)
    travlines, ctrllines = _getTravCtrlLines(whizzFile)
    plandlines, unplandlines = _getPlannedLines(whizzFile)

    with h5py.File(filename, 'r') as f:
        whizzHeader = list(f.keys())[0]
    print("\n", whizzHeader)
    print("\nOdd Lines:")
    util._print_wrappedlist(oddlines)
    print("\nEven Lines:")
    util._print_wrappedlist(evenlines)
    print("\nTraverse Lines:")
    util._print_wrappedlist(travlines)
    print("\nControl Lines:")
    util._print_wrappedlist(ctrllines)


def reportChannels(whizzFile, channel='', verbose=False):
    """
    Prints a short summary of the channel names in a HDF5 Whizz file.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    channel : string, optional
        If provided, only the details of that channel are given, otherwise all.
    verbose : Bool, optional
        If True, print all details, otherwise just the names. Default False.
    Returns
    -------
    None.

    """
    filename = str(whizzFile)
        
    with h5py.File(filename, 'r') as f:
        
        whizzHeader = list(f.keys())[0]
        g = f[whizzHeader]        
        gLines = g['Lines']
        lineNames = list(gLines.keys())
        numLines = len(lineNames)
        lineGroups = list(gLines.values())
        if channel == '':
            channelNames = list(lineGroups[0].keys())
        else:
            channelNames = [channel]
        numChannels = len(channelNames)
        
        print(whizzHeader)
        print(f'\n{numChannels} channels:\n')

        if verbose or channel != '':
            print(f'\033[1m  {"channel":<20} {"units":<14} {"description"}\033[0m')
            print('--------------------------------------------------')
            linegroup = gLines[lineNames[0]]
            for channel in channelNames:
                # myChanGroup = gLines[lineNames[0]][channel]
                my_units = rd.getChannelAttrs(linegroup, channel, myattribute='Units')
                my_description = rd.getChannelAttrs(linegroup, channel, myattribute='Description')
                print(f'  {channel:<20} {my_units:<14} {my_description}')
                # if my_units == '':
                #     print(f'  {channel:<20} {my_units:<14} {my_description}')
                #     ylabelstr = f'{channel}'
                # else:
                #     print(f'  {channel:<20} {myChanGroup.attrs["Units"]:<14} {myChanGroup.attrs["Description"]}')
                #     ylabelstr = f'{channel} [{my_units}]'
                # first, check Units attr exists!!!! TBD
                #             my_units = rd.getChannelAttrs(g[lines[0]], channel)

                

        else:
            print(f'\n{numChannels} channels:')
            util._print_wrappedlist(channelNames)
            # if channel != '':
            #     myChanGroup = gLines[lineNames[0]][channel]
            #     chanAttrs = list(myChanGroup.attrs)
            #     print(f'\nChannel {myChanGroup}')
            #     for attribute in chanAttrs:
            #         print(f'    {attribute}: {myChanGroup.attrs[attribute]}')


def reportFlights(whizzFile, flightChannel='', lines=[], detailed=False):
    """
    Prints a summary of the flight numbers in a HDF5 Whizz file.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    flightChannel : String, optional
        The name of the channel containing the flight numbers. The default is '' (get the channel name from attributes).
    lines : String Array, optional
        The array of line numbers, each formatted as a string, to report. The default is [] and all lines.
    detailed : Bool, optional
        If true, report the line numbers flown in each flight. The default is False.

    Returns
    -------
    None.
                gg = gLines[line]
                    gg.attrs['Flight'] = this_flight

            latitude = f[groupName]['CoordinateFrame'].attrs['LatitudeChannel']

    """
    filename = str(whizzFile)
        
    with h5py.File(filename, 'r') as f:
        
        whizzHeader = list(f.keys())[0]
        g = f[whizzHeader]
        gAttributeNames = list(g.attrs)
        gLines = g['Lines']

        # The data are stored with flights belonging to lines; we invert this relationship
        flight_dict = {}
        if lines == []:
            lines = list(gLines.keys())
        numLines = len(lines)
        for line in lines:
            if flightChannel != '':
                this_flight = gLines[line][flightChannel][0]
            else:
                if whizzAttrExists(gLines[line], 'Flight'):
                    this_flight = gLines[line].attrs['Flight']
                elif whizzAttrExists(gLines[line], 'flight'):
                    this_flight = gLines[line].attrs['flight']
                else:
                    print('No flight data found in data file so no report possible.')
                    return
            if this_flight in flight_dict:
                flight_dict[this_flight].append(line)
            else:
                flight_dict[this_flight] = [line]
        sorted_keys = sorted(flight_dict.keys())
        sorted_flights = {key:flight_dict[key] for key in sorted_keys}

        print(whizzHeader)
        for attribute in gAttributeNames:
            print(f'    {attribute:.20}: {g.attrs[attribute]}')

        print(f'\n{len(sorted_flights.keys())} flights including {numLines} lines.')

        flightreport = '\nFlights\n'
        len_strline = 0
        if detailed:
            for flight in sorted_flights:
                flightreport += f'    {flight:.0f}\n      '
                len_strline = 6
                for line in sorted(sorted_flights[flight]):
                    flightreport += f'L{line} '
                    len_strline += len(f'L{line} ')
                    if len_strline > 72:
                        flightreport += '\n      '
                        len_strline = 6
                flightreport += '\n'
        else:
            flightreport += '    '
            len_strline = 4
            for flight in sorted_flights:
                flightreport += f'{flight:.0f} '
                len_strline += len(f'{flight:.0f} ')
                if len_strline > 72:
                    flightreport += '\n    '
                    len_strline = 4
            flightreport += '\n'
        print(flightreport)


def reportSampling(whizzFile, timeChannel='', xChannel='', yChannel=''):
    """
    Prints a summary of the sample rate of the data in a HDF5 Whizz file.
    The sample rate in time are usually constant but the minimum, maximum and
    mean sample rates are reported for time and horizontal distance.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    timeChannel : String, optional
        The name of the channel containing the time data.  The default is '' which
        causes the channel stored in the CoordinateFram groups TimeChannel attribute 
        to be used.
    xChannel : String, optional
        The name of the channel containing the x data.  The default is '' which
        causes the channel stored in the CoordinateFram groups XChannel attribute 
        to be used.
    yChannel : String, optional
        The name of the channel containing the y data.  The default is '' which
        causes the channel stored in the CoordinateFram groups YChannel attribute 
        to be used.

    Returns
    -------
    None.
    """
    filename = str(whizzFile)
        
    with h5py.File(filename, 'r') as f:
        
        whizzHeader = list(f.keys())[0]
        g = f[whizzHeader]
        gAttributeNames = list(g.attrs)
        gLines = g['Lines']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if timeChannel == '':
            timeChannel = f[groupName]['CoordinateFrame'].attrs['TimeChannel']


        time_deltas = []
        x_deltas = []
        y_deltas = []
        lines = list(gLines.keys())
        numLines = len(lines)
        for line in lines:
            time_deltas = np.append(time_deltas, np.diff(rd.getLineData(gLines[line], timeChannel)))#gLines[line][timeChannel])))
            x_deltas = np.append(x_deltas, np.diff(rd.getLineData(gLines[line], xChannel))) #np.array(gLines[line][xChannel])))
            y_deltas = np.append(y_deltas, np.diff(rd.getLineData(gLines[line], yChannel))) #np.array(gLines[line][yChannel])))
        mean_dt = np.mean(time_deltas)
        min_dt = np.min(time_deltas)
        max_dt = np.max(time_deltas)
        std_dt = np.std(time_deltas)
        dd = util._distance(x_deltas, y_deltas)
        mean_dd = np.mean(dd)
        min_dd = np.min(dd)
        max_dd = np.max(dd)
        std_dd = np.std(dd)

        print(whizzHeader)
        for attribute in gAttributeNames:
            print(f'    {attribute:.20}: {g.attrs[attribute]}')

        print(f'\nSample time and distance statistics')
        print(f'  Min   = {min_dt:.3f} s, {min_dd:.1f} m')
        print(f'  Max   = {max_dt:.3f} s, {max_dd:.1f} m')
        print(f'  Mean  = {mean_dt:.3f} s, {mean_dd:.1f} m')
        print(f'  Stdev = {std_dt:.3g} s, {std_dd:.1g} m')


def _distanceFlown(whizzFile, x = '', y = '', lines=[]):
    """
    The total distance flown on the survey.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    x : String, optional
        The name of the channel of X positions (in metres). The default is '' which
        causes the channel stored in the CoordinateFram groups XChannel attribute 
        to be used.
    y : String, optional
        The name of the channel of Y positions (in metres). The default is '' which
        causes the channel stored in the CoordinateFram groups YChannel attribute 
        to be used.
    lines : array of strings, optional
        An array of line identifiers whose total distance will be returned.
        Default all lines in whizzFile.

    Returns
    -------
    count : Integer
        The total number of lines in the file.
    Float
        The total distance in km flown over all lines in the survey.

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        
        if x == '':
            if whizzAttrExists(f[groupName]['CoordinateFrame'], 'XChannel'):
                x = f[groupName]['CoordinateFrame'].attrs['XChannel']
            else:
                return None
        if y == '':
            if whizzAttrExists(f[groupName]['CoordinateFrame'], 'YChannel'):
                y = f[groupName]['CoordinateFrame'].attrs['YChannel']
            else:
                return None

        lineDistance = 0.0
        count = 0
    
        if lines == []:
            lines = list(g.keys())
        for line in lines:
            xPos = rd.getLineData(g[line], x)
            yPos = rd.getLineData(g[line], y)
            lineDistance += _lineLength(xPos, yPos)
            count += 1
            
    print(f'{count} lines: total distance flown [km] = {lineDistance/1000.0:,.1f}')
    return (count, lineDistance/1000.0)
    

def _lineLength(x, y):
    """
    Calculates the point-to-point length of a line defined by discrete (x,y) points.

    Parameters
    ----------
    x (numpy 1D array): a 1D array of x positions (units metres).

    y (numpy 1D array): a 1D array of y positions (units metres).

    Returns
    -------
    lenOfLine (Float) : the length of the line in metres.

    """
    lenOfLine = 0.0
    for ii in range(0, len(x)-1):
        lenOfLine += util.e2norm(x[ii] - x[ii+1], y[ii] - y[ii+1])
    return lenOfLine


def whizzAttrExists(group, my_attr):
    """
    True if `my_attr` is an attribute of the HDF5 `group`.
    """
    return my_attr in group.attrs.keys()









