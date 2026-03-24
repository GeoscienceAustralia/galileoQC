#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Retrieve data from a `geoWhizz` file.

Author: Mark Helm Dransfield

Created: 2023

License: CC BY-SA
"""

import numpy as np
import h5py

import galileoQC.config as config

groupName = config.groupName
projectName = config.projectName

def getWhizzData(whizzFile, line, channel):
    """
    Returns a numpy array containg the specified channel of
    data for the given line.

    Parameters
    ----------
    whizzFile : String or pathlib Path

        Name of a HDF5 Whizz file, including path and extension.

    line : String

        A flightline, e.g. '1000110.0'.

    channel : String

        The (case-sensitive) name of a channel in the database, e.g. 'EASTING'.

    Returns
    -------
    my_data : numpy array

        The requested data.

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        my_data = getLineData(g[line], channel)
            
    return my_data


def getLineData(linegroup, channel):
    """
    Returns a numpy array containg the specified channel of
    data for the line in the given linegroup.

    Parameters
    ----------
    linegroup : HDF5 Group

        A flight-line group.

    channel : String

        The name of a channel in the database, e.g. 'EASTING'.

    Returns
    -------
    my_data : numpy array

        The requested data.

    """
    my_data = np.array([])
    for datachannel in linegroup.items():
        if datachannel[0].upper() == channel.upper():
            # print(f'datachannel {datachannel[0]}; channel {channel}')
            my_data = np.array(linegroup[datachannel[0]])
    if np.array(my_data.size) == 0:
        print(f'ERROR - channel "{channel}" not found.')

    return my_data


def getChannelAttrs(linegroup, channel, myattribute='Units'):
    """
    Returns the requested attribute for the specified channel of
    data for the given line.

    Parameters
    ----------
    linegroup : HDF5 Group

        A flight-line group.

    channel : String

        The name of a channel in the database, e.g. 'EASTING'.

    myattribute : String, optional

        The name of the desired attribute, default 'Units'.

    Returns
    -------
    attr_value : String
        The `Units` attribute for channel, empty string if `Units` was not found.

    """
    for datachannel in linegroup.items():
        if datachannel[0].upper() == channel.upper():
            # print(f'datachannel {datachannel[0]}; channel {channel}')
            myChanGroup = linegroup[datachannel[0]]
            chanAttrs = list(myChanGroup.attrs)
            if myattribute in chanAttrs:
                attr_value = myChanGroup.attrs[myattribute]
                return attr_value
            break
                
    return ''


def getLineXChannel(whizzFile, line, x, channel):
    """
    Returns 1D numpy arrays of x and channel in line from geoWhizz filename. The
    inputs are just two channels in the data and the names 'x' and 'channel' do
    not carry intrinsic meaning.

    Intended use for analysis and plotting.

    Parameters
    ----------
    whizzFile : String or pathlib Path

        Name of a HDF5 Whizz file, including path and extension.

    line : String

        A flightline, e.g. '1000110.0'.

    x : String

        The name of the x variable.

    channel : String

        The name of the channel.

    Returns
    -------
    xData : numpy 1D array

        A numpy array of data, float32 or float64.

    yData : numpy 1D array

        A numpy array of data, float32 or float64.

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        xData = getLineData(g[line], x)
        yData = getLineData(g[line], channel)
    return xData, yData


def getChannels(whizzFile):
    """
    Returns an array of the channel names in a geoWhizz file.

    Parameters
    ----------
    whizzFile : String or pathlib Path

        Name of a geoWhizz file, including path and extension.

    Returns
    -------
    List of strings. The channel names.

    """
    filename = str(whizzFile)
        
    with h5py.File(filename, 'r') as f:
        
        whizzHeader = list(f.keys())[0]
        g = f[whizzHeader]        
        gLines = g['Lines']
        lineGroups = list(gLines.values())
        channelNames = list(lineGroups[0].keys())

    return channelNames


def getLines(whizzFile):
    """
    Returns an array of the line names in a geoWhizz file.

    Parameters
    ----------
    whizzFile : String or pathlib Path

        Name of a geoWhizz file, including path and extension.

    Returns
    -------
    List of strings. The Line numbers.

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        lines = list(g.keys())
    return lines


 