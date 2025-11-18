#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions to add or update channel metadata.

Author: Mark Helm Dransfield

Created: 2023

License: CC BY-SA
"""

# import necessary modules
import h5py
import pegasusQC.config as config

groupName = config.groupName
projectName = config.projectName


def updateChannelAttributes(whizzFile, channel, name='', units='', alias='', description='', chan_precision=-1):
    """
    Updates the channel attributes for all lines in the geoWhizz HDF5 file. For any
    attribute, the default is to not change its value.

    Parameters
    ----------
    whizzFile : String or pathlib Path

        Name of a HDF5 Whizz file, including path and extension.

    channel : String

        The name of the channel whose attributes are to be changed.

    name : String, optional

        The new name. The default is ''.

    units : String, optional

        The new units. The default is ''.

    alias : String, optional

        The new alias. The default is ''.

    description : String, optional

        The new description. The default is ''.

    chan_precision : Integer, optional

        The new chan_precision (number of places after decimal point). The default is -1 (unknown).

    Returns
    -------
    None.

    """
    
    filename = str(whizzFile)

    with h5py.File(filename, 'r+') as f:
        # create all the data structure ready for the datasets
        g = f[groupName]['Lines']
        changed = False
        
        for line in g.keys():
            dd = g[line][channel]
            if name != '':
                dd.attrs['Name'] = name
                changed = True
            if units != '':
                if units in ['gu', 'µm/s/s', 'um/s/s', 'um/s2']:
                    units = 'µm/s/s'
                dd.attrs['Units'] = units
                changed = True
            if alias != '':
                dd.attrs['Alias'] = alias
                changed = True
            if description != '':
                dd.attrs['Description'] = description
                changed = True
            if chan_precision > -1:
                dd.attrs['chan_precision'] = chan_precision
                changed = True
        if changed:
            print(f'Changed channel attribute(s) for {channel} in {whizzFile.name}.')
        else:
            print(f'No changes made.')
    return


