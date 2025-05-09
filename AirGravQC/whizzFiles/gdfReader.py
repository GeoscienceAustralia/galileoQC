#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Import ASEG-GDF2 data and write to `geoWhizz` format.
"""
import numpy as np
import h5py
from pathlib import Path

# aseg_gdf2 uses pandas which is producing annoying Future Warnings, so this to suppress them.
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import aseg_gdf2 as aseg

import AirGravQC.config as config

groupName = config.groupName
projectName = config.projectName


def asegToHDF(gdf_datfile, whizzFile='', lineChannel='LINE', flightChannel='FLIGHT', dateChannel='DATE', omitChannels=[]):
    '''
    Reads the data from the ASEG-GDF2 survey file and writes it to a new Whizz
    HDF5 survey file. Uses the aseg_gdf2 package by Kent Inverarity at:
        https://github.com/kinverarity1/aseg_gdf2

    Parameters
    ----------
    gdf_datfile : pathlib.PosixPath
        The pathlib Path to the input ASEG-GDF2 DAT file.
    whizzFile : pathlib.PosixPath, optional
        The pathlib Path to the output Whizz HDF5 file. The default is '' and
        the output file is then the same as the input file with the extension
        changed to '.hdf5'.
    lineChannel : String, optional
        The name of the ASEG GDF2 channel containing the line numbers, defaults to 'LINE'.
    flightChannel : String, optional
        The name of the ASEG GDF2 channel containing the flight numbers, defaults to 'FLIGHT'.
    dateChannel : String, optional
        The name of the ASEG GDF2 channel containing the dates, defaults to 'DATE'.
    omitChannels : [String]
        An array of channel or field names to omit from the saved geoWhizz HDF5 file.

    Returns
    -------
    None.

    '''

    if whizzFile == '':
        whizzFile = gdf_datfile.with_suffix('.hdf5')
    
    # First, check DFN file is ASCII
    if not _ascii_check(gdf_datfile):
        return
    
    # open GDF, pull out channel names, the units, and the description
    gdf = aseg.read(str(gdf_datfile), engine="dask")#, method='fixed-widths')
    df = gdf.df()
    
    channelNames = gdf.field_names()
    channelsOut = channelNames # these are the channels we will save
        
    # next, identify the column containing line numbers etc
    # we won't save these columns but do want the first value for each line
    # so re-get the indices.

    # now if index not found - error/warning; else pop off the channelName list
    lineIdx = _index_containing_substring(channelsOut, lineChannel)
    if lineIdx < 0:
        print(f'ERROR - column index for {lineChannel} not found.\n')
        return
    else:
        channelsOut.pop(lineIdx)
    foundChannels = lineChannel
        
    haveFlights = False
    flightIdx = _index_containing_substring(channelNames, flightChannel)
    if flightIdx < 0:
        print('WARNING - no flight channel found.\n')
    else:
        channelsOut.pop(flightIdx)
        haveFlights = True
        foundChannels += ', ' + flightChannel
        
    haveDates = False
    dateIdx = _index_containing_substring(channelNames, dateChannel)
    if dateIdx < 0:
        print('WARNING - no date channel found.\n')
    else:
        channelsOut.pop(dateIdx)
        haveDates = True
        foundChannels += ', ' + dateChannel
        
    print(f'Key channels for line attributes found:\n  {foundChannels}.\n')
    
        
    if not omitChannels == []:
        for channel in omitChannels:
            tempIdx = _index_containing_substring(channelsOut, channel.lower())
            if tempIdx < 0:
                print(f'WARNING - {channel} to omit not found.\n')
            else:
                print(f'Omitted {channelsOut[tempIdx]} in column {tempIdx}.\n')
                channelsOut.pop(tempIdx)
            
    print(f'{len(channelsOut)} channels to be written to geoWhizz file: ')
    print(channelsOut)
    
    # now, iterate through GDF collecting line numbers and sizes
    print('\nGetting line numbers ...')
    lineNumbers = df[lineChannel].unique().compute()
    numLines = len(lineNumbers)

    print('... done. Ready to create HDF file and transfer data into it.\n')

    with h5py.File(str(whizzFile), 'w') as f:
        print(f'Write to geoWhizz file:')
        # create all the data structure ready for the datasets
        g = f.create_group(groupName)
        g.attrs['ProjectName'] = projectName
        
        gCoord = g.create_group('CoordinateFrame')
        gLines = g.create_group('Lines')
        
        # create all the line groups
        for current_line in lineNumbers:
            
            # create a line group and metadata
            gg = gLines.create_group(f'{current_line}')
            gg.attrs['LineNumber'] = current_line
            line_data = df[df[lineChannel] == current_line].compute()
            # Get the flight number and date for line attributes.
            if haveFlights:
                gg.attrs['Flight'] = line_data[flightChannel].values[0]
            if haveDates:
                gg.attrs['Date_Local'] = line_data[dateChannel].values[0]

            # for each line group, create the DataSets with attributes
            for channelName in channelsOut:
                my_data = np.array(line_data[channelName].values)                
                dd = gg.create_dataset(channelName, data=my_data, compression="gzip", compression_opts=4) #, dtype='float64'
                dd.attrs['Name'] = channelName
                dd.attrs['Alias'] = channelName
                fieldDef = gdf.get_field_definition(channelName)
                dd.attrs['Units'] = fieldDef['unit'].split(':')[0]
                dd.attrs['Description'] = fieldDef['long_name']
            gg.attrs['NumberOfFids'] = my_data.size

    return


def _ascii_check(gdf_datfile):
    """
    Returns True if all characters in `textfile` are ASCII, otherwise False.

    Parameters
    ----------
    gdf_datfile : Path
        The Path to the ASEG-GDF2 data file.

    Returns
    -------
    None.

    """
    allASCII = True
    
    dfn_file = gdf_datfile.with_suffix('.DFN')
    if not dfn_file.is_file():
        dfn_file = gdf_datfile.with_suffix('.dfn')
    if not dfn_file.is_file():
        print(f'ERROR. DFN file {dfn_file} not found')
        return False
    
    with open(dfn_file, mode='r', encoding='ascii') as file:
        try:
            lines = file.readlines()
        except:
            print(f'ERROR. Fails ASEG-GDF2 standard: non-ASCII character(s) in {dfn_file}')
            allASCII = False

    return allASCII


def _index_containing_substring(the_list, substring):
    '''
    Returns the index to the substring in the_list; -1 if not found.

    Parameters
    ----------
    the_list : String
        The string to be searched.
    substring : String
        The substring whose position in the_list is desired.

    Returns
    -------
    Integer
        The index to the start of the substring, -1 if not found.

    '''
    if len(substring) < 1:
        return -1
    lowlist = [name.upper() for name in the_list]
    for i, s in enumerate(lowlist):
        if substring.upper() == s:
              return i
    return -1



