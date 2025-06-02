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


# channel indexing wrong!
# flights and dates not attributes
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
        return None
    
    # open GDF, pull out channel names, the units, and the description
    gdf = aseg.read(str(gdf_datfile), engine="dask", method='fixed-widths')
    df = gdf.df()
    
    channelsOut, channelindices, haveFlights, haveDates = _getDesiredChannels(gdf, lineChannel, flightChannel, dateChannel, omitChannels)
    if channelsOut is None:
        return None

# get line numbers and numlines [and numrecords[lines]]
    
    lineNumbers, numLines, numRecs = _getlinenumbers(df, lineChannel)

# create and open whizzfile
# set projectname, create coordframe and lines groups
    with h5py.File(str(whizzFile), 'w') as f:
        print(f'Writing to geoWhizz file: {str(whizzFile)}')
        # create all the data structure ready for the datasets
        g = f.create_group(groupName)
        g.attrs['ProjectName'] = projectName
        
        gCoord = g.create_group('CoordinateFrame')
        gLines = g.create_group('Lines')

        with open(gdf_datfile, 'r') as datfid:
            numrecords = numRecs
            chans = gdf.field_names() # not quite right
            num_readchans = len(chans) - len(omitChannels)#omit_chans)
            names = gdf.record_types.df().name.values
            widths = gdf.record_types.df().width.values
            pytypes = gdf.record_types.df().inferred_dtype.values

            for lidx, current_line in enumerate(lineNumbers):

                # Get the data for this line for allchannels from .DAT
                mydata = np.zeros((numrecords[lidx], num_readchans))
                readchans = [None] * num_readchans
                reccount = -1
                for recidx in range(0, numrecords[lidx]):
                    myrec = datfid.readline()
                    reccount += 1
                    start = 0
                    end = 0
                    arrayidx = 0
                    
                    typefailure_report = ''
                    typefailure_count = 0

                    for channum, chan in enumerate(chans):
                        for idx, name in enumerate(names):
                            if name == chan:
                                mywidth = int(widths[idx])
                                mytype = pytypes[idx]
                        
                        end += mywidth
                        extract = myrec[start:end]
                        start = end
                        if not chan in channelsOut:#omitChannels:#omit_chans:
                            continue
                        try:
                            if mytype is float:
                                mydata[reccount, arrayidx] = float(extract)
                                readchans[arrayidx] = chan
                                arrayidx += 1
                            elif mytype is int:
                                mydata[reccount, arrayidx] = int(extract)
                                readchans[arrayidx] = chan
                                arrayidx += 1
                            elif mytype is str:
                                mydata[reccount, arrayidx] = 0.0
                                readchans[arrayidx] = chan
                                arrayidx += 1
                            else:
                                typefailure_count += 1
                                typefailure_report += f'else - {myrec[0:59]} |{extract}| {reccount}, {channum}\n'
                                break
                        except:
                            typefailure_count += 1
                            typefailure_report += f'except - {myrec[0:59]} |{extract}| {reccount}, {channum}\n'
                            break
                    if typefailure_count > 1:
                        break               

                # create a line group with metadata in whizzfile, then ...
                gg = gLines.create_group(f'{current_line}')
                gg.attrs['LineNumber'] = current_line
                gg.attrs['NumberOfFids'] = numrecords[lidx]#my_data[:,1].size
                # if haveFlights:
                #     gg.attrs['Flight'] = line_data[flightChannel].values[0]
                # if haveDates:
                #     gg.attrs['Date_Local'] = line_data[dateChannel].values[0]

                # ... create the desired DataSets with attributes
                for channum, channelName in enumerate(channelsOut):
                    if channelName == readchans[channum]:
                        dd = gg.create_dataset(channelName, data=mydata[:,channum], compression="gzip", compression_opts=4) #, dtype='float64'
                        dd.attrs['Name'] = channelName
                        dd.attrs['Alias'] = channelName
                        fieldDef = gdf.get_field_definition(channelName)
                        dd.attrs['Units'] = fieldDef['unit'].split(':')[0]
                        dd.attrs['Description'] = fieldDef['long_name']
                    else:
                        print(f'ERROR - {channelName} != {readchans[channum]} in line {current_line}')
                        print(readchans)
                        print(channelsOut)
                        print(omitChannels)
                        return None

    print('Complete.')
    return whizzFile


def _getlinenumbers(df, lineChannel):
    """
    Returns the line numbers, number of lines, and number
    of records in each line, from a pandas dataframe read
    from an ASEG-GDF2 file.
    """
    print('\nGetting line numbers ...')
    lineNumbers = df[lineChannel].unique().compute()
    numLines = len(lineNumbers)
    gby = df.groupby(df[lineChannel])
    numRecs = gby.count().compute().values[:,0]
    return lineNumbers, numLines, numRecs



def _getDesiredChannels(gdf, lineChannel, flightChannel, dateChannel, omitChannels):

    channelNames = gdf.field_names()
    channelsOut = channelNames # these are the channels we will save
    channelindices = list(range(0,len(channelsOut)))

    # next, identify the column containing line numbers etc
    # we won't save these columns but do want the first value for each line
    # so re-get the indices.

    # now if index not found - error/warning; else pop off the channelName list
    lineIdx = _index_containing_substring(channelsOut, lineChannel)
    if lineIdx < 0:
        print(f'ERROR - column index for line channel {lineChannel} not found.\n')
        return None, None, None, None
    else:
        channelsOut.pop(lineIdx)
        channelindices.pop(lineIdx)
    foundChannels = lineChannel
        
    haveFlights = False
    flightIdx = _index_containing_substring(channelNames, flightChannel)
    if flightIdx < 0:
        print('WARNING - no flight channel found.\n')
    else:
        channelsOut.pop(flightIdx)
        channelindices.pop(flightIdx)
        haveFlights = True
        foundChannels += ', ' + flightChannel
        
    haveDates = False
    dateIdx = _index_containing_substring(channelNames, dateChannel)
    if dateIdx < 0:
        print('WARNING - no date channel found.\n')
    else:
        channelsOut.pop(dateIdx)
        channelindices.pop(dateIdx)
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
                channelindices.pop(tempIdx)
            
    print(f'{len(channelsOut)} channels to be written to geoWhizz file: ')
    print(channelsOut)
    print(channelindices, haveFlights, haveDates)

    return channelsOut, channelindices, haveFlights, haveDates


def asegToHDF_old(gdf_datfile, whizzFile='', lineChannel='LINE', flightChannel='FLIGHT', dateChannel='DATE', omitChannels=[]):
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
    gdf = aseg.read(str(gdf_datfile), engine="dask", method='fixed-widths')
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

#

