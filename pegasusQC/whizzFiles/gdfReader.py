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

import pegasusQC.config as config
from pegasusQC.whizzFiles.updateLineAttributes import updateLineAttributes

groupName = config.groupName
projectName = config.projectName


# channel indexing wrong!
# flights and dates not attributes
def asegToHDF(gdf_datfile, whizzFile='', lineChannel='LINE', flightChannel='', dateChannel='', omitChannels=[], verbose=False):
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
        The name of the ASEG GDF2 channel containing the flight numbers, defaults to ''.
    dateChannel : String, optional
        The name of the ASEG GDF2 channel containing the dates, defaults to ''.
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
    
    channelNames, channelindices, haveFlights, flightIdx, haveDates, dateIdx = _getDesiredChannels(gdf, lineChannel, flightChannel, dateChannel, omitChannels)
    if channelNames is None:
        return None

    for channelName in channelNames:
        if channelName is None:
            print(f'ERROR - missing channel name in {channelNames}.')
            return

    chans = gdf.field_names()
    num_readchans = len(chans) - len(omitChannels)
    names = [gdf.get_field_definition(chans[i])['name'] for i in range(0, len(chans))]
    widths = [gdf.get_field_definition(chans[i])['width'] for i in range(0, len(chans))]
    pytypes = [gdf.get_field_definition(chans[i])['inferred_dtype'] for i in range(0, len(chans))]
    nulls = [gdf.get_field_definition(chans[i])['null'] for i in range(0, len(chans))]

    # Calculate the start and end position in each file
    # record for each channel.
    starts = np.zeros(len(widths))
    starts[1:] = np.cumsum(widths[0:-1])
    starts = np.asarray(starts, dtype=int)
    ends = np.asarray(starts + np.array(widths), dtype=int)
    if verbose:
        print(f'\nwidths: {widths}')
        print(f'starts: {starts}')
        print(f'ends: {ends}')
        print(f'types: {pytypes}\n')

    # create and open whizzfile
    # set projectname, create coordframe and lines groups
    with h5py.File(str(whizzFile), 'w') as wizfid:
        print(f'Writing to geoWhizz file: {str(whizzFile)}')
        # create all the data structure ready for the datasets
        g = wizfid.create_group(groupName)
        g.attrs['ProjectName'] = projectName
        
        gCoord = g.create_group('CoordinateFrame')
        gLines = g.create_group('Lines')
        flightlines = set() # for uniqueness check.

        with open(gdf_datfile, 'r') as datfid:

            first_record = True
            for record in datfid:
                record_list = extract_from_record(record, starts, ends)
                line_in_rec = extract_line(record_list, gdf, lineChannel)
                if first_record:
                    flightlines.add(line_in_rec)
                    current_line = line_in_rec
                    record_lists = []
                    first_record = False
                    record_lists.append(record_list)
                    count = 1
                elif current_line == line_in_rec:
                    record_lists.append(record_list)
                    count += 1
                else:
                    # finalise previous flight-line
                    print(f'L {current_line}... count = {count}')
                    gLine = save_rec_lists(gLines, current_line, channelNames, pytypes, nulls, record_lists, gdf)
                    if haveFlights:
                        gLine.attrs['FlightNumber'] = int(record_lists[0][flightIdx])
                    if haveDates:
                        gLine.attrs['Date'] = float(record_lists[0][dateIdx])
                    
                    # check next flight-line is unique
                    if line_in_rec in flightlines:
                        print(f'ERROR - duplicate flight-line {line_in_rec} in {gdf_datfile}.')
                        return
                    flightlines.add(line_in_rec)
                    # prepare for next flight-line
                    current_line = line_in_rec
                    record_lists = []
                    record_lists.append(record_list)
                    count = 1

            # finalise last flight-line
            print(f'L {current_line}... count = {count}')
            gLine = save_rec_lists(gLines, current_line, channelNames, pytypes, nulls, record_lists, gdf)
            if haveFlights:
                gLine.attrs['FlightNumber'] = int(record_lists[0][flightIdx])
            if haveDates:
                gLine.attrs['Date'] = float(record_lists[0][dateIdx])
    
    updateLineAttributes(whizzFile, flight_chan=flightChannel, date_chan=dateChannel, verbose=False)
    print('Complete.')
    return whizzFile


def extract_from_record(myrec, starts, ends):
    """
    Returns a list of strings, each string being a data element from the record.
    """
    return [myrec[starts[i]:ends[i]] for i in range(0,len(starts))]


def extract_line(record_list, gdf, lineChannelName):
    """
    Returns the line number from the record list as a float.
    """
    lineIdx = _getDesiredChannel(gdf, lineChannelName)
    fieldDef = gdf.get_field_definition(lineChannelName)
    if fieldDef['inferred_dtype'] is int:
        return int(record_list[lineIdx])
    else:
        return float(record_list[lineIdx])


def save_rec_lists(wizLines, current_line, channelNames, pytypes, nulls, record_lists, gdf):
    """
    Saves a full flight-line of data in the array of string arrays
    into the wizfid lines group. Returns the flight-line group.
    """
    # print(f'Processing line {current_line}\n')

    # create a line group with metadata in whizzfile, then ...
    gg = wizLines.create_group(f'{current_line}')
    gg.attrs['LineNumber'] = current_line
    gg.attrs['NumberOfFids'] = len(record_lists)

    # create a float np array of size = size(record_lists)
    mydata = record_list_to_float(record_lists, pytypes, nulls)

    # convert str to float and store in mydata

    # ... create the desired DataSets with attributes
    for chanIdx, channelName in enumerate(channelNames):
        chanstr = str(channelName)
        fieldDef = gdf.get_field_definition(chanstr)
        if fieldDef['inferred_dtype'] is int:
            dd = gg.create_dataset(chanstr, data=mydata[:,chanIdx].astype(int), compression="gzip", compression_opts=4)
        else:
            dd = gg.create_dataset(chanstr, data=mydata[:,chanIdx], compression="gzip", compression_opts=4)
        dd.attrs['Name'] = chanstr
        dd.attrs['Alias'] = chanstr
        dd.attrs['Units'] = fieldDef['unit']#.split(':')[0]
        dd.attrs['Description'] = fieldDef['long_name']

    return gg


def record_list_to_float(record_lists, pytypes, nulls):
    """
    Convert the 2D string list into a 2D numpy float array.
    """
    num_chans = len(record_lists[0])
    num_recs = len(record_lists)
    mydata = np.zeros((num_recs, num_chans))

    typefailure_report = ''
    typefailure_count = 0

    for idx, record_list in enumerate(record_lists):
        for jdx, item in enumerate(record_list):
            try:
                if pytypes[jdx] is float:
                    mydata[idx, jdx] = float(item)
                elif pytypes[jdx] is int:
                    mydata[idx, jdx] = int(item)
                elif pytypes[jdx] is str:
                    mydata[idx, jdx] = 0.0
                else:
                    typefailure_count += 1
                    typefailure_report += f'else - |{item}| {idx}, {jdx}\n'
                    break
            except:
                typefailure_count += 1
                typefailure_report += f'except - |{item}| {idx}, {jdx}\n'
                break
            if type(nulls[jdx]) is str:
                if nulls[jdx] in item:
                    mydata[idx, jdx] = np.nan

        if typefailure_count > 1:
            print(typefailure_report)
            break


    return mydata


def _getDesiredChannel(gdf, channelName):

    allChannelNames = gdf.field_names()
    channelindices = list(range(0,len(allChannelNames)))

    # now if index not found - error/warning; else pop off the channelName list
    chanIdx = _index_containing_substring(allChannelNames, channelName)
    if chanIdx < 0:
        print(f'ERROR - column index for channel {channelName} not found.\n')
        print(allChannelNames)
        return None

    # print(f'ASEG-GDF2 channel {channelName} found at index {chanIdx}.\n')
    return chanIdx


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
        print(channelsOut)
        return None, None, None, None, None, None
    # else:
        # channelsOut.pop(lineIdx)
        # channelindices.pop(lineIdx)
    foundChannels = lineChannel + f' at {lineIdx}'
        
    haveFlights = False
    flightIdx = _index_containing_substring(channelNames, flightChannel)
    if flightIdx < 0:
        print('WARNING - no flight channel found.\n')
    else:
        # channelsOut.pop(flightIdx)
        # channelindices.pop(flightIdx)
        haveFlights = True
        foundChannels += ', ' + flightChannel + f' at {flightIdx}'
        
    haveDates = False
    dateIdx = _index_containing_substring(channelNames, dateChannel)
    if dateIdx < 0:
        print('WARNING - no date channel found.\n')
    else:
        # channelsOut.pop(dateIdx)
        # channelindices.pop(dateIdx)
        haveDates = True
        foundChannels += ', ' + dateChannel + f' at {dateIdx}' 
        
    print(f'Key channels for linegroup attributes found:\n  {foundChannels}.\n')
    
        
    if not omitChannels == []:
        for channel in omitChannels:
            tempIdx = _index_containing_substring(channelsOut, channel.lower())
            if tempIdx < 0:
                print(f'WARNING - {channel} to omit not found.\n')
            else:
                print(f'Omitted {channelsOut[tempIdx]} in column {tempIdx}.\n')
                # channelsOut.pop(tempIdx)
                # channelindices.pop(tempIdx)
            
    print(f'{len(channelsOut)} channels to be written to geoWhizz file: ')
    print(channelsOut)
    # print(channelindices, haveFlights, haveDates)

    return channelsOut, channelindices, haveFlights, flightIdx, haveDates, dateIdx


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






