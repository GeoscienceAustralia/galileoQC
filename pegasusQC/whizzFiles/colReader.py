#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Import columnar table ASCII data and write to `geoWhizz` format.
Author: Mark Helm Dransfield
Created: 2023
License: CC BY-SA
"""

import numpy as np
import h5py
from pathlib import Path

import pegasusQC.config as config
from pegasusQC.whizzFiles.updateLineAttributes import updateLineAttributes

groupName = config.groupName
projectName = config.projectName


def colsToHDF(col_datfile, lineColumn, dateColumn, latColumn, lonColumn, dumstr='*', whizzFile=''):
    '''
    Reads the data from the ASCII survey file and writes it to a new Whizz
    HDF5 survey file.

    Parameters
    ----------
    col_datfile : pathlib.PurePath
        The pathlib Path to the input ASCII data file.
    whizzFile : pathlib.PurePath, optional
        The pathlib Path to the output Whizz HDF5 file. The default is '' and
        the output file is then the same as the input file with the extension
        changed to '.hdf5'.
    lineColumn : Integer
        The column number containing the line numbers.

    Returns
    -------
    whizzFile : pathlib.PurePath
        The name of the created whizz file.

    '''
    if whizzFile == '':
        whizzFile = col_datfile.with_suffix('.hdf5')
    
    # First, check data file is ASCII
    if not _ascii_check(col_datfile):
        return None
    
    # open file and pull out channel names
    channelNames = _getChannels(col_datfile)
    if channelNames is None:
        return None

    for channelName in channelNames:
        if channelName is None:
            print(f'ERROR - missing channel name in {channelNames}.')
            return

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

        with open(col_datfile, 'r') as datfid:

            first_record = True
            for record in datfid:
                if record[0] == "/":
                    continue
                record_list = record.split()
                line_in_rec = record_list[lineColumn]
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
                    gLine = save_rec_lists(gLines, current_line, channelNames, record_lists, dumstr)
                    
                    # check next flight-line is unique
                    if line_in_rec in flightlines:
                        print(f'ERROR - duplicate flight-line {line_in_rec} in {col_datfile}.')
                        return
                    flightlines.add(line_in_rec)
                    # prepare for next flight-line
                    current_line = line_in_rec
                    record_lists = []
                    record_lists.append(record_list)
                    count = 1

            # finalise last flight-line
            print(f'L {current_line}... count = {count}')
            gLine = save_rec_lists(gLines, current_line, channelNames, record_lists, dumstr)
    
    print('Complete.')
    return whizzFile


def save_rec_lists(wizLines, current_line, channelNames, record_lists, dumstr):
    """
    Saves a full flight-line of data in the array of string arrays
    into the wizfid lines group. Returns the flight-line group.
    """
    # print(f'Processing line {current_line}\n')

    # create a line group with metadata in whizzfile, then ...
    gg = wizLines.create_group(f'{current_line}')
    gg.attrs['LineNumber'] = float(current_line)
    gg.attrs['NumberOfFids'] = len(record_lists)

    # create a float np array of size = size(record_lists)
    mydata = np.array(record_lists) #record_list_to_float(record_lists, pytypes, nulls)

    # convert str to float and store in mydata

    # ... create the desired DataSets with attributes
    for chanIdx, channelName in enumerate(channelNames):
        column_data = mydata[:,chanIdx]
        column_data = np.array([s.replace(dumstr, 'nan') for s in column_data])
        chanstr = str(channelName)
        if chanstr == "DATE":
            column_data = _decimal_year(column_data, 'nan')
        if chanstr == "LATITUDE":
            column_data = _decimal_degrees(column_data, 'nan')
        if chanstr == "LONGITUDE":
            column_data = _decimal_degrees(column_data, 'nan')
        dd = gg.create_dataset(chanstr, data=column_data.astype(float), compression="gzip", compression_opts=4)
        dd.attrs['Name'] = chanstr
        dd.attrs['Alias'] = chanstr

    return gg


def _decimal_degrees(angles, dumstr):
    """
    Given an array of strings, each of the form '-10.12.45.12234' implying
    -10 deg, 12 min, 45.12234 sec, returns a Numpy array of the equivalent
    values in decimal degrees.
    """
    decdeg = np.zeros((len(angles),))
    for i, angle in enumerate(angles):
        if angle == dumstr:
            decdeg[i] = np.nan
            continue
        parts = angle.split('.')
        if len(parts) == 2:
            decdeg[i] = float(angle)
        elif float(parts[0]) < 0:
            decdeg[i] = float(parts[0]) - (float(parts[1]) / 60.0 + float(parts[2]) / 3600.0 + float("0." + parts[3]) / 3600.0)
        else:
            decdeg[i] = float(parts[0]) + float(parts[1]) / 60.0 + float(parts[2]) / 3600.0 + float("0." + parts[3]) / 3600.0 
    return decdeg


def _decimal_year(dates, dumstr):
    """
    Given an array of strings, each of the form '2013/03/09' implying
    the date 9th March 2013, returns a Numpy array of the equivalent
    values in decimaly years.
    """
    dayspermonth = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    cumdays = np.cumsum(dayspermonth)
    lenyear = cumdays[-1]
    decyear = np.zeros((len(dates),))
    for i, date in enumerate(dates):
        if date == dumstr:
            decyear[i] = np.nan
            continue
        if '/' in date:
            parts = date.split('/')
            year = float(parts[0])
            month = int(parts[1])
            day = float(parts[2])
            days = cumdays[month - 2] + day
            if _leap(year) and int(month) > 2:
                days += 1
                lenyear += 1
            decyear[i] = year + days / lenyear
        else:
            decyear[i] = float(date)
    return decyear


def _leap(year):
    """
    Given a year, return True if it is a leap year, else False.
    """
    if (year % 400 == 0) and (year % 100 == 0):
        return True
    elif (year % 4 ==0) and (year % 100 != 0):
        return True
    else:
        return False


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
            break

    return mydata


def _getChannels(datafile):

    with open(datafile, 'r') as datfid:
        record = datfid.readline()

    if record[0] == "/":
        subrecs = record[1:].split()

    return subrecs


def _ascii_check(col_datfile):
    """
    Returns True if all characters in `col_datfile` are ASCII, otherwise False.

    Parameters
    ----------
    col_datfile : Path
        The Path to the data file.

    Returns
    -------
    None.

    """
    allASCII = True
        
    with open(col_datfile, mode='r', encoding='ascii') as file:
        try:
            lines = file.readlines()
        except:
            print(f'ERROR. non-ASCII character(s) in {col_datfile}')
            allASCII = False

    return allASCII






