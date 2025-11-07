#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Write the contents of an `XYZ` file to a `geoWhizz` file.
"""
import numpy as np
import h5py
from pathlib import Path
import pathlib
import filebrowser as fb

import pegasusQC.utility.utility as util
import pegasusQC.config as config

groupName = config.groupName
projectName = config.projectName


def xyzToHDF(xyzFilePath = '', hdfFilePath = '', projectName = '', verbose=False):
    """
    Read in a Geosoft XYZ file and write the contents to a Whizz HDF5 file.

    Parameters
    ----------
    xyzFilePath : pathlib Path
        The pathlib Path to the Geosoft XYZ file.
    hdfFilePath : pathlib Path, optional
        The pathlib Path to the Whizz HDF5 file to be created. The default is ''
        leaving the Whizz file to have the same path as the input file but with
        a ".HDF5" extension.
    projectName : String, optional
        The name of the survey or project. The default is ''.
    line_scale : float, optional
        All line numbers in the XYZ file are multiplied by line_scale before
        writing to the HDF5 file.
    verbose : Bool, optional
        If False (the default) the output is reduced.

    Returns
    -------
    None.

    """
        
    if xyzFilePath == '':
        xyzFilePath = fb.get_grid_filename()

    if isinstance(xyzFilePath, pathlib.Path):
        xyzFileStr = str(xyzFilePath)
    elif isinstance(xyzFilePath, str):
        xyzFileStr = xyzFilePath
        xyzFilePath = Path(xyzFileStr)
    else:
        print('Error - type of xyzFilePath not recognised. Must be Path or String')
        return
    if hdfFilePath == '':
        hdfFilePath = xyzFilePath.with_suffix('.hdf5')
        hdfFileStr = str(hdfFilePath)
    elif isinstance(hdfFilePath, pathlib.Path):
        hdfFileStr = str(hdfFilePath)
    elif isinstance(hdfFilePath, str):
        hdfFileStr = hdfFilePath
        hdfFilePath = Path(hdfFileStr)
    else:
        print('Error - type of hdfFilePath not recognised. Must be Path or String')
        return

    # access the data via Geosoft XYZ
    geosoftXYZ = readXYZ(xyzFileStr)
    if geosoftXYZ is None:
        print('xyzToHDF stopped. No geosoftXYZ data.')
        return hdfFilePath
    desiredFieldNames = [*geosoftXYZ[0]]
    num_lines = len(geosoftXYZ)
    lines = np.zeros((num_lines,))

    # get the line numbers from geosoftXYZ, checking each is unique.
    uniquelines = set()
    for count in range(0, num_lines):
        nextline = geosoftXYZ[count]['line_number']
        if nextline in uniquelines:
            print(f'ERROR - duplicate flight-line {nextline} in {xyzFileStr}.')
            return
        uniquelines.add(nextline)
        lines[count] = nextline
        
    print('Creating: ', hdfFileStr)
        
    with h5py.File(hdfFilePath, 'w') as f: # better data file
        # copy over survey metadata
        g = f.create_group(groupName)
        g.attrs['ProjectName'] = projectName
        gCoord = g.create_group('CoordinateFrame')
        gLines = g.create_group('Lines')
        
        # create first line group and metadata
        gg = gLines.create_group(f'{lines[0]:.3f}')
        gg.attrs['LineNumber'] = lines[0]

        # put this lines data in dataset; create the next line group and metadata
        for lineCount in range(0, len(lines)-1): #for aLine in lines:
            if verbose:
                print('  About to add data for line ', lines[lineCount], ' Lcount ', lineCount+1, '/', len(lines))
            # put all the data from last block into data set in current group
            # assume first field is line number
            for count in range(1, len(desiredFieldNames)):
                #print('Field name ', desiredFieldNames[count])
                dataArray = geosoftXYZ[lineCount][desiredFieldNames[count]]
                dd = gg.create_dataset(desiredFieldNames[count], \
                                       data = dataArray, dtype='float64', \
                                           compression="gzip", compression_opts=4)
                dd.attrs['Name'] = desiredFieldNames[count]
            
            #print(gg, lineCount)
            
            # create next line group and metadata
            gg = gLines.create_group(f'{lines[lineCount + 1]:.3f}')
            gg.attrs['LineNumber'] = lines[lineCount + 1]
        
        # put data in last line subgroup's dataset
        if verbose:
            print('  About to add data for line ', lines[len(lines)-1], ' Lcount ', len(lines), '/', len(lines), '\n')
        for count in range(1, len(desiredFieldNames)):
           dataArray = geosoftXYZ[len(lines)-1][desiredFieldNames[count]]
           #print(dataArray.shape)
           dd = gg.create_dataset(desiredFieldNames[count], \
                                  data = dataArray, dtype='float64', \
                                      compression="gzip", compression_opts=4)
           dd.attrs['Name'] = desiredFieldNames[count]
        
    return hdfFilePath


def read_xyz_header(filename):
    """
    Read in a Geosoft XYZ file and write the contents to a Whizz HDF5 file.

    Parameters
    ----------
    xyzFilePath : String
        The name of the Geosoft XYZ file.

    Returns
    -------
    Dictionary containing channel names and the first data record.

    """

    num_lines = 0
    num_head_recs = 0
    num_channels = 0
    need_first_data_rec = True

    # count number of header records, number of flight lines in file and number of channels
    with open(filename, 'r') as fid:
        for file_line in fid:
            if file_line[0] == '/':
                num_head_recs += 1
            elif file_line.startswith('Line') or file_line.startswith('Tie'):
                num_lines += 1
            elif num_lines == 1 and need_first_data_rec and not('*' in file_line):
                first_data_rec = file_line.split()
                field_precisions = np.zeros((len(first_data_rec),), dtype=int)
                for ii in range(0, len(first_data_rec)):
                    bits = first_data_rec[ii].split('.')
                    if bits[0] == first_data_rec[ii]:
                        field_precisions[ii] = 0
                    else:
                        field_precisions[ii] = len(bits[1])
                num_channels = len(first_data_rec)
                need_first_data_rec = False
            else:
                first_data_rec = file_line.split()
                continue
    fid.close()
    print(f'Found {num_head_recs} header records')
    print(f'Found {num_lines} lines')
    print(f'Found {num_channels} fields')
    print(f'Channel precisions (number of decimal places):')
    util._print_wrappedlist(f'{field_precisions}')
    if need_first_data_rec:
        print('Could not find a record without * - so no precisions calculated.')
        with open(filename, 'r') as fid:
            for file_line in fid:
                if file_line[0] == '/':
                    print(file_line)
                elif file_line.lstrip().upper().startswith('LINE') or file_line.lstrip().upper().startswith('TIE'):
                    print(file_line)
                elif need_first_data_rec:
                    #print(file_line)
                    need_first_data_rec = False
                else:
                    break
        fid.close()
        return
    
    # get channel names
    header_rec = 0
    with open(filename, 'r') as fid:
        for file_line in fid:
            temp_names = file_line[1:].split()
            header_rec += 1
            if len(temp_names) == num_channels:
                channelnames = temp_names
                for ii in range(0, len(channelnames)):
                    print(f'{channelnames[ii]} - precision {field_precisions[ii]}')
                break
            if header_rec > num_head_recs:
                print(f"Error - can't find header record with {num_channels} channel names.")
                break

    fid.close()
    first_dict = dict(zip(channelnames, first_data_rec))
    return first_dict

    
def readXYZ(filename, verbose=False):
    """
    Open a Geosoft XYZ file, read in the contents and return a dictionary
    array where each element corresponds to one flight line and contains the
    line number, flight number, date, and data for each channel.
    
    """

    num_lines = 0
    num_head_recs = 0
    num_channels = 0
    need_first_data_rec = True

    # count number of header records, number of flight lines in file and number of channels
    with open(filename, 'r') as fid:
        print(f'Accessing XYZ data in {filename}.\nFirst few records are:')
        for file_line in fid:
            if num_lines == 0:
                # Always useful to see the first few records in the file.
                util._print_wrappedlist(file_line)
                # print(file_line)
            if file_line.lstrip().startswith('/'):
                num_head_recs += 1
            elif file_line.lstrip().upper().startswith('LINE') or file_line.lstrip().upper().startswith('TIE'):
                num_lines += 1
            elif num_lines == 1 and need_first_data_rec:
                first_data_rec = file_line.split()
                field_precisions = np.zeros((len(first_data_rec),), dtype=int)
                for ii in range(0, len(first_data_rec)):
                    bits = first_data_rec[ii].split('.')
                    if bits[0] == first_data_rec[ii]:
                        field_precisions[ii] = 0
                    else:
                        field_precisions[ii] = len(bits[1])
                num_channels = len(first_data_rec)
                need_first_data_rec = False
            elif num_lines > 0:
                continue
            else:
                print(f'\n  Found {num_head_recs} header records')
                print(f'  Found {num_lines} lines')
                print(f'  Found {num_channels} fields\n')
                print(f'Channel precisions (number of decimal places):')
                util._print_wrappedlist(f'{field_precisions}')
                print('ERROR - no "LINE" records found in file, may not be a Geosoft XYZ file.')
                return None
    fid.close()
    print(f'\n  Found {num_head_recs} header records')
    print(f'  Found {num_lines} lines')
    print(f'  Found {num_channels} fields\n')
    print(f'Channel precisions (number of decimal places):')
    util._print_wrappedlist(f'{field_precisions}')
    
    # get channel names
    header_rec = 0
    with open(filename, 'r') as fid:
        for file_line in fid:
            if not file_line.lstrip().startswith('/'):
                continue
            file_line_trimmed = file_line.lstrip()[1:] # trim leading spaces and '/'
            temp_names = file_line_trimmed.split()
            header_rec += 1
            if len(temp_names) == num_channels:
                channelnames = temp_names
                if verbose:
                    print('  Found fields (channels):')
                    for ii in range(0, len(channelnames)):
                        print(f'    {channelnames[ii]} - precision {field_precisions[ii]}')
                break
            if header_rec > num_head_recs:
                print(f"Error - can't find header record with {num_channels} channel names.")
                break
    fid.close()
    
    # initialise counter for the number of fids per line and storage for lines, flights and dates
    num_fids = np.zeros((num_lines,), dtype=int)
    line_nos = np.zeros((num_lines,))
    line_ctr = 0

    # get lines, flights and dates
    with open(filename, 'r') as fid:
        for file_line in fid:
            if file_line.lstrip().startswith('/'):
                if not file_line.lstrip().startswith('//'):  # already got channel names
                    continue
                elif file_line.lstrip().upper().startswith('//FLIGHT'): # we will get flight number from channel
                    continue
                elif file_line.lstrip().upper().startswith('//DATE'): # we will get date from channel
                    continue
            elif file_line.lstrip().upper().startswith('LINE') or file_line.lstrip().upper().startswith('TIE'):
                current_line = float(file_line.split()[1])
                line_nos[line_ctr] = f'{current_line:.3f}' # FIX THIS TODO
                line_ctr += 1
            elif line_ctr > 0:
                num_fids[line_ctr-1] += 1
            else:
                print('ERROR - no "LINE" records found in file, may not be a Geosoft XYZ file.')
                return None
    fid.close()

    # initialise dictionary list
    geosoftXYZ = [{'line_number': line_nos[0]}]
    for jj in range(num_channels):
        mynans = np.zeros(((num_fids[0]),))
        mynans[mynans == 0] = np.nan
        geosoftXYZ[0][channelnames[jj]] = mynans

    for ii in range(1, num_lines):
        temp = {'line_number': line_nos[ii]}
        for jj in range(num_channels):
            mynans = np.zeros(((num_fids[ii]),))
            mynans[mynans == 0] = np.nan
            temp[channelnames[jj]] = mynans

        geosoftXYZ.append(temp)

    # line_ctr = 0 # for each survey line
    line_ctr = 0

    with open(filename, 'r') as fid:
        for file_line in fid:
            if file_line.lstrip().startswith('/'):
                continue
            elif file_line.lstrip().upper().startswith('LINE') or file_line.lstrip().upper().startswith('TIE'):
                current_line = float(file_line.split()[1])
                geosoftXYZ[line_ctr]['line_number'] = f'{current_line:.3f}' # TODO FIX THIS
                fid_ctr = 0
                line_ctr += 1
            #elif file_line.count('*') != 0:
            #    print('Skip record for dummy')
            #    continue
            else:
                line_str = file_line.split()
                if len(line_str) != num_channels:
                    if verbose:
                        print(f'Skipping file record, size of data ({len(line_str)}) mis-match to {num_channels} channels')
                    continue
                for ii in range(num_channels):
                    try:
                        geosoftXYZ[line_ctr-1][channelnames[ii]][fid_ctr] = float(line_str[ii])
                    except:
                        geosoftXYZ[line_ctr-1][channelnames[ii]][fid_ctr] = float('nan')

                fid_ctr += 1

    fid.close()

    return geosoftXYZ


