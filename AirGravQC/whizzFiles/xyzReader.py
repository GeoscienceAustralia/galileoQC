import numpy as np
import h5py
# import matplotlib.pyplot as plt
from pathlib import Path
import pathlib
# from scipy.interpolate import CloughTocher2DInterpolator
# from scipy import interpolate
import filebrowser as fb
# import aseg_gdf2 as aseg

# import AirGravQC.gridFiles.gridfiles as grd
# import AirGravQC.qc.qualityAnalysis as qc
import AirGravQC.utility.utility as util
import AirGravQC.config as config

groupName = config.groupName
projectName = config.projectName


def xyzToHDF(xyzFilePath = '', hdfFilePath = '', projectName = '', verbose=False):
    """
    Read in a Geosoft XYZ file and write the contents to a Whizz HDF5 file.

    Parameters
    ----------
    xyzFilePath : pathlib.PosixPath
        The pathlib Path to the Geosoft XYZ file.
    hdfFilePath : pathlib.PosixPath, optional
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
    desiredFieldNames = [*geosoftXYZ[0]]
    num_lines = len(geosoftXYZ)
    lines = np.zeros((num_lines,))

    for count in range(0, num_lines):
        lines[count] = geosoftXYZ[count]['line_number']
        
    print('Creating: ', hdfFileStr)
        
    with h5py.File(hdfFilePath, 'w') as f: # better data file
        # copy over survey metadata
        g = f.create_group(groupName)
        g.attrs['ProjectName'] = projectName
        gCoord = g.create_group('CoordinateFrame')
        gLines = g.create_group('Lines')
        
        # create first line group and metadata
        gg = gLines.create_group(f'{lines[0]}')
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
    
    **TODO**:
    
        1. remove reliance on Flight, Date and Line fields if possible.

        2. replace the entire routine with a xyzToHdf() function because this function
        will run out of memory if the XYZ file is large.
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
                print('  ', file_line)
            if file_line[0] == '/':
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
            else:
                continue
    fid.close()
    print(f'\n  Found {num_head_recs} header records')
    print(f'  Found {num_lines} lines')
    print(f'  Found {num_channels} fields\n')
    
    # get channel names
    header_rec = 0
    with open(filename, 'r') as fid:
        for file_line in fid:
            temp_names = file_line[1:].split()
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
    
    # rename key channels TODO : is this necessary now that whizz Files have 'XChannel' etc?
    # for jj in range(num_channels):
    #     if channelnames[jj] == 'x' or channelnames[jj] == 'X' or channelnames[jj] == 'easting' or channelnames[jj] == 'Easting':
    #         channelnames[jj] = 'X'
    #     if channelnames[jj] == 'y' or channelnames[jj] == 'Y' or channelnames[jj] == 'northing' or channelnames[jj] == 'Northing':
    #         channelnames[jj] = 'Y'

    # initialise counter for the number of fids per line and storage for lines, flights and dates
    num_fids = np.zeros((num_lines,), dtype=int)
    # flight_nos = np.zeros((num_lines,), dtype=int)
    line_nos = np.zeros((num_lines,))
    # flight_dates = [[1, 1, 1980] for k in range(num_lines)]
    line_ctr = 0
    # flight_no = 0
    # flight_date = [1, 1, 1980]

    # get lines, flights and dates
    with open(filename, 'r') as fid:
        for file_line in fid:
            if file_line[0] == '/':
                if file_line[1] != '/':  # already got channel names
                    continue
                elif file_line.lstrip().upper().startswith('//FLIGHT'):
                    flight_no = int(file_line.split()[1])
                elif file_line.lstrip().upper().startswith('//DATE'):
                    test = file_line.split()[1]
                    y, m, d = test.split('/')
                    flight_date = [int(d), int(m), int(y)]
            elif file_line.lstrip().upper().startswith('LINE') or file_line.lstrip().upper().startswith('TIE'):
                current_line = float(file_line.split()[1])
                line_nos[line_ctr] = f'{current_line:.3f}'
 #               flight_nos[line_ctr] = flight_no
  #              flight_dates[line_ctr] = flight_date
                line_ctr += 1
    #        elif file_line.count('*') != 0:
    #            continue  # print('Skip record for dummy')
            else:
                num_fids[line_ctr-1] += 1
    fid.close()

    # initialise dictionary list
    geosoftXYZ = [{'line_number': line_nos[0]}]
                   # 'flight_number': flight_nos[0],
                   # 'date_flown': flight_dates[0]
    for jj in range(num_channels):
        mynans = np.zeros(((num_fids[0]),))
        mynans[mynans == 0] = np.nan
        geosoftXYZ[0][channelnames[jj]] = mynans

    for ii in range(1, num_lines):
        temp = {'line_number': line_nos[ii]}
                # 'flight_number': flight_nos[ii],
                # 'date_flown': flight_dates[ii]}
        for jj in range(num_channels):
            mynans = np.zeros(((num_fids[ii]),))
            mynans[mynans == 0] = np.nan
            temp[channelnames[jj]] = mynans

        geosoftXYZ.append(temp)

    # line_ctr = 0 # for each survey line
    line_ctr = 0

    with open(filename, 'r') as fid:
        for file_line in fid:
            if file_line[0] == '/':
                continue
            elif file_line.lstrip().upper().startswith('LINE') or file_line.lstrip().upper().startswith('TIE'):
                current_line = float(file_line.split()[1])
                geosoftXYZ[line_ctr]['line_number'] = f'{current_line:.3f}'#file_line.split()[1]
                fid_ctr = 0
                line_ctr += 1
            #elif file_line.count('*') != 0:
            #    print('Skip record for dummy')
            #    continue
            else:
                line_str = file_line.split()
                if len(line_str) != num_channels:
                    print(line_ctr, fid_ctr)
                    print(line_str)
                    print('Number of channel names mis-match to size of data ', len(line_str))
                    break
                for ii in range(num_channels):
                    if line_str[ii] == '*':
                        line_str[ii] = 'nan'
                    # if line_str[ii] contains "/" then it is a date
                    #     decode date to decimal year string.
                    geosoftXYZ[line_ctr-1][channelnames[ii]][fid_ctr] = \
                        float(line_str[ii])

                fid_ctr += 1

    fid.close()

    return geosoftXYZ


