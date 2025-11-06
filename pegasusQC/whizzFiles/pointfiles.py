#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions to add or update data or metadata.
"""
# import necessary modules

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy import interpolate
import xarray
import pathlib

# import pegasusQC.gridFiles.gridfiles as grd
import pegasusQC.gridFiles.read_ers as grd
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.utility.utility as util
import pegasusQC.config as config

groupName = config.groupName
projectName = config.projectName


def _translate_date(decimal_year):
    """
    Translate decimal year into some useful format TBD. STUB does nothing.

    Parameters
    ----------
    decimal_year : Float
        The year.

    Returns
    -------
    decimal_year.

    """
    return decimal_year


def _trim_ends(x0, x1, y1):
    """
    Remove elements from x1 (and corresponding elements from y1)
    outside the range [x0[0], x0[-1]].

    Parameters
    ----------
    x0 : 1D numpy float array
        The independent variable (plan) of the inputs.
    x1 : 1D numpy float array
        The independent variable to be trimmed to within the range of x0.
    y1 : 1D numpy float array
        The dependent variable to interpolate onto.

    Returns
    -------
    newtime : 1D numpy float array
        The values of dataIn interpolated onto timeOut.
    newdata : 1D numpy float array
        To be kept synchronised with timeOut and returned.

    """
    max_allowed = np.max([x0[0], x0[-1]])
    min_allowed = np.min([x0[0], x0[-1]])
    tmptime = x1[x1 > min_allowed]
    tmpdata = y1[x1 > min_allowed]
    newtime = tmptime[tmptime < max_allowed]
    newdata = tmpdata[tmptime < max_allowed]

    return newtime, newdata


def _trim_monotonic(data_in, sync=[]):
    """
    Force data_in to be monotonic by removing samples; keep
    sync aligned by removing corresponding samples there 
    as well.

    Parameters
    ----------
    data_in : 1D numpy float array
        The independent variable (plan) of the inputs.
    sync : 1D numpy float array, optional
        Data to be kept synchronised with data_in.

    Returns
    -------
    data_trim : 1D numpy float array
        The values of data_in after monotonic trimming.
    sync_out : 1D numpy float array
        The sync data after removing samples corresponding to samples removed from data_in.

    """
    b = np.diff(data_in)
    # print(len(b[b > 0]), len(b[b < 0]))
    if len(b[b > 0]) < len(b[b < 0]):
        # if it is mostly decreasing, reverse, ...
        data_temp = data_in[::-1]
        if sync.size != 0:
            sync_temp = sync[::-1]
    else:
        data_temp = data_in
        if len(sync) != 0:
            sync_temp = sync
    # Now mostly increasing, keep only the increases
    b = np.concatenate(([0],np.diff(data_temp)))
    data_trim = data_temp[b > 0]
    if len(sync) != 0:
        sync_out = sync_temp[b > 0]
    else:
        sync_out = []
    return data_trim, sync_out


def updateProject(whizzFile, projectName='', blockID='', acquirer='', acquirerProjectID='', reportName=''):
    """
    Change any of the project attributes in the HDF5 Whizz file. Typically most of this information
    is not available when the Whizz file is created and this routine is used to add it when it is
    available.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    projectName : String, optional
        The name of the survey. The default is ''.
    blockID : String, optional
        The name of the survey block (useful if there is more than 1 block in a survey). The default is ''.
    acquirer : String, optional
        The name of the company that acquired the survey data. The default is ''.
    acquirerProjectID : String, optional
        The project or job number used by the acquirer to identify the survey. The default is ''.
    reportName : String, optional
        The title of the logistics & processing report. The default is ''.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r+') as f:
        # create all the data structure ready for the datasets
        g = f[groupName]
        if projectName != '':
            g.attrs['ProjectName'] = projectName
            print(f'Setting ProjectName = {projectName} for {whizzFile.name}.')
        if blockID != '':
            g.attrs['BlockID'] = blockID
            print(f'Setting BlockID = {blockID} for {whizzFile.name}.')
        if acquirer != '':
            g.attrs['Acquirer'] = acquirer
            print(f'Setting Acquirer = {acquirer} for {whizzFile.name}.')
        if acquirerProjectID != '':
            g.attrs['AcquirerProjectID'] = acquirerProjectID
            print(f'Setting AcquirerProjectID = {acquirerProjectID} for {whizzFile.name}.')
        if reportName != '':
            g.attrs['ReportName'] = reportName
            print(f'Setting ReportName = {reportName} for {whizzFile.name}.')
         
    
def updateCoordFrame(whizzFile, lat='', lon='', geoDatum='', alt='', htDatum='', x='', y='', projection='', utmz='', time='', timeDatum='', fid=''):
    """
    Change any of the attributes of the coordinate frame for the survey. This includes the names of the x,y,z,t coordinate fields
    and their respective datums. Typically these attributes are not entered to the Whizz file at creation and must be added later,
    by this routine.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    lat : String, optional
        The name of the latitude field. The default is ''.
    lon : String, optional
        The name of the longitude field. The default is ''.
    geoDatum : String, optional
        The geographic datum (eg 'WGS84', 'GDA94'). The default is ''.
    alt : String, optional
        The name of the altitude field. The default is ''.
    htDatum : String, optional
        The height datum (eg 'AHD71'). The default is ''.
    x : String, optional
        The name of the x or easting field. The default is ''.
    y : String, optional
        The name of the y or northing field. The default is ''.
    projection : String, optional
        The map projection (eg 'UTM'). The default is ''.
    utmz : String, optional
        If the projection is 'UTM', then the UTM zone (eg '51S'). The default is ''.
    time : String, optional
        The name of the time field. The default is ''.
    timeDatum : String, optional
        The time datum (eg 'UTC', 'local', 'GPS'. The default is ''.
    fid : String, optional
        The name of the fiducial field. The default is ''.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r+') as f:
        # create all the data structure ready for the datasets
        g = f[groupName]
        cf = g['CoordinateFrame']
        changed = False
        if lat != '':
            cf.attrs['LatitudeChannel'] = lat
            changed = True
        if lon != '':
            cf.attrs['LongitudeChannel'] = lon
            changed = True
        if geoDatum != '':
            cf.attrs['GeoDatum'] = geoDatum
            changed = True
        if alt != '':
            cf.attrs['AltitudeChannel'] = alt
            changed = True
        if htDatum != '':
            cf.attrs['HeightDatum'] = htDatum
            changed = True
        if x != '':
            cf.attrs['XChannel'] = x
            changed = True
        if y != '':
            cf.attrs['YChannel'] = y
            changed = True
        if projection != '':
            cf.attrs['Projection'] = projection
            changed = True
        if utmz != '':
            cf.attrs['UTMZone'] = utmz
            changed = True
        if time != '':
            cf.attrs['TimeChannel'] = time
            changed = True
        if timeDatum != '':
            cf.attrs['TimeDatum'] = timeDatum
            changed = True
        if fid != '':
            cf.attrs['FidChannel'] = fid
            changed = True
        if changed:
            print(f'Changed CoordFrame attribute(s) for {whizzFile.name}.')
 
    
def renameChannels(nchannels, chanNames):
    """
    Any channel named 'x' is renamed 'X'.

    Redundant function since the use of `XChannel` metadata.

    Parameters
    ----------
    nchannels : Int
        The number of channels.
    chanNames : [String]
        Array of channel names.

    Returns
    -------
    chanNames : [String]
        Array of renamed channel names.

    """
    
    for ii in range(0, nchannels): # TODO: expand this to cover a large range of possibilities
        if chanNames[ii] == 'x':
            chanNames[ii] = 'X'
            
    return chanNames

    
def addWhizzToWhizz(inputWhizzFile, outputWhizzFile, lines=[]):
    """
    Adds all the ['Lines'] sub-groups (with all they contain) from `inputWhizzFile` to the
    ['Lines'] group in `outputWhizzFile`

    TODO - Need to translate channel names
    """
    infile = str(inputWhizzFile)
    outFile = str(outputWhizzFile)
    if True:
        print("Can't be used until channel name translation is in place.")
        return
    
    with h5py.File(outFile, 'r+') as outf:
        outLines = outf[groupName]['Lines']

        with h5py.File(infile, 'r+') as inf:
            if lines == []:
                inLines = inf[groupName]['Lines']
                numLines = len(inLines.items())
                lines = inLines.values()

            # message = ''
            # num_lines_exceeded = 0
            # num_lines_unplanned = 0

            for line in lines:
                lineNumber = line.attrs['LineNumber']
                print(line, lineNumber)
                print(outLines['1001.0'])
                inf.copy(line, outLines)


def addLineToWhizz(hdf5FileId, databaseName, lineNo, lineType, plannedX = [], plannedY = [], plannedZ = [], distanceUnits = ''):
    """
    WARNING - A ONE-OFF

    Adds a line sub-group to a database sub-group of the project.

    Parameters
    ----------
    hdf5FileId : TYPE
        DESCRIPTION.
    databaseName : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    project = hdf5FileId["project"]
    database = project[databaseName]
    lineName = f'{lineType}{lineNo}'
    line = database.create_group(lineName)
    line.attrs['lineNumber'] = lineNo
    line.attrs['lineType'] = lineType
    line.attrs['qcChecked'] = False
    line.attrs['qcPassed'] = False
    
    if plannedX.count == plannedY.count == plannedZ.count:
        if plannedX.count > 1:
            xData = line.create_dataset('xPlan', data = plannedX, compression="gzip", compression_opts=4)
            xData.attrs['Units'] = distanceUnits
            xData.attrs['Precision'] = 0.01
            xData.attrs['Description'] = 'xPlan'
            yData = line.create_dataset('yPlan', data = plannedX, compression="gzip", compression_opts=4)
            yData.attrs['Units'] = distanceUnits
            yData.attrs['Precision'] = 0.01
            yData.attrs['Description'] = 'yPlan'
            zData = line.create_dataset('zPlan', data = plannedX, compression="gzip", compression_opts=4)
            zData.attrs['Units'] = distanceUnits
            zData.attrs['Precision'] = 0.01
            zData.attrs['Description'] = 'zPlan'
    return        


def doyToISO8601(doy):
    """
    Converts day-of-year in the format yyyymmdd into an ISO 8601 datetime string. No timezone is recorded.

    Parameters
    ----------
    doy : Int
        Day-of-year encoded as yyyymmdd and stored as Int.

    Returns
    -------
    Str
        The date as a naive (to timezone) ISO 8601 string.

    """
    return datetime.datetime.strptime(str(doy), '%Y%m%d').isoformat()


