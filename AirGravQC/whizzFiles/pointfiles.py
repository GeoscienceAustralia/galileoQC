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

import AirGravQC.gridFiles.gridfiles as grd
import AirGravQC.utility.utility as util
import AirGravQC.config as config

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
    
    # In the hdfFile, for all lines, update the attributes of the named channel if,
    # and only if, the passed attribute is not empty.
    
    filename = str(whizzFile)

    with h5py.File(filename, 'r+') as f:
        # create all the data structure ready for the datasets
        g = f[groupName]['Lines']
        changed = False
        
        for line in g.keys():
            dd = g[line][channel] # ToDo: make this case insensitive
            if name != '':
                dd.attrs['Name'] = name
                changed = True
            if units != '':
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
    return


def interpolateGridOntoLine(gridPath, hdfPath, lines=[]):
    """
    Interpolates the data in a grid file onto a new channel, with the same name
    as the gridPath stem, in the geoWhizz HDF5 file. Fills empty channel samples
    with nans.

    Parameters
    ----------
    gridPath : String or pathlib Path
        Name of a ERMapper file, including path and extension.
    hdfFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.

    Returns
    -------
    None.

    """
    gridFile = str(gridPath)
    hdfFile = str(hdfPath)

    # retrieve the grid information
    eg, ng, zg, datum, projection = grd.read_ers_image(gridPath)
    ng0 = np.min(ng)
    ngd = np.abs(ng[1] - ng[0])
    eg0 = np.min(eg)
    egd = np.abs(eg[1] - eg[0])
    zg = zg[::-1, :]
    newChannelName = gridPath.stem
    
    print('\nGrid file read for channel ', newChannelName)
    print(ng0, eg0, ngd, egd)
        
    with h5py.File(hdfFile, 'r+') as f:
        g = f[groupName]['Lines']
        lineGroups = list(g.values())
        channelNames = list(lineGroups[0].keys())
        x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        y = f[groupName]['CoordinateFrame'].attrs['YChannel']
        t = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
        height = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        if lines == []:
            lines = list(g.keys())
    #TODO: if newChannelName is in channelNames, append '_01' and check again.
    #TODO: handle x or y or zg = nan  
        for line in lines:
            lineNo = line
            lineText = 'Line ' + lineNo
            em = getLineData(g[line], x)
            nm = getLineData(g[line], y)
            time = getLineData(g[line], t)
            zm = np.empty(em.shape)
            zm[:] = np.nan
            dist = np.empty(em.shape)
            dist[:] = np.nan
            sampleSpacing = util.e2norm(em[1]-em[0], nm[1]-nm[0])
            print(lineText)
            print('  North min, max ', np.min(nm), np.max(nm), '  East min, max ', np.min(em), np.max(em))
            
            # estimate by weighted average about twice per grid cell
            subRate = int(max(1, int(min(ngd, egd)) / (2 * sampleSpacing)))
            ttemp = time[0:len(zm):subRate]
            ztemp = np.empty(ttemp.shape)
            ztemp[:] = np.nan
            kz = 0
            for kk in range(0, len(zm), subRate):
                ii = int(np.rint((nm[kk] - ng0) / ngd))
                jj = int(np.rint((em[kk] - eg0) / egd))
                #print('ii, jj ', ii, jj)
                
                if zg[ii, jj] == np.nan:
                    print('Error - nan in zGrid')
                    break
                    continue
                
                nCell = 0
                x = np.zeros((9,))
                y = np.zeros((9,))
                z = np.zeros((9,))
                x[:] = np.nan
                y[:] = np.nan
                z[:] = np.nan

                for icount in range(ii-1, ii+1):
                    if icount < 0:
                        continue
                    for jcount in range(jj-1, jj+1):
                        if jcount < 0:
                            continue
                        if zg[icount, jcount] == np.nan:
                            print('Error - nan in zGrid')
                            break
                            continue
                        x[nCell] = ng[icount]
                        y[nCell] = eg[jcount]
                        z[nCell] = zg[icount, jcount] # SWAP???
                        nCell += 1

                ztemp[kz] = _weightedAverage(x, y, z, em[kk], nm[kk])
                #print(x, y, z, em[kk], nm[kk], kz, ztemp[kz])                        
                kz += 1
                                        
            print(lineText, ' sampled ', kz, 'of ', len(zm), '. Max value = ', np.nanmax(ztemp))
            # Now interpolate through gaps by cubic spline
            (zm, _) = interpolateLine(ttemp, ztemp, time)
            
            # and store in line sub-group
            if newChannelName in g[line]:
                print('already there')
            else:
                dd = g[line].create_dataset(newChannelName, data = zm, compression="gzip", compression_opts=4)
                dd.attrs['Name'] = newChannelName
           
            fail, numExc = util._failsDeviation(zm - getLineData(g[line], height), 20.0, 13)
            if fail and numExc > 13:
 #           if np.abs(np.max(getLineData(g[line], 'altitude') - zm)) > 20.0:
                plotTime = time - time[0]
                fig = plt.figure()
                ax1 = fig.add_subplot(2,1,1)
                ax1.plot(plotTime, getLineData(g[line], 'altitude'), plotTime, zm)
                ax1.set_xlim(plotTime[0], plotTime[-1])
                plotTitle = f'Line {line}, altitude and planned drape height'
                plt.title(plotTitle, fontsize = 8)
                plt.grid(True)
                for label in ax1.get_xticklabels(): label.set_fontsize(6)
                for label in ax1.get_yticklabels(): label.set_fontsize(6)
                
                ax2 = fig.add_subplot(2,1,2)
                ax2.plot(plotTime, getLineData(g[line], 'altitude') - zm)
                ax2.set_xlim(plotTime[0], plotTime[-1])
                plotTitle = f'Difference'
                plt.title(plotTitle, fontsize = 8)
                plt.grid(True)
                for label in ax2.get_xticklabels(): label.set_fontsize(6)
                for label in ax2.get_yticklabels(): label.set_fontsize(6)
            
                plt.show()
            
    print("TODO")
    
    # find index of nearest grid cell to sample point at xs
    # where xo is the origin of the grid and dx is the cell size.
    # ii = np.rint((xs - xo) / dx).astype(int)
    return

    
def interpolateLine(timeIn, dataIn, timeOut, spare=[], plot_flag=False):
    """
    Interpolates dataIn, sampled at timeIn, onto the samples timeOut.
    These three input arrays are pre-processed to ensure that timeIn
    and timeOut are monotonically increasing, whilst keeping dataIn
    synchronised with timeIn, and spareOut synchronised with timeOut.

    Parameters
    ----------
    timeIn : 1D numpy float array
        The independent variable of the inputs to be interpolated.
    dataIn : 1D numpy float array
        The input dependent variable to be interpolated.
    timeOut : 1D numpy float array
        The independent variable to interpolate onto.
    spare : 1D numpy float array
        To be kept synchronised with timeOut and returned.
    plot_flag : Bool, optional
        If True, plot the np.diff() of dataIn and timeIn.
        Default False.

    Returns
    -------
    out : 1D numpy float array
        The values of dataIn interpolated onto timeOut.
    spareOut : 1D numpy float array
        To be kept synchronised with timeOut and returned.

    """
    
    timeIn = timeIn[np.logical_not(np.isnan(dataIn))]
    dataIn = dataIn[np.logical_not(np.isnan(dataIn))]
    spareOut = spare
    # print(np.min(np.diff(timeIn)), np.max(np.diff(timeIn)))
    # model = interpolate.InterpolatedUnivariateSpline(timeIn, dataIn)
    # return model(timeOut)

    min_length = 100
    if timeIn.size < min_length or dataIn.size < min_length or timeOut.size < min_length:
        out = np.zeros(timeOut.shape)
        out[:] = np.nan
        return out, spare
    if len(spare) < min_length:
        spare = []

    (timeIn_trim, dataIn_trim) = _trim_monotonic(timeIn, sync=dataIn)
    (timeOut_trim, spare_trim) = _trim_monotonic(timeOut, sync=spare)
    if len(spare) > min_length:
        spareOut = spare_trim

    if plot_flag:
        print(f'Len t: {len(timeIn)}; Len d: {len(dataIn)}')
        print(f'Shapes: dataIn {dataIn.shape}, {np.diff(dataIn).shape}; timeIn {timeIn.shape}')
        plt.plot(timeIn, dataIn, 'r', timeIn[1:], np.diff(dataIn), 'b', timeIn[1:], np.diff(timeIn), 'g')
        plt.show()
    spl = interpolate.splrep(timeIn_trim, dataIn_trim, k=3, s=0)
    # out = interpolate.splev(timeOut, spl)
    out = interpolate.splev(timeOut_trim, spl)
    return out, spareOut


def _trim_monotonic(data_in, sync=[]):
    """
    Force input to be monotonic by removing samples; keep
    sync aligned by removing corresponding samples there 
    as well.

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


def _invDist(x, y):
    """
    Returns the inverse of `distance(x,y)`
    """
    return 1.0 / util.e2norm(x, y)


def _weightedAverage(x, y, z, xo, yo):
    x = x[np.logical_not(np.isnan(z))]
    y = y[np.logical_not(np.isnan(z))]
    z = z[np.logical_not(np.isnan(z))]
    dx = x - xo
    dy = y - yo
    weight = np.zeros(x.shape)
    for ii in range(0, len(x)):
        weight[ii] = _invDist(dx[ii], dy[ii])
    
    return np.average(z, weights = weight)


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

