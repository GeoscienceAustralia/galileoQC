#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 16:43:31 2020

@author: markdransfield

A library of geophysical QA/QC functions for line data.

References:
    C. Aiken, J. M. Brozena, B. Coakley, D. Dater, G. Flanagan, R. Forsberg, 
    J. Kellogg, R. Kucks, X. Li, A. Mainville, R. Morin, M. Pilkington, 
    D. Plouff, D. Roman, J. Urrutia-Fucugauchi, M. V ́eronneau, and D. Winester. 
    New standards for reducing gravity data: The North American gravity 
    database. Geophysics, 70(4):J25, 2005.
    
    C. Jekeli, Theoretical fundamentals of airborne gradiometry. In Airborne 
    Gravity for Geodesy Summer School, 23-27 May, 2016.
"""

"""
Simplistic XYZ_TO_HDF5 converter for airborne gravity data

(c) Mark Dransfield 19 Jul 2020

"""
# import necessary modules
import numpy as np
import h5py
import matplotlib.pyplot as plt
#from scipy.signal import butter, lfilter
from pathlib import Path
import pathlib
import aseg_gdf2 as aseg
from scipy.interpolate import CloughTocher2DInterpolator
from scipy import interpolate
from src import gridfiles as erm
from src import qualityAnalysis as qa
from src import config
import filebrowser as fb

groupName = config.groupName
projectName = config.projectName

def getWhizzData(whizzFile, line, channel):
    '''
    Returns a numpy array containg the specified channel of
    data for the given line.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.

    Returns
    -------
    my_data : numpy array
        DESCRIPTION.

    '''
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        my_data = np.array(g[line][channel])
            
    return my_data


def updateLineAttributes(whizzFile, line_type='', line='', planned_line=0):
    '''
    For each line group, use the line_type field to set the associated planned
    line number, the segment number, and the reflight number for the line.

    Parameters
    ----------
    whizzFile : TYPE
        DESCRIPTION.
    line_type : TYPE, optional
        DESCRIPTION. The default is ''.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)

    with h5py.File(filename, 'r+') as f:
        # all the lines are in g:
        g = f[groupName]['Lines']

        if line_type == 'Xcal_nsw' or line_type == 'Xcal_can' or line_type == 'SGL_GA':
            print(f'\nSetting Line attributes for {whizzFile.name} according to the {line_type} scheme.')
            print('  {:<14} {:<14} {:<14} {:<14} '.format('Line No.','Plan Line No.', 'Segment No.', 'Re-flight No.'))
            for line in g:
                gg = g[line]
                current_line = gg.attrs['LineNumber']
                if line_type == 'Xcal_nsw':
                    if current_line < 8999999.0:
                        gg.attrs['PlannedLine'] = np.floor(current_line / 100.0) * 10.0
                        gg.attrs['Segment'] = 0
                        gg.attrs['ReflightNumber'] = int(current_line - np.floor(current_line / 10.0) * 10)
                    else:
                        gg.attrs['PlannedLine'] = np.floor(current_line / 10000.0)
                        gg.attrs['Segment'] = 0
                        gg.attrs['ReflightNumber'] = int(current_line - np.floor(current_line / 10000.0) * 10000)
                elif line_type == 'Xcal_can':
                    if current_line < 8999999.0:
                        gg.attrs['PlannedLine'] = np.floor(current_line / 10.0) * 10.0
                        gg.attrs['Segment'] = 0
                        gg.attrs['ReflightNumber'] = int(current_line - np.floor(current_line / 10.0) * 10)
                    else:
                        gg.attrs['PlannedLine'] = np.floor(current_line / 10000.0)
                        gg.attrs['Segment'] = 0
                        gg.attrs['ReflightNumber'] = int(current_line - np.floor(current_line / 10000.0) * 10000)
                elif line_type == 'SGL_GA':
                    if current_line < 7000:
                        gg.attrs['PlannedLine'] = np.floor(current_line * 10.0) / 10.0
                        current_segment = int(np.round(10 * (current_line - np.floor(current_line))))
                        gg.attrs['Segment'] = current_segment
                        gg.attrs['ReflightNumber'] = int(np.round(100 * (current_line - np.floor(current_line)))
                                                         - 10 * current_segment)
                    else:
                        gg.attrs['PlannedLine'] = np.floor(current_line)
                        gg.attrs['Segment'] = 0
                        gg.attrs['ReflightNumber'] = int(100 * (current_line - np.floor(current_line)))
                print('  {:<14} {:<14} {:<14} {:<14} '.\
                    format(current_line, gg.attrs['PlannedLine'], gg.attrs['Segment'], gg.attrs['ReflightNumber']))
        elif line != '' and planned_line != '':
            gg = g[line]
            gg.attrs['PlannedLine'] = planned_line
        else:
            print('NO ACTION TAKEN - line_type was neither "Xcal_nsw" or "SGL_GA" and planned_line = 0.')
    return


def updateChannelAttributes(whizzFile, channel, name='', units='', alias='', description='', chan_precision=-1):
    '''
    Updates the channel attributes for all lines in the geoWhizz HDF5 file. For any
    attribute, the default is to not change its value.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
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

    '''
    
    # In the hdfFile, for all lines, update the attributes of the named channel if,
    # and only if, the passed attribute is not empty.
    
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


def _distanceFlown(whizzFile, x = '', y = '', lines=[]):
    '''
    The total distance flown on the survey.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    x : String, optional
        The name of the channel of X positions (in metres). The default is 'X'.
    y : String, optional
        The name of the channel of Y positions (in metres). The default is 'Y'.
    lines : array of strings, optional
        An array of line identifiers whose total distance will be returned.
        Default all lines in whizzFile.

    Returns
    -------
    count : Integer
        The total number of lines in the file.
    Float
        The total distance in km flown over all lines in the survey.

    '''
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']

        lineDistance = 0.0
        count = 0
    
        if lines == []:
            lines = list(g.keys())
        for line in lines:
            xPos = np.array(g[line][x])
            yPos = np.array(g[line][y])
            lineDistance += _lineLength(xPos, yPos)
            count += 1
            
    print(f'{count} lines: total distance flown [km] = {lineDistance/1000.0:,.1f}')
    return (count, lineDistance/1000.0)
    

def _lineLength(x, y):
    """
    Calculates the point-to-point length of a line defined by discrete (x,y) points.

    Parameters
    ----------
    x (numpy 1D array): a 1D array of x positions (units metres).

    y (numpy 1D array): a 1D array of y positions (units metres).

    Returns
    -------
    lenOfLine (Float) : the length of the line in metres.

    """
    lenOfLine = 0.0
    for ii in range(0, len(x)-1):
        lenOfLine += distance(x[ii] - x[ii+1], y[ii] - y[ii+1])
    return lenOfLine


def displayGridded(filename, x, y, channel):
    """
    Given survey lines in a geoWhizz HDF5 file, each with x,y position fields and
    a value field z, grid z and display as a located shaded image. Just a quick
    and dirty routine.

    Parameters
    ----------
    filename (String) : the name of a geoWhizz HDF5 file.

    x (String) : The name of the field of X values (in metres). The default is 'X'.

    y (String) : The name of the field of Y values (in metres). The default is 'Y'.

    channel (String) : The name of the field of Z values (in metres). The default is 'X'.

    Returns
    -------
    None.

    """
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        em = np.zeros((0,))
        nm = np.zeros((0,))
        zm = np.zeros((0,))

        for line in g.keys():
            em = np.append(em, np.array(g[line][x]), axis=0)
            nm = np.append(nm, np.array(g[line][y]), axis=0)
            zm = np.append(zm, np.array(g[line][channel]), axis=0)

    numcells = (max(em) - min(em))/1000.0
    E = np.linspace(min(em), max(em), numcells.astype(np.int))
    numcells = (max(nm) - min(nm))/1000.0
    N = np.linspace(min(nm), max(nm), numcells.astype(np.int))
    E, N = np.meshgrid(E, N)  # 2D grid for interpolation
    interp = CloughTocher2DInterpolator(list(zip(em, nm)), zm)
    Z = interp(E, N)
    
    
    erm.displayShadedImage(E, N, Z, "hello")


def InterpolateGridOntoLine(gridPath, hdfPath, lines=[]):
    '''
    Interpolates the data in a grid file onto a new channel, with the same name
    as the gridPath stem, in the geoWhizz HDF5 file. Fills empty channel samples
    with nans.

    Parameters
    ----------
    gridPath : String or pathlib.PosixPath
        Name of a ERMapper file, including path and extension.
    hdfFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.

    Returns
    -------
    None.

    '''
    gridFile = str(gridPath)
    hdfFile = str(hdfPath)

    # retrieve the grid information
    eg, ng, zg, datum, projection = erm.read_ers_image(gridPath)
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
            em = np.array(g[line][x])
            nm = np.array(g[line][y])
            time = np.array(g[line][t])
            zm = np.empty(em.shape)
            zm[:] = np.nan
            dist = np.empty(em.shape)
            dist[:] = np.nan
            sampleSpacing = distance(em[1]-em[0], nm[1]-nm[0])
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

                ztemp[kz] = weightedAverage(x, y, z, em[kk], nm[kk])
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
           
            fail, numExc = qa.failsDeviation(zm - g[line][height], 20.0, 13)
            if fail and numExc > 13:
 #           if np.abs(np.max(g[line]['altitude'] - zm)) > 20.0:
                plotTime = time - time[0]
                fig = plt.figure()
                ax1 = fig.add_subplot(2,1,1)
                ax1.plot(plotTime, g[line]['altitude'], plotTime, zm)
                ax1.set_xlim(plotTime[0], plotTime[-1])
                plotTitle = f'Line {line}, altitude and planned drape height'
                plt.title(plotTitle, fontsize = 8)
                plt.grid(True)
                for label in ax1.get_xticklabels(): label.set_fontsize(6)
                for label in ax1.get_yticklabels(): label.set_fontsize(6)
                
                ax2 = fig.add_subplot(2,1,2)
                ax2.plot(plotTime, g[line]['altitude'] - zm)
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
    '''
    Interpolates dataIn, sampled at timeIn, onto the samples timeOut.
    These three input arrays are pre-processed to ensure that timeIn
    and timeOut are monotonically increasing, whilst keeping dataIn
    synchronised with timeIn, and spareOut synchronised with timeOut.

    Parameters
    ----------
    timeIn : 1D numpy float array
        DESCRIPTION.
    dataIn : 1D numpy float array
        DESCRIPTION.
    timeOut : 1D numpy float array
        DESCRIPTION.
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

    '''
    
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
        print(f'')
        plt.plot(timeIn[1:], np.diff(dataIn), 'b', timeIn[1:], np.diff(timeIn), 'g')
        plt.show()
    spl = interpolate.splrep(timeIn_trim, dataIn_trim, k=3, s=0)
    # out = interpolate.splev(timeOut, spl)
    out = interpolate.splev(timeOut_trim, spl)
    return out, spareOut


def _trim_monotonic(data_in, sync=[]):
    '''
    Force input to be monotonic by removing samples; keep
    sync aligned by removing corresponding samples there 
    as well.

    '''
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


def distance(x, y):
    return np.sqrt(x * x + y * y)


def invDist(x, y):
    return 1.0 / distance(x, y)


def weightedAverage(x, y, z, xo, yo):
    x = x[np.logical_not(np.isnan(z))]
    y = y[np.logical_not(np.isnan(z))]
    z = z[np.logical_not(np.isnan(z))]
    dx = x - xo
    dy = y - yo
    weight = np.zeros(x.shape)
    for ii in range(0, len(x)):
        weight[ii] = invDist(dx[ii], dy[ii])
    
    return np.average(z, weights = weight)


def getLineXChannel(whizzFile, line, x, channel):
    '''
    Returns 1D numpy arrays of x and channel in line from geoWhizz filename. The
    inputs are just two channels in the data and the names 'x' and 'channel' do
    not carry intrinsic meaning.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    line : String
        A flightline, e.g. '1000110.0'.
    line : TYPE
        DESCRIPTION.
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

    '''
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        xData = np.array(g[line][x])
        yData = np.array(g[line][channel])
    return xData, yData

    
def plotChannelLines(whizzFile, channel, flightLines, x, xOffset=True):
    '''
    For the given channel in the geoWhizz file filename, plot the values for all
    specified flightlines versus x.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    channel : String
        The name of the channel or field to plot.
    flightLines : [String]
        A list of flightlines, e.g. ['1000110.0', '1000210.0', '1000310.0'] .
    x : String
        The name of the independent variable for the plot.
    xOffset : Bool, optional
        If True, the x data will be offset to start at zero. The default is True.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        plotTitle = f'{projName} : {channel}'
        xDel = 0.0
        
        for line in flightLines:
            yData = np.array(g[line][channel])
            xData = np.array(g[line][x])
            if xOffset and line == flightLines[0]:
                xDel = xData[0]
            xData = xData - xDel
            myPlot, = ax.plot(xData, yData, color='blue', lw=0.3)
            plotTitle += f', L{line}'
            
    plt.xlabel(x, fontsize = 6)
    plt.ylabel(channel, fontsize = 6)
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    plt.show()


def reportWhizz(whizzFile, line='', channel=''):
    '''
    Prints a short summary of the data in a HDF5 Whizz file.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    line : String, optional
        The line number, formatted as a string, to report in detail. The default is '' and no line details.
    channel : String, optional
        The channel or field name to report in detail. The default is '' and no channel details.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)
        
    with h5py.File(filename, 'r') as f:
        
        whizzHeader = list(f.keys())[0]
        g = f[whizzHeader]
        
        gAttributeNames = list(g.attrs)
        
        gCoord = g['CoordinateFrame']
        coordAttrs = list(gCoord.attrs)
        
        gLines = g['Lines']
        lineNames = list(gLines.keys())
        numLines = len(lineNames)
        lineGroups = list(gLines.values())
        channelNames = list(lineGroups[0].keys())
        numChannels = len(channelNames)
        
        print(whizzHeader)
        for attribute in gAttributeNames:
            print(f'    {attribute:.20}: {g.attrs[attribute]}')
        print('\nCoordinates')
        for attribute in coordAttrs:
            print(f'    {attribute}: {gCoord.attrs[attribute]}')
        if line != '':
            myLineGroup = gLines[line]
            lineAttrs = list(myLineGroup.attrs)
            print(f'\nLine {myLineGroup}')
            for attribute in lineAttrs:
                print(f'    {attribute}: {myLineGroup.attrs[attribute]}')
        else:
            _ = _distanceFlown(whizzFile)
            print(f'\n{numLines} lines:\n', lineNames)

        print(f'\n{numChannels} channels:\n', channelNames)
        if channel != '':
            if line == '':
                line = lineNames[0]
                
            myChanGroup = gLines[line][channel]
            chanAttrs = list(myChanGroup.attrs)
            print(f'\nChannel {myChanGroup}')
            for attribute in chanAttrs:
                print(f'    {attribute}: {myChanGroup.attrs[attribute]}')


def reportFlights(whizzFile, flightChannel='FLIGHT', lines=[], detailed=False):
    '''
    Prints a summary of the flight numbers in a HDF5 Whizz file.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    flightChannel : String, optional
        The name of the channel containing the flight numbers. The default is 'FLIGHT'.
    lines : String Array, optional
        The array of line numbers, each formatted as a string, to report. The default is [] and all lines.
    detailed : Bool, optional
        If true, report the line numbers flown in each flight. The default is False.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)
        
    with h5py.File(filename, 'r') as f:
        
        whizzHeader = list(f.keys())[0]
        g = f[whizzHeader]
        gAttributeNames = list(g.attrs)
        gLines = g['Lines']

        # The data are stored with flights belonging to lines; we invert this relationship
        flight_dict = {}
        if lines == []:
            lines = list(gLines.keys())
        numLines = len(lines)
        for line in lines:
            this_flight = gLines[line][flightChannel][0]
            if this_flight in flight_dict:
                flight_dict[this_flight].append(line)
            else:
                flight_dict[this_flight] = [line]
        sorted_keys = sorted(flight_dict.keys())
        sorted_flights = {key:flight_dict[key] for key in sorted_keys}

        print(whizzHeader)
        for attribute in gAttributeNames:
            print(f'    {attribute:.20}: {g.attrs[attribute]}')

        print(f'\n{len(sorted_flights.keys())} flights over {numLines} lines.')

        print("\nFlights")
        for flight in sorted_flights:
            print(f'    {flight:.0f}')
            if detailed:
                print(f'        ', end = '')
                for line in sorted(sorted_flights[flight]):
                    print(f'L{line}', end = ' ')
                print('')


def reportSampling(whizzFile, timeChannel='', xChannel='', yChannel=''):
    filename = str(whizzFile)
        
    with h5py.File(filename, 'r') as f:
        
        whizzHeader = list(f.keys())[0]
        g = f[whizzHeader]
        gAttributeNames = list(g.attrs)
        gLines = g['Lines']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if timeChannel == '':
            timeChannel = f[groupName]['CoordinateFrame'].attrs['TimeChannel']


        time_deltas = []
        x_deltas = []
        y_deltas = []
        lines = list(gLines.keys())
        numLines = len(lines)
        for line in lines:
            time_deltas = np.append(time_deltas, np.diff(np.array(gLines[line][timeChannel])))
            x_deltas = np.append(x_deltas, np.diff(np.array(gLines[line][xChannel])))
            y_deltas = np.append(y_deltas, np.diff(np.array(gLines[line][yChannel])))
        mean_dt = np.mean(time_deltas)
        min_dt = np.min(time_deltas)
        max_dt = np.max(time_deltas)
        dd = qa._distance(x_deltas, y_deltas)
        mean_dd = np.mean(dd)
        min_dd = np.min(dd)
        max_dd = np.max(dd)

        print(whizzHeader)
        for attribute in gAttributeNames:
            print(f'    {attribute:.20}: {g.attrs[attribute]}')

        print(f'\nSample time and distance statistics')
        print(f'  Min  = {min_dt:.3f} s, {min_dd:.1f} m')
        print(f'  Max  = {max_dt:.3f} s, {max_dd:.1f} m')
        print(f'  Mean = {mean_dt:.3f} s, {mean_dd:.1f} m')


def updateProject(whizzFile, projectName='', blockID='', acquirer='', acquirerProjectID='', reportName=''):
    '''
    Change any of the project attributes in the HDF5 Whizz file. Typically most of this information
    is not available when the Whizz file is created and this routine is used to add it when it is
    available.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
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

    '''
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
    '''
    Change any of the attributes of the coordinate frame for the survey. This includes the names of the x,y,z,t coordinate fields
    and their respective datums. Typically these attributes are not entered to the Whizz file at creation and must be added later,
    by this routine.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
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

    '''
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

def asegReportChannels(datFilePath):
    '''
    Prints out the indices to the first channel name containg each of ['line',
    'flight', 'date', 'zone'] and then a list of all channel or field names in
    the given ASEG-GDF2 data. Useful when you want to know which channels to
    omit when creating a geoWhizz HDF5 file (asegToHDF()) or when you want to 
    know the correct name of a particular set of channels to report their
    values (asegReportFirst()).

    Parameters
    ----------
    datFilePath : pathlib.PosixPath
        The pathlib Path to the input ASEG-GDF2 file.

    Returns
    -------
    None.

    '''
    # open GDF, pull out channel names, the units, and the description
    gdf = aseg.read(str(datFilePath))
    channelNames = gdf.field_names()
    channelsOut = channelNames # these are the channels we will save
        
    # next, identify the column containing line numbers etc
    # we won't save these columns but do want the first value for each line
    # so re-get the indices.
    lineIdx = index_containing_substring(channelsOut, 'line')
    flightIdx = index_containing_substring(channelNames, 'flight')
    dateIdx = index_containing_substring(channelNames, 'date')
    zoneIdx = index_containing_substring(channelNames, 'zone')
    print(f' Indices - line {lineIdx}, flight {flightIdx}, date {dateIdx}, zone {zoneIdx}')

    print(' Channels: ')
    print(channelNames)

    return    


def asegReportFirst(datFilePath, channels):
    '''
    When converting an ASEG-GDF2 file to geoWhizz HDF5, it is useful to know the
    values of certain key parameters, stored as fields in the ASEG-GDF2 flat
    ASCII data table. This routine prints out the value for the first record
    of the ASEG-GDF2 file for each specified channel in channels. (See asegReportChannels())
    to generate a list of all channel names).

    Parameters
    ----------
    datFilePath : pathlib.PosixPath
        The pathlib Path to the input ASEG-GDF2 file.
    channels :  : [String]
        An array of channel or field names to report on.

    Returns
    -------
    None.

    '''
    # open GDF, pull out channel names, the units, and the description
    gdf = aseg.read(str(datFilePath))
    channelNames = gdf.field_names()
    chunk = gdf.df_chunked(chunksize=1).get_chunk()
    dataArray = np.array(chunk.values)
     
    print('First Record of selected channels')
    for channel in channels:
        chanIdx = index_containing_substring(channelNames, channel.lower())
        print(f'    {channel}: {dataArray[0, chanIdx]}')
        
    return    


def asegToHDF(datFilePath, outputHdf = '', omitChannels=['RT'], dontSave=False):
    '''
    Reads the data from the ASEG-GDF2 survey file and writes it to a new Whizz
    HDF5 survey file. Uses the aseg_gdf2 package by Kent Inverarity at:
        https://github.com/kinverarity1/aseg_gdf2

    Parameters
    ----------
    datFilePath : pathlib.PosixPath
        The pathlib Path to the input ASEG-GDF2 file.
    outputHdf : pathlib.PosixPath, optional
        The pathlib Path to the output Whizz HDF5 file. The default is '' and
        the output file is then the same as the input file with the extension
        changed to '.hdf5'.
    omitChannels : [String]
        An array of channel or field names to omit from the saved geoWhizz HDF5 file.
    dontSave : Bool
        If true, the data are not saved to an output HDF5 file.

    Returns
    -------
    None.

    '''

    if outputHdf == '':
        outputHdf = datFilePath.with_suffix('.hdf5')
    
    # open GDF, pull out channel names, the units, and the description
    gdf = aseg.read(str(datFilePath), method='fixed-widths')
    channelNames = gdf.field_names()
    channelsOut = channelNames # these are the channels we will save
        
    # next, identify the column containing line numbers etc
    # we won't save these columns but do want the first value for each line
    # so re-get the indices.
    lineIdx = index_containing_substring(channelsOut, 'line')
    flightIdx = index_containing_substring(channelNames, 'flight')
    dateIdx = index_containing_substring(channelNames, 'date')
    zoneIdx = index_containing_substring(channelNames, 'zone')
    print(f' Indices - line {lineIdx}, flight {flightIdx}, date {dateIdx}, zone {zoneIdx}')

    # now if index not found - error/warning; else pop off the channelName list
    if omitChannels == []:
        if lineIdx < 0:
            print('ERROR - lineIdx not found')
            return
        else:
            channelsOut.pop(lineIdx)
        tempIdx = index_containing_substring(channelsOut, 'flight')
        if flightIdx < 0:
            print('WARNING - flightIdx not found')
        else:
            channelsOut.pop(tempIdx)
        tempIdx = index_containing_substring(channelsOut, 'date')
        if dateIdx < 0:
            print('WARNING - dateIdx not found')
        else:
            channelsOut.pop(tempIdx)
        tempIdx = index_containing_substring(channelsOut, 'zone')
        if zoneIdx < 0:
            print('WARNING - zoneIdx not found')
        else:
            channelsOut.pop(tempIdx)
    else:
        for channel in omitChannels:
            tempIdx = index_containing_substring(channelsOut, channel.lower())
            if tempIdx < 0:
                print(f'WARNING - {channel} to omit not found')
            else:
                channelsOut.pop(tempIdx)
            
    print('channels out: ')
    print(channelsOut)
    
    # first, some initialisations
    chunkSize = 10 #000
    lastLineNumber = 0
    numLines = 0
    tempLineSizes = []
    lineNumbers = []
    flightNumbers = []
    dates = []
    zones = []
    numRecordsInLine = 0
    
    # now, iterate through GDF collecting line numbers and sizes
    print('\nGetting line numbers and lengths ...')
    for chunk in gdf.df_chunked(chunksize=chunkSize):
        dataArray = np.array(chunk.values)
        linesInChunk = dataArray[:, lineIdx]
        for i, v in enumerate(linesInChunk):
            # if the line number is new, create a new line sub-group and populate its attributes
            if v != lastLineNumber:
                tempLineSizes.append(numRecordsInLine)
                lineNumbers.append(v)
                if flightIdx >= 0:
                    flightNumbers.append(dataArray[i, flightIdx])
                if dateIdx >= 0:
                    dates.append(dataArray[i, dateIdx])
                if zoneIdx >= 0:
                    zones.append(dataArray[i, zoneIdx])
                numRecordsInLine = 1
                numLines += 1
                lastLineNumber = v

            else:
                numRecordsInLine += 1

    print('... done. Ready to create HDF file and transfer data into it.\n')

    tempLineSizes.append(numRecordsInLine)
    tempStarts = np.array(tempLineSizes)
    lineStarts = np.array(tempStarts.cumsum())
    lineSizes = np.array(tempLineSizes[1:])
    # maxLineLen = np.max(lineSizes) # longest line size is the size of the temp array for storing a line of data
    
    if len(lineSizes) != numLines:
        print(f'lineSizes.count = {len(lineSizes)} != numLines = {numLines}')

    with h5py.File(str(outputHdf), 'w') as f:
        # create all the data structure ready for the datasets
        g = f.create_group(groupName)
        g.attrs['ProjectName'] = projectName
        
        gCoord = g.create_group('CoordinateFrame')
        gLines = g.create_group('Lines')
        if zoneIdx >= 0:
            gCoord.attrs['UTMZone'] = zones[0]
        
        for ii in range(0, len(lineNumbers)):
            
            # create a line group and metadata
            gg = gLines.create_group(f'{lineNumbers[ii]}')
            gg.attrs['LineNumber'] = lineNumbers[ii]
            if flightIdx >= 0:
                gg.attrs['Flight'] = flightNumbers[ii]
            if dateIdx >= 0:
                gg.attrs['Date_Local'] = dates[ii]
            gg.attrs['NumberOfFids'] = lineSizes[ii]

            # create the DataSets with attributes
            for channelName in channelsOut:
                dd = gg.create_dataset(channelName, (lineSizes[ii],), dtype='float64', compression="gzip", compression_opts=4)
                dd.attrs['Name'] = channelName
                dd.attrs['Alias'] = channelName
                fieldDef = gdf.get_field_definition(channelName)
                dd.attrs['Units'] = cleanUnits(fieldDef['unit'])
                dd.attrs['Description'] = fieldDef['long_name']

        # now add the data sets, extracting from chunks as they come
        numChunk = 0
        lineIndex = 0
        jStart = 0
        iStart = 0
        printSum = False
        
        print('\nTransferring data ...')
        for chunk in gdf.df_chunked(chunksize=chunkSize):
            if lineIndex >= len(lineNumbers):
                break
            nextLineEnd = lineStarts[lineIndex+1] - numChunk * chunkSize
            blockLength = np.min([chunkSize, nextLineEnd])
            
            # transfer a block of data into each line within the chunk & prepare for next line
            while nextLineEnd < chunkSize:
                if printSum:
                    print(f'is={iStart}, js={jStart}, len={blockLength}, chunk no.={numChunk}')
                    print(f'Line: {lineNumbers[lineIndex]}, ch10= {channelsOut[10]}')
                    print(chunk.to_string())
                    
                transferBlock(lineNumbers[lineIndex], channelsOut, blockLength, iStart, chunk, jStart, gLines)
                lineIndex += 1
                if lineIndex > 7:
                    printSum = False
                if lineIndex >= len(lineNumbers):
                    break
                
                nextLineEnd = lineStarts[lineIndex+1] - numChunk * chunkSize
                iStart += blockLength
                jStart = 0
                blockLength = np.min([chunkSize, nextLineEnd]) - iStart
            
            if nextLineEnd == chunkSize:
                if printSum:
                    print(f'is={iStart}, js={jStart}, len={blockLength}, chunk no.={numChunk}')
                transferBlock(lineNumbers[lineIndex], channelsOut, blockLength, iStart, chunk, jStart, gLines)
                lineIndex += 1
                if lineIndex > 7:
                    printSum = False
                if lineIndex >= len(lineNumbers):
                    break
                iStart = 0
                jStart = 0
                
            elif nextLineEnd > chunkSize:
                if printSum:
                    print(f'is={iStart}, js={jStart}, len={blockLength}, chunk no.={numChunk}')
                transferBlock(lineNumbers[lineIndex], channelsOut, blockLength, iStart, chunk, jStart, gLines)
                iStart = 0
                if blockLength == chunkSize:
                    jStart += blockLength
                else:
                    jStart = blockLength
                
            numChunk += 1

    return


def cleanUnits(fieldDef):
    '''
    ASEG-GDF2 files generated by Atlas (proprietary geophysical processing software
    from Fugro Airborne, CGG Airborne and XCalibur) does not obey the ASEG-GDF2
    standard. One issue is the use of an extra colon where the standard requires
    a comma, making the units field incorrect. This routine takes a units fieldDef
    and fixes the error by removing the substring from the extra colon. 

    Parameters
    ----------
    fieldDef : String
        The field units from the dictionary fieldDef['units'].

    Returns
    -------
    String
        The field units after correction.

    '''
    return fieldDef.split(':')[0]


def transferBlock(line, channelsOut, length, iS, chunk, jS, gLines):
    '''
    A 'block' of data is all the data in the chunk that is part of the line. This
    routine copies that data from chunk to its correct position in gLines.

    Parameters
    ----------
    line : TYPE
        The survey line to which the data belongs.
    channelsOut : [String]
        The list of all the channels for which the data is to be extracted.
    length : Integer
        The number of records to be transferred to the data block array.
    iS : Integer
        The index in the chunk to the 1st record to be transferred.
    chunk : pandas data chunk from ASEG-GDF2
        A chunk of data from the ASEG-GDF2 data file.
    jS : Integer
        The index of the first record in the HDF datasets to where the data is
        to be transferred.
    gLines : HDF Group
        The geoWhizz HDF group containg the array of 'Lines' subGroups.

    Returns
    -------
    None.

    '''
    lineStr = f'{line}'
    # float_formatter = lambda x: "%.2f" % x
    for channelName in channelsOut:
        try:
            gLines[lineStr][channelName][jS:jS+length] = np.array(chunk[channelName][iS:iS+length])
        except:
            print(channelName)
            print(chunk[channelName][iS:iS+length])
                
        # if channelName == 'time':
        #     print(f'Data Copy - first 3 time values:')#' {chunk[channelName][iS:iS + 3]}')
        #     with np.printoptions(formatter={'float_kind':float_formatter}):
        #         print(f'{np.array(chunk[channelName][iS:iS + 3])}')
        #         print(f'{np.array(gLines[lineStr][channelName][jS:jS + 3])}')
    return
    
    
def index_containing_substring(the_list, substring):
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
    lowlist = [name.lower() for name in the_list]
    for i, s in enumerate(lowlist):
        if substring in s:
              return i
    return -1
    '''
    Reads the data from the ASEG-GDF2 survey file and returns an xarray.
    Uses the aseg_gdf2 package by Kent Inverarity at:
        https://github.com/kinverarity1/aseg_gdf2

    Parameters
    ----------
    datFilePath : pathlib.PosixPath
        The pathlib Path to the input ASEG-GDF2 file.
    omitChannels : [String]
        An array of channel or field names to omit from the saved geoWhizz HDF5 file.

    Returns
    -------
    None.

    '''

    # open GDF, pull out channel names, the units, and the description
    gdf = aseg.read(str(datFilePath))
    channelNames = gdf.field_names()
    channelsOut = channelNames # these are the channels we will save
    
    # now if index not found - error/warning; else pop off the channelName list
    if ~omitChannels == []:
        for channel in omitChannels:
            tempIdx = index_containing_substring(channelsOut, channel.lower())
            if tempIdx < 0:
                print(f'WARNING - {channel} to omit not found')
            else:
                channelsOut.pop(tempIdx)
            
    print('channels out: ')
    print(channelsOut)
    
    # first, some initialisations
    chunkSize = 10000
    
    print('... done. Ready to create HDF file and transfer data into it.\n')
    
    with h5py.File(str(outputHdf), 'w') as f:
        # create all the data structure ready for the datasets
        g = f.create_group(groupName)
        g.attrs['ProjectName'] = projectName
        
        gCoord = g.create_group('CoordinateFrame')
        gLines = g.create_group('Lines')
        
        # now add the data sets, extracting from chunks as they come
        numChunk = 0
        jStart = 0
        iStart = 0
        printSum = True
        
        print('\nTransferring data ...')
        for chunk in gdf.df_chunked(chunksize=chunkSize):
            if lineIndex >= len(lineNumbers):
                break
            nextLineEnd = lineStarts[lineIndex+1] - numChunk * chunkSize
            blockLength = np.min([chunkSize, nextLineEnd])
            
            # transfer a block of data into each line within the chunk & prepare for next line
            while nextLineEnd < chunkSize:
                if printSum:
                    print(f'is={iStart}, js={jStart}, len={blockLength}, chunk no.={numChunk}')
                transferBlock(lineNumbers[lineIndex], channelsOut, blockLength, iStart, chunk, jStart, gLines)
                lineIndex += 1
                if lineIndex > 7:
                    printSum = False
                if lineIndex >= len(lineNumbers):
                    break
                
                nextLineEnd = lineStarts[lineIndex+1] - numChunk * chunkSize
                iStart += blockLength
                jStart = 0
                blockLength = np.min([chunkSize, nextLineEnd]) - iStart
            
            if nextLineEnd == chunkSize:
                if printSum:
                    print(f'is={iStart}, js={jStart}, len={blockLength}, chunk no.={numChunk}')
                transferBlock(lineNumbers[lineIndex], channelsOut, blockLength, iStart, chunk, jStart, gLines)
                lineIndex += 1
                if lineIndex > 7:
                    printSum = False
                if lineIndex >= len(lineNumbers):
                    break
                iStart = 0
                jStart = 0
                
            elif nextLineEnd > chunkSize:
                if printSum:
                    print(f'is={iStart}, js={jStart}, len={blockLength}, chunk no.={numChunk}')
                transferBlock(lineNumbers[lineIndex], channelsOut, blockLength, iStart, chunk, jStart, gLines)
                iStart = 0
                if blockLength == chunkSize:
                    jStart += blockLength
                else:
                    jStart = blockLength
                
            numChunk += 1
    return


def xyzToHDF(xyzFilePath = '', hdfFileName = '', projectName = ''):
    '''
    Read in a Geosoft XYZ file and write the contents to a Whizz HDF5 file.

    Parameters
    ----------
    xyzFilePath : pathlib.PosixPath
        The pathlib Path to the Geosoft XYZ file.
    hdfFileName : pathlib.PosixPath, optional
        The pathlib Path to the Whizz HDF5 file to be created. The default is ''
        leaving the Whizz file to have the same path as the input file but with
        a ".HDF5" extension.
    projectName : String, optional
        The name of the survey or project. The default is ''.
    line_scale : float, optional
        All line numbers in the XYZ file are multiplied by line_scale before
        writing to the HDF5 file.

    Returns
    -------
    None.

    '''
        
    if xyzFilePath == '':
        xyzFilePath = fb.get_grid_filename()

    if xyzFilePath.__class__ == pathlib.PosixPath:

        xyzFile = str(xyzFilePath)
        if hdfFileName == '':
            hdfFileName = xyzFilePath.with_suffix('.hdf5')

    # access the data via Geosoft XYZ
    geosoftXYZ = readXYZ(xyzFile)
    desiredFieldNames = [*geosoftXYZ[0]]
    num_lines = len(geosoftXYZ)
    lines = np.zeros((num_lines,))

    for count in range(0, num_lines):
        lines[count] = geosoftXYZ[count]['line_number']
        
    print('Creating: ', str(hdfFileName))
        
    with h5py.File(hdfFileName, 'w') as f: # better data file
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
                # dd.attrs['Units'] = ... , eg metres, eotvos, etc would be nice TODO
            
            #print(gg, lineCount)
            
            # create next line group and metadata
            gg = gLines.create_group(f'{lines[lineCount + 1]}')
            gg.attrs['LineNumber'] = lines[lineCount + 1]
        
        # put data in last line subgroup's dataset
        print('  About to add data for line ', lines[len(lines)-1], ' Lcount ', len(lines), '/', len(lines), '\n')
        for count in range(1, len(desiredFieldNames)):
           dataArray = geosoftXYZ[len(lines)-1][desiredFieldNames[count]]
           #print(dataArray.shape)
           dd = gg.create_dataset(desiredFieldNames[count], \
                                  data = dataArray, dtype='float64', \
                                      compression="gzip", compression_opts=4)
           dd.attrs['Name'] = desiredFieldNames[count]
           # dd.attrs['Units'] = ... , eg metres, eotvos, etc would be nice TODO
        
    return hdfFileName


def read_xyz_header(filename):
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

    
def readXYZ(filename):
    """ Open a Geosoft XYZ file, read in the contents and return a dictionary
    array where each element corresponds to one flight line and contains the
    line number, flight number, date, and data for each channel.
    
    TODO: remove reliance on Flight, Date and Line fields if possible.
    TODO: replace the entire routine with a xyzToHdf() function because this function
            will run out of memory if the XYZ file is large.

    Mark Dransfield
    Jul 2020
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
                print('  Found fields (channels):')
                for ii in range(0, len(channelnames)):
                    print(f'    {channelnames[ii]} - precision {field_precisions[ii]}')
                break
            if header_rec > num_head_recs:
                print(f"Error - can't find header record with {num_channels} channel names.")
                break
    fid.close()
    
    # rename key channels TODO : is this necessary now that whizz Files have 'XChannel' etc?
    for jj in range(num_channels):
        if channelnames[jj] == 'x' or channelnames[jj] == 'X' or channelnames[jj] == 'easting' or channelnames[jj] == 'Easting':
            channelnames[jj] = 'X'
        if channelnames[jj] == 'y' or channelnames[jj] == 'Y' or channelnames[jj] == 'northing' or channelnames[jj] == 'Northing':
            channelnames[jj] = 'Y'

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
                line_nos[line_ctr] = f'{current_line:.2f}'
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
                geosoftXYZ[line_ctr]['line_number'] = f'{current_line:.2f}'#file_line.split()[1]
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
                    geosoftXYZ[line_ctr-1][channelnames[ii]][fid_ctr] = \
                        float(line_str[ii])

                fid_ctr += 1

    fid.close()

    return geosoftXYZ


def gdb2xyz(filename):
    import subprocess as sbp


    inputGdb = "E:/2020_salamander_bell/deliverables/lawin/processed/databases/databases/lawin_airftg_deliverable.gdb"
    outputGdb = "E:/2020_Salamander_Bell/Deliverables/Lawin/Processed/Databases/Databases/test.adb"
    maxChannels = 250
    maxLines = 1000
    maxBlobs = 1900                 
    crdict = dict(program = "ExportGDB",
                  version = "1.0.0",
                  inputgdb = inputGdb,
                  outputgdb = outputGdb,
                  outputmode = "Append",
                  gdbcompression = "1",
                  gdbMaxChannels = maxChannels,
                  gdbMaxLines = maxLines,
                  gdbMaxBlobs = maxBlobs,
                  channelIoCsv = "",
                  outputChannelPrefix = "",
                  outputChannelSuffix = "",
                  outputFidRangeOption = "0",
                  maskChannel = "",
                  lffFile = "",
                  dataResample = "No",
                  resampledBy = "Fid",
                  newFidIncr = "1",
                  distChannel = "",                                                                                                       
                  newDistIncr = 50)                                                                                                    
    
    crString = 'CONTROL_BEGIN\n'
    crString += f'Program              = {crdict["program"]}\n'
    crString += f'Version              = {crdict["version"]}\n'
    crString += f'InputGDB             = {crdict["inputgdb"]}\n'
    crString += f'OutputGDB            = {crdict["outputgdb"]}\n'
    crString += f'OutputMode           = {crdict["outputmode"]}\n'
    crString += f'GDBCompression       = {crdict["gdbcompression"]}\n'
    crString += f'GDBMaxChannels       = {crdict["gdbMaxChannels"]}\n'
    crString += f'GDBMaxLines          = {crdict["gdbMaxChannels"]}\n'
    crString += f'GDBMaxBlobs          = {crdict["gdbMaxBlobs"]}\n'
    crString += f'ChannelIO_CSV        = {crdict["channelIoCsv"]}\n'
    crString += f'OutChanPrefix        = {crdict["outputChannelPrefix"]}\n'
    crString += f'OutChanSuffix        = {crdict["outputChannelSuffix"]}\n'
    crString += f'OutputFidRangeOption = {crdict["outputFidRangeOption"]}\n'
    crString += f'MaskChannel          = {crdict["maskChannel"]}\n'
    crString += f'LFFFile              = {crdict["lffFile"]}\n'
    crString += f'DataResample         = {crdict["dataResample"]}\n'
    crString += f'ResampledBy          = {crdict["resampledBy"]}\n'
    crString += f'NewFidIncr           = {crdict["newFidIncr"]}\n'
    crString += f'DistChannel          = {crdict["distChannel"]}\n'
    crString += f'NewDistIncr          = {crdict["distChannel"]}\n'
    crString += 'CONTROL_END\n'
    
    # test run:
    myPath = "E:/2020_Salamander_Bell/Deliverables/Lawin/Processed/Databases/Databases"
    crFile = myPath + "/gdb2adb.cr"
    
    with open(crFile, 'wt') as fid:
         #  if fid == -1
         #      fprintf('Error opening control file (%s)\n', crFile)
         #      error(message)
         #  end
         # now print out control file 
         fid.write(crString)
         
    fid.close()
    
    myCmd = "C:/Atlas/Atlas.exe " + crFile
    sbp.check_output(myCmd, shell=True).decode()


def renameChannels(nchannels, chanNames):
    '''
    

    Parameters
    ----------
    nchannels : TYPE
        DESCRIPTION.
    chanNames : TYPE
        DESCRIPTION.

    Returns
    -------
    chanNames : TYPE
        DESCRIPTION.

    '''
    
    for ii in range(0, nchannels): # TODO: expand this to cover a large range of possibilities
        if chanNames[ii] == 'x':
            chanNames[ii] = 'X'
            
    return chanNames

    
def addLineToWhizz(hdf5FileId, databaseName, lineNo, lineType, plannedX = [], plannedY = [], plannedZ = [], distanceUnits = ''):
    '''
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

    '''
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


