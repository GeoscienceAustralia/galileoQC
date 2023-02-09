# -3.402823466385289e+38
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 14 20:28:20 2021

@author: markdransfield

Each 'public' function is a "check" function which checks the data against
some QC requirement. Each `check` function should:
    1.  Report what it is doing.
    2.  Report each found error in the data against line number, specifying
        the fiducial, time or position where it occurred.
    3.  Report the number of lines that failed.
    Example:

Velocities not known - will calculate from positions
Nominal speed 67.5; allowed 54.0 : 81.0 for < 13.0 seconds.

25 lines had some short exceedance(s).
0 lines failed for exceedance > allowed distance.

TODO:
1. Trap instances where a channel or line name is not found as a dataset or group and suggest a check on spelling and case.
"""


# import necessary modules
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter
#import pathlib as Path
from src import gridfiles as erm
from src import config
from src import pointfiles as mhd
from src import whizzPlot as wpl

groupName = config.groupName

def diffGroundGrid(whizzFile, whizzChannel, whizzLine, gridPath, plot_title='Ground & Airborne Comparison'):
    '''
    Samples data from the grid file onto the line from whizzFile
    and displays those data and the whizzChannel on a plot.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    whizzChannel : String
        Name of the channel of line data to be compared.
    whizzLine : String
        The identifier (eg '10030.0') for the line from which the
        whizzChannel data are taken for comparison.
    gridPath : Path
        The path (PathLib) to the ERS gridfile from which data is
        to be extracted for comparison.
    plot_title : String, optional
        A tiotle for the plot. Default is 'Ground & Airborne Comparison'.

    Returns
    -------
    None.

    '''

    # retrieve the grid information
    eg, ng, zg, datum, projection = erm.read_ers_image(gridPath)
    ngmin = np.min(ng)
    ngmax = np.max(ng)
    ngd = np.abs(ng[1] - ng[0])
    egmin = np.min(eg)
    egmax = np.max(eg)
    egd = np.abs(eg[1] - eg[0])
    zg = zg[::-1, :]
    newChannelName = gridPath.stem
    
    print('\nGrid file read for channel ', newChannelName)
    print('  SW Corner = ',ngmin, egmin, '  NE Corner = ',ngmax, egmax, '. Spacings = ', ngd, egd)

    # retrieve the line information
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        y = f[groupName]['CoordinateFrame'].attrs['YChannel']
        em = np.array(g[whizzLine][x])
        nm = np.array(g[whizzLine][y])
        zm = np.array(g[whizzLine][whizzChannel])

    # Some useful summary info, particularly if the line is not inside
    # the grid area.
    lineNo = whizzLine
    lineText = 'Line ' + lineNo
    print(lineText)
    print('  North min, max ', np.min(nm), np.max(nm), '  East min, max ', np.min(em), np.max(em))

    # Check that the endpoints of the line are both within the grid area.
    i_float, i_disp = divmod(np.nanmin(nm) - ngmin, ngd)
    j_float, j_disp = divmod(np.nanmin(em) - egmin, egd)
    if i_float < 0 or j_float < 0 or i_float > len(ng) or j_float > len(eg):
        print('ERROR - the flight line is not wholly within the area of the grid')
        return
    i_float, i_disp = divmod(np.nanmax(nm) - ngmin, ngd)
    j_float, j_disp = divmod(np.nanmax(em) - egmin, egd)
    if i_float < 0 or j_float < 0 or i_float > len(ng) or j_float > len(eg):
        print('ERROR - the flight line is not wholly within the area of the grid')
        return

    # A point on the line is close to the projection of a grid cell centre if
    # it is closer than eps (half the sample spacing along the line).
    sample_spacing = _distance(em[1] - em[0], nm[1] - nm[0])
    eps = sample_spacing / 2.0

    # initialise the arrays to store the sampled data
    num_grid_samples = int(float(len(em)) * np.sqrt(2) * sample_spacing / max(ngd, egd))
    ns = np.zeros((num_grid_samples,))
    es = np.zeros((num_grid_samples,))
    zs = np.zeros((num_grid_samples,))
    diffs = np.zeros((num_grid_samples,))
    m = 0

    # traverse along the line extracting grid data as it occurs
    for k in range(0, len(em)):
        # find the closest grid cell centres (i,j), (i+1,j), (i,j+1), (i+1,j+1)
        i_float, i_disp = divmod(nm[k] - ngmin, ngd)
        j_float, j_disp = divmod(em[k] - egmin, egd)

        if i_float < 0 or j_float < 0 or i_float > len(ng) or j_float > len(eg):
            print('ERROR - the flight line is not wholly within the area of the grid')
            return

        # if the nearest grid cell centre is on the line, use the grid data
        if i_disp < eps and j_disp < eps:
            i0 = int(i_float)
            j0 = int(j_float)
            ns[m] = nm[k]
            es[m] = em[k]
            zs[m] = zg[i0,j0]
            diffs[m] = zs[m] - zm[k]
            m += 1
        # if the nearest grid cell centre is on the same northing, interpolate the grid data over easting
        elif i_disp < eps:
            i0 = int(i_float)
            j0 = int(j_float)
            j1 = j0 + 1
            ns[m] = nm[k]
            es[m] = em[k]
            zs[m] = zg[i0,j0] + (zg[i0,j1] - zg[i0,j0]) * (em[k] - eg[j0]) / (eg[j1] - eg[j0])
            diffs[m] = zs[m] - zm[k]
            m += 1
        # similarly ... for same easting, interpolate over northing
        elif j_disp < eps:
            i0 = int(i_float)
            i1 = i0 + 1
            j0 = int(j_float)
            ns[m] = nm[k]
            es[m] = em[k]
            zs[m] = zg[i0,j0] + (zg[i1,j0] - zg[i0,j0]) * (nm[k] - ng[i0]) / (ng[i1] - ng[i0])
            diffs[m] = zs[m] - zm[k]
            m += 1

    # clean up a couple of things before plotting
    if (np.ptp(em) > np.ptp(nm)):
        plotx_s = es[:m-1]
        plotx_m = em
    else:
        plotx_s = ns[:m-1]
        plotx_m = nm
    plotz_s = zs[:m-1]
    plotz_m = zm * 10.0
    # results ...
    print(f'Stdev(diff) = {np.nanstd(diffs):.1f}')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(plotx_s, plotz_s, 'b', label=newChannelName)
    ax.plot(plotx_m, plotz_m, 'g', label=lineText + ' ' + whizzChannel)
    ax.set_title(plot_title)
    ax.set_xlabel('Easting [m]')
    ax.set_ylabel('Bouguer Gravity [um/s/s]')
    ax.legend()
    ax.grid()


def _distance(x, y):
    return np.sqrt(x * x + y * y)


def checkGNSS(whizzFile, num_sats, pdop, vdop, hdop, nsats_min=4, max_pdop=6, max_vdop=4, max_hdop=4):
    '''
    Checks that the data in a whizzFile meets the requirements for the minimum
    number of satellites, and maximum PDOP, VDOP and HDOP.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    num_sats : String
        Name of the channel containing the number of satellites visible for
        each measurement.
    pdop : String
        Name of the channel containing the PDOP for each measurement.
    vdop : String
        Name of the channel containing the VDOP for each measurement.
    hdop : String
        Name of the channel containing the HDOP for each measurement.
    nsats_min : Integer, optional
        The minimum number of satellites required, default 4.
    max_pdop : Integer, optional
        The maximum PDOP allowed, default 6
    max_vdop : Integer, optional
        The maximum VDOP allowed, default 6
    max_hdop : Integer, optional
        The maximum HDOP allowed, default 6

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']

        error_count = 0
        for line in g.keys():
            min_nsats_data = np.nanmin(np.array(g[line][num_sats]))
            max_pdop_data = np.nanmin(np.array(g[line][pdop]))
            max_vdop_data = np.nanmin(np.array(g[line][vdop]))
            max_hdop_data = np.nanmin(np.array(g[line][hdop]))
            if min_nsats_data < nsats_min:
                print(f'Line {line} fail: min sats = {min_nsats_data:.0f} < {nsats_min}')
                error_count += 1
            if max_pdop_data > max_pdop:
                print(f'Line {line} fail: max pdop = {max_pdop_data:.0f} > {max_pdop}')
                error_count += 1
            if max_vdop_data > max_vdop:
                print(f'Line {line} fail: max vdop = {max_vdop_data:.0f} > {max_vdop}')
                error_count += 1
            if max_hdop_data > max_hdop:
                print(f'Line {line} fail: max hdop = {max_hdop_data:.0f} > {max_hdop}')
                error_count += 1
        print(f'In {projName}, checked num sats, PDOP, VDOP and HDOP. Found {error_count} errors.')
        

def checkLatCorr(whizzFile, latCorr, latitude=''):
    '''
    Subtracts the latitude correction in the data file from one calculated using
    Hinze et al (2005) and the latitude data in the data file. The min, max, mean
    and standard deviation of the difference is calculated for each line and
    presented in an extended whisker plot.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    latCorr : String
        The name of the geoWhizz field or channel containing the latitude correction
        (sometimes called "normal gravity").
    latitude : String, optional
        The name of the geoWhizz field or channel containing the latitude. The
        default is to read the latitude field name from the Coordinate Frame.

    Returns
    -------
    None.

    '''         
    
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if latitude == '':
            latitude = f[groupName]['CoordinateFrame'].attrs['LatitudeChannel']

        flightLine = list(g.keys())[0]
        corr_units = g[flightLine][latCorr].attrs['Units']
        if corr_units == 'mGal':
            unit_scale = 10.0
        elif corr_units == 'gu' or corr_units == 'µm/s/s':
            unit_scale = 1.0
        else:
            print('ERROR - correction units not recognised, expected mGal or µm/s/s (gu)')
            return

        numLines = len(g.items())
        diffMin = np.zeros((numLines,))
        diffMax = np.zeros((numLines,))
        diffMean = np.zeros((numLines,))
        diffStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            lat_data = np.array(g[line][latitude])
            cor_data = np.array(g[line][latCorr])
            if line == '8474.0':
                fig = plt.figure()
                ax = fig.add_subplot(3,1,1)
                ax.plot(cor_data * 10)
                plt.title('lat corr in data')
                ax = fig.add_subplot(3,1,2)
                ax.plot(lat_data)
                plt.title('Latitude')
                ax = fig.add_subplot(3,1,3)
                ax.plot(_normalGravity(lat_data))
                plt.title('My Estimate')
                
            err_data = cor_data * unit_scale + _normalGravity(lat_data)
            diffMin[count] = np.min(err_data)
            diffMax[count] = np.max(err_data)
            diffMean[count] = np.mean(err_data)
            diffStd[count] = np.std(err_data)
            lineNo[count] = line
            count += 1
        
        fig = plt.figure()
        fig.suptitle(f'check Latitude Correction - {projName}', fontsize=10)
        fig.subplots_adjust(top=0.85)
        
        ax = fig.add_subplot(1,1,1)
        ax.plot(lineNo, diffMin, 'bo', mfc='w')
        ax.plot(lineNo, diffMax, 'bo', mfc='w')
        ax.errorbar(lineNo, diffMean, diffStd, capsize=3, marker='s', c='blue', ls='None')
        plt.ylabel('Latitude Correction Error um/s/s', fontsize = 6)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        plt.show()
    return
    

def checkEotvosCorr(whizzFile, eotCorr, latitude='', x='', y='', GRS80_height='', time='', east_vel='', north_vel=''):
    '''
    Subtracts the eotvos correction in the data file from one calculated using
    Hinze et al (2005) and the latitude and position data in the data file.
    The min, max, mean and standard deviation of the difference is calculated
    for each line and presented in an extended whisker plot.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    eotCorr : String
        The name of the geoWhizz field or channel containing the eotvos correction.
    latitude : String, optional
        The name of the geoWhizz field or channel containing the latitude. The
        default is to read the latitude field name from the Coordinate Frame.
    x : String, optional
        The name of the geoWhizz field or channel containing the x position. The
        default is to read the xChannel field name from the Coordinate Frame.
    y : String, optional
        The name of the geoWhizz field or channel containing the y position. The
        default is to read the yChannel field name from the Coordinate Frame.
    GRS80_height : String, optional
        The name of the geoWhizz field or channel containing the GRS80_height. The
        default is to read the altitudeChannel field name from the Coordinate Frame.
    time : String, optional
        The name of the geoWhizz field or channel containing the time. The
        default is to read the timeChannel field name from the Coordinate Frame.
    east_vel : String, optional
        The name of the geoWhizz field or channel containing the velocity in the
        east direction. The default ('') is to calculate this from the x and time
        channels.
    north_vel : String, optional
        The name of the geoWhizz field or channel containing the velocity in the
        north direction. The default ('') is to calculate this from the y and time
        channels.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if latitude == '':
            latitude = f[groupName]['CoordinateFrame'].attrs['LatitudeChannel']
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if GRS80_height == '':
            GRS80_height = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        if time == '':
            time = f[groupName]['CoordinateFrame'].attrs['TimeChannel']            

        flightLine = list(g.keys())[0]
        corr_units = g[flightLine][eotCorr].attrs['Units']
        if corr_units == 'mGal':
            unit_scale = 10.0
        elif corr_units == 'gu' or corr_units == 'µm/s/s':
            unit_scale = 1.0
        else:
            print('ERROR - correction units not recognised, expected mGal or µm/s/s (gu)')
            return

        numLines = len(g.items())
        diffMin = np.zeros((numLines,))
        diffMax = np.zeros((numLines,))
        diffMean = np.zeros((numLines,))
        diffStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            lat_data = np.array(g[line][latitude])
            x_data = np.array(g[line][x])
            y_data = np.array(g[line][y])
            ht_data = np.array(g[line][GRS80_height])
            time_data = np.array(g[line][time])
            cor_data = np.array(g[line][eotCorr])
            if (east_vel == '')  | (north_vel == ''):
                (n_speed, e_speed) = _calc_speed(x_data, y_data, time_data)
            else:
                n_speed = np.array(g[line][north_vel])
                e_speed = np.array(g[line][east_vel])
            cal_data = _eotvosCorrection(e_speed, n_speed, lat_data, ht_data)
            err_data = cor_data * unit_scale + cal_data
            diffMin[count] = np.min(err_data)
            diffMax[count] = np.max(err_data)
            diffMean[count] = np.mean(err_data)
            diffStd[count] = np.std(err_data)
            lineNo[count] = line
            
            if np.abs(diffMin[count] - diffMax[count]) > 10.0:
                fig = plt.figure()
                fig.suptitle(f'{projName} L{lineNo[count]}', fontsize=10)
                fig.subplots_adjust(top=0.85)
                ax = fig.add_subplot(2,1,1)
                ax.plot(x_data, cor_data*unit_scale, x_data, -cal_data)
                ax = fig.add_subplot(2,1,2)
                ax.plot(x_data, err_data)
                
            count += 1
        
        fig = plt.figure()
        fig.suptitle(f'check Eotvos Correction - {projName}', fontsize=10)
        fig.subplots_adjust(top=0.85)
        
        ax = fig.add_subplot(1,1,1)
        ax.plot(lineNo, diffMin, 'bo', mfc='w')
        ax.plot(lineNo, diffMax, 'bo', mfc='w')
        ax.errorbar(lineNo, diffMean, diffStd, capsize=3, marker='s', c='blue', ls='None')
        plt.ylabel('Eotvos Correction Error [um/s/s]', fontsize = 6)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        plt.show()
    return


def _calc_speed(e, n, t):
    '''
    Returns the (e_speed, n_speed) velocity given the (e, n) positions as
    a function of t (time). Uses a simple numpy.diff approach.
 
    Parameters
    ----------
    e : Numpy 1D array
        The easting data.
    n : Numpy 1D array
        The northing data.
    t : Numpy 1D array
        The time data.

    Returns
    -------
    e_speed : Numpy 1D array
        The velocity in the east direction.
    n_speed : Numpy 1D array
        The velocity in the north direction.

    '''
    n_speed = np.diff(n) / np.diff(t)
    last = n_speed[-1]
    n_speed = np.append(n_speed, last)
    e_speed = np.diff(e) / np.diff(t)
    last = e_speed[-1]
    e_speed = np.append(e_speed, last)
    return (e_speed, n_speed)
    

def checkAtmosEffect(whizzFile, atmosCorr, GRS80_height=''):
    '''
    Subtracts the atomspheric correction in the data file from one calculated using
    Hinze et al (2005) and the height data in the data file.
    The min, max, mean and standard deviation of the difference is calculated
    for each line and presented in an extended whisker plot.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    atmosCorr : String
        The name of the geoWhizz field or channel containing the atmospheric correction.
    GRS80_height : String, optional
        The name of the geoWhizz field or channel containing the GRS80_height. The
        default is to read the altitudeChannel field name from the Coordinate Frame.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if GRS80_height == '':
            GRS80_height = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']

        flightLine = list(g.keys())[0]
        corr_units = g[flightLine][atmosCorr].attrs['Units']
        if corr_units == 'mGal':
            unit_scale = 10.0
        elif corr_units == 'gu' or corr_units == 'µm/s/s':
            unit_scale = 1.0
        else:
            print('ERROR - correction units not recognised, expected mGal or µm/s/s (gu)')
            return

        numLines = len(g.items())
        diffMin = np.zeros((numLines,))
        diffMax = np.zeros((numLines,))
        diffMean = np.zeros((numLines,))
        diffStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            ht_data = np.array(g[line][GRS80_height])
            cor_data = np.array(g[line][atmosCorr])
           
            cal_data = _atmosEffect(ht_data)
            err_data = cor_data * unit_scale - cal_data 
            diffMin[count] = np.min(err_data)
            diffMax[count] = np.max(err_data)
            diffMean[count] = np.mean(err_data)
            diffStd[count] = np.std(err_data)
            lineNo[count] = line
            count += 1
            
        fig = plt.figure()
        fig.suptitle(f'check Atmospheric Correction - {projName}', fontsize=10)
        fig.subplots_adjust(top=0.85)
        
        ax = fig.add_subplot(1,1,1)
        ax.plot(lineNo, diffMin, 'bo', mfc='w')
        ax.plot(lineNo, diffMax, 'bo', mfc='w')
        ax.errorbar(lineNo, diffMean, diffStd, capsize=3, marker='s', c='blue', ls='None')
        plt.ylabel('Atmospheric Correction Error [um/s/s]', fontsize = 6)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        plt.show()
    return

    
def checkHeading(whizzFile, nominalHeading, x='', y='', tolerance=10.0):
    '''
    Checks heading in degrees is within +/- tolerance (in degrees) of nominal (in degrees).

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    x : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    y : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)
    upper_limit = nominalHeading + tolerance
    lower_limit = nominalHeading - tolerance
    print(f'limits: {lower_limit}, {upper_limit}')

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']

        for line in g.keys():
            dx = np.diff(np.array(g[line][x]))
            dy = np.diff(np.array(g[line][y]))
            heading = np.arctan(dx / dy) * 180.0 / np.pi
            mean_heading = np.mean(heading)
            exceeds = any(h > upper_limit for h in heading)
            exceeds = exceeds or any(h < lower_limit for h in heading)
            
            if exceeds:
                print(f'Line {line}: heading limit exceeded. Mean {mean_heading}')
                fig = plt.figure()
                fig.suptitle(f'Heading Check Line {line}', fontsize=10)
                fig.subplots_adjust(top=0.85)
                
                ax = fig.add_subplot(1,1,1)
                ax.plot(np.array(g[line][x])[1:], heading, 'b', mfc='w')
                plt.ylabel('Estimated heading [deg]', fontsize = 6)
                plt.grid(True)
                for label in ax.get_xticklabels(): label.set_fontsize(6)
                for label in ax.get_yticklabels(): label.set_fontsize(6)
                plt.show()
    return
    

def checkFreeAirCorr(whizzFile, faCorr, latitude='', GRS80_height=''):
    '''
    Subtracts the free-air correction in the data file from one calculated using
    Hinze et al (2005) and the latitude and height data in the data file.
    The min, max, mean and standard deviation of the difference is calculated
    for each line and presented in an extended whisker plot.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    faCorr : String
        The name of the geoWhizz field or channel containing the free-air correction.
    latitude : String, optional
        The name of the geoWhizz field or channel containing the latitude. The
        default is to read the latitude field name from the Coordinate Frame.
    GRS80_height : String, optional
        The name of the geoWhizz field or channel containing the GRS80_height. The
        default is to read the altitudeChannel field name from the Coordinate Frame.

    Returns
    -------
    None.

    '''
    
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if latitude == '':
            latitude = f[groupName]['CoordinateFrame'].attrs['LatitudeChannel']
        if GRS80_height == '':
            GRS80_height = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']

        flightLine = list(g.keys())[0]
        corr_units = g[flightLine][faCorr].attrs['Units']
        if corr_units == 'mGal':
            unit_scale = 10.0
        elif corr_units == 'gu' or corr_units == 'µm/s/s':
            unit_scale = 1.0
        else:
            print('ERROR - correction units not recognised, expected mGal or µm/s/s (gu)')
            return

        numLines = len(g.items())
        diffMin = np.zeros((numLines,))
        diffMax = np.zeros((numLines,))
        diffMean = np.zeros((numLines,))
        diffStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            lat_data = np.array(g[line][latitude])
            ht_data = np.array(g[line][GRS80_height])
            cor_data = np.array(g[line][faCorr])
           
            cal_data = _freeAirCorrection(ht_data, lat_data)
            err_data = cor_data * unit_scale + cal_data 
            diffMin[count] = np.min(err_data)
            diffMax[count] = np.max(err_data)
            diffMean[count] = np.mean(err_data)
            diffStd[count] = np.std(err_data)
            lineNo[count] = line
            count += 1
            
            if line == '8334.0':
                fig = plt.figure()
                fig.suptitle(f'{projName}', fontsize=10)
                fig.subplots_adjust(top=0.85)
                ax = fig.add_subplot(3,1,1)
                ax.plot(cor_data*1e6)
                ax = fig.add_subplot(3,1,2)
                ax.plot(cal_data)
                ax = fig.add_subplot(3,1,3)
                ax.plot(err_data)
        
        fig = plt.figure()
        fig.suptitle(f'check Free Air Correction - {projName}', fontsize=10)
        fig.subplots_adjust(top=0.85)
        
        ax = fig.add_subplot(1,1,1)
        ax.plot(lineNo, diffMin, 'bo', mfc='w')
        ax.plot(lineNo, diffMax, 'bo', mfc='w')
        ax.errorbar(lineNo, diffMean, diffStd, capsize=3, marker='s', c='blue', ls='None')
        plt.ylabel('Free-air Correction Error [um/s/s]', fontsize = 6)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        plt.show()
    return
    

def _eotvosCorrection(eSpeed, nSpeed, latitude, height=0):
    """
    Calculates the Eotvos correction for moving-base gravimetry. Uses the exact
    equation (see for example Jekeli (2016), slide 51)

    Parameters
    ----------
    eSpeed (Float) : the aircraft speed in the East direction in m/s/s.
    
    nSpeed (Float) : the aircraft speed in the North direction in m/s/s.
    
    latitude (Float) : latitude in degrees (N pos, S neg).
    
    height (Float) : the height in metres above the GRS80 ellipsoid; optional, default 0.0.

    Returns
    -------
    eotvos (Float) : the eotvos correction for gravity in um/s/s.

    """
    
    radius = 6378137.0
    ellipticity = 0.0818191908426
    ellSquared = ellipticity * ellipticity
    angularVelocity = 7.2921158553E-5
    cosLat = np.cos(latitude * np.pi / 180)
    sinLat = np.sin(latitude * np.pi / 180)
    ellSinLatSquared = ellSquared * sinLat * sinLat
    
    angTerm = -2.0 * angularVelocity * eSpeed * cosLat
    eastTerm = - (eSpeed * eSpeed) / (height + radius / (np.sqrt(1 - ellSinLatSquared)))
    northTerm = - (nSpeed * nSpeed) / (height + (radius * (1 - ellSquared)) / ((1 - ellSinLatSquared) * np.sqrt(1 - ellSinLatSquared)))
    
    eotvos = 1.0E6 * (angTerm + eastTerm + northTerm)
    return eotvos
                    

def _normalGravity(latitude):
    """
    Returns the ellipsoid theoretical gravity in µm/s/s of Hinze et al (2005), eqn 2,
    based on the GRS80 ellipsoid.    

    Parameters
    ----------
    latitude (Float) : latitude in degrees (N pos, S neg).
    

    Returns
    -------
    latitude correction (Float) : ellipsoid theoretical gravity in um/s/s.

    """
    
    gNormal = 9780326.7715
    k = 0.001931851353
    eSquared = 0.0066943800229
    sinLatSquared = np.sin(latitude * np.pi / 180.0)
    sinLatSquared = sinLatSquared * sinLatSquared
    
    
    ratio = (1 + k * sinLatSquared) / np.sqrt(1 - eSquared * sinLatSquared)
    return gNormal * (ratio)# - 1.0)


def _atmosEffect(height):
    """
    Returns the atmospheric gravity correction in um/s/s of Aiken et al (2005), eqn 3.

    Parameters
    ----------
    height (Float) : the height in metres above the GRS80 ellipsoid; optional, default 0.0.


    Returns
    -------
    atmospheric correction (Float) : in m/s/s.

    """
    return 10.0 * (0.874 - 9.9E-5 * height + 3.56E-9 * height * height)


def _freeAirCorrection(height, latitude):
    """
    Returns the free-air gravity correction in um/s/s of Hinze et al (2005), eqn 3.
    
    Parameters
    ----------
    height (Float) : the height in metres above the GRS80 ellipsoid; optional, default 0.0.

    latitude (Float) : latitude in degrees (N pos, S neg).

    Returns
    -------
    free-air gravity correction in um/s/s.

    """
    
    sinLatSquared = np.sin(latitude * np.pi / 180.0)
    sinLatSquared = sinLatSquared * sinLatSquared
    freeAir = -(0.3087691 - 0.0004398 * sinLatSquared) * height
    freeAir += 7.2125E-8 * height * height
    # Hinze et al work in mGal but we want um/s/s so x 10
    return 10 * freeAir    
    
                   
def diurnal(filename, lines = [], name = 'Basemag', rangeLimit = 5.0, nSamples = 3000, showPlot=False):
    """
    QC Standard: May not exceed rangeLimit nT over nSamples.

    From each individual observation:
    
    Examine nSamples / 2 into the past and nSamples / 2 into the future.
    Interpolate linearly.
    Measure the deviation from the interpolated line.
    The deviation may not be larger than rangeLimit.

    Parameters
    ----------
    name : String
        Name of the channel containg 'basemag' data.
    basemag : Numpy array
        the basemag data to be checked.
    rangeLimit : Float, optional
        The largest allowable deviation from a straight chord (usu nT). The default is 5.0.
    nSamples : Int, optional
        The number of samples of 'basemag' for the windowing for the QC check 
        (at 10 Hz, 3000 samples = 5 min) . The default is 3000.
    showPlot : Bool, optional
        if showPlot and ~ok: plot the results.

    Returns
    -------
    ok : Bool
        True if 'basemag' passed the test.
    report : String
        "" if ok == True, else an error description.

    """
    finalStatus = True
    report = ''
    with h5py.File(filename, 'r') as f:
        g = f[groupName]
        if lines == []:
            lines = g.keys()
        
        for line in lines:
            basemag = np.array(g[line][name])
            lineStatus = True
            if nSamples > len(basemag):
                nSamples = len(basemag)
            if (nSamples % 2) == 0:
                nSamples = nSamples - 1
            nSam = (nSamples - 1) // 2
            for ii in range(nSam, len(basemag)-nSam):
                localData = basemag[ii-nSam:ii+nSam]
                localSlope = (localData[-1] - localData[0]) / nSamples
                # deviation = np.zeros(localData.shape)
                deviation = localData - localSlope * range(0, len(localData)) - localData[0]
                # for jj in range(0, len(localData)):
                #     deviation[jj] = localData[jj] - localData[0] - localSlope * jj
                extremum = np.max(deviation) if np.max(deviation) > -np.min(deviation) else -np.min(deviation)
                if extremum > rangeLimit:
                    lineStatus = False
                    finalStatus = False
                    print(f'\n  Line {line} Diurnal for {name} reaches {extremum:.1f}, exceeding {rangeLimit:.1f} - FAIL')
                    report += f'\n  Line {line} Diurnal for {name} reaches {extremum:.1f}, exceeding {rangeLimit:.1f} - FAIL'
                    break
                
            print(lineStatus)   
            if lineStatus == False and showPlot == True:
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                ax.plot(basemag, lw=0.5)
                plotTitle = f'Line {line} Channel {name}: reaches {extremum:.1f} at {ii}, exceeding {rangeLimit:.1f} - FAIL'
                plt.title(plotTitle, fontsize = 8)
                plt.grid(True)
                for label in ax.get_xticklabels(): label.set_fontsize(6)
                for label in ax.get_yticklabels(): label.set_fontsize(6)
                fig.tight_layout()
                plt.show()
                
    return finalStatus, report
            

def checkPhase(filename, channel1, channel2):
    """
    For every survey line in a geoWhizz HDF5 file, given two fields, calculate
    the phase shift (in number of samples) required to maximise the correlation
    between the two fields

    Parameters
    ----------
    filename (String) : the name of a geoWhizz HDF5 file.

    channel1 (String) : The name of the first field.
    
    channel2 (String) : The name of the second field.
    
    Returns
    -------
    None.

    """
    from scipy.signal import correlate
    with h5py.File(filename, 'r') as f:
        g = f[groupName]
        numLines = len(g.items())
        offsets = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            A = np.array(g[line][channel1])
            B = np.array(g[line][channel2])
            nsamples = A.size

            # regularize datasets by subtracting mean and dividing by s.d.
            A -= A.mean(); A /= A.std()
            B -= B.mean(); B /= B.std()

            # Put in an artificial time shift between the two datasets
            #time_shift = 20
            #A = numpy.roll(A, time_shift)

            # Find cross-correlation
            xcorr = correlate(A, B)
            
            # delta time array to match xcorr
            dt = np.arange(1-nsamples, nsamples)
            recovered_time_shift = dt[xcorr.argmax()]
            offsets[count] = recovered_time_shift
            count += 1

            time = np.arange(dt[0], dt[-1], 0.1)
            # Now interpolate through gaps by cubic spline
            (xcorrInt, _) = mhd.interpolateLine(dt, xcorr, time)
            recovered_time_shift2 = time[xcorrInt.argmax()]
            print(f'Line {line}: Recovered time shift = {recovered_time_shift2:.1f}')
            
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(time, xcorrInt)
        plt.ylabel('Correlation [arbitrary units]', fontsize = 6)
        plt.xlabel('FID difference [s]', fontsize = 6)
        plotTitle = f'Line {line}: Correlation of {channel1} v {channel2}'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        plt.show()
        print(f'Offset MEAN = {np.mean(offsets):.2f}; STD = {np.std(offsets):.3f}')
        
        
def ilsNoiseVturb(whizzFile, diagComponent1, diagComponent2, diagComponent3, vertvelocity):
    '''
    For a Bell Air-FTG. For each line, reports the standard deviation of the in-line sums,
    and plots these as a scatter plot against the standard deviation of the vertical
    acceleration (estimated as the time difference of the vertical velocity).

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    diagComponent1 : String
        The name of the channel containing the first tensor diagonal (in-line) component.
    diagComponent2 : String
        The name of the channel containing the second tensor diagonal (in-line) component.
    diagComponent3 : String
        The name of the channel containing the third tensor diagonal (in-line) component.
    vertvelocity : String
        The name of the channel containing the vertical velocity field.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        numLines = len(g.items())
        accStd = np.zeros((numLines,))
        ilsStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            data = g[line][vertvelocity]
            accel = np.diff(data, n = 1)
            accStd[count] = np.std(accel)

            data1 = np.array(g[line][diagComponent1])
            data2 = np.array(g[line][diagComponent2])
            data3 = np.array(g[line][diagComponent3])
            ilsStd[count] = np.std(_inLineSum(data1, data2, data3))
            print(line, ' ', ilsStd[count])
            lineNo[count] = line
            count += 1
        
        x = {'label': 'Vertical acceleration [m/s/s]', 'data': accStd}
        y = {'label': 'Inline Sum [E]', 'data': ilsStd}
        wpl.plotxy(y, x, plotTitle = 'In-line Sum Noise versus Turbulence', xOffset=False, plot_symbol='+')


def diffGravVturb(whizzFile, turbulence, aD, bD, error_spec=5.0, low_cut=0.001, measX='', measY=''):
    '''
    For a Falcon AGG. For each line, reports the gD difference noise,
    stdev(aD-bD)/2, and plots this as a scatter plot against the mean turbulence.
    All lines with difference noise greater than `error_spec` are reported.
    
    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    turbulence : String
        The name of the channel containing the turbulence field.
    aD : String
        The name of the channel containing the A complement gravity data.
    bD : String
        The name of the channel containing the B complement gravity data.
    error_spec : Float, optional
        The value above which the difference noise is excessive and is
        reported. Default 5.0.
    low_cut : Float, optional
        The low frequency cut-off frequency (in Hz) for the band-pass filtering
        applied before differencing. Default 0.001 (ie 1 mHz).

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)

    num_lines_failed = 0
    report = 'diffGravVturb() estimates the noise in (A+B)/2 as the stdev(A-B)/2.\n'
    period = 1.0 / low_cut
    wavelength = 60.0 * period / 1000.0 # km
    report += f'The input data are band-pass filtered at [{period}, 1.0] sec or [{wavelength}, 0.06] km at 60m/s.\n'
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        numLines = len(g.items())
        turbMean = np.zeros((numLines,))
        errmean = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        line_strs = list(g.keys())
        flightLine = line_strs[0]
        turb_units = g[flightLine][turbulence].attrs['Units']
        err_units = g[flightLine][aD].attrs['Units']

        if measX == '':
            measX = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if measY == '':
            measY = f[groupName]['CoordinateFrame'].attrs['YChannel']


        failed_lines = 0
        for line in g.keys():
            xM = np.array(g[line][measX])
            yM = np.array(g[line][measY])
            line_length = _displacement2(xM[0], xM[-1], yM[0], yM[-1]) / 1000.0

            turb = np.array(g[line][turbulence])
            A_d = np.array(g[line][aD])
            B_d = np.array(g[line][bD])
            idx = np.where(~np.isnan(A_d + B_d))
            Ad = butter_bandpass_filter(A_d[idx], low_cut, 1.0, 8.0, order=3)
            Bd = butter_bandpass_filter(B_d[idx], low_cut, 1.0, 8.0, order=3)
            turb_clean = turb[idx]
            err_data = (Ad - Bd)/2.0
            err_data = err_data - np.mean(err_data)

            turbMean[count] = np.mean(turb_clean)
            errmean[count] = np.std(err_data)
            
            if errmean[count] > error_spec:
                report += f'{line} (length {line_length:6.1f} km) fails with noise {errmean[count]:.1f} > {error_spec}.\n'
                num_lines_failed += 1

            count += 1

        fig = plt.figure()
        fig.suptitle(f'Noise vs Turbulence - {projName}', fontsize=12)
        fig.subplots_adjust(top=0.85)
        ax = fig.add_subplot(1,1,1)
        neplot, = ax.plot(turbMean, errmean, 'bo', label='$E_{D}$')
        for ii in range(0, len(turbMean)):
            ax.text(turbMean[ii], errmean[ii], f'{line_strs[ii]}', fontsize=8)
        plt.ylabel(f'Difference Noise [{err_units}]', fontsize = 8)
        plt.xlabel(f'Turbulence [{turb_units}]', fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(8)
        for label in ax.get_yticklabels(): label.set_fontsize(8)
        fig.tight_layout()
        plt.show()
        report = f'{num_lines_failed} failed of {count} lines total.' + report
    print(report)
        

def diffNoiseVturb(whizzFile, turbulence, aNE='', aUV='', bNE='', bUV='', eNE='', eUV='', error_spec = 5.0):
    '''
    For a Falcon AGG. For each line, reports the mean NE and UV difference noise,
    and plots these as a scatter plot against the mean turbulence.
    
    TODO - fix algorithm, calculation not quite right

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    turbulence : String
        The name of the channel containing the turbulence field.
    aNE : String
        The name of the channel containing the A_NE field. If '', then the
        eNE and eUV channels are used, otherwise, the aNE, aUV, bNE
        and bUV channels are used.
    aUV : String
        The name of the channel containing the A_UV field. If '', then the
        eNE and eUV channels are used, otherwise, the aNE, aUV, bNE
        and bUV channels are used.
    bNE : String
        The name of the channel containing the B_NE field. If '', then the
        eNE and eUV channels are used, otherwise, the aNE, aUV, bNE
        and bUV channels are used.
    bUV : String
        The name of the channel containing the B_UV field. If '', then the
        eNE and eUV channels are used, otherwise, the aNE, aUV, bNE
        and bUV channels are used.
    eNE : String, optional
        The name of the channel containing the NE difference error. If '',
        then the aNE, aUV, bNE and bUV channels are used.
    eUV : String
        The name of the channel containing the UV difference error. If '',
        then the aNE, aUV, bNE and bUV channels are used.
    error_spec : Float
        The value above which the difference noise is excessive and is reported.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)
    if aNE=='' and aUV=='' and bNE=='' and bUV=='':
        need_calc = False
    elif eNE=='' and eUV=='':
        need_calc = True
    else:
        print('ERROR - must specify either all four raw channels or both error channels')
        return

    report = ''
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        numLines = len(g.items())
        turbMean = np.zeros((numLines,))
        errNEmean = np.zeros((numLines,))
        errUVmean = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0
        flightLine = list(g.keys())[0]


        turb_units = g[flightLine][turbulence].attrs['Units']
        if need_calc:
            err_units = g[flightLine][aNE].attrs['Units']
        else:
            err_units = g[flightLine][eNE].attrs['Units']

        failed_lines = 0
        for line in g.keys():
            turb = np.array(g[line][turbulence])
            if need_calc:
                A_ne = np.array(g[line][aNE])
                A_uv = np.array(g[line][aUV])
                B_ne = np.array(g[line][bNE])
                B_uv = np.array(g[line][bUV])
                idx = np.where(~np.isnan(A_ne + A_uv + B_ne + B_uv))
                Ane = A_ne[idx]
                Auv = A_uv[idx]
                Bne = B_ne[idx]
                Buv = B_uv[idx]
                turb_clean = turb[idx]
    
                errNE_data = (Ane - Bne)/np.sqrt(8)
                errUV_data = (Auv - Buv)/np.sqrt(8)
            else:
                E_ne = np.array(g[line][eNE])
                E_uv = np.array(g[line][eUV])
                idx = np.where(~np.isnan(E_ne + E_uv))
                errNE_data = E_ne[idx]
                errUV_data = E_uv[idx]
                turb_clean = turb[idx]

            turbMean[count] = np.mean(turb_clean)
            errNEmean[count] = np.std(errNE_data)
            errUVmean[count] = np.std(errUV_data)
            
            if errNEmean[count] > error_spec or errUVmean[count] > error_spec:
                report += f'{line} fails with noise > {error_spec}.\n'
            count += 1

        fig = plt.figure()
        fig.suptitle(f'Noise vs Turbulence - {projName}', fontsize=12)
        fig.subplots_adjust(top=0.85)
        ax = fig.add_subplot(1,1,1)
        ax.vlines(turbMean, errNEmean, errUVmean, 'k', lw=0.3)
        neplot, = ax.plot(turbMean, errNEmean, 'bo', label='$E_{ne}$')
        uvplot, = ax.plot(turbMean, errUVmean, 'go', label='$E_{uv}$')
        plt.ylabel(f'Difference Noise [{err_units}]', fontsize = 8)
        plt.xlabel(f'Turbulence [{turb_units}]', fontsize = 8)
        plt.grid(True)
        ax.legend(handles=[neplot, uvplot])
        for label in ax.get_xticklabels(): label.set_fontsize(8)
        for label in ax.get_yticklabels(): label.set_fontsize(8)
        fig.tight_layout()
        plt.show()
    print(report)
        

def _inLineSum(il1, il2, il3, fs=1.0, lowcut=0.03, highcut=0.1):
    '''
    Calculates the filtered in-line sum of the in-line components.

    Parameters
    ----------
    il1 : numpy 1D array
        The first in-line component data.
    il2 : numpy 1D array
        The second in-line component data.
    il3 : numpy 1D array
        The third in-line component data.
    fs : Float
        The sample frequency.
    lowcut : Float
        The low-pass frequency in Hz of the filter.
    highcut : Float
        The high-pass frequency in Hz of the filter.

    Returns
    -------
    numpy 1D array
        The filtered in-line sum.

    '''
    fs = 1.0
    lowcut = 0.03
    highcut = 0.1
    order = 6
    ils = (il1 + il2 + il3) / np.sqrt(3.0)
    ils = ils - np.mean(ils)
    return butter_bandpass_filter(ils, lowcut, highcut, fs, order = order)

    
def _FTGTransform(il1, il2, il3, cr1, cr2, cr3):
    '''
    Returns Txx, Txy, Txz, Tyy, Tyz, Tzz from an algebraic transform of
    the three inline and three cross components of the FTG gradiometer.
    Assumes, I believe, that the FTG is oriented so that the horizontal
    projection of the spin axis of GGI3 is north (the y-axis). From equation (5)
    in J. Brewster, Comparison of gravity gradiometer designs using the 3D
    sensitivity function. In SEG International Exposition and 86th Annual
    Meeting, 2016.

    Parameters
    ----------
    il1 : numpy 1D array
        The first in-line component data.
    il2 : numpy 1D array
        The second in-line component data.
    il3 : numpy 1D array
        The third in-line component data.
    cr1 : numpy 1D array
        The first cross component data.
    cr2 : numpy 1D array
        The second cross component data.
    cr3 : numpy 1D array
        The third cross component data.

    Returns
    -------
    Txx : numpy 1D array
        The xx component of the gravity gradient in ENU coordinates.
    Txy : numpy 1D array
        The xy component of the gravity gradient in ENU coordinates.
    Txz : numpy 1D array
        The xz component of the gravity gradient in ENU coordinates.
    Tyy : numpy 1D array
        The yy component of the gravity gradient in ENU coordinates.
    Tyz : numpy 1D array
        The yz component of the gravity gradient in ENU coordinates.
    Tzz : numpy 1D array
        The zz component of the gravity gradient in ENU coordinates.

    '''
    sq2 = np.sqrt(2.0)
    sq3 = np.sqrt(3.0)
    Txx = 2.0 / 3.0 * ( il3 - sq2 * cr1)
    Txy = 0.75/np.sqrt(2) * (np.sqrt(2/3) * il3 - np.sqrt(1.5) * (cr1 - cr2))
    Txz = np.sqrt(6)/(2+3*np.sqrt(2)) * (cr1 - cr2 + 2 * il3)
    Tyy = sq2 / 3.0 * cr1 - (sq2 * cr2 + cr3) / sq3  - il3 / 3.0
    Tyz = 3/np.sqrt(2) * (2/3 * (cr3 - il1 + il2) - (cr1 + cr2)/3)
    Tzz = sq2 / 3.0 * cr1 + (sq2 * cr2 + cr3) / sq3  - il3 / 3.0
    
    return Txx, Txy, Txz, Tyy, Tyz, Tzz


def eigenPlot(whizzFile, lines = [], noiselimit=30.0):
    '''
    Reports the noise for each line in lines from filename (a whizz file)
    which exceeds noiseLimit. Here the noise is calculated by `_FTGeigen` as
    the Frobenius norm.
    TODO : Why is the function name misleading?

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    lines : String list, optional.
        The line numbers to be checked. Default is all lines in the whizzFile.
    noiselimit : Float
        The maximum allowable in-line noise on a line.

    Returns
    -------
    None

    '''
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if lines == []:
            lines = g.keys()
        
        for line in lines:
            i1 = np.array(g[line]['Inline1_raw'])
            i2 = np.array(g[line]['Inline2_raw'])
            i3 = np.array(g[line]['Inline3_raw'])
            c1 = np.array(g[line]['Cross1_raw'])
            c2 = np.array(g[line]['Cross2_raw'])
            c3 = np.array(g[line]['Cross3_raw'])
            (Gxx, Gxy, Gxz, Gyy, Gyz, Gzz) = _FTGTransform(i1, i2, i3, c1, c2, c3)
            Txx = butter_bandpass_filter(Gxx, 0.1, 0.49, 1.0, order = 6)
            Txy = butter_bandpass_filter(Gxy, 0.1, 0.49, 1.0, order = 6)
            Txz = butter_bandpass_filter(Gxz, 0.1, 0.49, 1.0, order = 6)
            Tyy = butter_bandpass_filter(Gxy, 0.1, 0.49, 1.0, order = 6)
            Tyz = butter_bandpass_filter(Gyz, 0.1, 0.49, 1.0, order = 6)
            Tzz = butter_bandpass_filter(Gzz, 0.1, 0.49, 1.0, order = 6)
            noise = _FTGeigen(Txx, Txy, Txz, Tyy, Tyz, Tzz, line, noiselimit)
            print(f'Check line {line}. Noise = {noise:.1f}')
    

def checkRawAGG(whizzFile, ane, auv, bne, buv, turb, time='', lines=[], noiseLimit=10.0):
    '''
    Reports the high frequency noise for each line in lines from filename (a whizz file)
    which exceeds noiseLimit. See QC Report for Lawin for details of method.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    aNE : String
        The name of the channel containing the A_NE field.
    aUV : String
        The name of the channel containing the A_UV field.
    bNE : String
        The name of the channel containing the B_NE field.
    bUV : String
        The name of the channel containing the B_UV field.
    turbulence : String
        The name of the channel containing the turbulence field.
    lines : String list, optional
        The line numbers to be checked. Default is all lines in the whizzFile.
    noiselimit : Float, optional
        The maximum allowable high frequency noise on a line. Default=10.0.

    Returns
    -------
    None

    '''
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if lines == []:
            lines = g.keys()
        
        for line in lines:
            reportStr = f'Line {line} Noise: '
            if time == '':
                time = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
            time_data = np.array(g[line][time])
            time_data = time_data - time_data[0]
            fs = 1.0 / abs((time_data[1] - time_data[0]))
            turb_data = np.array(g[line][turb])
            for channel in [ane, auv, bne, buv]:
                data = np.array(g[line][channel])
                noSlope = np.zeros((len(data),))
                filtered = np.zeros((len(data),))
                myStd = np.zeros((len(data)-50,))
        
                stData = np.mean(data[:10])
                enData = np.mean(data[-10:])
                slope = (enData - stData) / len(data)
                for ii in range(0, len(data)):
                    noSlope[ii] = data[ii] - stData - ii * slope
                        
                filtered = butter_bandpass_filter(noSlope, 0.15,
                                                  3.6, fs, order = 6)
                
                for ii in range(0, len(data)-50):
                    myStd[ii] = np.std(filtered[ii:ii+50])
                    
                reportStr += f'{channel} = {np.max(myStd):.1f}; '
                
                if np.max(myStd) > noiseLimit:
                    fig = plt.figure()
                    ax1 = fig.add_subplot(4,1,1)
                    ax1.plot(time_data, data, time_data, noSlope, lw=0.5)
                    ax1.set_xlim(time_data[0], time_data[-1])
                    plotTitle = f'Line {line}, Channel {channel}: orig and slope removed'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax1.get_xticklabels(): label.set_fontsize(6)
                    for label in ax1.get_yticklabels(): label.set_fontsize(6)
                    ax2 = fig.add_subplot(4,1,2)
                    ax2.plot(time_data, filtered, lw=0.5)
                    ax2.set_xlim(time_data[0], time_data[-1])
                    ax2.set_ylim(-10.0, 10.0)
                    plotTitle = f'filtered [0.15,3.6]; fs={fs:.1f}'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax2.get_xticklabels(): label.set_fontsize(6)
                    for label in ax2.get_yticklabels(): label.set_fontsize(6)
                    ax3 = fig.add_subplot(4,1,3)
                    sTime = time_data[25:25+len(myStd)]
                    ax3.plot(sTime, myStd, lw=0.5)
                    ax3.set_xlim(time_data[0], time_data[-1])
                    ax3.set_ylim(0.0, 15.0)
                    plotTitle = 'rolling stdev'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax3.get_xticklabels(): label.set_fontsize(6)
                    for label in ax3.get_yticklabels(): label.set_fontsize(6)
                    ax4 = fig.add_subplot(4,1,4)
                    ax4.plot(time_data, turb_data, lw=0.5)
                    ax4.set_xlim(time_data[0], time_data[-1])
                    ax4.set_ylim(0.0, 2.0)
                    plotTitle = 'turb'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax4.get_xticklabels(): label.set_fontsize(6)
                    for label in ax4.get_yticklabels(): label.set_fontsize(6)
                    fig.tight_layout()
                    plt.show()
            print(reportStr)


def checkRawFTG(whizzFile, lines=[], noiseLimit=50):
    '''
    Reports the high frequency noise for each line in lines from filename (a whizz file)
    which exceeds noiseLimit. See QC Report for Lawin for details of method.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    lines : String list, optional
        The line numbers to be checked. Default is all lines in the whizzFile.
    noiselimit : Float
        The maximum allowable high frequency noise on a line.

    Returns
    -------
    None

    '''
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if lines == []:
            lines = g.keys()
        
        for line in lines:
            reportStr = f'Line {line} Noise: '
            time = g[line]['time']
            time = time - time[0]
            vv = g[line]['vertvelocity']
            turb = np.diff(vv, n = 1)
            for channel in ['Cross1_raw', 'Cross2_raw', 'Cross3_raw', 'Inline1_raw', 'Inline2_raw', 'Inline3_raw']:
                data = np.array(g[line][channel])
                noSlope = np.zeros((len(data),))
                filtered = np.zeros((len(data),))
                myStd = np.zeros((len(data)-50,))
        
                stData = np.mean(data[:10])
                enData = np.mean(data[-10:])
                slope = (enData - stData) / len(data)
                for ii in range(0, len(data)):
                    noSlope[ii] = data[ii] - stData - ii * slope
                        
                filtered = butter_bandpass_filter(noSlope, 0.1, 0.48, 1, order = 6)
                
                for ii in range(0, len(data)-50):
                    myStd[ii] = np.std(filtered[ii:ii+50])
                    
                reportStr += f'{channel} = {np.max(myStd):.1f}; '
                
                if np.max(myStd) > noiseLimit:
                    fig = plt.figure()
                    ax1 = fig.add_subplot(4,1,1)
                    ax1.plot(time, data, time, noSlope, lw=0.5)
                    ax1.set_xlim(time[0], time[-1])
                    plotTitle = f'Line {line}, Channel {channel}: orig and slope removed'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax1.get_xticklabels(): label.set_fontsize(6)
                    for label in ax1.get_yticklabels(): label.set_fontsize(6)
                    ax2 = fig.add_subplot(4,1,2)
                    ax2.plot(time, filtered, lw=0.5)
                    ax2.set_xlim(time[0], time[-1])
                    ax2.set_ylim(-100.0, 100.0)
                    plotTitle = 'filtered'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax2.get_xticklabels(): label.set_fontsize(6)
                    for label in ax2.get_yticklabels(): label.set_fontsize(6)
                    ax3 = fig.add_subplot(4,1,3)
                    sTime = time[25:25+len(myStd)]
                    ax3.plot(sTime, myStd, lw=0.5)
                    ax3.set_xlim(time[0], time[-1])
                    ax3.set_ylim(0.0, 50.0)
                    plotTitle = 'rolling stdev'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax3.get_xticklabels(): label.set_fontsize(6)
                    for label in ax3.get_yticklabels(): label.set_fontsize(6)
                    ax4 = fig.add_subplot(4,1,4)
                    ax4.plot(time[1:], turb, lw=0.5)
                    ax4.set_xlim(time[0], time[-1])
                    ax4.set_ylim(0.0, 2.0)
                    plotTitle = 'turb'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax4.get_xticklabels(): label.set_fontsize(6)
                    for label in ax4.get_yticklabels(): label.set_fontsize(6)
                    fig.tight_layout()
                    plt.show()
            print(reportStr)

    
def _FTGeigen(Txx, Txy, Txz, Tyy, Tyz, Tzz, line = "", noiselimit=30.0):
    '''
    Returns Txx, Txy, Txz, Tyy, Tyz, Tzz from an algebraic transform of
    the three inline and three cross components of the FTG gradiometer.
    Assumes, I believe, that the FTG is oriented so that the horizontal
    projection of the spin axis of GGI3 is north (the y-axis). From equation (5)
    in J. Brewster, Comparison of gravity gradiometer designs using the 3D
    sensitivity function. In SEG International Exposition and 86th Annual
    Meeting, 2016.

    Parameters
    ----------
    Txx : numpy 1D array
        The xx component of the gravity gradient in ENU coordinates.
    Txy : numpy 1D array
        The xy component of the gravity gradient in ENU coordinates.
    Txz : numpy 1D array
        The xz component of the gravity gradient in ENU coordinates.
    Tyy : numpy 1D array
        The yy component of the gravity gradient in ENU coordinates.
    Tyz : numpy 1D array
        The yz component of the gravity gradient in ENU coordinates.
    Tzz : numpy 1D array
        The zz component of the gravity gradient in ENU coordinates.
    line : String, optional.
        The line number containing the above data, just used for
        reporting. Default ''.
    noiselimit : Float, optional.
        The maximum allowed standard deviation of the Frobenius norm.
        Default 30.0 E.

    Returns
    -------
    Float
        The standard deviation of the Frobenius norm of the gravity gradient.

    '''
    numSamples = len(Txx)
    trace = np.zeros((numSamples,))
    det = np.zeros((numSamples,))
    I2 = np.zeros((numSamples,))
    frob = np.zeros((numSamples,))
    
    for ii in range(0, numSamples):
        a = np.array([[Txx[ii], Txy[ii], Txz[ii]], [Txy[ii], Tyy[ii], Tyz[ii]], [Txz[ii], Tyz[ii], Tzz[ii]]])
        w, v = np.linalg.eig(a)
        trace[ii] = w.sum()
        det[ii] = np.cbrt(w[0] * w[1] * w[2])
        I2[ii] = w[0] * w[1] - (w[0] + w[1]) * (w[0] + w[1])
        frob[ii] = np.linalg.norm(a, 'fro')

    if np.std(frob) > noiselimit:
        myTitle = 'Trace for Line ' + line
        fig = plt.figure(figsize=(5,8))
        ax1 = fig.add_subplot(4,1,1)
        ax1.plot(trace, '.', ms=2)
        plt.ylabel('Trace', fontsize = 6)
        plotTitle = myTitle
        plt.title(plotTitle, fontsize = 6)
        plt.grid(True)
        for label in ax1.get_xticklabels(): label.set_fontsize(6)
        for label in ax1.get_yticklabels(): label.set_fontsize(6)

        ax2 = fig.add_subplot(4,1,2)
        ax2.plot(det, '.', ms=2)
        plt.ylabel('Det', fontsize = 6)
        plotTitle = 'Det'
        plt.title(plotTitle, fontsize = 6)
        #ax.set_ylim(0.0, 0.0)
        plt.grid(True)
        for label in ax2.get_xticklabels(): label.set_fontsize(6)
        for label in ax2.get_yticklabels(): label.set_fontsize(6)
        
        ax3 = fig.add_subplot(4,1,3)
        ax3.plot(-np.sqrt(-I2), '.', ms=2)
        plt.ylabel('I2', fontsize = 6)
        plotTitle = 'I2'
        plt.title(plotTitle, fontsize = 6)
        #ax.set_ylim(0.0, 0.0)
        plt.grid(True)
        for label in ax3.get_xticklabels(): label.set_fontsize(6)
        for label in ax3.get_yticklabels(): label.set_fontsize(6)
        
        ax4 = fig.add_subplot(4,1,4)
        ax4.plot(frob, '.', ms=2)
        plt.ylabel('Frob Norm', fontsize = 6)
        plotTitle = f'Frob Norm: std = {np.std(frob):.1f}'
        plt.title(plotTitle, fontsize = 6)
        ax4.set_ylim(0.0, 100.0)
        plt.grid(True)
        for label in ax4.get_xticklabels(): label.set_fontsize(6)
        for label in ax4.get_yticklabels(): label.set_fontsize(6)
        
        fig.tight_layout()
        plt.show()
    return np.std(frob)
        
        
def calcDrift(whizzFile, time, gradient):
    '''
    NOT USED
    '''
    filename = str(whizzFile)
    from sklearn.linear_model import LinearRegression
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        
        for line in g.keys():
            t = np.array(g[line][time]).reshape((-1,1))
            gamma = g[line][gradient]
            
            model = LinearRegression().fit(t, gamma)
            r_sq = model.score(t, gamma)
            print('coefficient of determination:', r_sq, 'intercept:', model.intercept_, 'slope:', model.coef_)


def checkLineLengths(whizzFile, min_len=50.0, measX='', measY=''):
    '''
    Checks that all lines in whizzFile are at least min_len km long.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    min_len : TYPE, optional
        The minimum allowed line length in km. The default is 50.0.
    measX : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    measY : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.

    Returns
    -------
    None.

    '''
    measFile = str(whizzFile)
    
    with h5py.File(measFile, 'r') as f:
        gMeas = f[groupName]['Lines']
        if measX == '':
            measX = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if measY == '':
            measY = f[groupName]['CoordinateFrame'].attrs['YChannel']
        
        num_failed_lines = 0
        for line in gMeas.keys():
            xM = np.array(gMeas[line][measX])
            yM = np.array(gMeas[line][measY])
            line_length = _displacement2(xM[0], xM[-1], yM[0], yM[-1])
            if line_length < min_len * 1000.0:
                num_failed_lines += 1
                print(f'Line {line} length = {line_length:.1f} less than allowed min {min_len*1000.0:.1f}')
        print(f'Number failed lines = {num_failed_lines}')
            

def checkOverlaps(whizzFile, min_overlap = 7.6, lines = [], verbose=False, plot_flag=False):
    '''
    For every line in the file whizzFile, calculate the overlap with each other
    line that has the same prefix. Plot a map. Report overlaps.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    min_overlap : Float, optional
        The minimum overlap distance in km, default 7.6 km.
    lines : String list, optional.
        The line numbers to be checked. Default is all lines in the whizzFile.
    verbose : Bool, optional
        If True, report status of all overlaps, else only report errors. Default False.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line.

    Returns
    -------
    None.

    '''
    from matplotlib.ticker import StrMethodFormatter
    measFile = str(whizzFile)
    if plot_flag:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    
    with h5py.File(measFile, 'r') as f:
        gMeas = f[groupName]['Lines']
        east = f[groupName]['CoordinateFrame'].attrs['XChannel']
        nrth = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if lines == []:
            lines = list(gMeas.keys())
        num_coinc_lines = 0
        report = ''

        for idx1 in range(0, len(lines)): #line1 in lines: #gMeas.keys():
            line1 = lines[idx1]
            line1_plan = gMeas[line1].attrs['PlannedLine']

            for idx2 in range (idx1 + 1, len(lines)):#line2 in gMeas.keys():
                line2 = lines[idx2]
                line2_plan = gMeas[line2].attrs['PlannedLine']
                # if the second line isn't the first line but has the same planned line no.
                if line1 != line2 and line1_plan == line2_plan:
                    num_coinc_lines += 1
                    # extract positions
                    n1 = np.array(gMeas[line1][nrth])
                    e1 = np.array(gMeas[line1][east])
                    n2 = np.array(gMeas[line2][nrth])
                    e2 = np.array(gMeas[line2][east])
                    # get line direction in radians
                    dirn = np.arctan((e1[-1] - e1[0])/(n1[-1] - n1[0]))
                    
                    # make life easier by transforming to a 2D problem
                    n0 = n1[0]
                    e0 = e1[0]
                    (y1, x1) = _rotateCoords(e1 - e0, n1 - n0, -dirn)
                    (y2, x2) = _rotateCoords(e2 - e0, n2 - n0, -dirn)

                    # find ends of lines
                    min1 = np.nanmin(x1)
                    max1 = np.nanmax(x1)
                    min2 = np.nanmin(x2)
                    max2 = np.nanmax(x2)
                    overlap = 0.0
                    whole_line = 0
                    
                    # find overlap length if any
                    if (min2 < min1 and min1 < max2) and not (min2 < max1 and max1 < max2):
                        overlap = abs(max2 - min1)
                    elif (min2 < max1 and max1 < max2) and not (min2 < min1 and min1 < max2):
                        overlap = abs(min2 - max1)
                    elif (min2 < max1 and max1 < max2) and (min2 < min1 and min1 < max2):
                        overlap = abs(max1 - min1)
                        whole_line = 1
                    elif (min1 < max2 and max2 < max1) and (min1 < min2 and min2 < max1):
                        overlap = abs(max2 - min2)
                        whole_line = 2

                    if plot_flag:
                        plotline1, = ax.plot(e1, n1, color='blue', lw=0.2)
                    
                    if whole_line == 1:
                        if verbose:
                            report += f'  OK: All of line {line1} in {line2}: length {overlap:.0f} m.\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='green', lw=0.2)
                    elif whole_line == 2:
                        if verbose:
                            report += f'  OK: All of line {line2} in {line1}: length {overlap:.0f} m.\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='green', lw=0.2)
                    elif overlap < 1:
                        if verbose:
                            report += f'  OK: Non-overlapping lines {line1}, {line2}: overlap {overlap:.0f} m < {min_overlap * 1000}.\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='green', lw=0.2)
                    elif overlap < min_overlap * 1000.0:
                        report += f'\n  ERROR: Repeat line {line1}, {line2}: overlap {overlap:.0f} m < {min_overlap * 1000}.\n\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='red', lw=0.2)
                    else:
                        if verbose:
                            report += f'  OK: Repeat line {line1}, {line2}: overlap by {overlap:.0f} m.\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='green', lw=0.2)

    print(f'{num_coinc_lines} coincident lines found.')
    print(report)
    if num_coinc_lines > 0 and plot_flag:
        ax.set_aspect('equal')
        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        plt.xlabel('X [m]', fontsize = 10)
        plt.ylabel('Y [m]', fontsize = 10)
        plt.suptitle('Overlap Map', fontsize = 12)
        plt.title('[1st line blue, accepted line green, error line red]', fontsize = 10)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(10)
        for label in ax.get_yticklabels(): label.set_fontsize(10)
        plt.show()


def checkXYPlan(planPath, measPath, planX='', planY='', measX='', measY='', allowance=200.0, maxCounter=14, maxDistance=0, known='', plot_flag=False):
    '''
    Reports exceedances of actual horizontal position from planned horizontal
    positions for an airborne survey Whizz database.
    The positions (planX, planY) of the start and end of each planned survey line
    are read from planPath. The measured positions (measX, measY) are read from
    measPath and the perpendicular distance of each from the planned line is
    calculated. If this distance exceeds allowance for maxCounter consecutive
    fiducial, or maxDistance consecutive metres, then an out-of-specification
    exceedance is reported for that line.

    Parameters
    ----------
    planPath : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, with the survey
        positions plan.
    planX : String, optional
        The name of the geoWhizz field or channel containing the planned x position. The
        default is to read the xChannel field name from the Coordinate Frame.
    planY : String, optional
        The name of the geoWhizz field or channel containing the planned y position. The
        default is to read the yChannel field name from the Coordinate Frame.
    measPath : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, with the survey
        measured data.
    measX : String, optional
        The name of the geoWhizz field or channel containing the measured x position. The
        default is to read the xChannel field name from the Coordinate Frame.
    measY : String, optional
        The name of the geoWhizz field or channel containing the measured y position. The
        default is to read the yChannel field name from the Coordinate Frame.
    allowance : Float, optional
        The allowed horizontal distance for the measured line from the planned
        line. If any portion of a measured line is further than this from the
        planned position for more than the maximum allowed number of fids or the
        maximum allowed distance, then the line fails. The default is 200.0.
    maxCounter : Int, optional
        The maximum number of consecutive fids for which an exceedance
        greater than allowance is permitted. If 0, then the constraint is
        ignored. The default is 14.
    maxDistance : Float, optional
        The maximum number of consecutive metres for which an exceedance
        greater than allowance is permitted. If 0, then the constraint is
        ignored. The default is 0.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line.

    Returns
    -------
    None.

    '''
    start_x = 0.0
    end_x = 0.0
    start_y = 0.0
    end_y = 0.0

    exceedances_known = False
    this_exc_known = False
    number_known = 0

    planfile = str(planPath)
    measFile = str(measPath)
    
    with h5py.File(planfile, 'r') as fp:
        gPlan = fp[groupName]['Lines']
        if planX == '':
            planX = fp[groupName]['CoordinateFrame'].attrs['XChannel']
        if planY == '':
            planY = fp[groupName]['CoordinateFrame'].attrs['YChannel']

        with h5py.File(measFile, 'r+') as fm:
            gMeas = fm[groupName]['Lines']
            if measX == '':
                measX = fm[groupName]['CoordinateFrame'].attrs['XChannel']
            if measY == '':
                measY = fm[groupName]['CoordinateFrame'].attrs['YChannel']
            numLines = len(gMeas.items())

            message = ''
            num_lines_exceeded = 0
            total_num_excs = 0
            num_lines_unplanned = 0

            lines = list(gMeas.keys())
            for line in lines:
                planLine = f"{gMeas[line].attrs['PlannedLine']:.1f}"
                if planLine in gPlan:
                    xP = np.array(gPlan[planLine][planX])
                    yP = np.array(gPlan[planLine][planY])
                    xM = np.array(gMeas[line][measX])
                    yM = np.array(gMeas[line][measY])
                    max_deviation = 0.0

                    if known != '':
                        exceedances_known = True
                        exc_known = np.array(gMeas[line][known])

                    # rotate to line of x ~ 0 using line direction in radians
                    # if abs(yP[-1] - yP[0]) > epsilon:
                    dirn = np.arctan2((xP[-1] - xP[0]), (yP[-1] - yP[0]))
                    # else:
                    #     dirn = np.pi / 2.0    

                    [x, y] = _rotateCoords(xM - xP[0], yM - yP[0], -dirn)

                    num_fids_in_exceedance = 0
                    exceedance_in_line = False

                    for fid in range(0, len(x)):
                        # There is an exceedance ...
                        if np.abs(x[fid]) > allowance:
                            # If a new exceedance, then initialise variables;
                            if num_fids_in_exceedance == 0:
                                start_x = xM[fid]
                                start_y = yM[fid]
                                num_fids_in_exceedance = 1
                                max_deviation = abs(x[fid]) - allowance
                                this_exc_known = False

                            # Else increment and update on the current exceedance.
                            else:
                                num_fids_in_exceedance += 1
                                max_deviation = max(np.abs(x[fid]) - allowance, max_deviation)
                                if exceedances_known:
                                    if exc_known[fid] > 0:
                                        this_exc_known = True

                        else:
                            if num_fids_in_exceedance > 0: # the current exceedance has ended
                                end_x = xM[fid]
                                end_y = yM[fid]
                                len_exceedance = _dist([start_x, end_x], [start_y, end_y])[1]
                                if _exceedance_fail(num_fids_in_exceedance, len_exceedance, maxCounter, maxDistance):
                                    message += f'\nL {line} deviates more than {allowance} m for '
                                    message += f'{num_fids_in_exceedance} fids ({len_exceedance:.0f} m), max exceedance = {max_deviation:.0f} m.'
                                    message += f'\n  From ({start_x:.0f} E {start_y:.0f} N) to ({end_x:.0f} E {end_y:.0f} N).'
                                    total_num_excs += 1
                                    if this_exc_known:
                                        number_known += 1
                                        message += ' Known exceedance.'
                                        this_exc_known = False
                                    exceedance_in_line = True
                                num_fids_in_exceedance = 0
                else:
                    print(f'Line {line} / {planLine} not in plan.')
                    num_lines_unplanned += 1
                if exceedance_in_line:
                    num_lines_exceeded += 1
                    if plot_flag:
                        _plot_exceeding_line(x, y, allowance, line, planLine, dirn)

            message = f'\n{num_lines_exceeded} lines with horizontal exceedances.\n' + message # 5 DEC
            message = f'\n{total_num_excs} horizontal exceedances.\n' + message # 5 DEC
            message = f'\n{num_lines_unplanned} lines not in plan and not checked.\n' + message # 5 DEC
            message = f'\n{number_known} exceedances known in the database.\n' + message # 5 DEC
            print(message)
            if plot_flag:
                plt.show()


def _exceedance_fail(num_fids_in_exceedance, len_exceedance, maxCounter, maxDistance):
    '''
    Given a specification, maxCounter, on number of fids and, maxDistance,
    on distance, checks to see if either  num_fids_in_exceedance > maxCounter
    or len_exceedance > maxDistance. If either is True, returns True.

    Parameters
    ----------
    num_fids_in_exceedance : Integer
        DESCRIPTION.
    len_exceedance : Float
        DESCRIPTION.
    maxCounter : Integer
        DESCRIPTION.
    maxDistance : Float
        DESCRIPTION.

    Returns
    -------
    Bool.

    '''
    if maxCounter < 1 and maxDistance < 1:
        return False
    if maxCounter > 0:
        if num_fids_in_exceedance > maxCounter:
            return True
    if maxDistance > 0:
        if len_exceedance > maxDistance:
            return True
    return False


def _plot_exceeding_line(x, y, allowance, line, planLine, dirn):
    '''
    Plots a standard exceedance figure for checkXYPlan().

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    allowance : TYPE
        DESCRIPTION.
    line : TYPE
        DESCRIPTION.
    planLine : TYPE
        DESCRIPTION.
    dirn : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(x, y, 'b')
    ax.plot(-allowance * np.ones(x.shape), y, 'r')
    ax.plot(allowance * np.ones(x.shape), y, 'r')
    plt.xlabel('deviation from planned line [m]', fontsize = 9)
    plt.ylabel('distance along line', fontsize = 9)
    plt.title(f'Line {line} (plan: {planLine}); Bearing {dirn * 180 / np.pi:.1f} deg'\
        , fontsize = 10)
    plt.grid()
    return

   
def checkVertPlan(planPath, measPath, planX='', planY='', planZ='', measX='',
                  measY='', measZ='', allowance=30.0, maxCounter=13,
                  maxDistance=0.0, known='', plot_flag=False):
    '''
    Reports exceedances of actual vertical position from planned vertical positions
    for an airborne survey Whizz database.
    The positions (`planX`, `planY`, `planY`) of each planned survey line
    are read from `planPath`. The measured positions (`measX`, `measY`, `measZ`) are read from
    `measPath` and the vertical distance of each from the planned line is
    calculated. If this distance exceeds `allowance` for a distance greater than
    `maxDistance` m, then an out-of-specification exceedance is reported for that line.
    If `maxDistance` is less than 1.0, then the test is instead against `maxCounter`
    consecutive positions. The default for `maxCounter` is 13.

    Parameters
    ----------
    planPath : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, of survey plan.
    planX : String, optional
        The name of the geoWhizz field or channel containing the planned x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    planY : String, optional
        The name of the geoWhizz field or channel containing the planned y positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    measPath : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, of measured data.
    measX : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    measY : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    allowance : Float, optional
        DESCRIPTION. The default is 30.0.
    maxCounter : Int, optional
        The maximum number of consecutive fids for which an exceedance
        greater than allowance is permitted. The default is 13.
    maxDistance : Float, optional
        The maximum number of consecutive metres for which an exceedance
        greater than allowance is permitted. If 0, then `maxCounter` is
        used instead of `maxDistance`. The default is 0.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line.

    Returns
    -------
    None.

    TODO: If maxCounter unspecified, and maxDistance is specified, test against maxDistance.
          Use _exceedance_fail() to do this.

    '''
    planfile = str(planPath)
    measFile = str(measPath)
    if maxDistance > 1.0:
        counting = False
    elif maxCounter > 0:
        counting = True
    else:
        print(f'ERROR. Require maxDistance > 1.0 m or maxCounter > 0')
        print(f'  maxDistance = {maxDistance}, maxCounter = {maxCounter}.')
    

    exceedances_known = False
    this_exc_known = False
    number_exc_known = 0

    with h5py.File(planfile, 'r') as fp:
        gPlan = fp[groupName]['Lines']
        if planX == '':
            planX = fp[groupName]['CoordinateFrame'].attrs['XChannel']
        if planY == '':
            planY = fp[groupName]['CoordinateFrame'].attrs['YChannel']
        if planZ == '':
            planZ = fp[groupName]['CoordinateFrame'].attrs['AltitudeChannel']

        with h5py.File(measFile, 'r') as fm:
            gMeas = fm[groupName]['Lines']
            if measX == '':
                measX = fm[groupName]['CoordinateFrame'].attrs['XChannel']
            if measY == '':
                measY = fm[groupName]['CoordinateFrame'].attrs['YChannel']
            if measZ == '':
                measZ = fm[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
            numErrors = int(0)
            numErrLines = 0
            num_lines_unplanned = 0
            report = ''

            for line in gMeas.keys():
                line_flagged = False
                planLine = f"{gMeas[line].attrs['PlannedLine']:.1f}"
                if planLine in gPlan: # 5 DEC
                    xP = np.array(gPlan[planLine][planX])
                    yP = np.array(gPlan[planLine][planY])
                    zP = np.array(gPlan[planLine][planZ])
                    xM = np.array(gMeas[line][measX])
                    yM = np.array(gMeas[line][measY])
                    zM = np.array(gMeas[line][measZ])

                    if known != '':
                        exceedances_known = True
                        exc_known = np.array(gMeas[line][known])
                    
                    # make life easier by transforming to a 2D problem
                    dirn = np.arctan((xP[-1] - xP[0])/(yP[-1] - yP[0]))
                    (xm, ym) = _rotateCoords(xM - xP[0], yM - yP[0], -dirn)
                    (xp0, yp0) = _rotateCoords(xP - xP[0], yP - yP[0], -dirn)
                    xp = xp0
                    yp = yp0
                    zp = zP
                    
                    # interpolate (ym, zM) onto (yp, zmp)
                    # This section dodgy - TODO : work out what the signs should be and use them
                    # TODO : better code to clean out the input data
                    if line == '9100002.0' or line == '9100003.0' or line == '9100100.0':
                        decreasing = yp0[-1] < yp0[0]
                        if decreasing:
                            jj = 0
                            for ii in range(len(yp0), 1, -1):
                                if yp0[ii] < yp0[ii-1]:
                                    xp[jj] = xp0[ii]
                                    yp[jj] = yp0[ii]
                                    zp[jj] = zP[ii]
                                    jj += 1
                            xp = xp[0:jj]
                            yp = yp[0:jj]
                            zp = zp[0:jj]
                        else:
                            jj = 0
                            for ii in range(1, len(yp0), 1):
                                if yp0[ii] > yp0[ii-1]:
                                    xp[jj] = xp0[ii]
                                    yp[jj] = yp0[ii]
                                    zp[jj] = zP[ii]
                                    jj += 1
                            xp = xp[0:jj]
                            yp = yp[0:jj]
                            zp = zp[0:jj]
                        (zmp, zM_trim) = mhd.interpolateLine(yp, zp, -ym, zM)
                    else:
                        (zmp, zM_trim) = mhd.interpolateLine(-yp, zp, -ym, zM, plot_flag=False)
                    # calculate the deviation vector
                    z_dev = zM_trim - zmp
                
                    # check vertical deviations
                    # initialise fiducial counter and error report for line, assume no exceedances
                    fid = int(0)
                    exceeding = False
                
                    for one_x in z_dev:
                        fid += 1
                        in_spec = abs(one_x) < allowance
                        
                        if not (exceeding or in_spec):
                            start_fid = fid
                            exceeding = True
                            exc_fids = 0
                            this_exc_known = False
                        elif exceeding and not in_spec:
                            exc_fids += 1
                            if exceedances_known:
                                if exc_known[fid] > 0:
                                    this_exc_known = True
                        elif exceeding and in_spec:
                            num_fids = exc_fids
                            exceeding = False
                            ex0 = xM[start_fid]
                            ex1 = xM[start_fid + num_fids]
                            ey0 = yM[start_fid]
                            ey1 = yM[start_fid + num_fids]
                            exc_dist = _displacement2(ex0, ex1, ey0, ey1)
                            if (counting and num_fids > maxCounter) or (not counting and exc_dist > maxDistance):
                                if not line_flagged:
                                    numErrLines += 1
                                    line_flagged = True
                                numErrors += 1
                                max_dev = np.nanmax(abs(z_dev[start_fid:start_fid + num_fids]))
                                report += f'\nL {line} deviates more than {allowance:.1f} m for'
                                report += f' {num_fids} fids ({exc_dist:.0f} m),'
                                report += f' max exceedance = {max_dev - allowance:.1f} m.'
                                report += f'\n  From ({ex0:.0f} E, {ey0:.0f} N) to ({ex1:.0f} E, {ey1:.0f} N).'
                                if this_exc_known:
                                    number_exc_known += 1
                                    report += ' Known exceedance.'
                                    this_exc_known = False
                    if plot_flag and line_flagged:
                        fig = plt.figure()
                        ax = fig.add_subplot(1,1,1)
                        ax.plot(z_dev, ym[1:], 'b')
                        ax.plot(-allowance * np.ones(z_dev.shape), ym[1:], 'r')
                        ax.plot(allowance * np.ones(z_dev.shape), ym[1:], 'r')
                        plt.xlabel('deviation from planned drape [m]')
                        plt.ylabel('distance along line')
                        plt.title(f'Line {line}')
                        plt.grid()
                else:
                    print(f'Line {line} not in plan.')
                    num_lines_unplanned += 1

        print(f'\n{num_lines_unplanned} lines not in plan and not checked.\n')
        print(f'Total number of exceedances = {numErrors} over {numErrLines} erroneous lines.')
        print(f'\n{number_exc_known} exceedances known in the database.\n')
        print(report)
        if plot_flag:
            plt.show()
                                       
   
def _rotateCoords(x, y, angle):
    '''
    Rotates Cartesian vectors x and y, by `angle` to xr and yr.

    Parameters
    ----------
    x : numpy Float 1D array
        DESCRIPTION.
    y : numpy Float 1D array
        DESCRIPTION.
    angle : Float
        Angle in radians by which coordinates are to be rotated.

    Returns
    -------
    xr : numpy Float 1D array
        DESCRIPTION.
    yr : numpy Float 1D array
        DESCRIPTION.

    '''
    xr = np.cos(angle) * x + np.sin(angle) * y
    yr = np.sin(angle) * x - np.cos(angle) * y
    return xr, yr
 
    
def checkClearance(whizzFile, altitude, terrain, nominalClearance, allowance=20.0, maxDistance=1000.0, xChannel='', yChannel='', only_low=False):
    '''
    Checks the data from Whizz HDF5 file for height exceedances against a specification
    requiring heights to be within a relative range (allowance) about a nominal
    value (nominalClearance) over a particular distance (maxDistance).
    
    The clearance = altitude - terrain, must be within +/- allowance of nominalClearance.
    For each line, the maximum absolute deviation of clearance from nominalClearance is
    reported. 
    If any samples along the line exceed the allowance, then profiles are plotted
    for visual analysis. No use is made of maxDistance (yet).
    The positions (xChannel and yChannel) are only used for plotting.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    altitude : String
        The name of the absolute altitude or height field in the Lines group
        of the Whizz HDF5 file.
    terrain : String
        The name of the absolute terrain or DTM height field in the Lines group
        of the Whizz HDF5 file.
    nominalClearance : Float
        The planned clearance between above the terrain in metres.
    allowance : Float, optional
        The absolute maximum deviation allowed from the planned clearance in
        metres. The default is 20.0.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    only_low : String, optional
        If True, then only check for exceedances too close to the ground ("safety
        check"). The default is False.

    Returns
    -------
    None.

    '''
    
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        num_failed_lines = 0

        for line in g.keys():
            alt = np.array(g[line][altitude])
            dtm = np.array(g[line][terrain])
            clearance = alt[()] - dtm[()]
            deviation = nominalClearance - clearance
            if only_low:
                maxDeviation = np.max(deviation)
            else:
                maxDeviation = np.max(abs(deviation))
            # print(f'Line {line}: Max deviation from {nominalClearance:.0f} m clearance = {maxDeviation:.0f} m.')
            if maxDeviation > allowance:
                num_failed_lines += 1
                print(f'Clearance deviation of {maxDeviation:.0f} m on line {line}')
                x = np.array(g[line][xChannel])
                y = np.array(g[line][yChannel])
                distance = _dist(x, y)
                fig = plt.figure()
                ax = fig.add_subplot(2,1,1)
                ax.plot(distance, alt)
                ax.plot(distance, dtm)
                plt.xlabel('distance along line [m]', fontsize = 6)
                plt.ylabel('height', fontsize = 6)
                if only_low:
                    plotTitle = projName + ': Minimum Clearance Check ' + line
                else:
                    plotTitle = projName + ': Absolute Clearance Check ' + line
                plt.title(plotTitle, fontsize = 8)
                plt.legend([altitude, terrain])
                plt.grid(True)
                for label in ax.get_xticklabels(): label.set_fontsize(6)
                for label in ax.get_yticklabels(): label.set_fontsize(6)
                ax2 = fig.add_subplot(2,1,2)
                ax2.plot(distance, deviation)
                plt.ylabel('deviation', fontsize = 6)
                plt.xlabel('distance along line [m]', fontsize = 6)
                plt.grid(True)
                for label in ax2.get_xticklabels(): label.set_fontsize(6)
                for label in ax2.get_yticklabels(): label.set_fontsize(6)
                plt.show()
        print(f'Number of failed lines = {num_failed_lines}.')
        

def checkDrape(whizzFile, altitude, drape, warningClearance = 20.0, xChannel = '', yChannel = ''):
    '''
    Checks actual altitude flown against planned drape in drape channel.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    altitude : TYPE
        The name of the geoWhizz field or channel containing the measured altitudes.
    drape : TYPE
        The name of the geoWhizz field or channel containing the measured drape heights.
    warningClearance : Float, optional
        DESCRIPTION. The default is 20.0.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        
        for line in g.keys():
            alt = g[line][altitude]
            drape = g[line][drape]
            clearance = np.abs(alt[()] - drape[()])
            maxDeviation = np.max(clearance)
            print(f'Drape {maxDeviation} on line {line}')
            if maxDeviation > warningClearance:
                print(f'Drape deviation of {maxDeviation} on line {line}')
                x = g[line][xChannel]
                y = g[line][yChannel]
                distance = _dist(x, y)
                fig = plt.figure()
                ax = fig.add_subplot(2,1,1)
                ax.plot(distance, alt)
                ax.plot(distance, drape)
                plt.xlabel('distance along line [m]', fontsize = 6)
                plt.ylabel('height', fontsize = 6)
                plotTitle = projName + ': Clearance Check ' + line
                plt.title(plotTitle, fontsize = 8)
                plt.legend([altitude, drape])
                plt.grid(True)
                for label in ax.get_xticklabels(): label.set_fontsize(6)
                for label in ax.get_yticklabels(): label.set_fontsize(6)
                ax2 = fig.add_subplot(2,1,2)
                ax2.plot(distance, clearance)
                plt.ylabel('deviation', fontsize = 6)
                plt.xlabel('distance along line [m]', fontsize = 6)
                plt.grid(True)
                for label in ax2.get_xticklabels(): label.set_fontsize(6)
                for label in ax2.get_yticklabels(): label.set_fontsize(6)
                plt.show()
        

def _dist(x, y):
    '''
    Returns a 1D numpy array where the i-th element is the distance from
    (x[0], y[0]) to (x[i], y[i]).

    Parameters
    ----------
    x : numpy Float array
        DESCRIPTION.
    y : numpy Float array
        DESCRIPTION.

    Returns
    -------
    numpy Float array
        DESCRIPTION.

    '''
    return  np.sqrt((x - x[0]) * (x - x[0]) + (y - y[0]) * (y - y[0]))

 
def checkVertAcc(whizzFile):
    '''
    Uses numpy.diff to estimate the vertical acceleration from the vertical velocity
    data (units of m/s) in the channel called `vertvelocity` in whizzFile. Plots the
    result (units of milli-g) for each survey line.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.

    Returns
    -------
    None

    '''
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']

        numLines = len(g.items())
        chStd = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            data = g[line]['vertvelocity']
            accel = np.diff(data, n = 1)
            chStd[count] = int(np.std(accel) * 100.0)
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            ax.plot(accel, '.', ms=2)
            plt.ylabel('Acceleration', fontsize = 6)
            plotTitle = groupName + ': Vert Accel [m/s/s] ' + str(chStd[count]) + 'mg ' + line
            plt.title(plotTitle, fontsize = 8)
            plt.grid(True)
            for label in ax.get_xticklabels(): label.set_fontsize(6)
            for label in ax.get_yticklabels(): label.set_fontsize(6)
            plt.show()
            count += 1


def checkVertAccStats(whizzFile):
    '''
    Uses numpy.diff to estimate the vertical acceleration from the vertical velocity
    data (units of m/s) in the channel called `vertvelocity` in whizzFile. Plots the
    min, max, mean and stdev (units of milli-g) for each survey line as a function of
    line.
    Assumes samples are 1 sec apart.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.

    Returns
    -------
    None

    '''
    outlierRatio = 6.0
    specAcc = 1.0
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]

        numLines = len(g.items())
        chMin = np.zeros((numLines,))
        chMax = np.zeros((numLines,))
        chMean = np.zeros((numLines,))
        chStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            data = g[line]['vertvelocity']
            accel = np.diff(data, n = 1)
            lineNo[count] = line
            chMin[count] = np.min(accel)
            chMax[count] = np.max(accel)
            chMean[count] = np.mean(accel)
            chStd[count] = np.std(accel)
            maxDevRatio = (chMax[count] - chMean[count]) / chStd[count]
            minDevRatio = (chMean[count] - chMin[count]) / chStd[count]
            if chStd[count] > specAcc:
                print(lineNo[count], ': 100 mill-g exceedance ', chStd[count])
            if maxDevRatio > outlierRatio:
                print(lineNo[count], ': max outlier ratio ', maxDevRatio)
            if minDevRatio > outlierRatio:
                print(lineNo[count], ': min outlier ratio', minDevRatio)
            count += 1

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(lineNo, chMin, 'bo', mfc='w')
    ax.plot(lineNo, chMax, 'bo', mfc='w')
    ax.errorbar(lineNo, chMean, chStd, capsize=3, marker='s', c='blue')
    plt.ylabel('Acceleration [m/s/s]', fontsize = 6)
    plotTitle = groupName + ': Acceleration Stats'
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    plt.show()


def checkInlineSum(whizzFile):
    '''
    Estimates the inline sum for each sample on each line in an FTG whizzFile.
    Plots the min, max, mean and stdev (units of eotvos) for each survey line
    as a function of line.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.

    Returns
    -------
    None

    '''
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']

        numLines = len(g.items())
        chMin = np.zeros((numLines,))
        chMax = np.zeros((numLines,))
        chMean = np.zeros((numLines,))
        chStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            data = g[line]['Inline1_raw']
            data1 = data[()]
            data = g[line]['Inline2_raw']
            data2 = data[()]
            data = g[line]['Inline3_raw']
            data3 = data[()]
            ils_BP = _inLineSum(data1, data2, data3)
            print(line, ' STD BPF = ', np.std(ils_BP))
            lineNo[count] = line
            chMin[count] = np.min(ils_BP)
            chMax[count] = np.max(ils_BP)
            chMean[count] = np.mean(ils_BP)
            chStd[count] = np.std(ils_BP)
            count += 1

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(lineNo, chMin, 'bo', mfc='w')
        ax.plot(lineNo, chMax, 'bo', mfc='w')
        ax.errorbar(lineNo, chMean, chStd, capsize=3, marker='s', c='blue')
        plt.ylabel('Inline Sum', fontsize = 6)
        plotTitle = groupName + ': ' + ' Inline Sum Stats'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        plt.show()


def _butter_bandpass(lowcut, highcut, fs, order=5):
    '''
    Calculates the parameters of a Butterworth bandpass filter.

    Parameters
    ----------
    lowcut : Float
        The lowpass cutoff frequency or wavenumber.
    highcut : Float
        The highpass cutoff frequency or wavenumber.
    fs : Float
        The sample frequency or wavenumber.
    order : Integer, optional
        The filter order. Default 5.

    Returns
    -------
    b : Float
        DESCRIPTION.
    a : Float
        DESCRIPTION.

    '''
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    '''
    Uses a Butterworth bandpass filter [`lowcut`, `highcut`] to filter data
    (sampled at `fs`).

    Parameters
    ----------
    data : Numpy 1D array
        The input data.
    lowcut : Float
        The lowpass cutoff frequency or wavenumber.
    highcut : Float
        The highpass cutoff frequency or wavenumber.
    fs : Float
        The sample frequency or wavenumber.
    order : Integer, optional
        The filter order. Default 5.

    Returns
    -------
    Numpy 1D array
        The output data.

    '''
    b, a = _butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def checkBasemag(whizzFile, basemag, peak = 0.5, nSamples = 3000):
    '''
    Checks the basemag channel in a whizzFile against the specification that
    the peak to peak variation over a set number of samples must not exceed
    some peak value.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    basemag : String
        The name of the channel in whizzFile containing the basemag data.
    peak : Float
        The maximum allowed peak to peak variation.
    nSamples : Integer
        The number of samples (moving window) over which the test is applied.

    Returns
    -------
    None

    '''
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']

        for line in g.keys():
            dataFail = False
            # plotTitle = line + ' ' + basemag + ' Peak-to-peak'
            data = np.array(g[line][basemag])
            data = data[np.logical_not(np.isnan(data))]
            if len(data) < 2:
                print(line, ' insufficient data')
                continue

            if _peakToPeak(data) < 2.0 * peak:
                print(line, ' passed easily.')
                continue
            
            if len(data) < nSamples:
                    if _peakToPeak(data) > 2.0 * peak:
                        dataFail = True
                        print(line, ' FAIL, peak to peak range = ', _peakToPeak(data))
                        continue
                    else:
                        print(line, ' passed on one segment.')
                        continue
            else:
                for ii in range(0, len(data) - nSamples):
                    if _peakToPeak(data[ii:ii+nSamples]) > 2.0 * peak:
                        dataFail = True
                        print(line, ' FAIL, peak to peak range = ', _peakToPeak(data[ii:ii+nSamples]))
                        break
                if dataFail == False:
                    print(line, ' passed on many segments.')
                    continue
    return
            

def _peakToPeak(data):
    return np.max(data) - np.min(data)
    

def checkDiurnal(whizzFile, basemag, rangeLimit = 5.0, nSamples = 3000, diff4Limit = 0.5):
    # TODO: add check for singleValueExceedance()
    
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        _reportAllDiurnal(g, basemag, rangeLimit=rangeLimit, nSamples=nSamples, diff4Limit=diff4Limit)
                    

def checkLineDiurnal(whizzFile, line, basemag, rangeLimit = 5.0, nSamples = 3000, diff4Limit = 0.5):
    # TODO: add check for singleValueExceedance()
    # NOT USED??
    
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]

        _reportLineDiurnal(g, line, basemag, rangeLimit, nSamples, diff4Limit)


def _reportAllDiurnal(group, basemag, rangeLimit = 5.0, nSamples = 3000, diff4Limit = 0.5):
    # check 4th difference noise.
    for line in group.keys():
        _reportLineDiurnal(group, line, basemag, rangeLimit, nSamples, diff4Limit)
        
        
def _reportLineDiurnal(group, line, basemag, rangeLimit = 5.0, nSamples = 3000, diff4Limit = 0.5):
    # get data, check we have sufficient.
    diurnalExceeded = False
    failedSample = 0
    bigExtremum = 0.0
    data = np.array(group[line][basemag])
    data = data[np.logical_not(np.isnan(data))]
        
    if nSamples > len(data):
        print(f'\n  Short line: {len(data)} < {nSamples}.')
        nSamples = len(data)
    if (nSamples % 2) == 0:
        nSamples = nSamples - 1
    nSam = (nSamples - 1) // 2

    for ii in range(nSam, len(data)-nSam):
        localData = data[ii-nSam:ii+nSam]
        localSlope = (localData[-1] - localData[0]) / nSamples
        # deviation = np.zeros(localData.shape)
        deviation = localData - localSlope * range(0, len(localData)) - localData[0]
        # for jj in range(0, len(localData)):
        #     deviation[jj] = localData[jj] - localData[0] - localSlope * jj
        extremum = np.max(deviation) if np.max(deviation) > -np.min(deviation) else -np.min(deviation)
        if extremum > rangeLimit:
            diurnalExceeded = True
            if extremum > bigExtremum:
                bigExtremum = extremum
                failedSample = ii
            
    if diurnalExceeded:
        print(f'\n  Diurnal for {basemag} at sample number {failedSample} diverges from chord by {bigExtremum:.2f}, exceeding {rangeLimit:.1f} - FAIL')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(data)
        plotTitle = f'Line {line} Channel {basemag}: reaches {bigExtremum:.2f} at {failedSample}, exceeding {rangeLimit} - FAIL'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        fig.tight_layout()
        plt.show()
                            

def checkTCDiff4(whizzFile, TCDiff4='', rawMag='', limit = 0.02, nSamples = 3000, plotAll = False):
    '''
    Checks the total magnetic field fourth difference channel in a whizzFile
    against the specification that the peak to peak variation over a set
    number of samples must not exceed some peak value. If `TCDiff4` is not
    available, then the rawMag channel may be used.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    TCDiff4 : String, optional
        The name of the channel in whizzFile containing the 4th difference mag data.
        Default is '', in which case rawMag is used.
    rawMag : String, optional
        The name of the channel in whizzFile containing the raw magnetic data data.
        Default is ''; if both rawMag and TCDiff4 are '', then an error is reported.
    limit : Float, optional
        The maximum allowed peak to peak variation. Default = 0.02
    nSamples : Integer, optional
        The number of samples (moving window) over which the test is applied.
        Default = 3000
    plotAll : Bool, optional
        If True, all plots are generated.

    Returns
    -------
    None

    '''
    
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        _groupMagDiff4(g, TCDiff4=TCDiff4, rawMag=rawMag, limit = limit, nSamples = nSamples, plotAll = plotAll)
        
        
def _groupMagDiff4(group, TCDiff4='', rawMag='', limit = 0.02, nSamples = 3000, plotAll = False):
    '''
    Checks the total magnetic field fourth difference channel in an HDF Group
    against the specification that the peak to peak variation over a set
    number of samples must not exceed some peak value. If `TCDiff4` is not
    available, then the rawMag channel may be used.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    TCDiff4 : String, optional
        The name of the channel in whizzFile containing the 4th difference mag data.
        Default is '', in which case rawMag is used.
    rawMag : String, optional
        The name of the channel in whizzFile containing the raw magnetic data data.
        Default is ''; if both rawMag and TCDiff4 are '', then an error is reported.
    limit : Float, optional
        The maximum allowed peak to peak variation. Default = 0.02
    nSamples : Integer, optional
        The number of samples (moving window) over which the test is applied.
        Default = 3000
    plotAll : Bool, optional
        If True, all plots are generated.

    Returns
    -------
    None

    '''
    # check 4th difference noise.
    for line in group.keys():
        if TCDiff4 == '':
            if rawMag == '':
                print('ERROR - no rawmag or 4th difference channel name supplied.')
            else:
                mag = np.array(group[line][rawMag])
                md4 = np.diff(mag, n=4)
                data = np.append(np.append(md4[0:2],md4),md4[-3:-1])
        else:
            data = np.array(group[line][TCDiff4])
        rangeTooHigh = False
        plotTitle = line + ' ' + TCDiff4 + ' Range'
        
        rangeTooHigh, numExceedances = failsDeviation(data, limit, nSamples)
        
        if rangeTooHigh:
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            ax.plot(data)
            plt.title(plotTitle)
            plt.grid(True)
            plt.show()
            print(line, ' ', numExceedances, ' > ', nSamples, ' - FAIL')
          
        elif numExceedances > 0:
            print(line, ' ', numExceedances, ' < ', nSamples, ' - PASS')
            if plotAll:
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                ax.plot(data)
                plt.title(plotTitle)
                plt.grid(True)
                plt.show()
            
        else:
            print(line, ' - PASS')


def failsDeviation(data, limit, nSamples):
    """
    Checks data, testing for more than `nSamples` consecutive instances where
    the data values are outside the range [-limit, limit].
    
    Parameters
    ----------
    data : Numpy 1D Float array
        The data to be checked.
    limit : Float
        The limit, assumed > 0.
    nSamples : Float
        Maximum allowed number of deviations allowed.

    Returns
    -------
    Bool
        True if the data failed.
    Float
        The number of exceedances (regardless of whether consecutive).
    
    """
    
    if np.max(data) < limit and np.min(data) > -limit:
        print('all ok')
        return False, 0
    
    indxSpikes = np.argwhere(np.logical_or(data > limit, data < -limit))
    numExceedances = indxSpikes.shape[0]
    if numExceedances < nSamples:
        print(f'only {numExceedances} exceedances')
        return False, numExceedances
    
    else:
        count = 1
        for jj in range(1, len(indxSpikes)):
            if indxSpikes[jj][0] == indxSpikes[jj-1][0] + 1:
                count += 1
            else:
                if count >= nSamples:
                    return True, count
                else:
                    count = 1
        if count >= nSamples:
            return True, count
        print(f'{numExceedances} but not consecutive, most {count}')
        return False, numExceedances
        

def failsPeak(data, limit):
    '''
    NOT USED
    '''
    if np.max(data) > set or np.min(data) < -limit:
        return False
    else:
        return True


def checkSpikes(whizzFile, fields = [], numStd = 8.0):
    """
    Checks for spikes in all the given fields of data.
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    fields : String List
        List of field names from the database to be checked.
    numStd : Float, optional
        maximum allowed number of standard deviations allowed. The default is 8.0.

    Returns
    -------
    None

    """
    
    filename = str(whizzFile)
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if fields == []:
            lineGroups = list(g.values())
            fields = list(lineGroups[0].keys())
        noSpikes = True
        report = ''
        for line in g.keys():
            for field in fields:
                deriv = np.diff(g[line][field], n = 1)
                deriv = deriv - np.mean(deriv)
                if len(deriv) > 10:
                    myStd = np.std(deriv)
                    extremum = np.max(deriv) if np.max(deriv) > -np.min(deriv) else -np.min(deriv)
                    if extremum > numStd * myStd:
                        report += f'\n  {line}; {field} Extremum: {extremum:.2f} > {(numStd * myStd):.2f} = {numStd:.2f} x STD of {myStd:.2f}'
                        noSpikes = False
                else:
                    report += f'\n  {line}; {field} data length {len(deriv)+1} is too short for analysis.'
                    noSpikes = False
    print(report)
    return


def checkConstantSlope(whizzFile, fields=[]):
    """
    Checks for constant slope (`np.diff`) in all the given fields of data.
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    fields : String List
        List of field names from the database to be checked.

    Returns
    -------
    None

    """
    filename = str(whizzFile)
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if fields == []:
            lineGroups = list(g.values())
            fields = list(lineGroups[0].keys())
        data_is_good = True
        report = ''
        for line in g.keys():
            for field in fields:
                deriv = np.diff(g[line][field], n = 1)
                mean_deriv = np.mean(deriv)
                deriv = deriv - mean_deriv
                if len(deriv) > 10:
                    extremum = np.max(deriv) if np.max(deriv) > -np.min(deriv) else -np.min(deriv)
                    if extremum > mean_deriv / 100.0:
                        report += f'\n  {line}; {field} Largest difference (= {extremum:.2f}) > 1% of mean difference (= {(mean_deriv / 100.0):.2f})'
                        data_is_good = False
                else:
                    report += f'\n  {line}; {field} data length {len(deriv)+1} is too short for analysis.'
                    data_is_good = False
    if data_is_good:
        report += 'All fields tested were either constant or of constant slope for all lines tested.'
    print(report)
    return


def checkSpeeds(whizzFile, xChannel='', yChannel='', tChannel='', vel_north='', vel_east='', nominalSpeed=60.0, 
    maxDuration=0.0, maxDistance=0.0, allowance=0.1, minSafeSpeed=42.0, plot_flag=False):
    '''
    Checks the data from Whizz HDF5 file for speed exceedances against a specification
    requiring ground speeds to be within a relative range (allowance) about a nominal
    value (nominalSpeed) over a particular distance (maxDistance).
    
    The positions (xChannel and yChannel) are assumed to be sampled uniformly in time
    and the first two time (tChannel) values for the first line in the file are 
    differenced to obtain the sampling interval. From this, the number of samples, N, 
    for the aircraft to travel maxDistance is calculated. The algorithm compares the
    actual distance flown in N samples with maxDistance. If the actual distance is
    greater than
        maxDistance + allowance * maxDistance
    or less than 
        maxDistance - allowance * maxDistance
    an error is printed.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    tChannel : String, optional
        The name of the geoWhizz field or channel containing the measured times. The
        default is to read the tChannel field name from the Coordinate Frame.
    nominalSpeed : Float, optional
        The specified ground speed in m/s. The default is 60.0.
    maxDuration : Float, optional
        The time in seconds over which the speed estimate is determined.
        The default is 13.3.
    maxDistance : Float, optional
        The distance in metres over which the speed estimate is determined.
        The default is 1000.0.
    allowance : Float, optional
        The magnitude, relative to the nominalSpeed, of the range of allowed
        speeds. The default is 0.1 (i.e. +/- 10% of nominal).
    minSafeSpeed : Float, optional
        The minimum allowed instantaneous safe speed in m/s. The default is 42.0.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line. Default False.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)
    if maxDistance > 1.0:
        check_by_time = False
    elif maxDuration > 1.0:
        check_by_time = True
    else:
        print(f'ERROR. Require maxDistance > 1.0 m or maxDuration > 1.0')
        print(f'  maxDistance = {maxDistance}, maxDuration = {maxDuration}.')
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if tChannel == '':
            tChannel = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
        title_str = ''
        project = f[groupName].attrs['ProjectName']
        block = f[groupName].attrs['BlockID']
        if project != '':
            title_str += f'{project}'
        if block != '':
            title_str += f' {block}'
        _reportSpeeds(g, maxDuration=maxDuration, maxDistance=maxDistance, xChannel=xChannel, yChannel=yChannel,
                tChannel=tChannel, vel_north=vel_north, vel_east=vel_east, nominalSpeed=nominalSpeed, 
                     allowance=allowance, minSafeSpeed=minSafeSpeed, title_str=title_str,
                     plot_flag=plot_flag)
        

def _reportSpeeds(group, maxDuration=0.0, maxDistance=0.0, xChannel='X', yChannel='Y',
                tChannel='time', vel_north='', vel_east='', nominalSpeed=60.0,
                allowance=0.1, minSafeSpeed=42.0, title_str='', plot_flag=False):
    '''
    Checks the data from Whizz Line group for speed exceedances against a specification
    requiring ground speeds to be within a relative range (`allowance`) about a nominal
    value (`nominalSpeed`) over a specified distance (`maxDistance`). If `maxDistance` is
    not provided, and `maxDuration` is provided, then the check is instead over the
    specified duration.
    
    If `vel_north` and `vel_east` are both provided, then the data in those channels
    is used to calculate speed along the line. If not, then speeds are calculated
    by differencing the positions (`xChannel` and `yChannel`) and dividing by the
    difference of the first two time (`tChannel`) values for the first line in the
    group. This assumes that data are sampled uniformly in time.
    
    From this, the number of samples, N, 
    for the aircraft to travel maxDistance is calculated. The algorithm compares the
    actual distance flown in N samples with maxDistance. If the actual distance is
    greater than
        maxDistance + allowance * maxDistance
    or less than 
        maxDistance - allowance * maxDistance
    an error is printed.

    Parameters
    ----------
    group : HDF5 Group
        The Whizz line group containing the survey line data.
    maxDuration : Float, optional
        The time in seconds over which the speed estimate is determined.
        The default is 0.0.
    maxDistance : Float, optional
        The distance in metres over which the speed estimate is determined.
        The default is 0.0.
    xChannel : String, optional
        The field in the line containing X positions. The default is 'X'.
    yChannel : String, optional
        The field in the line containing Y positions. The default is 'Y'.
    tChannel : String, optional
        The field in the line containing sample times. The default is 'time'.
    vel_north : String, optional
        The field in the line containing velocity north. The default is ''.
    vel_east : String, optional
        The field in the line containing velocity east. The default is ''.
    nominalSpeed : Float, optional
        The specified ground speed in m/s. The default is 60.0.
    allowance : Float, optional
        The magnitude, relative to the nominalSpeed, of the range of allowed
        speeds. The default is 0.1 (i.e. +/- 10% of nominal).
    minSafeSpeed : Float, optional
        The minimum allowed instantaneous safe speed in m/s. The default is 42.0.
    title_str : String, optional
        A title string for the plots. Default ''.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line. Default False.

    Returns
    -------
    None.

    '''
    calc_from_pos = False
    if vel_north == '' or vel_east == '':
        calc_from_pos = True
        print('Velocities not known - will calculate from positions')

    check_against_dist = False
    max_allowed_str = f'{maxDuration:.1f} s'
    if maxDuration < 1.0:
        check_against_dist = True
        max_allowed_str = f'{maxDistance:.1f} m'

    settings = f'Nominal ground speed {nominalSpeed:.1f} m/s; '
    settings += f'allowed {nominalSpeed-allowance*nominalSpeed:.1f} : '
    settings += f'{nominalSpeed+allowance*nominalSpeed:.1f} for < {max_allowed_str}.\n'
    print(settings)
    
    num_failed_lines = 0
    num_exceed_lines = 0
    total_num_lines = 0
    report = ''

    lines = list(group.keys())
    total_num_lines = len(lines)

    for line in lines:
        if title_str == '':
            plot_title = f'Line {line}'
        else:
            plot_title = f'{title_str}: Line {line}'

        x, y, dist, t, speed = _get_data(group[line], xChannel, yChannel, tChannel, vel_north, vel_east)

        if speed[speed < minSafeSpeed].size > 0:
            lineUnsafeSlow = True
            print(f'\n For at least one reading, the ground speed was < {minSafeSpeed} (might be unsafe).')
        
        rel_speed = (speed - nominalSpeed) / nominalSpeed

        speed_extreme = 0.0
        num_fids_in_exceedance = 0
        exceedance_in_line = False
        too_slow = False
        line_fails = False

        for fid in range(0, len(x)):
            # There is a speed exceedance ...
            if rel_speed[fid] > allowance or rel_speed[fid] < -allowance:
                # If a new exceedance, then initialise variables;
                if num_fids_in_exceedance == 0:
                    num_exceed_lines += 1
                    if rel_speed[fid] < -allowance:
                        too_slow = True
                    else:
                        too_slow = False
                    start_x = x[fid]
                    start_y = y[fid]
                    start_t = t[fid]
                    num_fids_in_exceedance = 1
                    speed_extreme = speed[fid]
                # Else increment and update on the current exceedance.
                else:
                    # check we haven't swapped from too slow to too fast or vice versa
                    if (too_slow and rel_speed[fid] > allowance) or ((not too_slow) and rel_speed[fid] < allowance):
                        print(f'WARNING: Exceedance reversed speed in one fid. NOT BELIEVABLE. At time={t[fid]:.3f}')
                    num_fids_in_exceedance += 1
                    if too_slow:
                        speed_extreme = min(speed[fid], speed_extreme)
                    else:
                        speed_extreme = max(speed[fid], speed_extreme)
            else:
                if num_fids_in_exceedance > 0: # the current exceedance has ended
                    end_x = x[fid]
                    end_y = y[fid]
                    end_t = t[fid]
                    dist_exceedance = _dist([start_x, end_x], [start_y, end_y])[1]
                    durn_exceedance = end_t - start_t
                    if too_slow:
                        speed_msg = "too slow"
                    else:
                        speed_msg = "too fast"
                    if _exceedance_fail(durn_exceedance, dist_exceedance, maxDuration, maxDistance):
                        report += f'\nL {line} {speed_msg} for {durn_exceedance} sec '
                        report += f'({dist_exceedance:.0f} m), peak exceedance = {speed_extreme:.0f} m/s.'
                        report += f'\n  From ({start_x:.0f} E {start_y:.0f} N) to ({end_x:.0f} E {end_y:.0f} N).'
                        exceedance_in_line = True
                    num_fids_in_exceedance = 0
                    too_slow = False
        if exceedance_in_line:
            num_failed_lines += 1
            if plot_flag:
                _plot_speed(t, speed, nominalSpeed * (1.0 - allowance), 
                   nominalSpeed * (1.0 + allowance), plot_title=plot_title)

    print(f'Checked {total_num_lines} lines and {num_exceed_lines} had some short exceedance(s).')
    print(f'{num_failed_lines} lines failed for exceedance > allowed.')
    print(report)
    if plot_flag:
        plt.show()


def _get_data(line_group, xChannel, yChannel, tChannel, vel_north, vel_east):
    '''
    '''
    x = np.array(line_group[xChannel])
    y = np.array(line_group[yChannel])
    t = np.array(line_group[tChannel])
    distance = np.sqrt(x * x + y * y)
    if vel_north == '' or vel_east == '':
        sampleTime = t[1] - t[0]
        xVel = np.diff(x) / sampleTime
        yVel = np.diff(y) / sampleTime
        temp = np.sqrt(xVel * xVel + yVel * yVel)
        speed = np.append(temp, np.mean(temp))
    else:
        xVel = np.array(line_group[vel_east])
        yVel = np.array(line_group[vel_north])
        speed = np.sqrt(xVel * xVel + yVel * yVel)

    return x, y, distance, t, speed


def _plot_speed(t, speed, min_speed=54, max_speed=66, plot_title=''):
    '''
    Plots the speed as a function of t and adds limit lines showing the allowed
    tolerance

    Parameters
    ----------
    t : TYPE
        Time (sec) data vector.
    speed : TYPE
        The speed data vector.
    min_speed : TYPE, optional
        DESCRIPTION. The default is 54.
    max_speed : TYPE, optional
        DESCRIPTION. The default is 66.
    plot_title : TYPE, optional
        DESCRIPTION. The default is ''.

    Returns
    -------
    None.

    '''
    fig = plt.figure()
    fig.suptitle(plot_title, fontsize=10)
    fig.subplots_adjust(top=0.85)
    
    ax = fig.add_subplot(1,1,1)
    ax.plot(t, speed, 'b', t, np.ones(t.size) * min_speed, 'r', t, np.ones(t.size) * max_speed, 'r', mfc='w')
    plt.xlabel('Time [s]', fontsize = 6)
    plt.ylabel('Speed [m/s]', fontsize = 6)
    plotTitle = 'Speed' + ' Stats'
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    
    
def _displacement2(x0, x1, y0, y1):
    '''
    The displacement distance from (xo, y0) to (x1, y1).

    Parameters
    ----------
    x0 : TYPE
        DESCRIPTION.
    x1 : TYPE
        DESCRIPTION.
    y0 : TYPE
        DESCRIPTION.
    y1 : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    return np.sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1))
            

def checkGaps(whizzFile):
    '''
    Checks every dataset for each channel and each survey line in filePath for
    gaps, and reports all gaps found.
    TODO: modify to allow gaps smaller than some minimum size.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        _reportGaps(g)
        
        
def _reportGaps(group):
    '''
    Checks every dataset for each channel and each survey line in the HDF5 group
    for gaps, and reports all gaps found.
    TODO: modify to allow gaps smaller than some minimum size.

    Parameters
    ----------
    group : HDF5 Whizz file 'Lines' group
        The group containing the survey line data.

    Returns
    -------
    None.

    '''
    lineGroups = list(group.values())
    channelNames = list(lineGroups[0].keys())
    num_channels = len(channelNames)
    num_lines_failed = 0
    total_num_lines = 0
    message = ''

    for line in group.keys():
        total_num_lines += 1
        gaps_on_line = False
        lineNo = line
        lineText = f'Line {lineNo}'
        for channel in channelNames:
            numberMissing = np.count_nonzero(np.isnan(group[line][channel]))
            if numberMissing > 0:
                lineText += f'\n    {channel}, nans: {numberMissing}'
                gaps_on_line = True
        if gaps_on_line:
            num_lines_failed += 1
            message += lineText + '\n'
    print(f'Checking for all gaps in all {num_channels} channels on all {total_num_lines} lines.')
    print(message)
    print(f'{num_lines_failed} lines failed.')


def eigen(G):
    '''
    NOT USED
    '''
    w, v = np.linalg.eig(G)
    
    print(w)
    print(v)


def checkRepeatLines(whizzFile, channel, flightLines=[], x='', z='', xOffset=True):
    '''
    For all lines in flightLines (assumed to be repeats), plot (x, channel) and
    report stats of differences to mean. This will require trimming to
    [minX, maxX] and interpolating to common x.

    TODO : optionally analyse d(channel)/dx as well as channel.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    channel : TYPE
        DESCRIPTION.
    flightLines : TYPE
        DESCRIPTION.
    x : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.
    xOffset : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None

    '''

    nSamples = 0
    minBigX = 1.0E12
    maxSmallX = -1.0E12
    
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if z == '':
            z = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        if flightLines == []:
            flightLines = list(g.keys())
        nLines = len(flightLines)
        baseLine = g[flightLines[0]].attrs['PlannedLine']
        #if the channel has an attribute 'Units'
        dd = g[flightLines[0]][channel]
        chan_y_label = channel
        if 'Units' in dd.attrs.keys():
            chan_y_units = dd.attrs['Units']
            chan_y_label += ' ' + chan_y_units
        ddz = g[flightLines[0]][z]
        if 'Units' in ddz.attrs.keys():
            chan_z_units = ddz.attrs['Units']

        # nSamples is the array width for data storage
        for line in flightLines:
            xs = np.array(g[line][x])
            nSamples = max(nSamples, xs.size)
            minBigX = min(max(xs), minBigX)
            maxSmallX = max(min(xs), maxSmallX)
            deltaX = np.abs(xs[1] - xs[0])
            # print(f'line {line}, minBigX {minBigX}, maxSmallX {maxSmallX}, deltaX {deltaX}, size {xs.size}')
            
        if minBigX < maxSmallX:
            return 0.0
        xBase = np.linspace(maxSmallX, minBigX, num=nSamples, endpoint=True)#np.arange(maxSmallX, minBigX, deltaX, 'float')
        print(f'{nLines} lines analysed, each with {nSamples} samples.')
        xData = np.empty((nLines, nSamples))
        xData[:] = np.nan
        yData = np.empty((nLines, nSamples))
        yData[:] = np.nan
        zData = np.empty((nLines, nSamples))
        zData[:] = np.nan
        lineCount = 0
        
        # read the data into the arrays
        for line in flightLines:
            xd = np.array(g[line][x])
            yd = np.array(g[line][channel])
            zd = np.array(g[line][z])
            # ensure ordered in increasing x
            if xd[1] < xd[0]:
                xd = xd[::-1]
                yd = yd[::-1]
                zd = zd[::-1]

            xStart = 0
            xEnd = xd.size - 1
            
            # trim data and store
            for xSample in range(0, xd.size):
                if xd[xSample] < (maxSmallX - deltaX / 2.0):
                    xStart = max(xSample, xStart)
                else:
                    break
            for xSample in range(xd.size-1, 0, -1):
                if xd[xSample] > (minBigX + deltaX / 2.0):
                    xEnd = min(xSample, xEnd)
                else:
                    break
                    
            # interpolate data
            (yOut, _) = mhd.interpolateLine(xd-xBase[0], yd, xBase-xBase[0])
            (zOut, _) = mhd.interpolateLine(xd-xBase[0], zd, xBase-xBase[0])

            vec_len = len(xBase)-1 # interpolateLine has lost a datapoint in outputs
            # print(f'line {line}, shapes: xBase {xBase.shape}, xData {xData.shape}')
            xData[lineCount, 0:vec_len] = xBase[1:]
            yData[lineCount, 0:vec_len] = yOut
            zData[lineCount, 0:vec_len] = zOut
            lineCount += 1
        
        xPlot = xBase
        if xOffset:
                xPlot = xPlot - xPlot[0]
        
    fig = plt.figure(figsize=(12,9))
    
    #plot the y data
    ax = fig.add_subplot(2,2,1)
    for line in range(0, nLines):
        y1 = yData[line,:]
        y1 = y1[np.logical_not(np.isnan(y1))]
        x1 = xData[line,:]
        x1 = x1[np.logical_not(np.isnan(x1))]
        myPlot, = ax.plot(x1, y1, lw=0.5, label=f'Line {flightLines[line]}')
            
    ax.legend(fontsize=8)
    plt.xlabel('x', fontsize = 10)
    plt.ylabel(chan_y_label, fontsize = 10)
    plotTitle = f'{baseLine:.0f} Repeat Lines {channel}'
    plt.title(plotTitle, fontsize = 12)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    
    #plot the z data
    ax = fig.add_subplot(2,2,2)
    for line in range(0, nLines):
        z1 = zData[line,:]
        z1 = z1[np.logical_not(np.isnan(z1))]
        x1 = xData[line,:]
        x1 = x1[np.logical_not(np.isnan(x1))]
        myPlot, = ax.plot(x1, z1, lw=0.5, label=f'Line {flightLines[line]}')
        ax.legend(fontsize=8)
            
    plt.xlabel('x', fontsize = 10)
    plt.ylabel(f'z {chan_z_units}', fontsize = 10)
    plotTitle = f'{baseLine:.0f} Repeat Lines {z}'
    plt.title(plotTitle, fontsize = 12)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    
    # plot the y differences and RMS
    ax = fig.add_subplot(2,2,3)
    yMean = np.mean(yData, axis=0)
    ySum = np.zeros(yMean.shape)
    for line in range(0, nLines):
        ySum = ySum + yData[line,:] - yMean
        yStd = np.nanstd(yData[line,:]-yMean)
        y1 = yData[line,:]-yMean
        y1 = y1[np.logical_not(np.isnan(y1))]
        x1 = xData[line,:]
        x1 = x1[np.logical_not(np.isnan(x1))]

        myPlot, = ax.plot(x1, y1, lw=0.5, label=f'RMS = {yStd:.1f}')
        ax.legend(fontsize=8)
        print(f'Line {flightLines[line]}: stdev({channel}) = {yStd:.2f} {chan_y_units}')
            
    plt.xlabel('x', fontsize = 10)
    plt.ylabel(chan_y_label, fontsize = 10)
    plotTitle = f'Differences to mean and RMS {channel}'
    plt.title(plotTitle, fontsize = 12)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    
    # plot the z differences and RMS
    ax = fig.add_subplot(2,2,4)
    zMean = np.mean(zData, axis=0)
    zSum = np.zeros(zMean.shape)
    for line in range(0, nLines):
        zSum = zSum + zData[line,:] - zMean
        zStd = np.nanstd(zData[line,:]-zMean)
        z1 = zData[line,:]-zMean
        z1 = z1[np.logical_not(np.isnan(z1))]
        x1 = xData[line,:]
        x1 = x1[np.logical_not(np.isnan(x1))]

        myPlot, = ax.plot(x1, z1, lw=0.5, label=f'RMS = {zStd:.2f}')
        ax.legend(fontsize=8)
        # ax.text(x1[0], z1[0], f'RMS = {zStd:.2f}', fontsize=8)
        print(f'Line {flightLines[line]}: stdev({z}) = {zStd:.1f} {chan_z_units}')
            
    plt.xlabel('x', fontsize = 10)
    plt.ylabel(f'z {chan_z_units}', fontsize = 10)
    plotTitle = f'Differences to mean and RMS {z}'
    plt.title(plotTitle, fontsize = 12)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    
    plt.subplots_adjust(0.1, 0.05, 0.9, 0.95)
    plt.show()
            
    return


def allChanStats(whizzFile, allChannels=[]):
    '''
    Generate statistical plots for the channels across all lines. The plots show
    the min, mean, max and stdev for each channel as a function of line number.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    allChannels : String List, optional
        A list of the channels or fields to plot. Default is all in whizzFile.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        gProject = f[groupName]
        g = gProject['Lines']
        lines = list(g.keys())
        if allChannels == []:
            lineGroups = list(g.values())
            allChannels = list(lineGroups[0].keys())
        numLines = len(g.items())
        
        for channel in allChannels:
            chMin = np.zeros((numLines,))
            chMax = np.zeros((numLines,))
            chMean = np.zeros((numLines,))
            chStd = np.zeros((numLines,))
            lineNo = np.zeros((numLines,))
            count = 0
            dd = g[lines[0]][channel]
            #if the channel has an attribute 'Units'
            chan_y_label = channel
            if 'Units' in dd.attrs.keys():
                chan_y_label += ' ' + dd.attrs['Units']
        
            for line in g.keys():
                if line != 'CoordinateFrame':
                    lineNo[count] = line
                    if np.sum(~np.isnan(g[line][channel])) > 3:
                        chMin[count] = np.nanmin(g[line][channel])
                        chMax[count] = np.nanmax(g[line][channel])
                        chMean[count] = np.nanmean(g[line][channel])
                        chStd[count] = np.nanstd(g[line][channel])
                    else:
                        chMin[count] = 0.0
                        chMax[count] = 0.0
                        chMean[count] = 0.0
                        chStd[count] = 0.0
                        print(f'Less than three real values in {lineNo[count]:.2f} for {channel}, no statistics.')
                    count += 1

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            ax.plot(lineNo, chMin, 'bo', mfc='w')
            ax.plot(lineNo, chMax, 'bo', mfc='w')
            ax.errorbar(lineNo, chMean, chStd, capsize=3, marker='s', c='blue')
            plt.ylabel(chan_y_label, fontsize = 6)
            plotTitle = groupName + ': ' + channel + ' Stats'
            plt.title(plotTitle, fontsize = 8)
            plt.grid(True)
            for label in ax.get_xticklabels(): label.set_fontsize(6)
            for label in ax.get_yticklabels(): label.set_fontsize(6)
            plt.show()
    return
    
    
def lineStats(whizzFile, channel):
    '''
    Generate a statistical plot for the channel across all lines. The plot shows
    the min, mean, max and stdev for the channel as a function of line number.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    channel : String
        The name of the channel or field to plot.

    Returns
    -------
    None.

    '''
    filename = str(whizzFile)
    outlierRatio = 5.0

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        
        numLines = len(g.items())
        chMin = np.zeros((numLines,))
        chMax = np.zeros((numLines,))
        chMean = np.zeros((numLines,))
        chStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0
                
        for line in g.keys():
            #acqDate = g[line].attrs['Date_Local']
            lineNo[count] = line
            chMin[count] = np.min(g[line][channel])
            chMax[count] = np.max(g[line][channel])
            chMean[count] = np.mean(g[line][channel])
            chStd[count] = np.std(g[line][channel])
            maxDevRatio = (chMax[count] - chMean[count]) / chStd[count]
            minDevRatio = (chMean[count] - chMin[count]) / chStd[count]
            if maxDevRatio > outlierRatio:
                print(lineNo[count], ': max outlier ratio ', maxDevRatio)
            if minDevRatio > outlierRatio:
                print(lineNo[count], ': min outlier ratio', minDevRatio)
            count += 1

    fig = plt.figure()
    fig.suptitle(f'{projName}', fontsize=10)
    #fig.text(0.8, 0.95, f'approxacq: {acqDate:.0f}', fontsize=7)
    fig.subplots_adjust(top=0.85)
    
    ax = fig.add_subplot(1,1,1)
    ax.plot(lineNo, chMin, 'bo', mfc='w')
    ax.plot(lineNo, chMax, 'bo', mfc='w')
    ax.errorbar(lineNo, chMean, chStd, capsize=3, marker='s', c='blue')
    plt.ylabel(channel, fontsize = 6)
    plotTitle = channel + ' Stats'
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    plt.show()
    

def psdChannelDiff(whizzFile, channel1, channel2, flightLines=[]):
    '''
    Plot the PSD (log-log Sqrt(Power) from welch method) of
    channel1 - channel2 in each flightLine. 

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    flightLines : String List, optional
        A list of flightline, e.g. ['1000110.0']. Default is all lines in whizzFile.
    channel1 : String
        The name of a channel.
    channel2 : String
        The name of a channel.

    Returns
    -------
    None.

    '''
    import scipy.signal as sig
    global mean_speed
    
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if flightLines == []:
            flightLines = list(g.keys())
        corr_units = g[flightLines[0]][channel1].attrs['Units']
        if not (g[flightLines[0]][channel1].attrs['Units'] == corr_units):
            print('Error: {channel1} and {channel2} do not have the same units.')
            return

        y_label = f'{channel1} - {channel2} [{corr_units}]'
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
            
        for line in flightLines:
            mean_speed = _mean_line_speed(f[groupName], line)
            f_sample = _time_frequency(f[groupName])
            data = np.array(g[line][channel1]) - np.array(g[line][channel2])

            freq, Pxx = sig.welch(data, nfft=4*4096, fs = f_sample)
            period = 1.0 / freq[1:]
            rootPwr = np.sqrt(Pxx[1:]) / freq[1:]
            maxPwr = np.max(rootPwr)
            #print(f'{line} - low-f limit: {rootPwr[0]:.2f}, max: {maxPwr:.2f}')
            plt.loglog(period, rootPwr, color='blue', lw=0.3)
            if rootPwr[0] > 3.0:
                ax.text(period[0], rootPwr[0], f'{line}', fontsize=6)
    
        plt.ylim([1,1E5])
        plt.xlabel('Period [s]', fontsize = 6)
        plt.ylabel(y_label, fontsize = 6)
        plotTitle = f'{projName} : {y_label}'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        secax = ax.secondary_xaxis('top',
            functions=(_period_to_dist, _dist_to_period))
        secax.set_xlabel('wavelength [m]', fontsize=6)
        for label in secax.get_xticklabels(): label.set_fontsize(6)
        plt.show()


def _period_to_dist(p):
    '''
    Converts period (sec) to distance by multiplying by the GLOBAL
    mean speed (m/s).
    For use in creating a secondary wavelength axis for `psdChannelDiff`.

    Parameters
    ----------
    p : Float
        The period in seconds.

    Returns
    -------
    Float. The distance in metres.

    '''
    global mean_speed
    return mean_speed * p


def _dist_to_period(x):
    '''
    Converts distance (metres) to period(sec) by dividing by the GLOBAL
    mean speed (m/s).
    For use in creating a secondary wavelength axis for `psdChannelDiff`.

    Parameters
    ----------
    x : Float
        The distance in metres.

    Returns
    -------
    Float. The period in seconds.

    '''
    global mean_speed
    return x / mean_speed


def _time_frequency(group):
    '''
    Returns the sample frequency of the data in the project group.
    Simply calculated as the inverse of the difference of the first
    two sample times in the first line. Relies on the TimeChannel
    attribute being set. 

    Parameters
    ----------
    group : HDF5 group
        Must be the project group (top level of hierarchy in whizz File).

    Returns
    -------
    Float. The sample frequency in Hz.

    '''
    flightLines = list(group['Lines'].keys())
    time = group['CoordinateFrame'].attrs['TimeChannel']
    t = np.array(group['Lines'][flightLines[0]][time])
    return 1.0 / np.abs(t[1] - t[0])


def _mean_line_speed(group, line):
    '''
    Returns the mean line speed of the data in the line (which is
    in the given project group). Instantaneous speeds are calculated
    for each sample along the line and their average value returned.
    Relies on the XChaneel, YChannel and TimeChannel attributes being set. 

    Parameters
    ----------
    group : HDF5 group
        Must be the project group (top level of hierarchy in whizz File).
    line : String
        The line identifier.

    Returns
    -------
    Float. The mean speed along the line.

    '''
    xChannel = group['CoordinateFrame'].attrs['XChannel']
    yChannel = group['CoordinateFrame'].attrs['YChannel']
    tChannel = group['CoordinateFrame'].attrs['TimeChannel']
    xPos = np.array(group['Lines'][line][xChannel])
    yPos = np.array(group['Lines'][line][yChannel])
    sampleTime = np.array(group['Lines'][line][tChannel])
    dt = np.gradient(sampleTime)
    xVel = np.gradient(xPos) / dt
    yVel = np.gradient(yPos) / dt
    sample_speed = np.sqrt(xVel * xVel + yVel * yVel)
    return np.mean(sample_speed)


def checkErsHeaders(folderPath='\.'):
    '''
    Compares all .ers files in the folder for certain key parameters and reports
    any that are different from the mode in any parameter value.

    Parameters
    ----------
    folderPath : Path, optional
        DESCRIPTION. The default is '\.'.

    Returns
    -------
    None.

    '''
    reportString = ''
    allOk = True
    fileOK = True
    
    # get a list of the ers file paths
    file_count = 0
    ersFiles = []
    folderFiles = folderPath.iterdir()
    for aFile in folderFiles:
        if aFile.is_file() and (aFile.suffix == '.ERS' or aFile.suffix == '.ers') and (aFile.name[0] != '.'):
            ersFiles.append(aFile)
            file_count += 1
    print(f'Found {file_count} .ers files ...')
    print(f'in: {str(folderPath)}')

    # get the header dict for each ers file
    ersDicts = []
    for aFile in ersFiles:
            [ncells, nrows, nbands, eastings, northings, nullcell,
                    precision, headerbytes, originalnullcell, byteorder, datum,
                    projection] = erm.read_ers_header(str(aFile))
            commonOK, reportStr = _commonErsHdrErrors(ncells, nrows, nbands, nullcell,
                    precision, headerbytes, originalnullcell, byteorder, datum,
                    projection)
            ersDicts.append([ncells, nrows, nbands, nullcell,
                    precision, headerbytes, originalnullcell, byteorder, datum,
                    projection])
            
    # compare the contents line-by-line
    print(f'Comparing ERS files against {ersFiles[0].name}.')
    firstDict = ersDicts[0]
    print(firstDict)
    for jj in range(1, len(ersDicts)):
        fileOK = True
        print(f'Checking file {ersFiles[jj].name}')
        for ii in range(0, len(firstDict)):
            if firstDict[ii] != ersDicts[jj][ii]:
                allOk = False
                fileOK = False
                print(f'  Different element {ii}. {firstDict[ii]} != {ersDicts[jj][ii]}')
        if fileOK:
            print('  Checked OK.')
    
    return


def _commonErsHdrErrors(ncells, nrows, nbands, nullcell, precision, headerbytes, originalnullcell,
    byteorder, datum, projection):
    '''
    Checks a set of ERS header fields for various common errors.
    TODO Check that there is a Units field.

    Parameters
    ----------
    ncells : Integer
        DESCRIPTION.
    nrows : Integer
        DESCRIPTION.
    nbands : Integer
        DESCRIPTION.
    nullcell : Integer
        DESCRIPTION.
    precision : Integer
        DESCRIPTION.
    headerbytes : Integer
        DESCRIPTION.
    originalnullcell : Integer
        DESCRIPTION.
    byteorder : Integer
        DESCRIPTION.
    datum : Integer
        DESCRIPTION.
    projection : Integer
        DESCRIPTION.

    Returns
    -------
    allOk : Bool
        True if no errors found, else False.
    reportStr : String
        A report listing all errors found.

    '''
    allOk = True
    reportStr = ''
    
    if projection == 'WGS84':
        reportStr += 'ERROR: projection cannot be WGS84'
        allOk = False
    
    return allOk, reportStr


def checkIntersectionHeights(whizzFile, controls, max_allowed_height=10.0):
    """
    Checks for .
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    max_allowed_height : Float, optional
        .

    Returns
    -------
    None

    """
    data_is_good = True
    report = ''
    filename = str(whizzFile)
    num_intersections_checked = 0
    num_failed_intersections = 0
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        zChannel = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        lines = list(g.keys())
        for linec in controls:
            x_ctrl = np.array(g[linec][xChannel])
            y_ctrl = np.array(g[linec][yChannel])
            z_ctrl = np.array(g[linec][zChannel])
            bear_ctrl = _calc_bearing(x_ctrl, y_ctrl)
            (y_ctrl, x_ctrl) = _rotateCoords(x_ctrl, y_ctrl, bear_ctrl)
            for linet in lines:
                if linet == linec:
                    continue
                x_trav = np.array(g[linet][xChannel])
                y_trav = np.array(g[linet][yChannel])
                z_trav = np.array(g[linet][zChannel])
                (y_trav, x_trav) = _rotateCoords(x_trav, y_trav, bear_ctrl)
                if _lines_cross(x_ctrl, y_ctrl, x_trav, y_trav):
                    # print(f'bearings: {bearingt}, {bearingt} -- {np.abs(np.cos(bearingc - bearingt))}')
                    dh = _intersection_height(x_trav, y_trav, z_trav, x_ctrl, y_ctrl, z_ctrl, bear_ctrl)
                    num_intersections_checked += 1
                    if dh > max_allowed_height:
                        num_failed_intersections += 1
                        # print(dh)
                        report += f'\n  {linet} : {linec} intersection height difference = {dh:.1f} > {max_allowed_height:.1f}.'
                        data_is_good = False
                # else:
                #     report += f'\n  {linet} : {linec} un-tested since not perpendicular.'
    if data_is_good:
        report += f'All {num_intersections_checked} intersection heights were less than {max_allowed_height:.1f}.'
    else:
        tmpstr = f'Of {num_intersections_checked} intersections checked'
        tmpstr += f', {num_failed_intersections} exceeded the '
        tmpstr += f'{max_allowed_height} m allowed height difference.\n'
        report = tmpstr + report
    print(report)
    return


def _intersection_height(x_trav, y_trav, z_trav, x_ctrl, y_ctrl, z_ctrl, bearingc):
    '''
    '''

    y = np.abs(y_trav - _mean_1std(y_ctrl))
    it_arr = np.where(y == y.min())
    it = it_arr[0]

    x = np.abs(x_ctrl - x_trav[it])
    ic_arr = np.where(x == x.min())
    ic = ic_arr[0]

    return np.abs(z_trav[it] - z_ctrl[ic])[0]


def _mean_1std(x):
    '''
    Calculate the mean of the values in x that fall in the range of +/- stdev(x).
    This is a simplistic "mean excluding outliers".
    '''
    mean1 = np.mean(x)
    std1 = np.std(x)
    idx = np.argwhere(np.logical_and(x > (mean1 - std1), x < (mean1 + std1)))
    return np.mean(x[idx])


def _lines_cross(x_ctrl, y_ctrl, x_trav, y_trav):
    '''
    Assumes x,y coordinates rotated so that y_ctrl is approximately constant (y-axis
    approximately parallel to the control line).
    '''

    # The traverse is 'north' or 'south' of the control line
    if y_trav.min() > y_ctrl.max() or y_trav.max() < y_ctrl.min():
        return False

    # We don't allow shallow angle crossovers (and won't count lines that were supposed to be parallel!).
    min_cosine = 0.1
    bear_ctrl = _calc_bearing(x_ctrl, y_ctrl)
    bear_trav = _calc_bearing(x_trav, y_trav)
    if np.abs(np.cos(bear_ctrl - bear_trav)) < min_cosine:
        return False

    # Find the index to the traverse sample whose y coordinate is closest to the mean control y coordinate
    y = np.abs(y_trav - _mean_1std(y_ctrl))
    it_arr = np.where(y == y.min())
    if np.size(it_arr):
        itc = it_arr[0]
    else:
        return False

    # Now can check if the traverse is 'east' or 'west' of the control line.
    if x_ctrl.min() > x_trav[itc] or x_ctrl.max() < x_trav[itc]:
        return False

    return True


def _calc_bearing(x, y):
    '''
    arctan(mean(diff(x) / mean(diff(y))))
    '''
    return np.arctan(_mean_1std(np.diff(x)) / _mean_1std(np.diff(y)))
