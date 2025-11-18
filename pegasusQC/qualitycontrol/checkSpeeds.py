#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check that aircraft ground speed is within specification.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter
from matplotlib.ticker import StrMethodFormatter
import matplotlib.ticker as tkr

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.gridFiles.read_ers as ers
import pegasusQC.gridFiles.gridfiles as grd
import pegasusQC.utility.utility as util
import pegasusQC.whizzPlots.whizzPlot as wpl

groupName = config.groupName


def checkSpeeds(whizzFile, lines=[], xChannel='', yChannel='', tChannel='', vel_north='', vel_east='', nominalSpeed=60.0, 
    maxDuration=0.0, maxDistance=0.0, allowance=0.1, allowed_range=[], minSafeSpeed=42.0, known='', plot_flag=False):
    """
    Checks the data from Whizz HDF5 file for speed exceedances against a specification
    requiring ground speeds to be within a relative range (allowance) about a nominal
    value (nominalSpeed) over a particular distance (maxDistance).
    
    The positions (xChannel and yChannel) are assumed to be sampled uniformly in time
    and the first two time (tChannel) values for the first line in the file are 
    differenced to obtain the sampling interval. From this, the number of samples, N, 
    for the aircraft to travel maxDistance is calculated. The algorithm compares the
    actual distance flown in N samples with maxDistance. If the actual distance is
    greater than

    *maxDistance + allowance * maxDistance*
    
    or less than 

    *maxDistance - allowance * maxDistance*
    
    an error is printed.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path

        The pathlib Path to the Whizz HDF5 file containing the survey line data.

    lines : Array{String}, optional

        Array of line numbers as strings. Default = [], meaning all lines are checked.

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
        speeds. The default is to use `allowed_range` or, failing that, 0.1 (i.e. +/- 10% of nominal).

    allowed_range : [Float], optional

        The minimum and maximum allowed speeds as an array, `[min_allowed, max_allowed]`. The
        default is to use `allowance`.

    minSafeSpeed : Float, optional

        The minimum allowed instantaneous safe speed in m/s. The default is 42.0.

    plot_flag : Bool, optional

        If True, plot exceedances for each failed line. Default False.

    Returns
    -------
    None.

    """
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
        _reportSpeeds(g, lines=lines, maxDuration=maxDuration, maxDistance=maxDistance, xChannel=xChannel, yChannel=yChannel,
                tChannel=tChannel, vel_north=vel_north, vel_east=vel_east, nominalSpeed=nominalSpeed, 
                     allowance=allowance, allowed_range=allowed_range, minSafeSpeed=minSafeSpeed, title_str=title_str, known=known,
                     plot_flag=plot_flag)
        

def _reportSpeeds(group, lines=[], maxDuration=0.0, maxDistance=0.0, xChannel='X', yChannel='Y',
                tChannel='time', vel_north='', vel_east='', nominalSpeed=60.0,
                allowance=0.1, allowed_range=[], minSafeSpeed=42.0, title_str='', known='', plot_flag=False):
    """
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

    lines : Array{String}, optional

        Array of line numbers as strings. Default = [], meaning all lines are checked.

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

    allowed_range : [Float], optional

        The minimum and maximum allowed speeds as an array, `[min_allowed, max_allowed]`. The
        default is to use `allowance`.

    minSafeSpeed : Float, optional

        The minimum allowed instantaneous safe speed in m/s. The default is 42.0.

    title_str : String, optional

        A title string for the plots. Default ''.

    plot_flag : Bool, optional

        If True, plot exceedances for each failed line. Default False.

    Returns
    -------
    None.

    """
    calc_from_pos = False
    if vel_north == '' or vel_east == '':
        calc_from_pos = True
        print('Velocities not known - will calculate from positions')

    check_against_dist = False
    max_allowed_str = f'{maxDuration:.1f} s'
    if maxDuration < 1.0:
        check_against_dist = True
        max_allowed_str = f'{maxDistance:.1f} m'

    num_failed_lines = 0
    num_failures = 0
    num_exceed_lines = 0
    total_num_lines = 0
    report = ''

    if lines == []:
        lines = list(group.keys())

    total_num_lines = len(lines)

    exceedances_known = False
    this_exc_known = False
    number_exc_known = 0

    if len(allowed_range) == 2:
        min_allowance = allowed_range[0]
        max_allowance = allowed_range[1]
    else:
        min_allowance = nominalSpeed * (1.0 - allowance)
        max_allowance = nominalSpeed * (1.0 + allowance)

    settings = f'Nominal ground speed {nominalSpeed:.1f} m/s; '
    settings += f'allowed {min_allowance:.1f} : '
    settings += f'{max_allowance:.1f} for < {max_allowed_str}.\n'
    print(settings)
    meanspeed = 0.0
    numfids = 0.0

    for line in lines:
        lineName = util._get_lineName(group[line])
        if title_str == '':
            plot_title = lineName
        else:
            plot_title = f'{title_str}: {lineName}'

        x, y, dist, t, speed = _get_speeddata(group[line], xChannel, yChannel, tChannel, vel_north, vel_east)
        cleaned_speed = speed[~np.isnan(speed)]
        meanspeed += np.nansum(cleaned_speed)
        numfids += len(cleaned_speed)

        if speed[speed < minSafeSpeed].size > 0:
            lineUnsafeSlow = True
            print(f' For at least one reading in L{lineName}, the ground speed was < {minSafeSpeed} (might be unsafe).')
            if plot_flag:
                wpl._plot_speed(dist, 'Distance from start of line [m]', speed, minSafeSpeed, max_allowance, plot_title=plot_title)
        
        speed_extreme = 0.0
        num_fids_in_exceedance = 0
        exceedance_in_line = False
        too_slow = False
        line_fails = False

        if known != '':
            exceedances_known = True
            exc_known = np.array(group[line][known])
            report_known = -1
        
        for fid in range(0, len(x)):
            # There is a speed exceedance ...
            if speed[fid] > max_allowance or speed[fid] < min_allowance:
                # If a new exceedance, then initialise variables;
                if num_fids_in_exceedance == 0:
                    num_exceed_lines += 1
                    if speed[fid] < min_allowance:
                        too_slow = True
                    else:
                        too_slow = False
                    start_x = x[fid]
                    start_y = y[fid]
                    start_t = t[fid]
                    num_fids_in_exceedance = 1
                    speed_extreme = speed[fid]
                    this_exc_known = False
                # Else increment and update on the current exceedance.
                else:
                    # check we haven't swapped from too slow to too fast or vice versa
                    if (too_slow and speed[fid] > max_allowance) or ((not too_slow) and speed[fid] < min_allowance):
                        print(f'WARNING: Exceedance reversed speed in one fid. NOT BELIEVABLE. At time={t[fid]:.3f}')
                    num_fids_in_exceedance += 1
                    if exceedances_known:
                        if exc_known[fid] > 0:
                            report_known = exc_known[fid]
                            this_exc_known = True
                    if too_slow:
                        speed_extreme = min(speed[fid], speed_extreme)
                    else:
                        speed_extreme = max(speed[fid], speed_extreme)
            else:
                if num_fids_in_exceedance > 0: # the current exceedance has ended
                    end_x = x[fid]
                    end_y = y[fid]
                    end_t = t[fid]
                    dist_exceedance = util._length([start_x, end_x], [start_y, end_y])[1]
                    durn_exceedance = end_t - start_t
                    if too_slow:
                        speed_msg = "too slow"
                    else:
                        speed_msg = "too fast"
                    if util._exceedance_fail(durn_exceedance, dist_exceedance, maxDuration, maxDistance):
                        report += f'\nL {lineName} {speed_msg} for {durn_exceedance:.1f} sec '
                        report += f'({dist_exceedance:.0f} m), peak exceedance = {speed_extreme:.0f} m/s.'
                        report += f'\n  From ({start_x:.0f} E {start_y:.0f} N) to ({end_x:.0f} E {end_y:.0f} N).'
                        exceedance_in_line = True
                        num_failures += 1
                        if this_exc_known:
                            number_exc_known += 1
                            report += f' Known exceedance {report_known}.'
                            this_exc_known = False
                    num_fids_in_exceedance = 0
                    too_slow = False
        if exceedance_in_line:
            num_failed_lines += 1
            if plot_flag:
                if check_against_dist:
                    wpl._plot_speed(dist, 'Distance from start of line [m]', speed, min_allowance, max_allowance, plot_title=plot_title)
                else:
                    wpl._plot_speed(t-t[0], 'Time from start of line [sec]', speed, min_allowance, max_allowance, plot_title=plot_title)

    print(f' Checked {total_num_lines} lines and {num_exceed_lines} had some short exceedance(s).')
    print(f' {num_failed_lines} lines failed for exceedance > allowed.')
    print(f' Total number of full exceedances = {num_failures}.')
    print(f'\n{number_exc_known} exceedances known in the database.')
    print(f'\nMean speed over all lines = {meanspeed / numfids:.1f}.\n')
    print(report)
    if plot_flag:
        plt.show()


def _get_speeddata(line_group, xChannel, yChannel, tChannel, vel_north, vel_east):
    """
    """
    x = rd.getLineData(line_group, xChannel) #np.array(line_group[xChannel])
    y = rd.getLineData(line_group, yChannel)
    t = rd.getLineData(line_group, tChannel)
    distance = util._length(x, y)
    if vel_north == '' or vel_east == '':
        sampleTime = t[1] - t[0]
        xVel = np.diff(x) / sampleTime
        yVel = np.diff(y) / sampleTime
        temp = np.sqrt(xVel * xVel + yVel * yVel)
        speed = np.append(temp, np.mean(temp))
    else:
        xVel = rd.getLineData(line_group, vel_east)
        yVel = rd.getLineData(line_group, vel_north)
        speed = np.sqrt(xVel * xVel + yVel * yVel)

    return x, y, distance, t, speed


