#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Update the line-spacing info for the CoordinateFrame group in a `geoWhizz` file.
Author: Mark Helm Dransfield
Created: 2023
License: CC BY-SA
"""

import numpy as np
from pathlib import Path
import h5py
import collections
from scipy import stats

import pegasusQC.utility.utility as util
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.config as config

groupName = config.groupName
projectName = config.projectName


def updateLineSpacing(whizzFile, trav_spacing=None, ctrl_spacing=None, x='', y='', lines=[]):
    """
    Writes the nominal traverse and control line-spacing as (float) attributes
    for the CoordinateFrame group in `whizzFile`. If these are not known, then
    uses the position data in flight-lines with 'LineVariety' metadata 'Traverse'
    or 'Control' to estimate the line spacing for each variety.
    
    Note that any existing line spacing metadata (attribute) will be overwritten.
    A check for existing metadata before using `calcLineSpacing` is recommended.
    
    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    trav_spacing : Float, optional
        The traverse line spacing. Default None is calculate the TraverseSpacing attribute.
    ctrl_spacing : Float, optional
        The control line spacing. Default None is calculate the ControlSpacing attribute.
    x : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    y : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    lines : Array{str}, optional
        Array of line numbers as strings. Default = [], and computation is over all lines.

    Returns
    -------
    trav_spacing, ctrl_spacing : list(float)
        The calculated traverse line and control line spacings.

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r+') as f:
        # Simple, explicit case.
        if not trav_spacing is None:
            f[groupName]['CoordinateFrame'].attrs['TraverseSpacing'] = trav_spacing
            print(f'Traverse spacing set to {trav_spacing:.0f}')
        if not ctrl_spacing is None:
            f[groupName]['CoordinateFrame'].attrs['ControlSpacing'] = ctrl_spacing
            print(f'Control spacing set to {ctrl_spacing:.0f}')
        if (not ctrl_spacing is None) or (not trav_spacing is None):
            return

        # More complex, calculation case.
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']
        g = f[groupName]['Lines']
        if lines == []:
            lines = list(g.keys())
        numLines = len(lines)
        
        # These arrays may be too big; so fill with NaNs and exploit nanmean etc.
        travs = np.full((numLines, 2), np.nan)
        numtravs = 0
        ctrls = np.full((numLines, 2), np.nan)
        numctrls = 0
        travs_sorted = []
        ctrls_sorted = []

        # pick a point in lines as origin
        gg = g[lines[0]]
        x_data = rd.getLineData(gg, x)
        y_data = rd.getLineData(gg, y)
        orig_x = x_data[0]
        orig_y = y_data[0]

        for i,line in enumerate(lines):
            # print(line)
            gg = g[line]
            if 'LineVariety' in gg.attrs:
                if gg.attrs['LineVariety'] == 'Traverse':
                    travs[i,0], travs[i,1], travs_rotated = _calcLineEqn(gg, x, y, orig_x, orig_y)
                    numtravs += 1
                elif gg.attrs['LineVariety'] == 'Control':
                    ctrls[i,0], ctrls[i,1], ctrls_rotated = _calcLineEqn(gg, x, y, orig_x, orig_y)
                    numctrls += 1
                else:
                    print(f'Warning - line {line} ignored. It is not a traverse or a control.')
                    print(f'    need to run updateLineAttributes with `line_type` set.')
            else:
                print(f'Warning - line {line} ignored since it has no LineVariety attribute.')
                print(f'    need to run updateLineAttributes with `line_type` set.')
                    
        if numtravs > 0:
            travs_dc_mean, travs_dc_stdv = _calcSpacingStats(travs[0:numtravs, 1])
            print(f'Mean projected traverse spacing = {travs_dc_mean:.1f} with std dev = {travs_dc_stdv:.1f}.')
            true_spacing = _calcTrueSpacing(travs_dc_mean, travs[0:numtravs,0], travs_rotated)

            if travs_dc_stdv > 0.5 * travs_dc_mean:
                print(f'Std dev too high and estimated line spacing uncertain. Not set.')
            else:
                f[groupName]['CoordinateFrame'].attrs['TraverseSpacing'] = true_spacing
                print(f'Traverse spacing set to {true_spacing:.0f}')
        else:
            trav_spacing = 0.0
            print('No traverse lines found; TraverseSpacing not set.')
        print("")
        if numctrls > 0:
            ctrls_dc_mean, ctrls_dc_stdv = _calcSpacingStats(ctrls[0:numctrls, 1])
            print(f'Mean projected control spacing = {ctrls_dc_mean:.1f} with std dev = {ctrls_dc_stdv:.1f}.')
            true_spacing = _calcTrueSpacing(ctrls_dc_mean, ctrls[0:numctrls,0], ctrls_rotated)

            if ctrls_dc_stdv > 0.5 * ctrls_dc_mean:
                print(f'Std dev too high and estimated line spacing uncertain. Not set.')
            else:
                f[groupName]['CoordinateFrame'].attrs['ControlSpacing'] = true_spacing
                print(f'Control spacing set to {true_spacing:.0f}')
        else:
            ctrl_spacing = 0.0
            print('No control lines found; ControlSpacing not set.')
    return
                    

def _calcLineEqn(gg, x, y, orig_x, orig_y):
    """
    Uses numpy's polyfit to estimate the slope and axis intersection
    of the flight-line based on its (x,y) coordinates. Built to be
    called by `updateLineSpacing`.

    Parameters
    ----------
    gg : HDF5 group
        The line group pointing to the flight-line data.
    x : numpy 1D float array
        The x coordinates of the points on the flight-line.
    y : numpy 1D float array
        The y coordinates of the points on the flight-line.
    orig_x : float
        Constant to subtract from all x data (to avoid very large numbers).
    orig_y : float
        Constant to subtract from all y data (to avoid very large numbers).

    Returns
    -------
    m, c, flag : (float, float, bool)
        The calculated slope and intersection. The flag is true if the
        x and y data were swapped before the calculation.

    """
    x_data = rd.getLineData(gg, x)
    y_data = rd.getLineData(gg, y)
    m, c = np.polyfit(x_data-orig_x, y_data-orig_y, 1)
    if abs(m) > 1.0:
        m, c = np.polyfit(y_data, x_data, 1)
        return m, c, True
    else:
        return m, c, False


def _calcSpacingStats(data):
    """
    Calculates the mean and standard deviation of the axis intersections
    ('c' values) returned by `_calcLineEqn`. The data have NaNs removed,
    then they are sorted, and differenced. These differences are estimates
    of the flight-line spacing. Very small distances are discarded as they
    are likely due to repeat flight-lines. The first and last 20-percentile
    are also discarded as outliers, likely resulting from poor flight-
    line following. The mean and standard deviation of the remaining
    differences are returned.

    Parameters
    ----------
    data : numpy 1D float array
        The axis intersections calculated by `_calcLineEqn`.

    Returns
    -------
    data_dc_mean, data_dc_stdv : list(float)
        The calculated mean and standard deviation.

    """
    # remove the nans
    cleaned_data = data[~np.isnan(data)]
    # print(cleaned_data)

    # sort the data
    data_sorted = np.sort(cleaned_data)
    # print(data_sorted)

    # diff the values of dim1, and take the mean and stdev
    data_dc = np.diff(data_sorted)

    # Lines closer than 25 m are likely repeats
    # It would be better to replace 25.0 with a value calculated from the data
    mask = data_dc > 25.0
    data_masked = data_dc[mask]

    # ignore outliers
    q1 = np.percentile(data_masked, 20)
    q3 = np.percentile(data_masked, 80)
    iqr = q3 - q1
    lower_fence = q1 - 1.5 * iqr
    upper_fence = q3 + 1.5 * iqr
    keep = np.where((data_masked > lower_fence) & (data_masked < upper_fence))
    data_dc_mean = np.nanmean(data_masked[keep])
    data_dc_stdv = np.nanstd(data_masked[keep])

    return data_dc_mean, data_dc_stdv


def _calcTrueSpacing(spacing, data, rotated):
    """
    Given the best estimate of the projected flight-line spacing and
    the list of the best-fit slopes of the flight-lines, calculate
    the best estimate of the true flight-line spacing.
    
    Parameters
    ----------
    spacing : float
        The estimated, projected flight-line spacing.
    data : numpy 1D float array
        The line slopes for each flight-line calculated by `_calcLineEqn`.
    rotated : bool
        Flag to indicate if  the x and y data were swapped before the
        slope was calculated (returned by `_calcLineEqn`).

    Returns
    -------
    true_spacing : float
        The calculated true flight-line spacing.

    """    # remove the nans
    m = data[~np.isnan(data)]
    # Get the best estimate of the nominal slope
    q1 = np.percentile(m, 20)
    q3 = np.percentile(m, 80)
    iqr = q3 - q1
    lower_fence = q1 - 1.5 * iqr
    upper_fence = q3 + 1.5 * iqr
    keep = np.where((m > lower_fence) & (m < upper_fence))
    m_mean = np.nanmean(m[keep])
    m_stdv = np.nanstd(m[keep])
    print(f'Best estimate of slope = {m_mean:.2g}; std dev = {m_stdv:.2g}.')
    if rotated:
        true_spacing = spacing * np.cos(np.arctan2(m_mean, 1.0))
    else:
        true_spacing = spacing * np.cos(np.arctan2(m_mean, 1.0))
    print(f'Calculated spacing = {true_spacing:.1f}.')
    true_spacing = round(true_spacing, -1)
    print(f'Final nominal spacing = {true_spacing:.1f}.')
    return true_spacing


