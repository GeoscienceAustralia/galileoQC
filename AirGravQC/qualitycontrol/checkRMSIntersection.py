#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check RMS differences at traverse - control line intersections.
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

import AirGravQC.config as config
import AirGravQC.utility.utility as util
from AirGravQC.qualitycontrol.checkIntersection import lines_by_variety

groupName = config.groupName
        

def checkRMSIntersection(whizzFile, controls=[], travs=[], xChannel='', yChannel='', zChannel='', max_allowed_deltaZ=10.0, divroot2=False, verbose=False):
    """
    Checks that, for each control line, the RMS difference between the values of the data in `zChannel` at the intersection
    with each traverse line is less than the maximum allowed delta.
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    controls : [String]
        A list of control flightlines, e.g. ['1000110.0', '1000210.0', '1000310.0']. Defaults
        to all flight lines with value < 1000.
    travs : [String], optional
        A list of traverse flightlines, e.g. ['1000110.0', '1000210.0', '1000310.0']. Defaults
        to all flight lines in the `whizzFile` with value > 999.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    zChannel : String, optional
        The name of the geoWhizz field or channel containing the data to be tested. The
        default is to read the AltitudeChannel field name from the Coordinate Frame.
    max_allowed_deltaZ : Float, optional
        The maximum allowed RMS difference in `zChannel` between the traverse and control lines
        at each intersection point along each control line. Defaults to 10.0.
    divroot2 : Bool, optional
        If True, the check is against the RMS difference divided by the square root of 2.
        Default False.
    verbose : Bool, optional
        If True, verbose reporting is given which is annoying if there are many errors.
        Default False.

    Returns
    -------
    None

    """
    print('checkRMSIntersection() has been replaced with checkIntersection(). Please use the latter function with `mode`="RMS" or "RMSroot2".')
    return
    data_is_good = True
    report = ''
    filename = str(whizzFile)
    num_travs_checked = 0
    num_failed_travs = 0
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if zChannel == '':
            zChannel = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        all_lines = list(g.keys())
        allcontrols, alltravs, anyunknowns = lines_by_variety(g, all_lines)
        print(f'{len(anyunknowns)} lines in database not attributed as traverse or as control.')
        if len(anyunknowns) > 0:
            print('The line type (control or traverse) can be set by `updateLineAttributes`.')
        if controls == []:
            # if any lines have line variety set:
            #     controls = all_lines_with_variety(g, 'controls')
            # else
            controls = allcontrols
        if travs == []:
            # if any lines have line variety set:
            #     travs = all_lines_with_variety(g, 'trav')
            # else
            lines = alltravs
        else:
            lines = travs
        #if the channel has an attribute 'Units'
        dd = g[lines[0]][zChannel]
        z_units = ''
        if 'Units' in dd.attrs.keys():
            z_units = dd.attrs['Units']

        for linetrav in lines:
            num_travs_checked += 1
            x_trv = np.array(g[linetrav][xChannel])
            y_trv = np.array(g[linetrav][yChannel])
            z_trv = np.array(g[linetrav][zChannel])

            bear_trv = util._calc_bearing(x_trv, y_trv)
            (y_trv1, x_trv1) = util._rotateCoords(x_trv-x_trv[0], y_trv-y_trv[0], -bear_trv)
            dh = np.array([])
            for linectrl in controls:
                # if linectrl == linetrav, then it is a control line, not a traverse.
                # TODO: compare the PlannedLine for linectrl and linetrav rather than the lines themselves,
                #       since we don't want to compare different segments of the same line!
                if linectrl == linetrav:
                    continue
                x_ctr = np.array(g[linectrl][xChannel])
                y_ctr = np.array(g[linectrl][yChannel])
                z_ctr = np.array(g[linectrl][zChannel])
                (y_ctr1, x_ctr1) = util._rotateCoords(x_ctr-x_trv[0], y_ctr-y_trv[0], -bear_trv)
                if _intersect(x_trv[0], y_trv[0], x_trv[-1], y_trv[-1], x_ctr[0], y_ctr[0], x_ctr[-1], y_ctr[-1]):
                    dh = np.append(dh, _intersection_height(x_ctr1, y_ctr1, z_ctr, x_trv1, y_trv1, z_trv, bear_trv))

            if len(lines) == 1 and verbose == True:
                print(dh)

            dh_rms = np.nanstd(dh)
            if divroot2:
                dh_rms = dh_rms / np.sqrt(2.0)
            if dh_rms > max_allowed_deltaZ:
                num_failed_travs += 1
                bc_deg = bear_trv * 180.0 / np.pi
                report += f'\n  {linetrav} [bearing={bc_deg:.1f}] intersection {zChannel} RMS difference = {dh_rms:.1f} > {max_allowed_deltaZ:.1f}'
                if z_units != '':
                    report += ' ' + z_units
                report += '.'
                data_is_good = False
    if data_is_good:
        report += f'All {num_travs_checked} traverse lines had {zChannel} RMS differences less than {max_allowed_deltaZ:.1f}'
        if z_units != '':
            report += ' ' + z_units
        report += '.'
    else:
        tmpstr = f'Of {num_travs_checked} traverse lines checked'
        tmpstr += f', {num_failed_travs} traverse lines exceeded the '
        tmpstr += f'{max_allowed_deltaZ} allowed {zChannel} RMS difference.\n'
        report = tmpstr + report
    print(report)
    return


def _ccw(x1, y1, x2, y2, x3, y3):
    return (y3-y1)*(x2-x1) > (y2-y1)*(x3-x1)


def _intersect(cx1, cy1, cx2, cy2, tx1, ty1, tx2, ty2):
    return _ccw(cx1, cy1, tx1, ty1, tx2, ty2) != _ccw(cx2, cy2, tx1, ty1, tx2, ty2) and _ccw(cx1, cy1, cx2, cy2, tx1, ty1) != _ccw(cx1, cy1,cx2, cy2, tx2, ty2)


def _intersection_height(x_trav, y_trav, z_trav, x_ctrl, y_ctrl, z_ctrl, bearingc):
    """
    """

    y = np.abs(y_trav - num_failed_travs(y_ctrl))
    it_arr = np.where(y == y.min())
    it = it_arr[0]
    if len(it) > 1:
        print('\nBug in _intersection_height - this intersection height cannot be calculated')
        print(it, type(it), len(it))
        return 0.0

    x = np.abs(x_ctrl - x_trav[it])
    ic_arr = np.where(x == x.min())
    ic = ic_arr[0]

    # print(it, z_trav[it], ic, z_ctrl[ic])
    return (z_trav[it] - z_ctrl[ic])[0]


def _lines_cross(x_ctrl, y_ctrl, x_trav, y_trav):
    """
    Assumes x,y coordinates rotated so that y_ctrl is approximately constant (y-axis
    approximately parallel to the control line).
    """

    # The traverse is 'north' or 'south' of the control line
    if y_trav.min() > y_ctrl.max() or y_trav.max() < y_ctrl.min():
        return False

    # We don't allow shallow angle crossovers (and won't count lines that were supposed to be parallel!).
    min_cosine = 0.1
    bear_ctrl = util._calc_bearing(x_ctrl, y_ctrl)
    bear_trav = util._calc_bearing(x_trav, y_trav)
    if np.abs(np.cos(bear_ctrl - bear_trav)) < min_cosine:
        return False

    # Find the index to the traverse sample whose y coordinate is closest to the mean control y coordinate
    y = np.abs(y_trav - num_failed_travs(y_ctrl))
    it_arr = np.where(y == y.min())
    if np.size(it_arr):
        itc = it_arr[0]
    else:
        return False

    # Now can check if the traverse is 'east' or 'west' of the control line.
    if x_ctrl.min() > x_trav[itc] or x_ctrl.max() < x_trav[itc]:
        return False

    return True



