#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check absolute differences at traverse - control line intersections.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

# import necessary modules
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr

import pegasusQC.config as config
import pegasusQC.utility.utility as util
import pegasusQC.whizzPlots.whizzPlot as wpl

groupName = config.groupName
        

def checkIntersection(whizzFile, controls=[], travs=[], xChannel='', yChannel='', zChannel='', max_allowed_deltaZ=10.0, mode='RMS', verbose=False, plot_flag=False):
    """
    Checks that the values of the data in `zChannel` at the intersection of traverse and control
    lines are different by no more than the maximum allowed delta. If `mode` = "value", then each
    intersection is assessed individually; if `mode` = "RMS", then the check is against the RMS
    difference along each traverse line, if `mode` = "RMSroot2", then the check along each traverse
    line is against the RMS difference divided by the square root of 2.
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path

        The pathlib Path to the Whizz HDF5 file containing the survey line data.

    controls : [String], optional

        A list of control flightlines, e.g. ['1000110.0', '1000210.0', '1000310.0']. Defaults
        to all flight lines with the `LineVariety` attribute set to "Control".

    travs : [String], optional

        A list of traverse flightlines, e.g. ['1000110.0', '1000210.0', '1000310.0']. Defaults
        to all flight lines with the `LineVariety` attribute set to "Traverse".

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

        The maximum allowed difference in `zChannel` between the traverse and control lines
        at each intersection point. Defaults to 10.0.

    mode : String, optional

        Must be one of "value", "RMS", or "RMSroot2".
        Default "RMS".

    verbose : Bool, optional

        If True, verbose reporting is given which is annoying if there are many errors.
        Default False.

    plot_flag : Bool, optional

        If True, plot a map of each pair of intersecting lines where the `zChannel`
        values differ by more than `max_allowed_deltaZ`. Default False.

    Returns
    -------
    None

    """
    if mode != 'RMS' and mode != 'RMSroot2' and mode != 'value':
        print(f'ERROR: mode is set to {mode} but must be one of "RMS", "RMSroot2", or "value".')
        return
    data_is_good = True
    report = ''
    filename = str(whizzFile)
    num_checked = 0
    num_failed = 0
    # num_intersections_checked = 0
    # num_failed_intersections = 0
    # num_travs_checked = 0
    # num_failed_travs = 0

    if plot_flag:
        plot_subtitle = 'Intersection analysis'

        fig = plt.figure(figsize=(8, 6))
        thou_format = tkr.FuncFormatter(util._space_thou)
        ax = fig.add_subplot(1,1,1)
        plotTitle = 'Flight Line Map: '
    
    with h5py.File(filename, 'r') as f:
        if plot_flag:
            plotTitle += wpl.make_plot_title(f[groupName])

        g = f[groupName]['Lines']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if zChannel == '':
            zChannel = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        if plot_flag:
            plot_subtitle += f': {zChannel}'
        all_lines = list(g.keys())
        allcontrols, alltravs, anyunknowns = lines_by_variety(g, all_lines)
        if len(anyunknowns) > 0:
            print(f'{len(anyunknowns)} lines in database not attributed as traverse or as control.')
            print('The line type (control or traverse) can be set by `updateLineAttributes`.')
        else:
            print(f'All {len(all_lines)} lines in database indexed as traverse or as control.')
        if controls == []:
            controls = allcontrols
        if travs == []:
            lines = alltravs
        else:
            lines = travs
        #if the channel has an attribute 'Units'
        dd = g[lines[0]][zChannel]
        z_units = ''
        if 'Units' in dd.attrs.keys():
            z_units = dd.attrs['Units']

        for line_trav in lines:
            x_trv = np.array(g[line_trav][xChannel])
            y_trv = np.array(g[line_trav][yChannel])
            z_trv = np.array(g[line_trav][zChannel])
            trav_name = util._get_lineName(g[line_trav])
            if plot_flag:
                travline, = ax.plot(x_trv, y_trv, color='blue', lw=0.6, alpha=0.7)

            bear_trv = util._calc_bearing(x_trv, y_trv)
            (y_trv1, x_trv1) = util._rotateCoords(x_trv-x_trv[0], y_trv-y_trv[0], -bear_trv)
            deltaheights = np.array([])
            for line_ctrl in controls:
                # if line_ctrl == line_trav, then it is a control line, not a traverse.
                # TODO: compare the PlannedLine for line_ctrl and line_trav rather than the lines themselves,
                #       since we don't want to compare different segments of the same line!
                if line_ctrl == line_trav:
                    continue
                x_ctl = np.array(g[line_ctrl][xChannel])
                y_ctl = np.array(g[line_ctrl][yChannel])
                z_ctl = np.array(g[line_ctrl][zChannel])
                ctrl_name = util._get_lineName(g[line_ctrl])
                if plot_flag:
                    ctrlline, = ax.plot(x_ctl, y_ctl, color='blue', lw=0.6, alpha=0.7)

                (y_ctl1, x_ctl1) = util._rotateCoords(x_ctl-x_trv[0], y_ctl-y_trv[0], -bear_trv)

                if _intersect(x_trv[0], y_trv[0], x_trv[-1], y_trv[-1], x_ctl[0], y_ctl[0], x_ctl[-1], y_ctl[-1]):
                    # print(f'bearings: {bearingt}, {bearingt} -- {np.abs(np.cos(bearingc - bearingt))}')
                    dh, ipoint = _intersection_height(x_ctl1, y_ctl1, z_ctl, x_trv1, y_trv1, z_trv)
                    deltaheights = np.append(deltaheights, dh)
                    if mode == "value":
                        abs_dh = np.abs(dh)
                        num_checked += 1
                        if abs_dh > max_allowed_deltaZ:
                            num_failed += 1
                            # print(dh)
                            bc_deg = bear_trv * 180.0 / np.pi
                            report += f'\n  {trav_name} : {ctrl_name} [track = {bc_deg:.1f} deg N] intersection {zChannel} difference = {abs_dh:.1f} > {max_allowed_deltaZ:.1f}'
                            if z_units != '':
                                report += ' ' + z_units
                            report += '.'
                            if plot_flag:
                                travpoint, = ax.plot(x_trv[ipoint], y_trv[ipoint], color='red', marker="o", ms=3, alpha=0.9)
                            data_is_good = False

            if mode == "RMS" or mode == "RMSroot2":
                mode_str = "RMS"
                dh_rms = np.nanstd(deltaheights)
                num_checked += 1
                if mode == "RMSroot2":
                    dh_rms = dh_rms / np.sqrt(2.0)
                    mode_str = "RMS/sqrt2"
                if dh_rms > max_allowed_deltaZ:
                    num_failed += 1
                    bc_deg = bear_trv * 180.0 / np.pi
                    report += f'\n  {trav_name} [bearing={bc_deg:.1f}] intersection {zChannel} {mode_str} difference = {dh_rms:.1f} > {max_allowed_deltaZ:.1f}'
                    if z_units != '':
                        report += ' ' + z_units
                    report += '.'
                    if plot_flag:
                        travline, = ax.plot(x_trv, y_trv, color='red', lw=0.8, alpha=0.9)
                    data_is_good = False

    if data_is_good:
        if mode == "RMS" or mode == "RMSroot2":
            report += f'All {num_checked} traverse lines had {zChannel} {mode_str} differences less than {max_allowed_deltaZ:.1f}'
        else:
            report += f'All {num_checked} intersection {zChannel} differences were less than {max_allowed_deltaZ:.1f}'
        if z_units != '':
            report += ' ' + z_units
        report += '.'
    else:
        if mode == "RMS" or mode == "RMSroot2":
            tmpstr = f'{num_checked} traverse lines checked,'
            tmpstr += f', {num_failed} exceeded the {max_allowed_deltaZ}'
            if z_units != '':
                tmpstr += ' ' + z_units
            tmpstr += f' allowed {mode_str} difference in {zChannel}.\n'
        else:
            tmpstr = f'{num_checked} intersections checked,'
            tmpstr += f', {num_failed} exceeded the {max_allowed_deltaZ}'
            if z_units != '':
                tmpstr += ' ' + z_units
            tmpstr += f' allowed absolute difference in {zChannel}.\n'
        report = tmpstr + report
    print(report)
    if plot_flag:
        ax.set_aspect('equal')
        ax.xaxis.set_major_formatter(thou_format)
        ax.yaxis.set_major_formatter(thou_format)
        plt.xlabel(f'{xChannel} [m]', fontsize = 10)
        plt.ylabel(f'{yChannel} [m]', fontsize = 10)
        plt.suptitle(plotTitle, fontsize = 12)
        plt.title(plot_subtitle, fontsize = 10)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(8)
        for label in ax.get_yticklabels(): label.set_fontsize(8)
        plt.show()

    return


def lines_by_variety(gLines, all_lines):
    """
    Returns lists of control lines, traverse lines and lines not attributed as either
    from the lines `all_lines` in group `gLines`.
    
    Parameters
    ----------
    gLines : HDF5 Group

        The geoWhizz lines group containing the survey line data.

    all_lines : [String]

        A list of all flightlines to search through.

    Returns
    -------
    controls : [String]

        A list of control flightlines, e.g. ['1000110.0', '1000210.0', '1000310.0'] with the
        `LineVariety` attribute set to "Control".

    traverses : [String]

        A list of traverse flightlines with the `LineVariety` attribute set to "Traverse".

    unknowns : [String]

        A list of flightlines with the `LineVariety` attribute not set to either "Control" or
        "Traverse".

    """
    controls = []
    traverses = []
    unknowns = []
    for line in all_lines:
        if gLines[line].attrs['LineVariety'] == "Control":
            controls.append(line)
        elif gLines[line].attrs['LineVariety'] == "Traverse":
            traverses.append(line)
        else:
            unknowns.append(line)
    return controls, traverses, unknowns


def _ccw(x1, y1, x2, y2, x3, y3):
    return (y3-y1)*(x2-x1) > (y2-y1)*(x3-x1)


def _intersect(cx1, cy1, cx2, cy2, tx1, ty1, tx2, ty2):
    return _ccw(cx1, cy1, tx1, ty1, tx2, ty2) != _ccw(cx2, cy2, tx1, ty1, tx2, ty2) and _ccw(cx1, cy1, cx2, cy2, tx1, ty1) != _ccw(cx1, cy1,cx2, cy2, tx2, ty2)


def _intersection_height(x_ctrl, y_ctrl, z_ctrl, x_trav, y_trav, z_trav):
    """
    Returns the difference in `z` values at the intersection of the traverse and control lines.

    We are given the horizontal positions (x_trav, y_trav) and (x_ctrl, y_ctrl) for sample locations
    along the traverse and control line respectively, and in coordinates where the nominal control
    line is parallel to the `x` axis, and the nominal traverse line is parallel to the `y` axis.
    We find the index, `it`, to the traverse position where `y_trav` is closest to the control line;
    this must be the closest traverse point to the intersection of the nominal lines. Then we find the
    index, ic, to the control position where `x_ctrl` is closest to the traverse line; similarly, this
    is the closest control point to the intersection. The difference in `z` is then 
         z_trav[it] - z_ctrl[ic]
    There is no interpolation done. Furthermore, these are not necessarily the nearest neighbours to
    each other, but rather the points nearest to the intersection of the nominal control and traverse
    lines.
    
    Parameters
    ----------
    x_trav : Numpy 1D array

        The `x` values of the traverse line.

    y_trav : Numpy 1D array

        The `y` values of the traverse line.

    z_trav : Numpy 1D array

        The `z` values of the traverse line.

    x_ctrl : Numpy 1D array

        The `x` values of the control line.

    y_ctrl : Numpy 1D array

        The `y` values of the control line.

    z_ctrl : Numpy 1D array

        The `z` values of the control line.

    Returns
    -------
    zdif : float

        the difference in `z` values at the intersection of the traverse and control lines.

    idx : float

        the index in the traverse arrays of the point closest to the intersection.

    """

    # print(np.std(y_ctrl), np.std(x_ctrl))
    # print(np.std(y_trav), np.std(x_trav))

    # The data are rotated onto coordinates so that the traverse line is nominally at y=0.
    # So the control line intersects at y=0; find the index to the point closest to y_ctrl = 0
    y = np.abs(y_ctrl)
    ic = np.where(y == y.min())[0]
    # then the index to the traverse point whose x value is closest to the control lines
    # x value at that intersection.
    x = np.abs(x_trav - x_ctrl[ic])
    it = np.where(x == x.min())[0]

    return (z_trav[it] - z_ctrl[ic])[0], it


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
    y = np.abs(y_trav - util._mean_1std(y_ctrl))
    it_arr = np.where(y == y.min())
    if np.size(it_arr):
        itc = it_arr[0]
    else:
        return False

    # Now can check if the traverse is 'east' or 'west' of the control line.
    if x_ctrl.min() > x_trav[itc] or x_ctrl.max() < x_trav[itc]:
        return False

    return True


def controls_lessthan_1000(all_lines):
    ctrl_strs = []
    trav_strs = []
    for line in all_lines:
        if float(line) < 999.99:
            ctrl_strs.append(line)
        else:
            trav_strs.append(line)
    return trav_strs, ctrl_strs


def controls_nineteen(all_lines):
    ctrl_strs = []
    trav_strs = []
    for line in all_lines:
        if line[:2] == "19":
            ctrl_strs.append(line)
        else:
            trav_strs.append(line)
    return trav_strs, ctrl_strs


