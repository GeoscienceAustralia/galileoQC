#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check GNSS statistics are within specification.
"""
import numpy as np
import h5py

import AirGravQC.config as config
import AirGravQC.whizzFiles.retrieveData as rd

groupName = config.groupName


def checkGNSS(whizzFile, num_sats, pdop, vdop, hdop, nsats_min=4, max_pdop=6, max_vdop=4, max_hdop=4, lines=[]):
    """
    Checks that the data in a whizzFile meets the requirements for the minimum
    number of satellites, and maximum PDOP, VDOP and HDOP.

    Parameters
    ----------
    whizzFile : String or pathlib Path
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
    lines : Array{String}, optional
        Array of line numbers as strings. Default = [], meaning all lines are checked.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    report = ''
    num_failed_lines = 0

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        xchan = f[groupName]['CoordinateFrame'].attrs['XChannel']
        ychan = f[groupName]['CoordinateFrame'].attrs['YChannel']

        error_count = 0

        if lines == []:
            lines = list(g.keys())

        for line in lines:
            x = rd.getLineData(g[line], xchan)
            y = rd.getLineData(g[line], ychan)

            nsats_data = rd.getLineData(g[line], num_sats)
            min_nsats_data = np.nanmin(nsats_data)
            if min_nsats_data < nsats_min:
                xmin_fail = np.nanmin(x[nsats_data < nsats_min])
                ymin_fail = np.nanmin(y[nsats_data < nsats_min])
                xmax_fail = np.nanmax(x[nsats_data < nsats_min])
                ymax_fail = np.nanmax(y[nsats_data < nsats_min])
                numfail = np.count_nonzero(nsats_data < nsats_min)
                report += f'Line {line} failed for {numfail} fids: min sats = {min_nsats_data:.0f} < {nsats_min}.'
                report += f' Easting [{xmin_fail:.0f}, {xmax_fail:.0f}],'
                report += f' Northing [{ymin_fail:.0f}, {ymax_fail:.0f}].\n'
                error_count += 1

            pdop_data = rd.getLineData(g[line], pdop)
            max_pdop_data = np.nanmax(pdop_data)
            if max_pdop_data > max_pdop:
                xmin_fail = np.nanmin(x[pdop_data > max_pdop])
                ymin_fail = np.nanmin(y[pdop_data > max_pdop])
                xmax_fail = np.nanmax(x[pdop_data > max_pdop])
                ymax_fail = np.nanmax(y[pdop_data > max_pdop])
                numfail = np.count_nonzero(pdop_data > max_pdop)
                report += f'Line {line} failed for {numfail} fids: max pdop = {max_pdop_data:.1f} > {max_pdop}.'
                report += f' Easting [{xmin_fail:.0f}, {xmax_fail:.0f}],'
                report += f' Northing [{ymin_fail:.0f}, {ymax_fail:.0f}].\n'
                error_count += 1

            vdop_data = rd.getLineData(g[line], vdop)
            max_vdop_data = np.nanmax(vdop_data)
            if max_vdop_data > max_vdop:
                xmin_fail = np.nanmin(x[vdop_data > max_vdop])
                ymin_fail = np.nanmin(y[vdop_data > max_vdop])
                xmax_fail = np.nanmax(x[vdop_data > max_vdop])
                ymax_fail = np.nanmax(y[vdop_data > max_vdop])
                numfail = np.count_nonzero(vdop_data > max_vdop)
                report += f'Line {line} failed for {numfail} fids: max vdop = {max_vdop_data:.1f} > {max_vdop}.'
                report += f' Easting [{xmin_fail:.0f}, {xmax_fail:.0f}],'
                report += f' Northing [{ymin_fail:.0f}, {ymax_fail:.0f}].\n'
                error_count += 1

            hdop_data = rd.getLineData(g[line], hdop)
            max_hdop_data = np.nanmax(hdop_data)
            if max_hdop_data > max_hdop:
                xmin_fail = np.nanmin(x[hdop_data > max_hdop])
                ymin_fail = np.nanmin(y[hdop_data > max_hdop])
                xmax_fail = np.nanmax(x[hdop_data > max_hdop])
                ymax_fail = np.nanmax(y[hdop_data > max_hdop])
                numfail = np.count_nonzero(hdop_data > max_hdop)
                report += f'Line {line} failed for {numfail} fids: max hdop = {max_hdop_data:.1f} > {max_hdop}.'
                report += f' Easting [{xmin_fail:.0f}, {xmax_fail:.0f}],'
                report += f' Northing [{ymin_fail:.0f}, {ymax_fail:.0f}].\n'
                error_count += 1

        if error_count == 1:
            errstr = f'Found 1 error.'
        elif error_count > 1:
            errstr = f'Found {error_count} errors.'
        else:
            errstr = f'Found 0 errors.'
        report = f'In {projName}, checked num sats, PDOP, VDOP and HDOP. {errstr}\n' + report
        print(report)

