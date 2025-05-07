#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Update the metadata for a `line group` in a `geoWhizz` file.
"""
import numpy as np
import h5py

import AirGravQC.config as config
import AirGravQC.utility.utility as util

groupName = config.groupName


def updateLineAttributes(whizzFile, planfile='', line_type='', line='', planned_line=0, flight_chan='', date_chan='', verbose=False):
    """
    For each line group, use the line_type field to set the associated planned
    line number, line variety, segment number, and reflight number for the line.

    If a planfile is given, then the line and planned_line parameters are ignored
    and the planned, segment and reflight attributes are set for all lines in
    whizzfile that also exist in planfile.

    If no planfile is given, and both line and planned_line have non-default values,
    then the planned attribute of line is set to planned_line. No segment or reflight
    attribute is set.

    Regardless of the last two paragraphs, if flight_chan and/or date_chan are non-default,
    then the flight and/or date attributes will be set for every line in whizzfile.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz observations file, including path and extension.
    planfile : String or pathlib Path, optional
        Name of a HDF5 Whizz plan file, including path and extension. Default is ''
        which results in no planned_line attributes being written to the observations
        whizz file.
    line_type : TYPE, optional
        Either 'Xcal_nsw' or 'SGL_GA' or 'SGL_NSW' or 'NRG' or 'Xcal_can'.
        The default is '' which causes the `line`s 'PlannedLine' attribute to be set
        to `planned_line`.
    line : String
        A flightline, e.g. '1000110.0'.
    planned_line : String, optional
        A flightline, e.g. '1000110.0' in a separate whizzFile containing the
        planned x,y,z locations of the survey flightlines. Defaults to '0'.
    flight_chan : String, optional
        The name of the data channel containing the flight numbers. If provided,
        the flight numbers will be stored as an attribute of each line group. The
        default is to not store them.
    date_chan : String, optional
        The name of the data channel containing the dates. If provided,
        the dates will be stored as an attribute of each line group. The
        default is to not store them.
    verbose : Bool, optional
        If False (the default) the output is reduced.

    Returns
    -------
    None.

    """

    if flight_chan != '':
        _updateFlightAttributes(whizzFile, flight_chan)
    if date_chan != '':
        _updateDateAttributes(whizzFile, date_chan)

    filename = str(whizzFile)
    if planfile == '':
        no_plan = True
    else:
        no_plan = False
        planfile_str = str(planfile)

    with h5py.File(filename, 'r+') as f:
        # all the lines are in g:
        g = f[groupName]['Lines']

        if line != '' and planned_line != 0:
            gg = g[line]
            gg.attrs['PlannedLine'] = planned_line
            print('One planned line updated.')
        elif no_plan:
            print('NO ACTION TAKEN ON LINE_TYPE - no plan file provided.')
        elif line_type in ['Xcal_nsw', 'Xcal_can', 'SGL_GA', 'SGL_NSW', 'NRG', 'SGL_GDF']:
            print(f'\nSetting Line attributes for {whizzFile.name} according to the {line_type} scheme.')
            print(f'Verifying planned line numbers against {planfile_str}.')

            with h5py.File(planfile_str, 'r') as pf:
                pfg = pf[groupName]['Lines']
                if verbose:
                    print('  {:<14} {:<14} {:<14} {:<14} '.format('Line No.','Plan Line No.', 'Segment No.', 'Re-flight No.'))
                for line in g:
                    gg = g[line]
                    current_line = gg.attrs['LineNumber']
                    gg.attrs['LineType'] = line_type
                    if line_type == 'Xcal_nsw':
                        if current_line < 8999999.0:
                            gg_planned_line = np.floor(current_line / 100.0) * 10.0
                            gg.attrs['Segment'] = 0
                            gg.attrs['ReflightNumber'] = int(current_line - np.floor(current_line / 10.0) * 10)
                        else:
                            gg_planned_line = np.floor(current_line / 10000.0)
                            gg.attrs['Segment'] = 0
                            gg.attrs['ReflightNumber'] = int(current_line - np.floor(current_line / 10000.0) * 10000)
                        seconddigit = str(current_line)[1]
                        if seconddigit == '9':
                            gg.attrs['LineVariety'] = "Control"
                        else:
                            gg.attrs['LineVariety'] = "Traverse"
                    elif line_type == 'Xcal_can':
                        if current_line < 8999999.0:
                            gg_planned_line = np.floor(current_line / 10.0) * 10.0
                            gg.attrs['Segment'] = 0
                            gg.attrs['ReflightNumber'] = int(current_line - np.floor(current_line / 10.0) * 10)
                        else:
                            gg_planned_line = np.floor(current_line / 10000.0)
                            gg.attrs['Segment'] = 0
                            gg.attrs['ReflightNumber'] = int(current_line - np.floor(current_line / 10000.0) * 10000)
                        seconddigit = str(current_line)[1]
                        if seconddigit == '9':
                            gg.attrs['LineVariety'] = "Control"
                        else:
                            gg.attrs['LineVariety'] = "Traverse"
                    elif line_type == 'NRG':
                        gg_planned_line = np.floor(current_line / 10.0) * 10.0
                        gg.attrs['Segment'] = 0
                        gg.attrs['ReflightNumber'] = int(current_line - np.floor(current_line / 10.0) * 10)
                        if current_line > 89999:
                            gg.attrs['LineVariety'] = "Control"
                        else:
                            gg.attrs['LineVariety'] = "Traverse"
                    elif line_type == 'SGL_GA':
                        if current_line < 7000:
                            gg_planned_line = np.floor(current_line * 10.0) / 10.0
                            current_segment = int(np.round(10 * (current_line - np.floor(current_line))))
                            gg.attrs['Segment'] = current_segment
                            gg.attrs['ReflightNumber'] = int(np.round(100 * (current_line - np.floor(current_line)))
                                                             - 10 * current_segment)
                        else:
                            gg_planned_line = np.floor(current_line)
                            gg.attrs['Segment'] = 0
                            gg.attrs['ReflightNumber'] = int(100 * (current_line - np.floor(current_line)))
                        if current_line < 999.9:
                            gg.attrs['LineVariety'] = "Control"
                        else:
                            gg.attrs['LineVariety'] = "Traverse"
                    elif line_type == 'SGL_NSW':
                        if current_line > 1000000:
                            gg_planned_line = np.floor(current_line / 100.0) / 10.0
                        else:
                            gg_planned_line = np.floor(current_line / 100.0) / 10.0
                        if current_line < 999.9:
                            gg.attrs['LineVariety'] = "Control"
                        else:
                            gg.attrs['LineVariety'] = "Traverse"
                    elif line_type == 'SGL_GDF':
                        if current_line < 7000:
                            gg_planned_line = np.floor(current_line / 1000.0)
                            current_segment = int(np.round(10 * (current_line - np.floor(current_line))))
                            gg.attrs['Segment'] = current_segment
                            gg.attrs['ReflightNumber'] = int(np.round(100 * (current_line - np.floor(current_line)))
                                                             - 10 * current_segment)
                        else:
                            gg_planned_line = np.floor(current_line / 1000.0)
                            gg.attrs['Segment'] = 0
                            gg.attrs['ReflightNumber'] = int(100 * (current_line - np.floor(current_line)))
                        if current_line < 999999.9:
                            gg.attrs['LineVariety'] = "Control"
                        else:
                            gg.attrs['LineVariety'] = "Traverse"

                    # check that gg_planned_line is in plan_file
                    lineInGroup, _ = util.getLineNumberInGroup(pfg, gg_planned_line)
                    if lineInGroup:
                        gg.attrs['PlannedLine'] = gg_planned_line
                    else:
                        print(f'Calculated planned line {gg_planned_line} not found in {planfile_str}.')
                    if verbose:
                        print('  {:<14} {:<14} {:<14} {:<14} '.\
                            format(current_line, gg.attrs['PlannedLine'], gg.attrs['Segment'], gg.attrs['ReflightNumber']))
                if verbose:
                    print('\n')
        else:
            print(f"NO ACTION TAKEN ON LINE_TYPE - line_type {line_type} not in ['Xcal_nsw', 'Xcal_can', 'SGL_GA', 'SGL_NSW', 'NRG'].")

    return



def _updateFlightAttributes(whizzFile, flight_chan):
    filename = str(whizzFile)
    print(f'Setting Line attributes for {filename} to include flight numbers from {flight_chan}.')

    with h5py.File(filename, 'r+') as f:
        g = f[groupName]['Lines']
        for line in g:
            gg = g[line]
            this_flight = gg[flight_chan][0]
            gg.attrs['Flight'] = this_flight


def _updateDateAttributes(whizzFile, date_chan):
    filename = str(whizzFile)
    print(f'Setting Line attributes for {filename} to include date values from {date_chan}.')

    with h5py.File(filename, 'r+') as f:
        g = f[groupName]['Lines']
        for line in g:
            gg = g[line]
            this_date = gg[date_chan][0]
            gg.attrs['Date'] = this_date
