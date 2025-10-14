#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Update the metadata for a `line group` in a `geoWhizz` file.
"""
import numpy as np
import h5py
import pathlib

import pegasusQC.config as config
import pegasusQC.utility.utility as util

groupName = config.groupName


def updateLineAttributes(whizzFile, planfiles=None, line_type='', line='', planned_line=0, flight_chan='', date_chan='', verbose=False):
    """
    For each line group, use the line_type field to set the associated planned
    line number, line variety, segment number, and reflight number for the line.

    If a planfile is given, then the line and planned_line parameters are ignored
    and the planned, segment and reflight attributes are set for all lines in
    whizzfile that also exist in planfiles.

    If no planfile is given, and both line and planned_line have non-default values,
    then the planned attribute of line is set to planned_line. No segment or reflight
    attribute is set.

    Regardless of the last two paragraphs, if flight_chan and/or date_chan are non-default,
    then the flight and/or date attributes will be set for every line in whizzfile.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz observations file, including path and extension.
    planfiles : Either String or pathlib Path of list of same, optional
        Name of HDF5 Whizz plan file(s), including path and extension. Default is None
        which results in no planned_line attributes being written to the observations
        whizz file.
    line_type : TYPE, optional
        Either 'Xcal_nsw' or 'SGL_GA' or 'SGL_NSW' or SGL_GDF' or 'NRG' or 'Xcal_can'
        or 'ARK'.
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
    line_types = ['Xcal_nsw', 'Xcal_can', 'SGL_GA', 'SGL_NSW', 'NRG', 'ARK', 'SGL_GDF', 'SGL_Kauring']

    if flight_chan != '':
        if  _channelExists(whizzFile, flight_chan):
            _updateFlightAttributes(whizzFile, flight_chan)
        else:
            print('NO ACTION TAKEN ON flight_chan - {flight_chan} not found.')

    if date_chan != '':
        if _channelExists(whizzFile, date_chan):
            _updateDateAttributes(whizzFile, date_chan)
        else:
            print('NO ACTION TAKEN ON date_chan - {date_chan} not found.')

    filename = str(whizzFile)
    if planfiles is None:
        no_plan = True
        print('NO ACTION TAKEN ON LINE_TYPE - no plan file provided.')
    elif type(planfiles) is list:
        for pfile in planfiles:
            if isinstance(pfile, pathlib.PurePath):
                pfile = str(pfile)
            elif not isinstance(planfiles[0], str):
                print("ERROR - the planfiles list needs to be a list of type String or Path but it is not.")
                return
        no_plan = False
    elif isinstance(planfiles, str):
        planfiles = [planfiles]
        no_plan = False
    elif isinstance(planfiles, pathlib.PurePath):
        planfiles = [str(planfiles)]
        no_plan = False
    else:
        print("ERROR - the input parameter needs to be a list of planfile names but it is not.")
        return

    with h5py.File(filename, 'r+') as f:
        # all the lines are in g:
        g = f[groupName]['Lines']

        if line != '' and planned_line != 0:
            gg = g[line]
            gg.attrs['PlannedLine'] = planned_line
            print(f'PlannedLine attribute of line {line} set to {planned_line}.')
            return
        if not line_type in line_types:
            print(f"NO ACTION TAKEN ON LINE_TYPE - line_type {line_type} not in {line_types}.")
            return

        print(f'\nSetting Line attributes for {whizzFile.name} according to the {line_type} scheme.')
        if not no_plan:
            print(f'Verifying planned line numbers against provided plan file(s).')

        report = ''
        if verbose:
            print('  {:<14} {:<14} {:<14} {:<14} '.format('Line No.','Plan Line No.', 'Segment No.', 'Re-flight No.'))
        for line in g:
            gg = g[line]
            current_line = gg.attrs['LineNumber']
            gg.attrs['LineType'] = line_type
            plan_line, segment_line, reflight_line, variety_line = _decode_linenumber(current_line, line_type)
            gg.attrs['Segment'] = segment_line
            gg.attrs['ReflightNumber'] = reflight_line
            gg.attrs['LineVariety'] = variety_line

            # check that plan_line is in plan_file
            if not no_plan:
                for planfile in planfiles:
                    planfile_str = str(planfile)
                    with h5py.File(planfile_str, 'r') as pf:
                        pfg = pf[groupName]['Lines']
                        lineInGroup = False
                        lineInGroup, _ = util.getLineNumberInGroup(pfg, plan_line)
                        if lineInGroup:
                            gg.attrs['PlannedLine'] = plan_line
                            break
                if not lineInGroup:
                    report += f'Line {current_line} - calculated planned line {plan_line} not found in any provided plan file.'
            if verbose:
                if no_plan or not lineInGroup:
                    print('  {:<14} {:<14} {:<14} {:<14} '.\
                        format(current_line, ' ', gg.attrs['Segment'], gg.attrs['ReflightNumber']))
                else:
                    print('  {:<14} {:<14} {:<14} {:<14} '.\
                        format(current_line, gg.attrs['PlannedLine'], gg.attrs['Segment'], gg.attrs['ReflightNumber']))
        if verbose:
            print('\n')

    return


def _decode_linenumber(current_line, line_type):
    if line_type == 'Xcal_nsw':
        if current_line < 8999999.0:
            plan_line = np.floor(current_line / 100.0) * 10.0
            segment_line = 0
            reflight_line = int(current_line - np.floor(current_line / 10.0) * 10)
        else:
            plan_line = np.floor(current_line / 10000.0)
            segment_line = 0
            reflight_line = int(current_line - np.floor(current_line / 10000.0) * 10000)
        seconddigit = str(current_line)[1]
        if seconddigit == '9':
            variety_line = "Control"
        else:
            variety_line = "Traverse"
    elif line_type == 'Xcal_can':
        if current_line < 8999999.0:
            plan_line = np.floor(current_line / 10.0) * 10.0
            segment_line = 0
            reflight_line = int(current_line - np.floor(current_line / 10.0) * 10)
        else:
            plan_line = np.floor(current_line / 10000.0)
            segment_line = 0
            reflight_line = int(current_line - np.floor(current_line / 10000.0) * 10000)
        seconddigit = str(current_line)[1]
        if seconddigit == '9':
            variety_line = "Control"
        else:
            variety_line = "Traverse"
    elif line_type == 'NRG':
        plan_line = np.floor(current_line / 10.0) * 10.0
        segment_line = 0
        reflight_line = int(current_line - np.floor(current_line / 10.0) * 10)
        if current_line > 89999:
            variety_line = "Control"
        else:
            variety_line = "Traverse"
    elif line_type == 'ARK':
        plan_line = np.floor(current_line / 10.0) * 10.0
        segment_line = 0
        reflight_line = int(current_line - np.floor(current_line / 10.0) * 10)
        if current_line > 8999:
            variety_line = "Control"
        else:
            variety_line = "Traverse"
    elif line_type == 'SGL_GA':
        if current_line < 7000:
            plan_line = np.floor(current_line * 10.0) / 10.0
            current_segment = int(np.round(10 * (current_line - np.floor(current_line))))
            segment_line = current_segment
            reflight_line = int(np.round(100 * (current_line - np.floor(current_line)))
                                             - 10 * current_segment)
        else:
            plan_line = np.floor(current_line)
            segment_line = 0
            reflight_line = int(100 * (current_line - np.floor(current_line)))
        if current_line < 999.9:
            variety_line = "Control"
        else:
            variety_line = "Traverse"
    elif line_type == 'SGL_NSW':
        if current_line > 1000000:
            plan_line = np.floor(current_line / 100.0) / 10.0
        else:
            plan_line = np.floor(current_line / 100.0) / 10.0
        if current_line < 999.9:
            variety_line = "Control"
        else:
            variety_line = "Traverse"
    elif line_type == 'SGL_GDF':
        plan_line = np.floor(current_line / 100.0) / 10.0
        segment_line = int((plan_line - np.floor(current_line / 1000.0)) * 10)
        reflight_line = int(current_line - plan_line * 1000.0)
        if current_line < 999999.9:
            variety_line = "Control"
        else:
            variety_line = "Traverse"
    elif line_type == 'SGL_Kauring':
        plan_line = np.floor(current_line / 100.0)
        remaining = (current_line - plan_line * 100.0)
        segment_line = np.floor(remaining / 10.0)
        reflight_line = int((remaining - 10.0 * segment_line))
        if current_line < 99999.9:
            variety_line = "Control"
        else:
            variety_line = "Traverse"

    return plan_line, segment_line, reflight_line, variety_line


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


def _channelExists(whizzFile, flight_chan):
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        lines = list(g.keys())
        for line in lines:
            gg = g[line]
            if flight_chan in gg.keys():
                return True
            else:
                return False
