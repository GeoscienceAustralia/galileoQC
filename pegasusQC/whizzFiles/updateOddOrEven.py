#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set line group attribute Parity for each line.
Author: Mark Helm Dransfield
Created: 2023
License: CC BY-SA
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from pathlib import Path
import xarray as xr
import netCDF4 as nc4
import filebrowser as fb
# import rioxarray
import h5py
import matplotlib.ticker as tkr
import collections

import pegasusQC.gridFiles.graphics as graphics
import pegasusQC.utility.utility as util
import pegasusQC.gridFiles.read_ers as ers
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.config as config
from pegasusQC.gridFiles.graphicsShaded import graphicsShaded
from pegasusQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from pegasusQC.gridFiles.xarray_to_grid import xarray_to_grid
from pegasusQC.gridFiles.xdImage import xdImage
import pegasusQC.gridFiles.gridutility as gut


groupName = config.groupName
projectName = config.projectName


def updateOddOrEven(whizzFile='', x='', y='', oddlines=[], evenlines=[], planFile='', verbose=False):
    """
    Saves the Parity (odd or even) attribute for each traverse line so that odd-even analysis can
    be easily performed.

    If either oddlines or evenlines contain a list of line numbers as text, then
    they will be flagged as odd or even respectively in `whizzFile`. If they are empty or left as default,
    and `whizzFile` is an empty string (the default) then the traverse lines in `planFile` will be
    analysed and set to either odd or even appropriately. If `whizzFile` and `planFile` are both provided,
    and are valid whizzFiles, then the traverse lines in both files will be set to odd or even according
    to the algorithm. If there are no oddlines or evenlines, and `planFile` is empty (the default),
    then an error is reported and no other action is taken.

    Parameters
    ----------
    whizzFile : Path or String, optional
        The Path to, or String name of, the whizz file in HDF5 format. Default is
        to not provide a planFile.
    x : String, optional
        The name of the x (easting) channel in `whizz_file`. Default is to use
        the `XChannel` attribute.
    y : String, optional
        The name of the y (northing) channel in `whizz_file`. Default is to use
        the `YChannel` attribute.
    oddlines : String Array, optional
        List of lines to have their parity set to True. Default is to calculate the
        parity based on analysis of the `planFile`. Only lines with the Variety attribute
        set to `Traverse` will have their parity set.
    evenlines : String Array, optional
        List of lines to have their parity set to False. Default is to calculate the
        parity based on analysis of the `planFile`. Only lines with the Variety attribute
        set to `Traverse` will have their parity set.
    planFile : Path or String, optional
        The Path to, or String name of, the plan file in HDF5 format. This is
        required if parities are to be set by calculation. Default is to not use
        calculation. Default is to not provide a planFile.
    verbose : Bool, optional
        The verbosity of output. Default False.

    Returns
    -------
    Nothing.

    """
    if len(oddlines) > 0 or len(evenlines) > 0:
        filename = str(whizzFile)
        with h5py.File(filename, 'r+') as f:
            lines_group = f[groupName]['Lines']
            for line in oddlines:
                linegroup = lines_group[line]
                if not (linegroup.attrs['LineVariety'] == 'Traverse'):
                    continue
                else:
                    linegroup.attrs['Parity'] = True
                    if verbose:
                        print(f'Line {line}: Parity = {linegroup.attrs["Parity"]}.')
            for line in evenlines:
                linegroup = lines_group[line]
                if not (linegroup.attrs['LineVariety'] == 'Traverse'):
                    continue
                else:
                    linegroup.attrs['Parity'] = False
                    if verbose:
                        print(f'Line {line}: Parity = {linegroup.attrs["Parity"]}.')
        return

    elif (not planFile == '') and whizzFile == '':
        print('Setting parities in plan file.')
        _set_plan_parities(planFile, x, y, verbose)
    elif not planFile == '' and not whizzFile == '':
        print('Setting parities in plan file.')
        _set_plan_parities(planFile, x, y, verbose)
        print('Setting parities in observations file.')
        _set_flown_parities(whizzFile, planFile, verbose)
    else:
        print('ERROR - no values provided for either planFile, or oddlines / evenlines; stopping.')
        return


def _set_plan_parities(planFile, x='', y='', verbose=False):
    filename = str(planFile)
    with h5py.File(filename, 'r+') as f:
        lines_group = f[groupName]['Lines']
        lines = list(lines_group.keys())
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if 'TravSpacing' in f[groupName]['CoordinateFrame'].attrs:
            trav_spacing = f[groupName]['CoordinateFrame'].attrs['TravSpacing']
        else:
            print('ERROR - TravSpacing not found in CoordinateFrame attributes.')
            print('    Try running updateLineSpacing() to set TravSpacing.')
            return

        trav_parity = True
        first = True

        for line in lines:
            linegroup = lines_group[line]
            if not (linegroup.attrs['LineVariety'] == 'Traverse'):
                continue

            if first:
                if _linelength(linegroup, x, y) < 100.0:
                    print(f'WARNING: line {line} is < 100 m long, skipped.')
                    continue
                first = False
                linegroup.attrs['Parity'] = trav_parity
                firstgroup = linegroup
                continue

            # if the number of traverse spacings is even
            nspacings = _approx_spacings(linegroup, firstgroup, trav_spacing, x, y)
            if nspacings == 0:
                if verbose:
                    print(f'Skipped line {line}: spacing to first line inconsistent with {trav_spacing}.')
                continue
            # if the number of traverse spacings is even, the parity is the same
            elif nspacings % 2 == 0:
                linegroup.attrs['Parity'] = trav_parity
            else:
                linegroup.attrs['Parity'] = not trav_parity
            if verbose:
                print(f'Line {line}: Parity = {linegroup.attrs["Parity"]}.')


def _set_flown_parities(whizzFile, planFile, verbose=False):
    whizzname = str(whizzFile)
    planname = str(planFile)

    with h5py.File(whizzname, 'r+') as fw:
        wlines_group = fw[groupName]['Lines']
        wlines = list(wlines_group.keys())

        with h5py.File(planname, 'r+') as fp:
            plines_group = fp[groupName]['Lines']
            plines = list(plines_group.keys())

            trav_parity = True
            first = True

            for line in wlines:
                linegroup = wlines_group[line]
                if not (linegroup.attrs['LineVariety'] == 'Traverse'):
                    continue

                all_ok, plinegroup = util.getPlannedLine(plines_group, linegroup)
                if all_ok:
                    if 'Parity' in plines_group[plinegroup].attrs:
                        linegroup.attrs['Parity'] = plines_group[plinegroup].attrs['Parity']
                    else:
                        print('No Parity in ', plines_group[plinegroup])
                        continue
                    if verbose:
                        print(f'Line {line} parity set to {linegroup.attrs["Parity"]}')


def _approx_spacings(linegroup, firstgroup, trav_spacing, x, y, tolerance=0.2):
    """
    Returns the integer number of times that the spacing between the flight-lines specified
    by `linegroup` and `firstgroup` are divisible by `trav_spacing`. Returns 0 if the
    remainder is greater than `tolerance`.
    """

    # find best straight line fit to flight-line locations
    a, b, c = _line_fit(firstgroup, x, y)
    if a == 0.0 and b == 0.0:
        print(f'ERROR - flight-line start and end are the same point!')
        print(f'Coefficients of straight line fit to flight-line = {a:.1f}, {b:.1f}, {c:.1f}.')
        return 0.0
    # print(f'{np.arctan(a/b)*180/np.pi:.2f} deg')

    dx = rd.getLineData(linegroup, x)
    dy = rd.getLineData(linegroup, y)

    separation = abs(a * dx[0] + b * dy[0] + c) / np.sqrt(a * a + b * b)

    spacings = round(separation / trav_spacing)
    rem_ratio = abs(separation - spacings * trav_spacing) / trav_spacing

    if rem_ratio > tolerance:
        print(f'WARNING - line {linegroup} separation is not a good multiple of the traverse spacing.')
        print(f'  Trav spacing {trav_spacing:.1f}; line separation {separation:.1f}; rem_ratio {rem_ratio:.2f}.')
        spacings = 0
    return spacings


def _line_fit(linegroup, x, y):
    """
    Given a flight-line group, with channels x and y, return the coefficients
    a, b, and c for the straight line ax + by + c = 0 that passes through the
    points near the start and end of the flight-line.
    """

    dx = rd.getLineData(linegroup, x)
    dy = rd.getLineData(linegroup, y)

    # Crude attempt to avoid points near ends of line where aircraft might
    # have started turning early.
    num_fids = len(dx)
    if num_fids > 100:
        i1 = int(len(dx) / 100.0)
        i2 = len(dx) - i1
        a = dy[i1] - dy[i2]
        b = dx[i2] - dx[i1]
        c = dx[i1] * dy[i2] - dx[i2] * dy[i1]
    else:
        a = dy[0] - dy[-1]
        b = dx[-1] - dx[0]
        c = dx[0] * dy[-1] - dx[-1] * dy[0]
    return a, b, c


def _linelength(linegroup, x, y):

    dx = rd.getLineData(linegroup, x)
    dy = rd.getLineData(linegroup, y)
    return util._displacement2(dx[0], dx[-1], dy[0], dy[-1])


