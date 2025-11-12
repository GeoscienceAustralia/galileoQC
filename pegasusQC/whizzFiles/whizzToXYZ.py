#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Write the contents of a `geoWhizz` file to `XYZ` format.
Author: Mark Helm Dransfield
Created: 2023
License: CC BY-SA
"""

import numpy as np
import h5py
from pathlib import Path
import pathlib
import filebrowser as fb

import pegasusQC.utility.utility as util
import pegasusQC.config as config

groupName = config.groupName
projectName = config.projectName


def whizzToXYZ(whizzfilePath, chans, wids, precisions, xyzFilePath='', verbose=False):
    """
    Read in a `geoWhizz` HDF5 file and write the contents to a `XYZ` file.

    Parameters
    ----------
    whizzfilePath : pathlib Path
        The pathlib Path to the Whizz HDF5 file to be read.
    chans : string array
        Array of the channels names (which must be in whizzfilePath) to be
        written to xyzFilePath.
    wids : number array
        Array of the widths (number of characters) of the formatted output channels.
    precisions : number array
        Array of the precisions (number after decimal point) of the formatted output channels.
    xyzFilePath : pathlib Path
        The pathlib Path to the output Geosoft XYZ file. The default is ''
        leaving the XYZ file to have the same path as the input file but with
        a ".xyz" extension.
    verbose : Bool, optional
        If False (the default) the output is reduced.

    Returns
    -------
    None.

    """

    # prepare file names
    if whizzfilePath == '':
        whizzfilePath = fb.get_grid_filename()
    if isinstance(whizzfilePath, pathlib.Path):
        whizzFileStr = str(whizzfilePath)
    elif isinstance(whizzfilePath, str):
        whizzFileStr = whizzfilePath
        whizzfilePath = Path(whizzFileStr)
    else:
        print('Error - type of whizzfilePath not recognised. Must be Path or String')
        return
    if xyzFilePath == '':
        xyzFilePath = whizzfilePath.with_suffix('.xyz')
        xyzFileStr = str(xyzFilePath)
    elif isinstance(xyzFilePath, pathlib.Path):
        xyzFileStr = str(xyzFilePath)
    elif isinstance(xyzFilePath, str):
        xyzFileStr = xyzFilePath
        xyzFilePath = Path(xyzFileStr)
    else:
        print('Error - type of xyzFilePath not recognised. Must be Path or String')
        return

    with h5py.File(whizzFileStr, 'r') as fwhizz:
        lines_group = fwhizz[groupName]['Lines']
        lines = list(lines_group.keys())
        firstlinegroup = lines_group[lines[0]]
        whizzchans = list(firstlinegroup.keys())
        allchanstr = ''
        for chan in chans:
            if chan in whizzchans:
                allchanstr.append(f'{chan}  ')
            else:
                print(f'ERROR: channel {chan} not in whizzfile.')
                return

        with open(filename, 'w') as fxyz:
            # write xyz header
            fxyz.write(f'/\n/\n/ {allchanstr}')

            for line in lines:
                fxyz.write(f' Line  {line}')

                for index in range(0, len(line[chan[0]].data)):
                    datastring = ''
                    for idx, chan in enumerate(chans):
                        chanstr = '{val:{wid}.{pr}f}'.format(wid=wids[idx], pr=precisions[idx], val=line[chan].data[index])
                        datastring.append(chanstr)
                    fxyz.write(datastring)















