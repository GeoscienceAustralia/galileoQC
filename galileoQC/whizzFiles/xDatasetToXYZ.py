#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Write an `xArray` DataSet to a `geoWhizz` file. A STUB - DON'T USE!
Author: Mark Helm Dransfield
Created: 2023
License: CC BY-SA
"""

import numpy as np
from pathlib import Path
import pathlib


def _xDatasetToXYZ(xd, chans, wids, precisions, xyzFilePath, line="LINE", verbose=False):
    """
     A STUB - DON'T USE! Read in a Whizz HDF5 file and write the contents to a Geosoft XYZ file.

    Parameters
    ----------
    xd : xArray Dataset
        Contains the data to be written..
    chans : string array
        Array of the channels names (which must be in `xd`) to be
        written to xyzFilePath.
    wids : number array
        Array of the widths (number of characters) of the formatted output channels.
    precisions : number array
        Array of the precisions (number after decimal point) of the formatted output channels.
    xyzFilePath : pathlib Path
        The pathlib Path to the output Geosoft XYZ file.
    line : string
        The name of the channel containing the line numbers.
    verbose : Bool, optional
        If False (the default) the output is reduced.

    Returns
    -------
    None.

    """

    # prepare file name
    if isinstance(xyzFilePath, pathlib.Path):
        xyzFileStr = str(xyzFilePath)
    elif isinstance(xyzFilePath, str):
        xyzFileStr = xyzFilePath
        xyzFilePath = Path(xyzFileStr)
    else:
        print('Error - type of xyzFilePath not recognised. Must be Path or String')
        return

    # get the unique line numbers from Dataset

    with open(filename, 'w') as fxyz:
        # write xyz header
        allchanstr = _getchanstr(xd, chans, omit='LINE')
        fxyz.write(f'/\n/\n/ {allchanstr}')

        currentline = 0.0
        fiducials = xd.fiducials.data
        for fid in fiducials:
            if line != currentline:
                fxyz.write(f' Line  {line}')
                currentline = line

            datastring = ''
            for idx, chan in enumerate(chans):
                value = xd[chan].sel(fiducials=slice(fid,fid)).values[0]
                chanstr = '{val:{wid}.{pr}f}'.format(wid=wids[idx], pr=precisions[idx], val=value)
                datastring.append(chanstr)
            fxyz.write(datastring)



def _getchanstr(xd, inputchans, omit='LINE'):
    """
    A string made of the variable names in the `xd` xArray Dataset is useful
    for writing data to XYZ format.

    Parameters
    ----------
    xd : xArray Dataset
        Contains the channel names to be checked against chans.
    inputchans : string array
        Array of the channels names (which must be in `xd`) to be
        written to `xyzFilePath`.
    omit : string
        The name of one channel to be omitted from the returned string.

    Returns
    -------
    chan_liststr : string
        A string of the channel names in inputchans (except for omit).

    """
    chan_liststr = ''
    for channame in inputchans:
        if channame == omit:
            continue
        if channame in list(xd.keys()):
            chan_liststr += channame + "  "
        else:
            print(f'ERROR - xDataset has {channame} channel but not found in input channel names.')
    return chan_liststr











