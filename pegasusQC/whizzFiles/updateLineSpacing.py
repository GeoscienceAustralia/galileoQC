#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Update the line-spacing info for the CoordinateFrame group in a `geoWhizz` file.
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


def updateLineSpacing(whizzFile, trav_spacing=0.0, ctrl_spacing=0.0):
    """
    Writes the mean traverse and mean control line-spacing as (float) attributes
    for the CoordinateFrame group in `whizzFile`.

    Parameters
    ----------
    whizzFile : Path or String
        The Path to, or String name of, the whizz file in HDF5 format.
    trav_spacing : Float, optional
        The traverse line spacing. Default (0.0) is to not set the TravSpacing attribute.
    ctrl_spacing : Float, optional
        The control line spacing. Default (0.0) is to not set the CtrlSpacing attribute.

    Returns
    -------
    Nothing.

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r+') as f:
        no_calc_needed = False
        if trav_spacing > 0.0:
            f[groupName]['CoordinateFrame'].attrs['TravSpacing'] = trav_spacing
        if ctrl_spacing > 0.0:
            f[groupName]['CoordinateFrame'].attrs['CtrlSpacing'] = ctrl_spacing



