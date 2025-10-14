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


def updateLineVariety(whizzFile, trav_lines=[], ctrl_lines=[]):
    """
    Set the 'LineVariety' metadata for each of the `trav_lines` to 'Traverse', 
    and for each of the `ctrl_lines` to 'Control'. Existing metadata are
    overwritten.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz observations file, including path and extension.
    trav_lines : Array{str}, optional
        Array of traverse line numbers as strings. Default = [], and no action taken.
    ctrl_lines : Array{str}, optional
        Array of control line numbers as strings. Default = [], and no action taken.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    trav_line_found = False
    ctrl_line_found = False

    with h5py.File(filename, 'r+') as f:
        # all the lines are in g:
        g = f[groupName]['Lines']
        lines = g.keys()

        for line in lines:
            gg = g[line]
            if line in trav_lines:
                gg.attrs['LineVariety'] = 'Traverse'
                trav_line_found = True
            elif line in ctrl_lines:
                gg.attrs['LineVariety'] = 'Control'
                ctrl_line_found = True

    if trav_line_found:
        print('"Traverse" line variety set for given flight-lines.')

    if ctrl_line_found:
        print('"Control" line variety set for given flight-lines.')

    return



