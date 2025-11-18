#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the flight-lines are of the required minimum length.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import h5py

import pegasusQC.config as config
import pegasusQC.utility.utility as util
import pegasusQC.whizzFiles.retrieveData as rd

groupName = config.groupName
    

def checkLineLengths(whizzFile, min_len=50.0, measX='', measY='', lines=[]):
    """
    Checks that all lines in whizzFile are at least min_len km long.

    Parameters
    ----------
    whizzFile : String or pathlib Path

        Name of a HDF5 Whizz file, including path and extension, to be checked.

    min_len : TYPE, optional

        The minimum allowed line length in km. The default is 50.0.

    measX : String, optional

        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.

    measY : String, optional

        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.

    lines : Array{String}, optional

        Array of line numbers as strings. Default = [], meaning all lines are checked.

    Returns
    -------
    None.

    """
    measFile = str(whizzFile)
    
    with h5py.File(measFile, 'r') as f:
        gMeas = f[groupName]['Lines']
        if measX == '':
            measX = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if measY == '':
            measY = f[groupName]['CoordinateFrame'].attrs['YChannel']
        
        num_failed_lines = 0
        
        if lines == []:
            lines = list(gMeas.keys())

        for line in lines:
            xM = rd.getLineData(gMeas[line], measX)
            yM = rd.getLineData(gMeas[line], measY)
            line_length = util._displacement2(xM[0], xM[-1], yM[0], yM[-1])
            if line_length < min_len * 1000.0:
                num_failed_lines += 1
                print(f'Line {line} length = {line_length:.1f} less than allowed min {min_len*1000.0:.1f}')
        print(f'Number failed lines = {num_failed_lines}')
