#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Image all ERS grids in directory.
"""
"""
Created on Sat Jul 18 16:43:31 2020

author: Mark Dransfield

"""

from pathlib import Path
import colorcet as cc
import AirGravQC.config as config
from AirGravQC.gridFiles.gridfiles import gridfile_to_xa
from AirGravQC import xdImage


groupName = config.groupName
projectName = config.projectName

def imageAllInDir(path_name):
    """
    Quick and dirty method to image all the ERS grids in a 
    directory for QC.

    Parameters
    ----------
    path_name : Path
        The Path where the ERS grid files are located.

    Returns
    -------
    None.

    """
    # get a list of the ers file paths
    file_count = 0
    ersFiles = []
    folderFiles = path_name.iterdir()
    for aFile in folderFiles:
        if aFile.is_file() and (aFile.suffix == '.ERS' or aFile.suffix == '.ers') and (aFile.name[0] != '.'):
            ersFiles.append(aFile)
            file_count += 1
    print(f'Found {file_count} .ers files ...')
    print(f'in: {str(path_name)}')

    for f in ersFiles:
        (dxt, _) = gridfile_to_xa(f)
        xdImage(dxt, str(f.name), colormap=cc.m_CET_R1)


