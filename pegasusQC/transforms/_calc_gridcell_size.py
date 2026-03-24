#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate a sensible gridding cell size.

Author: Mark Helm Dransfield

Created: Oct 2025

License: CC BY-SA
"""

import numpy as np
from pathlib import Path
import xarray as xr
import h5py

from pegasusQC.whizzFiles.updateLineSpacing import updateLineSpacing
import pegasusQC.config as config

groupName = config.groupName


def _calc_gridcell_size(whizzFile):
    """
    Returns a sensible grid cell size for gridding line-based data in xarray
    DataSet `ds`. Usefully provides a default for gridding functions.
    
    The algorithm uses the "TraverseLineSpacing", calling updateLineSpacing to
    set it if it does not exist in the whizzfile, and returns 1/4 of that value.
    
    Parameters
    ----------
    whizzFile : Path or String

        The Path to, or String name of, the whizz file in HDF5 format.

    Returns
    -------
    cell_size : float

        the recommended cell size for gridding.
        
    """

    filename = str(whizzFile)
    with h5py.File(filename, 'r+') as f:
        if 'TraverseSpacing' in f[groupName]['CoordinateFrame'].attrs:
            trav_spacing = f[groupName]['CoordinateFrame'].attrs['TraverseSpacing']
            cell_size = _nice_number(trav_spacing / 4.0)

        else:
            updateLineSpacing(whizzFile)
            if 'TraverseSpacing' in f[groupName]['CoordinateFrame'].attrs:
                trav_spacing = f[groupName]['CoordinateFrame'].attrs['TraverseSpacing']
                cell_size = _nice_number(trav_spacing / 4.0)
            else:
                print('ERROR - could not calculate cell size for gridding.')
                print('    TraverseSpacing not set in whizzFile, and value not estimatable.')
                print('    Recommend running updateLineSpacing() with explicit values.')
                print('    And/or re-running craig_transform with explicit cell_size.')
                return None
    path_substring = filename.split('/')
    local_filename = path_substring[-1].rsplit('.', 1)
    print(f'\nGridding {local_filename[0]} data with cell size = {cell_size}.')
    return cell_size
    

def _nice_number(number):
    """
    Given any float, returns a reasonable value to be used for gridding. Currently
    just rounded to one decimal place.
    
    Parameters
    ----------
    number : float

        The number to be transformed.

    Returns
    -------
     : float

        the transformed value.
        
    """
    return round(number, 1)