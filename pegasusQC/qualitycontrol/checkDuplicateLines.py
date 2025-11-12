#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check whizzFile for duplicate flight-lines.
Author: Mark Helm Dransfield
Created: ca 2023
License: CC BY-SA
"""

def checkDuplicateLines(whizzFile):
    """
    Check the flight-lines in whizzfile and report any duplicate group names.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r+') as f:
        # create all the data structure ready for the datasets
        g = f[groupName]['Lines']
        flightlines = list(g.keys())
    # `set` removes non-unique elements. 
    seen = set()

    # Check for duplicates
    unique = True
    for x in flightlines:
        if x in seen:
            unique = False
            print(f'WARNING - duplicate flight-line {x} in {filename}.')
            break
        seen.add(x)            
