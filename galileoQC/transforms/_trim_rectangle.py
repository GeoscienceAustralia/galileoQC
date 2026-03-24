#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trim a grid along each edge.

Author: Mark Helm Dransfield

Created: Oct 2025

License: CC BY-SA
"""

import numpy as np

def _trim_rectangle(grid):
    """
    Returns a smaller grid trimmed along each edge. Simply find the
    smallest and largest easting and northing that have a non-Nan
    value, and trim to the rectangle so defined.
    
    """
    
    lenx = grid.sizes['x']
    leny = grid.sizes['y']

    icount = 0
    jcount = 0
    for i in range(0, lenx-1):
        ith_slice = grid.isel(x=slice(i,i+1)).data
        if not np.isnan(ith_slice).all():
            break
        else:
            icount += 1
    for j in range(0, leny-1):
        jth_slice = grid.isel(y=slice(j,j+1)).data
        if not np.isnan(jth_slice).all():
            break
        else:
            jcount += 1

    iend = lenx
    jend = leny
    for i in range(lenx, 0, -1):
        ith_slice = grid.isel(x=slice(i-1,i)).data
        if not np.isnan(ith_slice).all():
            break
        else:
            iend += -1
    for j in range(leny, 0, -1):
        jth_slice = grid.isel(y=slice(j-1,j)).data
        if not np.isnan(jth_slice).all():
            break
        else:
            jend += -1
            
    return grid.isel(x=slice(icount, iend), y=slice(jcount, jend))
    