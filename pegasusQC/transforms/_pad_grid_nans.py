#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pad a grid with NaNs.

Author: Mark Helm Dransfield

Created: Oct 2025

License: CC BY-SA
"""

import numpy as np


def _pad_grid_nans(grid, pad_cells):
    """
    Adds `pad_cells` cells to every side of the rectangular `grid`,
    expanding it by (2 x `pad_cells`) in length and width. The additional
    cells all have data value NaN. The grid is also mean-corrected.
   
    Parameters
    ----------
    grid : xarray 2D DataArray

        Grid to be expanded.

    pad_cells : int, optional

        Number of cells to add to the grid along each of the four sides.
        The resultant grid is thus larger by `2 x pad_cells` in both directions.

    Returns
    -------
    xarray 2D DataArray : 

        Expanded, mean-corrected grid.
        
    """
    grid_mean = np.nanmean(grid.values)
    parr = grid.pad(x=(pad_cells, pad_cells), 
                    y=(pad_cells, pad_cells), 
                    keep_attrs=True,
                    mode="constant",
                    fill_value=np.nan,
                    stat_length=None)
    dx = grid.x.values[1] - grid.x.values[0]
    dy = grid.y.values[1] - grid.y.values[0]
    for i in range(0,pad_cells):
        parr.x.values[i] = parr.x.values[pad_cells] - (pad_cells - i) * dx
        parr.x.values[-1-i] = parr.x.values[-1-pad_cells] + (pad_cells - i) * dx
        parr.y.values[i] = parr.y.values[pad_cells] - (pad_cells - i) * dy
        parr.y.values[-1-i] = parr.y.values[-1-pad_cells] + (pad_cells - i) * dy

    # ... so need to mean-correct post-padding
    return parr - parr.mean()
