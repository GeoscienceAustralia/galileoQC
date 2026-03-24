#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pad a grid with NaNs.

Author: Mark Helm Dransfield

Created: Oct 2025

License: CC BY-SA
"""

import numpy as np


def _pad_grid_nans(grid, pad_cells, cell_size=None):
    """
    Adds `pad_cells` cells to every side of the rectangular `grid`,
    expanding it by (2 x `pad_cells`) in length and width. The additional
    cells all have data value NaN. The grid is also mean-corrected.
   
    Parameters
    ----------
    grid : xarray 2D DataArray

        Grid to be expanded.

    pad_cells : int

        Number of cells to add to the grid along each of the four sides.
        The resultant grid is thus larger by `2 x pad_cells` in both directions.

    cell_size : float, optional

        The width in metres of a single, square grid cell. Default None, in
        which case it is calculated from the grid coordinates. Preferably, it
        is passed since this avoids rounding errors.

    Returns
    -------
    xarray 2D DataArray : 

        Expanded, mean-corrected grid.
        
    """
    grid_mean = np.nanmean(grid.values)
    if cell_size is None:
        dx = grid.x.values[1] - grid.x.values[0]
        dy = grid.y.values[1] - grid.y.values[0]
    else:
        dx = cell_size
        dy = cell_size
    if (dx - dy) < 0.1:
        dy = dx
    e0 = grid.x.values[0] - dx * pad_cells
    n0 = grid.y.values[0] - dy * pad_cells
    ew = grid.x.values[-1] + dx * pad_cells
    nw = grid.y.values[-1] + dy * pad_cells
    nume = round((ew - e0) / dx) + 1
    numn = round((nw - n0) / dy) + 1
    e = np.linspace(e0, ew, nume, endpoint=True)
    n = np.linspace(n0, nw, numn, endpoint=True)
    parr = grid.pad(x=(pad_cells, pad_cells), 
                    y=(pad_cells, pad_cells), 
                    keep_attrs=True,
                    mode="constant",
                    fill_value=np.nan,
                    stat_length=None)
    parr = parr.assign_coords(x=e, y=n)   
    parr += -parr.mean()
    return parr
