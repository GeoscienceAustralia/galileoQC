#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pad two grids preparatory to FFT.
Author: Mark Helm Dransfield
Created: Oct 2025
License: CC BY-SA
"""

import numpy as np
import rioxarray # for the extension to load

from pegasusQC.transforms._pad_mean import _pad_mean
from pegasusQC.transforms._pad_regional import _pad_regional

def _pad_grids(Gne_grid, Guv_grid, pad_cells=0, mode='regional', regional_grid=None, firstorder=False):
    """
    `Gne_grid` and `Guv_grid` are mean-corrected and padded by `pad_cells` pixels in all four
    directions. If `pad_cells` < 1, or `mode` is unrecognised, or is "regional"
    but `regional_grid` is None, then `Gne_grid` and `Guv_grid` are mean-corrected but not padded.

    Parameters
    ----------
    Gne_grid : xarray 2D DataArray
        Gne grid to be expanded.
    Guv_grid : xarray 2D DataArray
        Guv grid to be expanded.
    pad_cells : int, optional
        Number of cells to add to the grid along each of the four sides.
        The resultant grid is thus larger by `2 x pad_cells` in both directions.
    mode : str, optional
        Plan is to have several modes but right now, this parameter is ignored
        and the xarray "mean" padding mode is used.
    regional_grid : xarray 2D DataArray, optional
        The regional grid. Required if `padding_mode` == "regional". Default None.
    firstorder : bool, optional
        If True, include first order Craig correction. Default False.

    Returns
    -------
    grid : xarray 2D DataArray
        the expanded, mean-corrected grid.
        
    """

    # No expansion: just mean-correct
    if pad_cells < 1:
        nanmask = np.isnan(Gne_grid.data)
        Gne_grid_padded = Gne_grid - Gne_grid.mean()
        Guv_grid_padded = Guv_grid - Guv_grid.mean(), nanmask
    
    # Here the mean-padding is along each vector, ...
    if mode == "mean":
        Gne_grid_padded, nanmask = _pad_mean(Gne_grid, pad_cells)
        Guv_grid_padded, nanmask = _pad_mean(Guv_grid, pad_cells)
        # nanmask = np.isnan(Gne_grid_padded.data)
    
    elif mode == "regional":
        Gne_grid_padded, Guv_grid_padded, nanmask = _pad_regional(Gne_grid,
            Guv_grid, pad_cells, regional_grid, firstorder=firstorder)

        # if nanmask is None:
        #     nanmask = np.isnan(Gne_grid_padded.data)
    else:
        print("\nWARNING in _pad_grid. The mode {mode} is unrecognised.")
        print("It should be one of 'mean' or 'regional'; no padding done.")
        nanmask = np.isnan(Gne_grid.data)
        Gne_grid_padded = Gne_grid - Gne_grid.mean()
        Guv_grid_padded = Guv_grid - Guv_grid.mean()

    return Gne_grid_padded, Guv_grid_padded, nanmask





