#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pad a grid with regional data.

Author: Mark Helm Dransfield

Created: Oct 2025

License: CC BY-SA
"""

import numpy as np

from pegasusQC.transforms._check_extension import _check_extension
from pegasusQC.transforms._pad_mean import _pad_mean
from pegasusQC.transforms.Kurvs_of_grav import Kurvs_of_grav
from pegasusQC.transforms._pad_grid_nans import _pad_grid_nans
from pegasusQC.transforms._grid_match import _grid_match
from pegasusQC.gridFiles.gridutility import report_gridStats
from pegasusQC.gridFiles.xdImage import xdImage
from pegasusQC.transforms._trim_rectangle import _trim_rectangle


def _pad_simple(local_grid, regional_grid, verbose=False):
    """
    `local_grid` is mean-corrected and padded to same area as `regional_grid`.

    If `local_grid` is of terrain-corrected data, then the `regional_data`
    should be too. Similarly for free-air data.
    
    The size of `regional_grid` is assumed appropriate for the work. Any checking
    of this assumption should be done prior to calling `_pad_simple.
    
    Parameters
    ----------
    local_grid : xarray 2D DataArray

        gravity grid to be expanded.

    regional_grid : xarray 2D DataArray

        Regional gravity grid

    verbose : Bool, optional

        If True, then prints out details. Default = False.

    Returns
    -------
    local_x_grid : xarray 2D DataArray

        the expanded grid.
        
    """
    """
    
    - NaN-pad `local_grid`
    - interpolate `regional_grid` onto locations of padded `local_grid`
    - constant-correct `regional_grid` so that its mean over coincident locations is zero
    - replace pad-zeros of padded `local_grid` with values from `regional_grid`
    - if any NaNs remain in padded `local_grid`
        - mean-pad
    
    """

    cell_size = local_grid.x.values[1] - local_grid.x.values[0]
    pad_xm = round((min(local_grid.x.data) - min(regional_grid.x.data)) / cell_size)
    pad_xp = round((max(regional_grid.x.data) - max(local_grid.x.data)) / cell_size)
    pad_ym = round((min(local_grid.y.data) - min(regional_grid.y.data)) / cell_size)
    pad_yp = round((max(regional_grid.y.data) - max(local_grid.y.data)) / cell_size)

    # Pad-extend local grid
    local_grid.data = local_grid.data - local_grid.mean().data
    parr = local_grid.pad(x=(pad_xm, pad_xp), 
                        y=(pad_ym, pad_yp), 
                        keep_attrs=True,
                        mode="constant",
                        fill_value=np.nan,
                        stat_length=None)
    for i in range(0, pad_xm):
        parr.x.values[i] = parr.x.values[pad_xm] - (pad_xm - i) * cell_size
    for i in range(0,pad_xp):
        parr.x.values[-1-i] = parr.x.values[-1-pad_xp] + (pad_xp - i) * cell_size
    for i in range(0,pad_ym):
        parr.y.values[i] = parr.y.values[pad_ym] - (pad_ym - i) * cell_size
    for i in range(0,pad_yp):
        parr.y.values[-1-i] = parr.y.values[-1-pad_yp] + (pad_yp - i) * cell_size

    # get the locations of NaNs in the local data
    data_mask = np.isnan(parr.data)
    data_unmask = np.logical_not(data_mask)

    # interpolate the regional grid onto the same locations
    reg_match = _grid_match(parr, regional_grid)
    reg_mean = reg_match.where(data_unmask, np.nan).mean()

    test = reg_match - parr - reg_mean - local_grid.mean()
    print('\n Grid stats: regional subtract local.')
    print(' Ideally, the mean should be close to 0 (abs value < 1),')
    print(' and the range should be roughly symmetrical.')
    if verbose:
        report_gridStats(test)
    xdImage(_trim_rectangle(test), 'Regional subtract local', hs=False)

    # replace nans with zero in local grids for sum
    local_arith = parr.fillna(0.0)

    # zero out regional data where there are not local NaNs
    reg_masked = data_mask * reg_match
    
    # add local and regional, replacing NaNs in local
    local_final = reg_masked + local_arith

    return local_final, reg_match, data_mask


def _scaling_doubtful(grid1, grid2):
    """
    If the ratio of the ranges (max - min) of the input grids is
    not nearly 1.0, then the scaling is doubtful.
    
    Parameters
    ----------
    grid1 : xarray 2D DataArray

        A grid for comparison.

    grid2 : xarray 2D DataArray

        A grid for comparison.

    Returns
    -------
    True if abs(log10(ratio of ranges)) > 0.5, else False.
        
    """
    range1 = np.abs(grid1.max().data.item() - grid1.min().data.item())
    range2 = np.abs(grid2.max().data.item() - grid2.min().data.item())
    if abs(np.log10(range1 / range2)) > 0.5:
        return True

    return False

    
    