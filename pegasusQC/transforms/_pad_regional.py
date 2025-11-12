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


def _pad_regional(Gne_grid, Guv_grid, pad_cells, regional_grid, firstorder=False):
    """
    `Gne_grid` and `Guv_grid` are mean-corrected and padded by `pad_cells`
    pixels in all four directions using values calculated from the `regional_grid`.

    If `Gne_grid` and `Guv_grid` are of terrain-corrected data, then the
    `regional_data` should be too. Similarly for free-air data.
    
    If the `regional_grid` is too small to pad the required number of cells,
    then padding will be applied to the extent of `regional_grid`, and then
    mean padding will be used to achieve the remaining padding.
    
    Parameters
    ----------
    Gne_grid : xarray 2D DataArray
        Gne grid to be expanded.
    Guv_grid : xarray 2D DataArray
        Guv grid to be expanded.
    pad_cells : int
        Number of cells to add to the grid along each of the four sides.
        The resultant grid is thus larger by `2 x pad_cells` in both directions.
    regional_grid : xarray 2D DataArray
        Regional (vertical) gravity grid
    firstorder : bool, optional
        If True, include first order Craig correction. Default False.

    Returns
    -------
    Gne_x_grid, Guv_x_grid : (xarray 2D DataArray, xarray 2D DataArray)
        the expanded grids.
        
    """
    """
    
    - mean-correct `grid`
    - NaN-pad `grid`
    - transform `regional_grid` to `component`, R_c
    - interpolate R_c onto locations of padded `grid`
    - constant-correct R_c so that its mean over coincident locations is zero
    - replace pad-zeros of padded `grid` with values from R_c
    - if any NaNs remain in padded `grid`
        - mean-pad
    
    """
    # Check regional grid is big enough
    if not _check_extension(Gne_grid, regional_grid, pad_cells):
        print(f'\nWARNING - Regional grid too small for padding by {pad_cells}.')
        print(f'    Using mean padding instead!')
        Gne_x_grid, data_mask_ne = _pad_mean(Gne_grid, pad_cells)        
        Guv_x_grid, data_mask_uv = _pad_mean(Guv_grid, pad_cells)        
        return Gne_x_grid, Guv_x_grid, data_mask_ne

    # Transform regional to differential curvatures
    gne_reg, guv_reg = Kurvs_of_grav(regional_grid, firstorder=firstorder)

    print('\nGrid Statistics')
    print('  Gne regional')
    report_gridStats(gne_reg)
    print('  Gne local')
    report_gridStats(Gne_grid)
    print('  Guv regional')
    report_gridStats(guv_reg)
    print('  Guv local')
    report_gridStats(Guv_grid)
    print('\n')

    if _scaling_doubtful(gne_reg, Gne_grid):
        print('\nWARNING - grid statistical ranges suggest units error.\n')

    # Pad-extend local grids
    Gne_padded = _pad_grid_nans(Gne_grid, pad_cells)
    Guv_padded = _pad_grid_nans(Guv_grid, pad_cells)

    # interpolate the regional grids onto the same locations
    gne_reg_match = _grid_match(Gne_padded, gne_reg)
    guv_reg_match = _grid_match(Guv_padded, guv_reg)

    # get the locations of NaNs in the local data
    data_mask_ne = np.isnan(Gne_padded.data)
    data_mask_uv = np.isnan(Guv_padded.data)

    # get the locations of not-NaNs in the local data
    # data_good_ne = not data_mask_ne
    # data_good_uv = not data_mask_uv

    test_ne = gne_reg_match - Gne_padded
    print('\n Grid stats Gne regional subtract local')
    report_gridStats(test_ne)
    xdImage(_trim_rectangle(test_ne), 'Gne regional subtract local (E)', hs=False)

    test_uv = guv_reg_match - Guv_padded
    print('\n Grid stats Guv regional subtract local')
    report_gridStats(test_uv)
    xdImage(_trim_rectangle(test_uv), 'Guv regional subtract local (E)', hs=False)

    # replace nans with zero in local grids for transform
    Gne_arith = Gne_padded.fillna(0.0)
    Guv_arith = Guv_padded.fillna(0.0)

    # zero out regional data where there are local NaNs
    gne_reg_masked = data_mask_ne * gne_reg_match
    guv_reg_masked = data_mask_uv * guv_reg_match
    
    # add local and regional, replacing NaNs in local
    tne = gne_reg_masked + Gne_arith
    tuv = guv_reg_masked + Guv_arith

    return tne, tuv, data_mask_ne


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

    
    