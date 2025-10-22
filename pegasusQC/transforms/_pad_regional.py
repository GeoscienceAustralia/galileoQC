import numpy as np

from pegasusQC.transforms._check_extension import _check_extension
from pegasusQC.transforms._pad_mean import _pad_mean
from pegasusQC.transforms.Kurvs_of_grav import Kurvs_of_grav
from pegasusQC.transforms._pad_grid_nans import _pad_grid_nans
from pegasusQC.transforms._grid_match import _grid_match


def _pad_regional(Gne_grid, Guv_grid, pad_cells, regional_grid):
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
        print(f'WARNING - Regional {regional_grid} too small for padding by {pad_cells}.')
        print(f'    Using mean padding instead!')
        Gne_x_grid = _pad_mean(Gne_grid, pad_cells)        
        Guv_x_grid = _pad_mean(Guv_grid, pad_cells)        
        return Gne_x_grid, Guv_x_grid, None       

    # Transform regional to differential curvatures
    gne_reg, guv_reg = Kurvs_of_grav(regional_grid)

    # Pad-extend local grids
    Gne_padded = _pad_grid_nans(Gne_grid, pad_cells)
    Guv_padded = _pad_grid_nans(Guv_grid, pad_cells)

    # interpolate the regional grids onto the same locations
    gne_reg_match = _grid_match(Gne_padded, gne_reg)
    guv_reg_match = _grid_match(Guv_padded, guv_reg)

    # get the locations of NaNs in the local data
    data_mask_ne = np.isnan(Gne_padded.data)
    # plt.imshow(data_mask_ne)
    data_mask_uv = np.isnan(Guv_padded.data)

    # replace nans with zero in local grids for transform
    Gne_arith = Gne_padded.fillna(0.0)
    Guv_arith = Guv_padded.fillna(0.0)

    # zero out regional data on local NaNs
    gne_reg_masked = data_mask_ne * gne_reg_match
    guv_reg_masked = data_mask_uv * guv_reg_match
    
    # add local and regional, replacing NaNs in local
    tne = gne_reg_masked + Gne_arith
    tuv = guv_reg_masked + Guv_arith

    return tne, tuv, data_mask_ne
    