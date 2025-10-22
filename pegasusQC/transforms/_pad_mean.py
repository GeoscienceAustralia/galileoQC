import numpy as np

from pegasusQC.transforms._trim_rectangle import _trim_rectangle

def _pad_mean(grid, pad_cells):
    """
    The `grid` is mean-corrected and padded by `pad_cells` pixels in all four
    directions.

    Parameters
    ----------
    grid : xarray 2D DataArray
        Grid to be expanded.
    pad_cells : int, optional
        Number of cells to add to the grid along each of the four sides.
        The resultant grid is thus larger by `2 x pad_cells` in both directions.

    Returns
    -------
    xarray 2D DataArray
        the expanded, mean-corrected grid.
        
    """
    # Trim edges of inputs to remove NaNs, 
    grid = _trim_rectangle(grid)

    # First, fill any holes in the rectangular grid array
    # print(Gne_grid)
    # rioxarray does this really well but it wants a CRS first
    grid.rio.write_crs(4283, inplace=True)
    grid.rio.write_nodata(np.nan, encoded=True, inplace=True)
    grid_fill = grid.rio.interpolate_na()

    grid_mean = np.nanmean(grid_fill.values)
    parr = grid_fill.pad(x=(pad_cells, pad_cells),
                         y=(pad_cells, pad_cells), 
                         keep_attrs=True,
                         mode="mean",
                         fill_value=grid_mean,
                         stat_length=None)
    dx = grid_fill.x.values[1] - grid_fill.x.values[0]
    dy = grid_fill.y.values[1] - grid_fill.y.values[0]
    for i in range(0,pad_cells):
        parr.x.values[i] = parr.x.values[pad_cells] - (pad_cells - i) * dx
        parr.x.values[-1-i] = parr.x.values[-1-pad_cells] + (pad_cells - i) * dx
        parr.y.values[i] = parr.y.values[pad_cells] - (pad_cells - i) * dy
        parr.y.values[-1-i] = parr.y.values[-1-pad_cells] + (pad_cells - i) * dy

    # ... so need to mean-correct post-padding
    return parr - parr.mean()
