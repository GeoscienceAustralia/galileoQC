

def _calc_padcells(grid):
    """
    Returns a sensible number of grid cells to be added as padding to
    a grid before discrete Fourier transforms.
    
    Smarter ideas may come. 
    
    Parameters
    ----------
    grid : xarray 2D DataArray
        Grid to be expanded.

    Returns
    -------
    pad_cells : int
        the number of grid cells to be added to each of the four
        sides of the grid.
        
    """
    pad_cells = max(grid.x.shape[0], grid.y.shape[0])
    
    return pad_cells
