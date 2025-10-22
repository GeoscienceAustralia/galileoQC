import numpy as np


def _check_extension(local, regional, pad_cells):
    """
    The `regional` grid must extend beyond the `local` grid
    for a distance of at least `pad_cells` x cell size in
    each direction.
   
    Parameters
    ----------
    local : xarray 2D DataArray
        Grid to be expanded.
    regional : xarray 2D DataArray
        The grid whose data will be used in the expansion.
    pad_cells : int, optional
        Number of cells to add to the grid along each of the four sides.
        The resultant grid is thus larger by `2 x pad_cells` in both directions.

    Returns
    -------
    boolean : True if regional is large enough, else False.
        
    """
    dx = np.abs(local.x[1] - local.x[0]) * pad_cells
    dy = np.abs(local.y[1] - local.y[0]) * pad_cells
    
    max_x =  max(local.x) + dx > max(regional.x)
    min_x =  min(local.x) - dx < min(regional.x)
    max_y =  max(local.y) + dy > max(regional.y)
    min_y =  min(local.y) - dy < min(regional.y)
    
    if max_x or min_x or max_y or min_y:
        print('Error - regional grid does not extend far enough for padding.')
        return False
    return True
