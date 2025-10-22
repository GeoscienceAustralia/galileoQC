def _calc_gridcell_size(ds):
    """
    Returns a sensible grid cell size for gridding line-based data in xarray
    DataSet `ds`. Usefully provides a default for gridding functions.
    
    The algorithm is described here: "return 500.0"
    Later we will use the "TraverseLineSpacing", if available in the whizzfile.
    
    Parameters
    ----------
    ds : xarray 1D DataSet
        Airborne survey line-based DataSet indexed by fiducial.

    Returns
    -------
    cell_size : float
        the recommended cell size for gridding.
        
    """
    cell_size = 500.0
    
    return cell_size
    