import numpy as np
import matplotlib.path as mpltPath


def report_gridStats(my_grid, mask_polygon=[]):
    """
    Reports the standard deviation of the z values in the portion of my_grid inside the mask_polygon.

    Parameters
    ----------
    my_grid : 2D xarray DataArray
        The data to be reported.
    mask_polygon : numpy 2D array, optional
        If the size of mask_polygon > 0, then my_grid will be masked to the area
        within the polygon defined by it.

    Returns
    -------
    None.

    """
    if np.array(mask_polygon).size > 0:
        my_grid = maskGridByPolygon(my_grid, mask_polygon, x_chan='x', y_chan='y')
    print(f'RMS of result = {my_grid.std().data.item():.2f} {my_grid.attrs["units"]}')


def maskGridByPolygon(my_grid, mask_polygon, x_chan='x', y_chan='y'):
    """
    Sets the z values of my_grid to NaN everywhere outside mask_polygon.

    Parameters
    ----------
    my_grid : 2D xarray DataArray
        The data to be reported.
    channel : String
        The channel or field name to analyse. Must exist in `whizz_file`.
    grid_space : Float
        The width of the grid cell to be used in gridding. Recommend: 1/5 - 1/4 line spacing.
    oddlines : Array of String, optional
        An array of line numbers that will constitute the odd lines. The default is NOT WORKING! ...to take the first,
        and then every second traverse thereafter.
    evenlines : Array of String, optional
        An array of line numbers that will constitute the even lines. The default is NOT WORKING! ...to take every
        second traverse (alternates to the oddlines).
    mask_polygon : numpy 2D array, optional
        If the size of mask_polygon > 0, then data_array will be masked to the area
        within the polygon defined by it.

    Returns
    -------
    None.

    """
    xlen = my_grid[x_chan].values.shape[0]
    ylen = my_grid[y_chan].values.shape[0]
    xx, yy = np.meshgrid(my_grid[x_chan].values, my_grid[y_chan].values)
    xflat = np.reshape(xx, (ylen * xlen, ))
    yflat = np.reshape(yy, (ylen * xlen, ))
    points = np.transpose(np.vstack((xflat, yflat)))
    path = mpltPath.Path(mask_polygon)
    inside2 = path.contains_points(points)
    mymask = np.reshape(inside2, (ylen, xlen))
    return my_grid.where(mymask)


def automask(xd_grid, xd_located, width):
    """
    Sets the z-value in xd_grid to NaN for all elements whose x position
    and y position are both further away than width from all points in
    xd_located.

    Parameters
    ----------
    xd_grid : 2D xarray DataArray
        The data to be masked.
    xd_located : 1D xarray DataArray
        The located data whose x,y locations decide the masking.
    width : Float
        The search distance in each of x and y used in the masking calculation.

    Returns
    -------
    None.


    for j in grid_pts:
        if value_j != np.nan:
            close_j = False
            for i in located_pts:
                if |x_i - x_j| < width and |y_i - y_j| < width:
                    close_j = True
                    break
            if not close_j:
                z_j = NaN
    return grid

    """