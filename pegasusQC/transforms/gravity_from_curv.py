import numpy as np
import xarray as xr

from pegasusQC.gridFiles.grid_to_xarray import gridfile_to_xa
from pegasusQC.gridFiles.xarray_to_grid import xarray_to_grid
from pegasusQC.transforms._calc_padcells import _calc_padcells
from pegasusQC.transforms._pad_grids import _pad_grids
from pegasusQC.transforms.grav_of_Kurvs import grav_of_Kurvs


def gravity_from_curv(Ane, Auv, cell_size, altitude=None,
                      unit_scale=1000.0, 
                      mask_polygon=None,
                      pad_cells=None,
                      padding_mode='regional',
                      regional_grid_file=None
):
    """
    Calculates the vertical gravity, gD, from the differential curvature
    gravity gradient components (in NED coordinates). The `Craig transform`
    is used.
    
    Parameters
    ----------
    Ane : xarray 1D DataSet
        Airborne survey line-based DataSet of the Gne component
        dimensioned by fiducial.
    Auv : xarray 1D DataSet
        Airborne survey line-based DataSet of the Guv component
        dimensioned by fiducial.
    cell_size : float
        the cell size for gridding.
    altitude : xarray 1D DataSet, optional
        Airborne survey line-based DataSet of the altitude dimensioned
        by fiducial. These must be referenced to a fixed datum (not
        the varying ground surface). Default None in which case the Gne
        and Guv data are assumed at constant altitude.
        If given, altitude is not used. TBD
    unit_scale : float, optional
        The resultant gravity data are multiplied by `unit_scale`.
        Usually the gradients are in eotvos, and the positions are
        in metres, giving gravity in nm/s/s. A `unit_scale' of 1000.0
        results in gravity units of um/s/s. Default 1000.0.
    mask_polygon : numpy 1D array, optional
        In order [min_x, max_x, min_y, max_y]. If the mask_polygon is given,
        then the output arrays will be masked to the area within the polygon
        defined by it. Default None.
    pad_cells : int, optional
        The number of grid cells to pad the data. Default None in which case
        an optimum number is calculated internally.
    padding_mode : String, optional
        The method to be used in filling padding around the grid before DFT.
        Choices are "mean" and "regional". Default is "regional" which requires
        that `regional_grid_file` is provided.
    regional_grid_file : String, optional
        Name of the ERS file containing the regional grid. Required if
        `padding_mode` == "regional". Default None.
        

    Returns
    -------
    gD : xarray 2D DataArray
        Grid of vertical gravity, gD.
    gD_err : xarray 2D DataArray
        Grid of the imaginary component output of the Craig transform.
        Useful for error estimation.

    """
    mode = padding_mode
    if mode == 'mean':
        regional_grid = None
        input_mask_pixels = 5
    elif mode == 'regional':
        reg_xa, _ = gridfile_to_xa(regional_grid_file)
        # Nicer if we use the same coordinates for everyone!
        regional_grid = xr.DataArray(data=reg_xa.data,
                                     dims=["y", "x"],
                                     coords={
                                         "y": reg_xa.N.values,
                                         "x": reg_xa.E.values
                                     })
        input_mask_pixels = 2
    elif mode == 'constant':
        regional_grid = None
        input_mask_pixels = 5
        mode = 'mean'
        print(f'WARNING - mode {mode} not available - using "mean".')
    elif mode == 'iterative':
        regional_grid = None
        input_mask_pixels = 5
        mode = 'mean'
        print(f'WARNING - mode {mode} not available - using "mean".')
    else:
        print(f'ERROR - mode {mode} not recognised.')
        return
    
    if pad_cells is None:
        pad_cells = _calc_padcells(Ane)
        
    # Find the coverage for each gridded data set
    e_chan = Ane.attrs['x_channel']
    n_chan = Ane.attrs['y_channel']
    ne_region = [
                    np.min(Ane[e_chan].values),
                    np.max(Ane[e_chan].values),
                    np.min(Ane[n_chan].values),
                    np.max(Ane[n_chan].values)
                ]
    uv_region = [
                    np.min(Auv[e_chan].values),
                    np.max(Auv[e_chan].values),
                    np.min(Auv[n_chan].values),
                    np.max(Auv[n_chan].values)
                ]

    # We are only interested in the data within the intersection of the regions.
    intersectregion = [
                        max(ne_region[0], uv_region[0]),
                        min(ne_region[1], uv_region[1]),
                        max(ne_region[2], uv_region[2]),
                        min(ne_region[3], uv_region[3]),
                        ]
    
   # Grid the data in common region
    Ane_grid, Ane_region = xarray_to_grid(Ane, cell_size, region=intersectregion,
                                             method='minc', mask_polygon=[], 
                                             mask_pixels=input_mask_pixels, 
                                             numneighbours=5, bdist=None, maxiters=100
                 )
    Auv_grid, Auv_region = xarray_to_grid(Auv, cell_size, region=intersectregion,
                                             method='minc',mask_polygon=[], 
                                             mask_pixels=input_mask_pixels, 
                                             numneighbours=5, bdist=None, maxiters=100
                 )

    # Pad and fill the data
    Ane_grid_pad, Auv_grid_pad, nan_mask = _pad_grids(Ane_grid, Auv_grid, pad_cells,
                              mode=mode, regional_grid=regional_grid)
    
    # Transform to gD
    if mask_polygon is None:
        mask_polygon = intersectregion
    gD_grid, gD_err = grav_of_Kurvs(Ane_grid_pad, Auv_grid_pad, mask_polygon=mask_polygon, nan_mask=nan_mask)

    return gD_grid/unit_scale, gD_err/unit_scale, Ane_grid
    