#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates the vertical gravity from the differential curvatures.

Author: Mark Helm Dransfield

Created: Oct 2025

License: CC BY-SA
"""

import numpy as np
import xarray as xr

from galileoQC.gridFiles.grid_to_xarray import gridfile_to_xa
from galileoQC.gridFiles.xarray_to_grid import xarray_to_grid
from galileoQC.transforms._calc_padcells import _calc_padcells
from galileoQC.transforms._pad_grids import _pad_grids
from galileoQC.gridFiles.xdImage import xdImage
from galileoQC.transforms.grav_of_Kurvs import grav_of_Kurvs
from galileoQC.gridFiles.gridutility import report_gridStats


def gravity_from_curv(Ane, Auv, cell_size,
                    gd_chan=None,
                    altitude=None,
                    result_units='um/s/s', 
                    survey_polygon=None,
                    pad_cells=None,
                    padding_mode='regional',
                    regional_grid=None,
                    firstorder=False,
                    plot_flag=False,
                    verbose=False
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

    gd_chan : String, optional

        The `long name` attribute of the gD output DataArray.
        Default None in which case, the code uses `gd_craig`.

    altitude : xarray 1D DataSet, optional

        Airborne survey line-based DataSet of the altitude dimensioned
        by fiducial. These must be referenced to a fixed datum (not
        the varying ground surface). Default None in which case the Gne
        and Guv data are assumed at constant altitude.
        If given, altitude is not used. TBD

    result_units : String, optional

        The gravity units of the final resultant grid. Must be either "mGal" or
        "gu" or "um/s/s". Default "um/s/s".

    survey_polygon : numpy 1D array, optional

        The polygon vertices of the survey boundary, as an array or sequence of (x,y)
        pairs, in either clockwise or counter-clockwise order around the boundary.
        For example, survey_polygon = [(0, 0), (1, 0), (1,1), (0,1)]. Final output will be
        trimmed to this polygon if provided. Default None.

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

    regional_grav_units : String, optional

        The gravity units of the regional grid. Must be either "mGal" or
        "gu" or "um/s/s". Required if `padding_mode` == "regional".
        Default "mGal".

    firstorder : bool, optional

        If True, include first order Craig correction. Default False.
    
    plot_flag : Bool, optional

        If True, images will be displayed of the grids at each stage of the processing. If False,
        then only the `local` input and output grids will be displayed. Default False.
        
    verbose : Bool, optional

        If True, then prints out details. Default = False.

    Returns
    -------
    gD : xarray 2D DataArray

        Grid of vertical gravity, gD.

    gD_err : xarray 2D DataArray

        Grid of the imaginary component output of the Craig transform.
        Useful for error estimation.

    """
    e_chan = Ane.attrs['x_channel']
    n_chan = Ane.attrs['y_channel']

    mode = padding_mode
    if mode == 'mean':
        regional_grid = None
        input_mask_pixels = 5
    elif mode == 'regional':
        input_mask_pixels = 3
    elif mode == 'constant':
        regional_grid = None
        input_mask_pixels = 5
        mode = 'mean'
        print(f'\nWARNING - mode {mode} not available - using "mean".')
    elif mode == 'iterative':
        regional_grid = None
        input_mask_pixels = 5
        mode = 'mean'
        print(f'\nWARNING - mode {mode} not available - using "mean".')
    else:
        print(f'\nERROR - mode {mode} not recognised.')
        return
    
    if pad_cells is None:
        pad_cells = _calc_padcells(Ane)
        
    # Find the coverage for each gridded data set
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
    print('\nInterpolating the local curvature gradients to grids.')
    Ane_grid, Ane_region = xarray_to_grid(
        Ane, cell_size, region=intersectregion,
        method='minc', mask_polygon=survey_polygon, 
        mask_pixels=input_mask_pixels, 
        numneighbours=5, bdist=None, maxiters=100
        )
    Auv_grid, Auv_region = xarray_to_grid(
        Auv, cell_size, region=intersectregion,
        method='minc',mask_polygon=survey_polygon, 
        mask_pixels=input_mask_pixels, 
        numneighbours=5, bdist=None, maxiters=100
        )
    if verbose:
        report_gridStats(Ane_grid, title="Ane_grid")
        report_gridStats(Auv_grid, title="Auv_grid")
    # Pad and fill the data
    print('\nPadding the local curvature gradient grids.')
    Ane_grid_pad, Auv_grid_pad, nan_mask = _pad_grids(
        Ane_grid, Auv_grid, pad_cells,
        mode=mode, regional_grid=regional_grid,
        firstorder=firstorder,
        verbose=verbose
        )
    if verbose:
        report_gridStats(Ane_grid_pad, title="Ane_grid_pad")
        report_gridStats(Auv_grid_pad, title="Auv_grid_pad")
    
    # Transform to gD
    print('\nTransforming the local curvature gradient grids to gravity.')
    if survey_polygon is None:
        survey_polygon = [
        (intersectregion[0], intersectregion[2]),
        (intersectregion[1], intersectregion[2]),
        (intersectregion[1], intersectregion[3]),
        (intersectregion[0], intersectregion[3])
        ]
    gD_grid, gD_err, gD_raw = grav_of_Kurvs(
        Ane_grid_pad, Auv_grid_pad,
        firstorder=firstorder, 
        survey_polygon=survey_polygon, 
        nan_mask=nan_mask,
        plot_flag=plot_flag,
        verbose=True
        )
    if plot_flag:
        xdImage(gD_grid, 'Post-gok: gD_grid (Em)', hs=False)

    # ..., scale to desired units
    # ... 10,000 Em per mGal
    if result_units in ('mGal', 'mGal'):
        gD_grid = gD_grid / 10000.0
        gD_grid.attrs['Units'] = 'mGal'
        gD_err = gD_err / 10000.0
        gD_err.attrs['Units'] = 'mGal'
        gD_raw = gD_raw / 10000.0
        gD_raw.attrs['Units'] = 'mGal'
    else:
        gD_grid = gD_grid / 1000.0
        gD_grid.attrs['Units'] = 'um/s/s'
        gD_err = gD_err / 1000.0
        gD_err.attrs['Units'] = 'um/s/s'
        gD_raw = gD_raw / 1000.0
        gD_raw.attrs['Units'] = 'um/s/s'

    # ..., and store attributes
    if gd_chan is None:
        gd_chan = 'gD_craig'
    imag_chan = f'{gd_chan}_imag'

    gD_grid.attrs['x_channel'] = e_chan
    gD_grid.attrs['y_channel'] = n_chan
    gD_grid.attrs['long_name'] = gd_chan
    gD_grid.attrs['title'] = 'gD via Craig transform'

    gD_err.attrs['x_channel'] = e_chan
    gD_err.attrs['y_channel'] = n_chan
    gD_err.attrs['long_name'] = imag_chan
    gD_err.attrs['title'] = 'gD imaginary part'

    gD_raw.attrs['x_channel'] = e_chan
    gD_raw.attrs['y_channel'] = n_chan
    gD_raw.attrs['long_name'] = 'gD_craig_raw'
    gD_raw.attrs['title'] = 'gD via Craig transform (unmasked)'

    return gD_grid, gD_err, gD_raw
    