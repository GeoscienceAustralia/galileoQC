#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Conform a local gravity grid to a regional grid, both read from file.

Author: Mark Helm Dransfield

Created: Feb 2026

License: CC BY-SA
"""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import xrft

from galileoQC.transforms.cos2_filter import cos_square_radial_filter
from galileoQC.transforms._pad_grid_nans import _pad_grid_nans
from galileoQC.transforms._trim_rectangle import _trim_rectangle
from galileoQC.transforms._grid_match import _grid_match
from galileoQC.gridFiles.xdImage import xdImage
from galileoQC.gridFiles.grid_to_xarray import gridfile_to_xa
from galileoQC.gridFiles.gridutility import report_gridStats


def conform_from_file(
    loc_file, reg_file, loc_units, reg_units, survey_polygon, 
    rim=16, low_lambda=None, hi_lambda=None, pad_cells=None, plot_flag=False
):
    """
    Conforms (Dransfield, 2008) the `local` grid in `loc_file` to the `regional` grid in `reg_file`.
    
    If `loc_file` and `reg_file` both equal the string "synth", then a synthetic model is
    constructed and conformed. This is for testing.
    
    Parameters
    ----------
    loc_file : Path

        The path to the grid file containing the `local` grid (.ers or .nc).
    
    reg_file : Path

        The path to the grid file containing the `regional` grid (.ers or .nc).

    loc_units : String

        The units ("mGal" or "ums2") of the local gravity data.
        
    reg_units : String

        The units ("mGal" or "ums2") of the regional gravity data.
    
    survey_polygon : List

        The polygon vertices of the survey boundary, as an array or list of (x,y)
        pairs, in either clockwise or counter-clockwise order around the boundary.
        For example, survey_polygon = [(0, 0), (1, 0), (1,1), (0,1)]. Final output will be
        trimmed to this polygon.
    
    rim : Int, optional

        The width in number of pixels for a narrow strip around the `local` grid which will
        be filled by interpolation to the padding area. A value between 8 and 16 is recommended.
        Default 16.
        
    low_lambda : Float, optional

        The short wavelength defining the cosine squared conforming filter. Defaults to a value
        depending on the size of the local grid.
    
    hi_lambda : Float, optional

        The long wavelength defining the cosine squared conforming filter. Defaults to 1.2 times
        `low_lambda`.
    
    pad_cells : Int, optional

        The number of grid cells to pad the `local` data for the Fourier transform. Defaults to
        a value dependent on `low_lambda`.
    
    plot_flag : Bool, optional

        If True, images will be displayed of the grids at each stage of the processing. If False,
        then only the `local` input and output grids will be displayed. Default False.
    

    RETURNS
    ----------
    conformed_local : xarray DataArray
    
        The `local` grid conformed to the long wavelengths of the `regional` grid. Output in units of um/s/s.
        
    NOTES
    -----
    ToDo:
    
    - allow survey_polygon to be None, and default to smallest enclosing cardinal rectangle;
    - use convex hull instead of rectangle at 'pad out local with NaNs' stage;
    - investigate the use of windows in the FFTs, and optimise;
    
    """
    stoperror = False

    # setup for the synthetic model
    if "synth" in str(loc_file) or "synth" in str(reg_file):
        local = harmonicSurvey(size=6000.)
        regional = harmonicSurvey(size=50000.)
        local.attrs['units'] = "um/s/s"
        regional.attrs['units'] = "um/s/s"
        if low_lambda == None or hi_lambda == None:
            low_lambda = 6000.
            hi_lambda = 9000.
        if pad_cells == None:
            cell_size = round((local.x.data[1] - local.x.data[0]) / 2.0, 1)
            pad_cells = int(np.round(low_lambda / cell_size, 0))
    
    # setup for data read from file
    else:
        local, _ = gridfile_to_xa(whizzFile=loc_file)
        regional, _ = gridfile_to_xa(whizzFile=reg_file)
        
        # Local - check units
        if loc_units in ["mGal", "mgal"]:
            local = local * 10.0
            local.attrs['units'] = "um/s/s"
        elif loc_units in ["gu", "ums2", "um/s/s"]:
            local.attrs['units'] = "um/s/s"
        else:
            print('ERROR - loc_units not one of "mGal", "ums2", "um/s/s"')
            stoperror = True
            
        # Regional - check units
        if reg_units in ["mGal", "mgal"]:
            regional = regional * 10.0
            regional.attrs['units'] = "um/s/s"
        elif reg_units in ["gu", "ums2", "um/s/s"]:
            regional.attrs['units'] = "um/s/s"
        else:
            print('ERROR - reg_units not one of "mGal", "ums2", "um/s/s"')
            stoperror = True
            
        if stoperror:
            return None
        
    return conform(local, regional, survey_polygon, rim, low_lambda, hi_lambda, pad_cells, plot_flag)


def conform(
    local_in, regional, survey_polygon, rim=16, low_lambda=None, hi_lambda=None, 
    pad_cells=None, original_mask=None, plot_flag=False
):
    """
    Conforms (Dransfield, 2008) the `local` grid in `loc_file` to the `regional` grid in `reg_file`.
    
    Parameters
    ----------
    local_in : xarray.core.dataarray.DataArray

        The `local` grid to be conformed.
    
    regional : xarray.core.dataarray.DataArray

        The regional grid used to conform.

    survey_polygon : List

        The polygon vertices of the survey boundary, as an array or list of (x,y)
        pairs, in either clockwise or counter-clockwise order around the boundary.
        For example, survey_polygon = [(0, 0), (1, 0), (1,1), (0,1)]. Final output will be
        trimmed to this polygon.
    
    rim : Int, optional

        The width in number of pixels for a narrow strip around the `local` grid which will
        be filled by interpolation to the padding area. A value between 8 and 16 is recommended.
        Default 16.
        
    low_lambda : Float, optional

        The short wavelength defining the cosine squared conforming filter. Defaults to a value
        depending on the size of the local grid.
    
    hi_lambda : Float, optional

        The long wavelength defining the cosine squared conforming filter. Defaults to 1.2 times
        `low_lambda`.
    
    pad_cells : Int, optional

        The number of grid cells to pad the `local` data for the Fourier transform. Defaults to
        a value dependent on `low_lambda`.

    original_mask : numpy mask array, optional

        Used to mask the final conformed gravity before completion. Default is None, in
        which case, `original_mask` is calculated from `local_in`.
    
    plot_flag : Bool, optional

        If True, images will be displayed of the grids at each stage of the processing. If False,
        then only the `local` input and output grids will be displayed. Default False.
    

    RETURNS
    ----------
    conformed_local : xarray.core.dataarray.DataArray
    
        The `local` grid conformed to the long wavelengths of the `regional` grid. Output in units of um/s/s.
        
    NOTES
    -----
    ToDo:
    
    - allow survey_polygon to be None, and default to smallest enclosing cardinal rectangle;
    - use convex hull instead of rectangle at 'pad out local with NaNs' stage;
    - investigate the use of windows in the FFTs, and optimise;
    
    """
    stoperror = False
    da_extremes = lambda da : f'; range = ({da.min().values:.0f}, {da.max().values:.0f})'
    print('\nConforming the local gravity to the regional.')
            
    # local grid check
    local_in, stoperror = _check_grid(local_in)
    if stoperror:
        print(f'ERROR in local grid.')
        return None

    # regional grid check
    regional, stoperror = _check_grid(regional)
    if stoperror:
        print(f'ERROR in regional grid.')
        return None

    low_lambda, hi_lambda, pad_cells = _set_defaults(low_lambda, hi_lambda, pad_cells, survey_polygon, local_in.attrs['cell_size'])
    print(f'Parameters: low_lambda = {low_lambda:.0f}, hi_lambda = {hi_lambda:.0f}, pad_cells = {pad_cells}')

    if plot_flag:
        xdImage(local_in, f'starting local grid {da_extremes(local_in)}', hs=True)
        xdImage(regional,f'starting regional grid {da_extremes(regional)}', hs=True)
    
    # This mask is so we can mask the final data back to the original extents after conforming
    if original_mask is None:
        original_mask = np.ma.masked_array(local_in, ~np.isnan(local_in)).mask
    if plot_flag:
        # xdImage(original_mask, f'original_mask', clipTo3Std=False)
        plt.imshow(original_mask, origin='lower')
        print(f'orig mask 1 shape = {original_mask.shape}')


    # Filling holes in the original local grid (may want to restrict this to smallest
    # enclosing rectangle or, better, a convex hull)
    local = local_in.copy()
    local.rio.write_crs(4283, inplace=True)
    local.rio.write_nodata(np.nan, encoded=True, inplace=True)
    local = local.rio.interpolate_na(method='nearest')
    if plot_flag:
        xdImage(local,f'local, interpolation-filled {da_extremes(local)}', hs=False)

    # now pad out local with NaNs
    loc_nans = _pad_grid_nans(local, pad_cells)
    if plot_flag:
        xdImage(loc_nans,f'local grid padded with NaNs {da_extremes(loc_nans)}', hs=False)
    
    # make the corresonding mask
    simple_mask = np.ma.masked_array(loc_nans, ~np.isnan(loc_nans))
    z = simple_mask.mask

    # extend the mask by the rim width
    w = rim
    N = z.shape[1]
    M = z.shape[0]
    newz = z.copy()
    for k in range(w):
        for j in range(k,N-k):
            for i in range(k,M-k):
                if z[i+k,j-k] or z[i+k,j] or z[i+k,j+k] or z[i,j-k] or z[i,j] or z[i,j+k] or z[i-k,j-k] or z[i-k,j] or z[i-k,j+k]:
                    newz[i,j] = True
    
    # pad with a linear ramp to 0.0 at the edge
    grid_mean = np.nanmean(local.values)
    padding = pad_cells + w
    dx = local.x.values[1] - local.x.values[0]
    dy = local.y.values[1] - local.y.values[0]
    if (dx - dy) < 0.1:
        dy = dx
    e0 = local.x.values[0] - dx * padding
    n0 = local.y.values[0] - dy * padding
    ew = local.x.values[-1] + dx * padding
    nw = local.y.values[-1] + dy * padding
    nume = round((ew - e0) / dx) + 1
    numn = round((nw - n0) / dy) + 1
    e = np.linspace(e0, ew, nume, endpoint=True)
    n = np.linspace(n0, nw, numn, endpoint=True)
    # ugly, but we need to ensure coordinates align, avoiding rounding error in `linspace`
    e[padding:-padding] = local.x.values
    n[padding:-padding] = local.y.values
    parr = local.pad(
        y=(padding, padding),
        x=(padding, padding),
        keep_attrs=True,
        mode='linear_ramp',
        fill_value=grid_mean,
        stat_length=None)
    parr = parr.assign_coords(y=n, x=e)
    if plot_flag:
        xdImage(parr,f'local, linear-ramp padded {da_extremes(parr)}', hs=True)
    
    # interpolate regional onto local (x,y) positions
    reg_match_pad = _grid_match(parr, regional)

    # (pad original_mask with `pad_cells` NaNs so that it matches the final dataarray shape)
    # original_mask = np.pad(original_mask, pad_cells, mode='constant', constant_values=False)
    # print(f'orig mask 2 shape = {original_mask.shape}')

    # smooth the padded grid in x, then y, directions
    test = parr.rolling(x=w, center=True).mean()
    test1 = test.rolling(y=w, center=True).mean()
    if plot_flag:
        xdImage(test1, f'padded, smoothed local grid {da_extremes(test1)}', hs=True)

    # replace original local grid, and the rim, with NaNs
    masked_test1 = test1[w:-w,w:-w].where(~newz)
    if plot_flag:
        xdImage(masked_test1, f'local grid, and rim removed {da_extremes(masked_test1)}', hs=True)

    # ensure regional covers same area as local
    reg_match = reg_match_pad[w:-w,w:-w]

    # put the original local grid back in place
    masked_test1[pad_cells:-pad_cells,pad_cells:-pad_cells] = local
    if plot_flag:
        xdImage(masked_test1,f'rimmed local grid padded with smoothed linear ramp {da_extremes(masked_test1)}', hs=True)

    # interpolate from original local data through w-rim to smoothed padding
    masked_test1.rio.write_crs(4283, inplace=True)
    masked_test1.rio.write_nodata(np.nan, encoded=True, inplace=True)
    masked_test1 = masked_test1.rio.interpolate_na(method='cubic')#,'linear', 'nearest', 'cubic'},)
    if plot_flag:
        xdImage(masked_test1, f'local grid padded with smoothed linear ramp (rim filled by interp) {da_extremes(masked_test1)}', hs=True)

    # set the wavenumbers for the filter
    k_low = 2 * np.pi / low_lambda
    k_hi = 2 * np.pi / hi_lambda
    
    # remove 2D linear function from reg_match
    reg_in = xrft.detrend(reg_match, ('x', 'y'), detrend_type='linear').transpose(transpose_coords=False)
    reg_in.attrs = reg_match.attrs
    reg_restore = reg_match - reg_in
    reg_restore.attrs = reg_match.attrs
    if plot_flag:
        xdImage(reg_match, f'regional grid interpolated to local coordinate space {da_extremes(reg_match)}', hs=True)
        xdImage(reg_in, f'regional grid minus linear trend {da_extremes(reg_in)}', hs=True)
        xdImage(reg_restore, f'linear trend removed from regional {da_extremes(reg_restore)}', hs=True)
    
    # low-pass reg_match (CHECK WINDOW!!)
    reg_fft = xrft.dft(reg_in.fillna(0.0), detrend=None, true_phase=True, true_amplitude=True)
    kx, ky = np.meshgrid(reg_fft.freq_x, reg_fft.freq_y)
    k = np.sqrt(kx * kx + ky * ky)
    print(f'Regional grid Wavenumber resolution = {k[0,1] - k[0,0]:.3g}')
    print(f'              Equivalent wavelength = {2.0 * np.pi / (k[0,1] - k[0,0]):.0f} m.')
    k_zero = max(k_low, k_hi)
    k_unity = min(k_low, k_hi)
    print(f'Filter 2dB = {(k_zero + k_unity) / 2.:.3g} per m.')
    print(f'Filter 2dB = {4. * np.pi / (k_zero + k_unity):.3g} m.')
    print(f'Low pass filtering at {k_unity:.3g}, {k_zero:.3g} per m.')
    lowfilter = cos_square_radial_filter(k, k_zero, k_unity)
    low_fft = reg_fft * lowfilter
    reg_out = xrft.idft(low_fft, window=True, true_phase=True, true_amplitude=True)
    reg_out.attrs = reg_match.attrs
    
    # restore 2D linear function to reg_match
    low_grav = reg_out.real# + reg_restore.real
    low_grav.values = low_grav.values + reg_restore.values
    low_grav.attrs = reg_match.attrs
    if plot_flag:
        xdImage(reg_out.real, f'low-pass filtered, detrended regional {da_extremes(reg_out.real)}', hs=True)
        xdImage(low_grav, f'low-pass filtered regional (trend restored) {da_extremes(low_grav)}', hs=True)
    
    # high-pass filter the prepared local data (CHECK WINDOW!!)
    loc_fft = xrft.fft(masked_test1.fillna(0.0), detrend='linear', window='hamming', true_phase=True, true_amplitude=True)#'hamming'
    kx, ky = np.meshgrid(loc_fft.freq_x, loc_fft.freq_y)
    k = np.sqrt(kx * kx + ky * ky)
    print(f'Local grid Wavenumber resolution = {k[0,0] - k[0,1]:.3g}')
    print(f'           Equivalent wavelength = {2.0 * np.pi / (k[0,0] - k[0,1]):.0f} m.')
    k_zero = min(k_low, k_hi)
    k_unity = max(k_low, k_hi)
    print(f'High pass filtering at {k_unity:.3g}, {k_zero:.3g} per m.')
    highfilter = cos_square_radial_filter(k, k_zero, k_unity)
    hi_fft = loc_fft * highfilter
    hi_grav = xrft.ifft(hi_fft, window=None, true_phase=True, true_amplitude=True)
    hi_grav.attrs = local_in.attrs
    if plot_flag:
        xdImage(hi_grav.real, f'high-pass filtered local grid {da_extremes(hi_grav.real)}', hs=True)
    
    # sum hi- and lo-pass filtered data
    sum_grav = low_grav.real# + reg_restore.real
    sum_grav.values = (sum_grav.values + hi_grav.real.values)
    if plot_flag:
        xdImage(local_in, f'local_in {da_extremes(local_in)}', clipTo3Std=False)
        xdImage(sum_grav, f'sum_grav {da_extremes(sum_grav)}', clipTo3Std=False)
        print(f'local_in shape = {local_in.shape}')
        print(f'sum_grav shape = {sum_grav.shape}')
    
    sum_clipped = sum_grav.sel(x=local_in.x, y=local_in.y, method="nearest")
    sum_clipped.attrs = local_in.attrs
    if plot_flag:
        report_gridStats(sum_clipped)
        xdImage(sum_clipped, f'sum_clipped {da_extremes(sum_clipped)}', clipTo3Std=False)


    # mask out data except where original data existed
    if plot_flag:
        plt.imshow(original_mask, origin='lower')
        print(f'orig mask 3 shape = {original_mask.shape}')
    sum_masked = sum_clipped.where(original_mask)
    if plot_flag:
        xdImage(sum_masked, f'sum_masked 2 {da_extremes(sum_masked)}', clipTo3Std=False)




    # mask back to the given survey_mask, and trim to smallest enclosing rectangle
    geometries = [
        {
            'type': 'Polygon',
            'coordinates': [survey_polygon]
        }
    ]
    sum_masked.rio.write_crs(4283, inplace=True)
    sum_masked.rio.write_nodata(np.nan, encoded=True, inplace=True)
    sum_clipped = _trim_rectangle(sum_masked.real.rio.clip(geometries))
    sum_clipped.attrs = local_in.attrs

    # show before and after grids
    if plot_flag:
        report_gridStats(sum_clipped)
        myfig, myax = plt.subplots(1,2, figsize=(12,6))

        local_in.plot(ax=myax[0])
        myax[0].set_xlabel(local_in.attrs['x_channel'])
        myax[0].set_ylabel(local_in.attrs['y_channel'])
        myax[0].set_title(f'starting local grid {da_extremes(local_in)}')

        sum_clipped.plot(ax=myax[1])
        myax[1].set_xlabel(sum_clipped.attrs['x_channel'])
        myax[1].set_ylabel(sum_clipped.attrs['y_channel'])
        myax[1].set_title(f'final conformed local grid {da_extremes(sum_clipped)}')
    return sum_clipped


def harmonicSurvey(size=6000.0, lambdas=None):
    """
    Creates an airborne survey grid from the sum of several sinusoidal functions.
    One-off function for testing conforming code on synthetic data.
    
    Parameters
    ----------
    size : Float, optional
    
        Length of each side of the square synthetic survey grid area. Default 6000.0.
        
    lambdas : List of Float, optional
    
        The wavelengths of the sinusoidals used to construct the model. Default None
        in which case, these wavelengths are generated randomly.

    Returns
    -------
    g : xarray DataArray

        The synthetic survey data.

    """
    cellsize = 25.0
    if lambdas == None:            
        lambdas = (5.0 * size - 10.0 * cellsize) * np.random.rand(20,) + 10.0 * cellsize
    # survey parameters - note that the loop order below is important.
    xstart = -size
    xend = size
    ystart = -size
    yend = size
    nx = round((xend - xstart) / cellsize)
    ny = round((yend - ystart) / cellsize)
    x = np.linspace(xstart, xend, nx)
    y = np.linspace(ystart, yend, ny)
    da = xr.DataArray(
        data=np.zeros((nx,ny)),
        dims=["y", "x"],
        coords={"y": y, "x": x},
    )
    X, Y = np.meshgrid(x, y)

    # Calculate the 2D sinusoid
    for wavelength in lambdas:
        k_x = 2 * np.pi / wavelength  # Wave number in x
        k_y = 2 * np.pi / wavelength
        A = 10.0
        phi_x = np.pi * np.random.rand()
        phi_y = np.pi * np.random.rand()
        da += A * np.sin(k_x * X + phi_x) * np.sin(k_y * Y + phi_y)

    return da


def _check_grid(grid):
    """
    Checks that the gravity `grid` has square grid cells, and units of um/s/s. Designed
    only for use by `conform`.
    """
    stoperror = False

    # check for square grid cells
    cell_size = grid.x.data[1] - grid.x.data[0]
    tst_cell_size = grid.y.data[1] - grid.y.data[0]
    if cell_size != tst_cell_size:
        print(f'ERROR grid does not have square cells: {cell_size} != {(grid.y.data[1] - grid.y.data[0])}')
        stoperror = True
    else:
        grid.attrs['cell_size'] = cell_size
        
    # local: check for valid units
    if not 'units' in grid.attrs.keys():
        print(f'ERROR grid does not have units.')
        stoperror = True
    elif not grid.attrs['units'] in ["gu", "ums2", "um/s/s"]:
        print(f'ERROR grid units {grid.attrs["units"]} not in ["gu", "ums2", "um/s/s"].')
        stoperror = True
    else:
        grid.attrs['units'] = "um/s/s"
        
    return grid, stoperror


def _set_defaults(low_lambda, hi_lambda, pad_cells, survey_polygon, cell_size):
    """
    Sets sensible defaults for filter wavelengths, and number of cells for grid
    padding, based on the survey size and the grid cell size.
    
    Designed only as a function for `conform`.
    """
    
    if (low_lambda is None) and (not hi_lambda is None):
        low_lambda = 0.8 * hi_lambda
    elif not (low_lambda is None) and (hi_lambda is None):
        hi_lambda = 1.2 * low_lambda
    elif (low_lambda is None) and (hi_lambda is None):
        minx = 1e10
        maxx = -1.e10
        miny = 1e10
        maxy = -1.e10
        for i in range(len(survey_polygon)):
            minx = min(minx, survey_polygon[i][0])
            maxx = max(maxx, survey_polygon[i][0])
            miny = min(miny, survey_polygon[i][1])
            maxy = max(maxy, survey_polygon[i][1])
        low_lambda = min((maxx - minx), (maxy - miny)) / 2.0
        hi_lambda = 1.2 * low_lambda
        
    if pad_cells is None:
        pad_cells = int(np.round(low_lambda / cell_size, 0))
    
    return low_lambda, hi_lambda, pad_cells


