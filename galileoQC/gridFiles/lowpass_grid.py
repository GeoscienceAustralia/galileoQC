"""
Steps:

1. fill holes in grid by nearest-neighbour interpolation
2. pad extension with zeros
3. set wavenumbers from wavelengths
4. remove 2D linear trend
5. FFT
6. calculate wave-vectors
7. apply cos2 filter
8. inverse FFT
9. mask out padding and filled holes
10. trim to smallest enclosing rectange
11. image result
12. report grid statistics
13. return resulting grid

"""
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

from galileoQC.transforms._pad_mean import _pad_mean
from galileoQC.transforms.conform import _set_defaults
from galileoQC.transforms.cos2_filter import cos_square_radial_filter
from galileoQC.transforms._pad_grid_nans import _pad_grid_nans
from galileoQC.transforms._trim_rectangle import _trim_rectangle
from galileoQC.transforms._grid_match import _grid_match
from galileoQC.gridFiles.xdImage import xdImage
from galileoQC.gridFiles.grid_to_xarray import gridfile_to_xa
from galileoQC.gridFiles.gridutility import report_gridStats


def lowpass_grid(agrid, low_lambda, filtername='cos2', plot_flag=False):
    """
    Low-pass filters `agrid`.
    
    Parameters
    ----------
    agrid : 2D xarray Dataset

        The grid to be filtered.
            
    low_lambda : Float

        The 0.5 gain wavelength of the filter.
    
    filtername : String, optional

        The name of the filter to be used. Currently only the default ('cos2')
        is available and any other choice is ignored.
    

    RETURNS
    ----------
    fgrid : 2D xarray DataArray
    
        The filtered grid.
        
    NOTES
    -----
    ToDo:
    
    - add other filter choices;
    
    """
    da_extremes = lambda da : f'; range = ({da.min().values:.0f}, {da.max().values:.0f})'

    # local grid check
    stoperror = False
    agrid, stoperror = _check_grid(agrid)
    if stoperror:
        print(f'ERROR in grid.')
        return None

    low_lambda, _, pad_cells = _set_defaults(low_lambda, None, None, None, agrid.attrs['cell_size'])

    # Filling holes in the original local grid (may want to restrict this to smallest
    # enclosing rectangle or, better, a convex hull)
    local = agrid.copy()
    local.rio.write_crs(4283, inplace=True)
    local.rio.write_nodata(np.nan, encoded=True, inplace=True)
    local = local.rio.interpolate_na(method='nearest')
    if plot_flag:
        xdImage(local,f'local, interpolation-filled {da_extremes(local)}', hs=False)

    # remove 2D linear function from hole-filled local
    loc_detrended = xrft.detrend(local, ('x', 'y'), detrend_type='linear').transpose(transpose_coords=False)
    loc_detrended.attrs = local.attrs
    trend_lin2 = local - loc_detrended
    trend_lin2.attrs = local.attrs
    if plot_flag:
        xdImage(local, f'hole-filled grid {da_extremes(local)}', hs=True)
        xdImage(loc_detrended, f'detrended grid {da_extremes(loc_detrended)}', hs=True)
        xdImage(trend_lin2, f'the linear trend removed from grid {da_extremes(trend_lin2)}', hs=True)

    # now pad out loc_detrended with NaNs
    loc_nans = _pad_grid_nans(loc_detrended, pad_cells)
    if plot_flag:
        xdImage(loc_nans,f'grid padded with NaNs {da_extremes(loc_nans)}', hs=False)
    
    # make the corresonding mask (Consider doing this before filling holes? TODO?)
    simple_mask = np.ma.masked_array(loc_nans, ~np.isnan(loc_nans))
    z = simple_mask.mask

	# now replace padding with mean (Consider a better padding method? TODO?)
    loc_padded, _ = _pad_mean(loc_detrended, pad_cells)

    # set the wavenumbers for the filter
    k_mid = 2 * np.pi / low_lambda
    k_lo = 0.8 * k_mid
    k_hi = 1.2 * k_mid

    
    # low-pass loc_padded (CHECK WINDOW!!)
    local_fft = xrft.dft(loc_padded.fillna(0.0), detrend=None, true_phase=True, true_amplitude=True)
    kx, ky = np.meshgrid(local_fft.freq_x, local_fft.freq_y)
    k = np.sqrt(kx * kx + ky * ky)
    report_gridStats(loc_padded, title='padded grid ready for FFT')

    k_zero = max(k_lo, k_hi)
    k_unity = min(k_lo, k_hi)
    print(f'Filter 2dB = {(k_zero + k_unity) / 2.:.3g} per m.')
    print(f'Filter 2dB = {4. * np.pi / (k_zero + k_unity):.3g} m.')
    print(f'Low pass filtering at {k_unity:.3g}, {k_zero:.3g} per m.')
    lowfilter = cos_square_radial_filter(k, k_zero, k_unity)
    local_low = local_fft * lowfilter
    local_out = xrft.idft(local_low, window=True, true_phase=True, true_amplitude=True)
    local_out.attrs = loc_padded.attrs

    local_clipped = local_out.sel(x=local.x, y=local.y, method="nearest")
    local_clipped.attrs = local.attrs

    # restore 2D linear function to local_out
    low_grav = local_clipped.real
    low_grav.values = low_grav.values + trend_lin2.values
    low_grav.attrs = local_clipped.attrs
    if plot_flag:
        report_gridStats(local_clipped, title='low-pass filtered grid')
        xdImage(local_out.real, f'low-pass filtered, detrended grid {da_extremes(local_out.real)}', hs=True)
        xdImage(low_grav, f'low-pass filtered grid (trend restored) {da_extremes(low_grav)}', hs=True)

    return low_grav


def _check_grid(grid):
    """
    Checks that the gravity `grid` has square grid cells. Designed
    only for use by `lowpass_grid`.
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
                
    return grid, stoperror



