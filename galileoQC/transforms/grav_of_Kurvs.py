#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Apply the Craig Transform to (`Guv` + i `Gne`).

Author: Mark Helm Dransfield

Created: Oct 2025

License: CC BY-SA
"""

import numpy as np
import xrft
import matplotlib.path as mpltPath

from galileoQC.gridFiles.xdImage import xdImage
from galileoQC.gridFiles.gridutility import report_gridStats
from galileoQC.transforms._trim_rectangle import _trim_rectangle
import galileoQC.gridFiles.gridutility as gut


def grav_of_Kurvs(Gne, Guv, firstorder=False, survey_polygon=None, nan_mask=None, plot_flag=False, verbose=False):
    """
    Apply the Craig Transform to (`Guv` + i `Gne`), mask the resulting grid
    to `mask_region` (NOT IMPLEMENTED),
    and return the real part and imaginary part as two grids.

    Calculation can be to zero-th or first order.
    
    - 2D FFT of Craig function
    - multiply by kernel
    - inverse 2D FFT

    Parameters
    ----------
    Gne : xarray 2D DataArray

        Grid of NE component data in eotvos.

    Guv : xarray 2D DataArray

        Grid of UV component data in eotvos.

    firstorder : bool, optional

        If True, include first order Craig correction. Default False.

    survey_polygon : numpy 1D array, optional

        The polygon vertices of the survey boundary, as an array or sequence of (x,y)
        pairs, in either clockwise or counter-clockwise order around the boundary.
        For example, survey_polygon = [(0, 0), (1, 0), (1,1), (0,1)]. Final output will be
        trimmed to this polygon if provided. Default None.

    nan_mask : numpy 2D mask array, optional

        blah blah. Default None.
    
    plot_flag : Bool, optional

        If True, images will be displayed of the grids at each stage of the processing. If False,
        then only the `local` input and output grids will be displayed. Default False.

    Returns
    -------
    gD_result : xarray 2D DataArray

        the calculated (real part of the) down component of gravity

    gD_err : xarray 2D DataArray

        the calculated imaginary part of the down component of gravity
        
    """
    
    # Form complex field ...
    kurv = Guv + 1j * Gne
    
    # ..., then take DFT of trimmed complex differential curvature.
    kurv_fft = xrft.dft(kurv, detrend='linear', window=True, true_phase=True, true_amplitude=True)
    
    # Next, form the wavenumbers ...
    kx, ky = np.meshgrid(kurv_fft.freq_x, kurv_fft.freq_y)
    k = np.sqrt(kx * kx + ky * ky)
    print(f'Wavenumber resolution = {k[0,1] - k[0,0]:.3g}')
    print(f'Equivalent wavelength = {2.0 * np.pi / (k[0,1] - k[0,0]):.3g} m.')
    
    # ..., find the indices of the wavenumber origin ... 
    ky_null = np.nonzero(kurv_fft.freq_y.values == 0)[0][0]
    kx_null = np.nonzero(kurv_fft.freq_x.values == 0)[0][0]
    
    # ..., and set k so that the kernel will be zero there.
    repl_zero = 1.0E-9
    kx[ky_null, kx_null] = repl_zero
    ky[ky_null, kx_null] = repl_zero
    k[ky_null, kx_null] = 0.0
    
    # Form and apply the Craig transform kernel ...
    craig = 2.0 * k / ((kx - 1j * ky) ** 2.0)
    if firstorder:
        tiltangle_deg = 2.0
        tiltangle_rad = tiltangle_deg / 180.0 * np.pi
        root2_theta = np.sqrt(2.0) * tiltangle_rad
        craig -= 2.0 * root2_theta * (1 + 1j) * (ky - kx) * (k ** 2.0) / ((kx - 1j * ky) ** 4.0)
        print("Making first order tilt correction to transform FROM gravity.")
    hD = craig * kurv_fft / 2.0 / np.pi
    
    # Then inverse transform, ...
    gD_tr = xrft.idft(hD, detrend='linear', window=True, true_phase=True, true_amplitude=True)
    
    # Separate real and imag parts now so masking works as desired.
    gD_im = gD_tr.imag
    gD_re = gD_tr.real
    if plot_flag:
        xdImage(gD_re, 'unmasked, untrimmed gD_grid (um/s/s)', hs=True)
        if verbose:
            report_gridStats(gD_re)

    # ... replace NaNs, ...
    if not nan_mask is None:
        gD_im = gD_im.where(~nan_mask, other=np.nan)
        gD_re = gD_re.where(~nan_mask, other=np.nan)
        if plot_flag:
            xdImage(gD_re, 'masked, untrimmed gD_grid (um/s/s)', hs=True)
            if verbose:
                print('    gD_re')
                report_gridStats(gD_re)
    # ... mask to survey boundary polygon, ...
    if not survey_polygon is None:
        x_chan = 'x'
        y_chan = 'y'
        gD_err = gut.maskGridByPolygon(gD_im, survey_polygon, x_chan=x_chan, y_chan=y_chan)
        x_chan = 'x'
        y_chan = 'y'
        gD_result = gut.maskGridByPolygon(gD_re, survey_polygon, x_chan=x_chan, y_chan=y_chan)
        if plot_flag:
            xdImage(gD_result, 'masked, trimmed gD_grid (um/s/s)', hs=True)
            if verbose:
                print(f'    gD_result, {gD_result.shape}')
                report_gridStats(gD_result)
    else:
        gD_err = gD_im
        gD_result = gD_re

    # ...,  and trim to smallest cardinal rectangle.
    gD_result_t = _trim_rectangle(gD_result)
    gD_err_t = _trim_rectangle(gD_err)
    if plot_flag:
        print(f'    gD_result_t, {gD_result_t.shape}')
        print(f'    gD_err_t, {gD_err_t.shape}')

    # We want to keep any infill for gD_raw so it helps conforming
    margin = 0.2 * np.abs(gD_result_t.x[1] - gD_result_t.x[0])
    min_x = gD_result_t.x.min().values - margin
    max_x = gD_result_t.x.max().values + margin
    min_y = gD_result_t.y.min().values - margin
    max_y = gD_result_t.y.max().values + margin
    test_poly = [
        [min_x, min_y],
        [max_x, min_y],
        [max_x, max_y],
        [min_x, max_y],
    ]
    gD_raw = gut.maskGridByPolygon(gD_tr.real, test_poly, x_chan=x_chan, y_chan=y_chan)
    gD_raw_t = _trim_rectangle(gD_raw)
    if plot_flag:
        print(f'    gD_raw, {gD_raw.shape}')
        print(f'    gD_raw_t, {gD_raw_t.shape}')

    # return _trim_rectangle(gD_result), _trim_rectangle(gD_err), _trim_rectangle(gD_raw)
    return gD_result_t, gD_err_t, gD_raw_t
