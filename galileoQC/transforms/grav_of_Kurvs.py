#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Apply the Craig Transform to (`Guv` + i `Gne`).

Author: Mark Helm Dransfield

Created: Oct 2025

License: CC BY-SA
"""

import numpy as np
import xarray as xr
import xrft
import matplotlib.path as mpltPath

from galileoQC.gridFiles.xdImage import xdImage
from galileoQC.gridFiles.gridutility import report_gridStats
from galileoQC.transforms._trim_rectangle import _trim_rectangle
import galileoQC.gridFiles.gridutility as gut


def grav_of_Kurvs(Gne, Guv, firstorder=False, survey_polygon=None, nan_mask=None, plot_flag=False, verbose=False):
    """
    Apply the Craig Transform to (`Guv` + i `Gne`), producing gravity and
    the gravity gradients, and return the resulting grids.

    Calculation can be to zero-th or first order.
    
    - 2D FFT of Craig function
    - multiply by kernel
    - inverse 2D FFT
    - mask and trim

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
        
    gD_raw : xarray 2D DataArray

        the calculated of the down component of gravity with filled NaNs retained
        
    g_nn : xarray 2D DataArray

        the transformed NN component of the gravity gradient tensor
        
    g_ne : xarray 2D DataArray

        the transformed NE component of the gravity gradient tensor
        
    g_nd : xarray 2D DataArray

        the transformed ND component of the gravity gradient tensor
        
    g_ee : xarray 2D DataArray

        the transformed EE component of the gravity gradient tensor
        
    g_ed : xarray 2D DataArray

        the transformed ED component of the gravity gradient tensor
        
    g_dd : xarray 2D DataArray

        the transformed DD component of the gravity gradient tensor
        
    """
    
    # Form complex field ...
    kurv = Guv + 1j * Gne
    
    # ..., then take DFT of trimmed complex differential curvature.
    kurv_fft = xrft.dft(kurv, detrend='linear', window=True, true_phase=True, true_amplitude=True)
    
    # Next, form the wavenumbers ...
    ke, kn = np.meshgrid(kurv_fft.freq_x, kurv_fft.freq_y)
    k = np.sqrt(ke * ke + kn * kn)
    print(f'Wavenumber resolution = {k[0,1] - k[0,0]:.3g}')
    print(f'Equivalent wavelength = {2.0 * np.pi / (k[0,1] - k[0,0]):.3g} m.')
    
    # ..., find the indices of the wavenumber origin ... 
    kn_null = np.nonzero(kurv_fft.freq_y.values == 0)[0][0]
    ke_null = np.nonzero(kurv_fft.freq_x.values == 0)[0][0]
    
    # ..., and set k so that the kernel will be zero there.
    repl_zero = 1.0E-9
    ke[kn_null, ke_null] = repl_zero
    kn[kn_null, ke_null] = repl_zero
    k[kn_null, ke_null] = 0.0
    
    # Form and apply the Craig transform kernel ... (**##!! tensor component transforms go here)
    curv_kernel = - 2.0 / ((kn + 1j * ke) ** 2.0)
    craig_d = k * curv_kernel
    craig_nn = - kn * kn * curv_kernel
    craig_ne = - ke * kn * curv_kernel
    craig_nd = 1j * kn * k  * curv_kernel
    craig_ee =  - ke * ke * curv_kernel
    craig_ed = 1j * ke * k * curv_kernel
    craig_dd = (k * k) * curv_kernel

    if firstorder:
        tiltangle_deg = 2.0
        tiltangle_rad = tiltangle_deg / 180.0 * np.pi
        root2_theta = np.sqrt(2.0) * tiltangle_rad
        craig_d -= 2.0 * root2_theta * (1 + 1j) * (kn - ke) * (k ** 2.0) / ((ke - 1j * kn) ** 4.0)
        print("Making first order tilt correction to transform FROM gravity.")
    hD = craig_d * kurv_fft / 2.0 / np.pi
    hnn = craig_nn * kurv_fft / 2.0 / np.pi
    hne = craig_ne * kurv_fft / 2.0 / np.pi
    hnd = craig_nd * kurv_fft / 2.0 / np.pi
    hee = craig_ee * kurv_fft / 2.0 / np.pi
    hed = craig_ed * kurv_fft / 2.0 / np.pi
    hdd = craig_dd * kurv_fft / 2.0 / np.pi
    
    # Then inverse transform, ...
    gD_tr = xrft.idft(hD, detrend='linear', window=True, true_phase=True, true_amplitude=True)
    g_nn_tr = xrft.idft(hnn, detrend='linear', window=True, true_phase=True, true_amplitude=True)
    g_ne_tr = xrft.idft(hne, detrend='linear', window=True, true_phase=True, true_amplitude=True)
    g_nd_tr = xrft.idft(hnd, detrend='linear', window=True, true_phase=True, true_amplitude=True)
    g_ee_tr = xrft.idft(hee, detrend='linear', window=True, true_phase=True, true_amplitude=True)
    g_ed_tr = xrft.idft(hed, detrend='linear', window=True, true_phase=True, true_amplitude=True)
    g_dd_tr = xrft.idft(hdd, detrend='linear', window=True, true_phase=True, true_amplitude=True)
    
    # Separate real and imag parts now so masking works as desired.
    gD_im = gD_tr.imag
    gD_re = gD_tr.real
    if plot_flag:
        ntemp = np.array(range(craig_d.real.shape[0]))
        etemp = np.array(range(craig_d.real.shape[1]))
        craigArray = xr.DataArray(craig_d.real, coords={"N": ntemp, "E": etemp}, dims=["N", "E"])
        xdImage(craigArray, 'Craig gD filter (Fourier real)', hs=False)
        xdImage(gD_re, 'unmasked, untrimmed gD_grid (um/s/s)', hs=True)
        if verbose:
            report_gridStats(gD_re, title='unmasked, untrimmed gD_grid (um/s/s)')

    g_nn_re = g_nn_tr.real
    g_ne_re = g_ne_tr.real
    g_nd_re = g_nd_tr.real
    g_ee_re = g_ee_tr.real
    g_ed_re = g_ed_tr.real
    g_dd_re = g_dd_tr.real
    
    # ... replace NaNs, ...
    if not nan_mask is None:
        gD_im = gD_im.where(~nan_mask, other=np.nan)
        gD_re = gD_re.where(~nan_mask, other=np.nan)
        if plot_flag:
            xdImage(gD_re, 'masked, untrimmed gD_grid (um/s/s)', hs=True)
            if verbose:
                report_gridStats(gD_re, title='masked, untrimmed gD_grid (um/s/s)')

        g_nn_re = g_nn_re.where(~nan_mask, other=np.nan)
        g_ne_re = g_ne_re.where(~nan_mask, other=np.nan)
        g_nd_re = g_nd_re.where(~nan_mask, other=np.nan)
        g_ee_re = g_ee_re.where(~nan_mask, other=np.nan)
        g_ed_re = g_ed_re.where(~nan_mask, other=np.nan)
        g_dd_re = g_dd_re.where(~nan_mask, other=np.nan)
        
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
                report_gridStats(gD_result, title='masked, trimmed gD_grid (um/s/s)')

        x_chan = 'x'
        y_chan = 'y'
        g_nn_result = gut.maskGridByPolygon(g_nn_re, survey_polygon, x_chan=x_chan, y_chan=y_chan)
        g_ne_result = gut.maskGridByPolygon(g_ne_re, survey_polygon, x_chan=x_chan, y_chan=y_chan)
        g_nd_result = gut.maskGridByPolygon(g_nd_re, survey_polygon, x_chan=x_chan, y_chan=y_chan)
        g_ee_result = gut.maskGridByPolygon(g_ee_re, survey_polygon, x_chan=x_chan, y_chan=y_chan)
        g_ed_result = gut.maskGridByPolygon(g_ed_re, survey_polygon, x_chan=x_chan, y_chan=y_chan)
        g_dd_result = gut.maskGridByPolygon(g_dd_re, survey_polygon, x_chan=x_chan, y_chan=y_chan)
    else:
        gD_err = gD_im
        gD_result = gD_re

        g_nn_result = g_nn_re
        g_ne_result = g_ne_re
        g_nd_result = g_nd_re
        g_ee_result = g_ee_re
        g_ed_result = g_ed_re
        g_dd_result = g_dd_re

    # ...,  and trim to smallest cardinal rectangle.
    gD_result_t = _trim_rectangle(gD_result)
    gD_err_t = _trim_rectangle(gD_err)
    if plot_flag:
        print(f'    gD_result_t, {gD_result_t.shape}')
        print(f'    gD_err_t, {gD_err_t.shape}')

    g_nn = _trim_rectangle(g_nn_result)
    g_ne = _trim_rectangle(g_ne_result)
    g_nd = _trim_rectangle(g_nd_result)
    g_ee = _trim_rectangle(g_ee_result)
    g_ed = _trim_rectangle(g_ed_result)
    g_dd = _trim_rectangle(g_dd_result)

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
    return gD_result_t, gD_err_t, gD_raw_t, g_nn, g_ne, g_nd, g_ee, g_ed, g_dd
