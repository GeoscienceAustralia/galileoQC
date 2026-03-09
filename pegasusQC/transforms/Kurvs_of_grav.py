#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates the differential curvatures of gravity.

Author: Mark Helm Dransfield

Created: Oct 2025

License: CC BY-SA
"""

from scipy.fftpack import fft2
import xrft
import numpy as np


def Kurvs_of_grav(gravity, firstorder=False, scale=1000.0):
    """
    Calculates the differential curvatures of gravity, via
    FFT.
   
    Parameters
    ----------
    gravity : xarray 2D DataArray

        Grid to be transformed.

    firstorder : bool, optional

        If True, include first order Craig correction. Default False.

    scale : float

        After the transform is complete, the output grids are
        multiplied by scale. Usually, the input is gravity in
        mGal or um/s/s and the desired outputs are gradients
        in E. For mGal, scale = 10000.0; for E, 1000.0.

    Returns
    -------
    gne, guv : list(xarray 2D DataArray)

        The differential curvatures of gravity.
        
    """
    if gravity.attrs['units'] == "mGal":
        scale = 10000.0
    else:
        scale = 1000.0

     # Take DFT of vertical gravity.
    grav_fft = xrft.dft(gravity, detrend='linear', window=True, true_phase=True, true_amplitude=True)
    
    # Next, form the wavenumbers ...
    kx, ky = np.meshgrid(grav_fft.freq_x, grav_fft.freq_y)
    k = np.sqrt(kx * kx + ky * ky)
    print(f'Wavenumber resolution = {kx[0,1] - kx[0,0]:.3g}')
    print(f'Equivalent wavelength = {2.0 * np.pi / (kx[0,1] - kx[0,0]):.3g} m.')
    
    # ..., find the indices of the wavenumber origin ... 
    ky_null = np.nonzero(grav_fft.freq_y.values == 0)[0][0]
    kx_null = np.nonzero(grav_fft.freq_x.values == 0)[0][0]
    
    # ..., and set k so that the kernel will be zero there.
    k[ky_null, kx_null] = 1.0E-9

    #
    potl_fft = grav_fft / k
    potl_fft[ky_null, kx_null] = 0.0
    
    # Form and apply the Craig transform kernel ...
    gne_fft = kx * ky * potl_fft
    guv_fft = (kx * kx - ky * ky) * potl_fft / 2.0

    # gne_fft = kx * ky / k
    if firstorder:
        tiltangle_deg = 2.0
        root8_theta = tiltangle_deg / 180.0 * np.pi / np.sqrt(2.0)
        gne_fft += root8_theta * (ky - kx) * grav_fft
        guv_fft -= root8_theta * (ky - kx) * grav_fft
        print("Making first order tilt correction to transform FROM gravity.")
    
    # Then inverse transform, ...
    gne = xrft.idft(gne_fft * 2.0 * np.pi, detrend='linear', window=True, true_phase=True, true_amplitude=True)
    guv = xrft.idft(guv_fft * 2.0 * np.pi, detrend='linear', window=True, true_phase=True, true_amplitude=True)

    # ..., scale to eotvos
    gne_real = gne.real * scale
    gne_real.attrs['units'] = 'eotvos'
    guv_real = guv.real * scale
    guv_real.attrs['units'] = 'eotvos'
    
    return gne_real, guv_real
