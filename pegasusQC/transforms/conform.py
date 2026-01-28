#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Conform survey grid to regional grid.

Author: Mark Helm Dransfield

Created: Nov 2025

License: CC BY-SA
"""

import numpy as np
import xrft
import matplotlib.pyplot as plt

from pegasusQC.transforms.cos2_filter import cos_square_radial_filter
import pegasusQC.gridFiles.gridutility as gut
from pegasusQC.transforms._grid_match import _grid_match
from pegasusQC.gridFiles.xdImage import xdImage

def conform(gD_grid, regional_grid, survey_polygon=None):
    """
    DFT regional
    low pass filter

    DFT local
    high pass filter

    add
    """
    # fig = plt.figure()

    # Ensure both grids are on same sampling
    reg_match = _grid_match(gD_grid, regional_grid)

    # Find the coverage for each gridded data set
    x_chan='x'
    y_chan='y'
    local_region = [
        np.min(gD_grid[x_chan].values),
        np.max(gD_grid[x_chan].values),
        np.min(gD_grid[y_chan].values),
        np.max(gD_grid[y_chan].values)
        ]
    reg_region = [
        np.min(reg_match[x_chan].values),
        np.max(reg_match[x_chan].values),
        np.min(reg_match[y_chan].values),
        np.max(reg_match[y_chan].values)
        ]

    # We are only interested in the data within the intersection of the regions.
    intersectregion = [
        (max(local_region[0], reg_region[0]), max(local_region[2], reg_region[2])),
        (min(local_region[1], reg_region[1]), max(local_region[2], reg_region[2])),
        (min(local_region[1], reg_region[1]), min(local_region[3], reg_region[3])),
        (max(local_region[0], reg_region[0]), min(local_region[3], reg_region[3])),
        (max(local_region[0], reg_region[0]), max(local_region[2], reg_region[2]))
        ]

    gD_grid = gut.maskGridByPolygon(gD_grid, intersectregion, x_chan=x_chan, y_chan=y_chan)
    reg_match = gut.maskGridByPolygon(reg_match, intersectregion, x_chan=x_chan, y_chan=y_chan)

    # low-pass
    reg_fft = xrft.dft(reg_match.fillna(0.0), detrend='constant', true_phase=True, true_amplitude=True)
    kx, ky = np.meshgrid(reg_fft.freq_x, reg_fft.freq_y)
    k = np.sqrt(kx * kx + ky * ky)
    print(f'Low pass Wavenumber resolution = {k[0,1] - k[0,0]:.3g}')
    print(f'Equivalent wavelength = {2.0 * np.pi / (k[0,1] - k[0,0]):.3g} m.')
    k_zero = 2 * np.pi / 40000.0
    k_unity = 2 * np.pi / 30000.0
    lowfilter = cos_square_radial_filter(k, k_zero, k_unity)
    low_fft = reg_fft * lowfilter
    low_grav = xrft.idft(low_fft * 2.0 * np.pi, detrend='linear', window=True, true_phase=True, true_amplitude=True)

    # high-pass
    loc_fft = xrft.dft(gD_grid.fillna(0.0), detrend='linear', window=True, true_phase=True, true_amplitude=True)
    kx, ky = np.meshgrid(loc_fft.freq_x, loc_fft.freq_y)
    k = np.sqrt(kx * kx + ky * ky)
    print(f'Low pass Wavenumber resolution = {k[0,1] - k[0,0]:.3g}')
    print(f'Equivalent wavelength = {2.0 * np.pi / (k[0,1] - k[0,0]):.3g} m.')
    k_zero = 2 * np.pi / 30000.0
    k_unity = 2 * np.pi / 40000.0
    highfilter = cos_square_radial_filter(k, k_zero, k_unity)
    hi_fft = loc_fft * highfilter
    hi_grav = xrft.idft(hi_fft * 2.0 * np.pi, detrend='linear', window=True, true_phase=True, true_amplitude=True)

    # print('lo ', low_grav)
    # print('hi ', hi_grav)
    
    # ax = fig.add_subplot(3,3,3)
    # ax.imshow(reg_fft.real)
    # ax.set_title("reg_fft")
    
    # ax = fig.add_subplot(3,3,4)
    # ax.imshow(low_fft.real)
    # ax.set_title("low_fft")
    
    sum_grav = low_grav.real + hi_grav.real
    # ax = fig.add_subplot(3,3,9)
    # ax.imshow(sum_grav)
    # ax.set_title("Sum")

    # ... mask to survey boundary polygon, ...
    if not survey_polygon is None:
        x_chan = 'x'
        y_chan = 'y'
        sum_grav = gut.maskGridByPolygon(sum_grav, survey_polygon, x_chan=x_chan, y_chan=y_chan)

    return sum_grav