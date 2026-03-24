#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a cos squared, harmonic-domain, 2D filter.

Author: Mark Helm Dransfield

Created: Nov 2025

License: CC BY-SA
"""

import numpy as np


def cos_square_radial_filter(k, k_zero, k_unity):
    """
    The cosine squared, radially symmetric filter. High-pass if k_zero < k_unity
    and low-pass if k_zero > k_unity.
    
    PARAMETERS
    ----------
    
    k : numpy 2D array
    
        radial wavenumbers
        
    k_zero :
    
        zero-gain wavenumber cut-off
        
    k_unity :
    
        unity-gain wavenumber cut-off
        
    RETURNS
    -------
    
    filter : numpy 2D array
    
        the filter to be used (same shape as input k)
        
    """
    cos2_filter = np.zeros(k.shape)
    
    # low-pass filter
    if k_zero > k_unity:
        cos2_filter[np.where(k < k_unity)] = 1.0
        cos2_filter[np.where(k > k_zero)] = 0.0
        cond1 = k >= k_unity
        cond2 = k <= k_zero
        cond = np.logical_and(cond1, cond2)
        cos2_filter[cond] = np.cos(np.pi / 2.0 * (k[cond] - k_unity) / (k_zero - k_unity)) ** 2.0
    elif k_zero < k_unity:
        cos2_filter[np.where(k > k_unity)] = 1.0
        cos2_filter[np.where(k < k_zero)] = 0.0
        cond1 = k <= k_unity
        cond2 = k >= k_zero
        cond = np.logical_and(cond1, cond2)
        cos2_filter[cond] = np.cos(np.pi / 2.0 * (k[cond] - k_unity) / (k_zero - k_unity)) ** 2.0
    else:
        print(f'ERROR - the filter wavenumbers are equal!')
                    
    return cos2_filter
