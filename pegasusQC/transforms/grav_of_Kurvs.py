import numpy as np
from pegasusQC.gridFiles.gridutility import report_gridStats
# import pegasusQC.gridFiles.gridutility as gut
from scipy.fftpack import fft2
import xrft


def grav_of_Kurvs(Gne, Guv, mask_polygon=None, nan_mask=None):
    """
    Apply the Craig Transform to (`Guv` + i `Gne`), mask the resulting grid
    to `mask_region` (NOT IMPLEMENTED),
    and return the real part and imaginary part as two grids.
    
    - 2D FFT of Craig function
    - multiply by kernel
    - inverse 2D FFT

    Parameters
    ----------
    Gne : xarray 2D DataArray
        Grid of NE component data in eotvos.
    Guv : xarray 2D DataArray
        Grid of UV component data in eotvos.
    mask_polygon : numpy 1D array, optional
        In order [min_x, max_x, min_y, max_y]. If the mask_polygon is given,
        then the output arrays will be masked to the area within the polygon
        defined by it. Default None.

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
    print(f'Wavenumber resolution = {kx[0,1] - kx[0,0]}')
    
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
    hD = craig * kurv_fft
    
    # Then inverse transform, ...
    gD_tr = xrft.idft(hD, detrend='linear', window=True, true_phase=True, true_amplitude=True)
    
    # Separate real and imag parts now so masking works as desired.
    gD_im = gD_tr.imag
    gD_re = gD_tr.real

    # ... replace NaNs and remove padding, ...
    if not nan_mask is None:
        gD_im_2 = gD_im.where(~nan_mask, other=np.nan)
        gD_re_2 = gD_re.where(~nan_mask, other=np.nan)
    else:
        gD_im_2 = gD_im
        gD_re_2 = gD_re
    # plt.imshow(gD_tr.real)
    # plt.imshow(gD_tr_2.real)
    if not mask_polygon is None:
        gD_err = gD_im_2.sel(x=slice(mask_polygon[0], mask_polygon[1]), y=slice(mask_polygon[2], mask_polygon[3]))
        gD_result = gD_re_2.sel(x=slice(mask_polygon[0], mask_polygon[1]), y=slice(mask_polygon[2], mask_polygon[3]))
    else:
        gD_err = gD_im_2
        gD_result = gD_re_2
    
    return gD_result, gD_err
