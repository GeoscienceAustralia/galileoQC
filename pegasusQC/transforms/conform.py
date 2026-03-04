import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import xrft
from pathlib import Path

from pegasusQC.transforms.cos2_filter import cos_square_radial_filter
from pegasusQC.transforms._pad_grid_nans import _pad_grid_nans
from pegasusQC.transforms._trim_rectangle import _trim_rectangle
from pegasusQC.transforms._grid_match import _grid_match
from pegasusQC.gridFiles.grid_to_xarray import gridfile_to_xa
from pegasusQC.gridFiles.xdImage import xdImage

def conform(loc_file, reg_file, loc_units, reg_units, pad_cells, mask_polygon, rim, low_lambda, hi_lambda, plot_flag=False):
    """
    Conforms (Dransfield, 2008) the `local` grid in `loc_file` to the `regional` grid in `reg_file`.
    
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
        
    pad_cells : Int

        The number of grid cells to pad the `local` data for the Fourier transform.
    
    mask_polygon : List

        The polygon vertices of the survey boundary, as an array or list of (x,y)
        pairs, in either clockwise or counter-clockwise order around the boundary.
        For example, mask_polygon = [(0, 0), (1, 0), (1,1), (0,1)]. Final output will be
        trimmed to this polygon.
    
    rim : Int

        The width in number of pixels for a narrow strip around the `local` grid which will
        be filled by interpolation to the padding area. A value between 8 and 16 is recommended.
    
    low_lambda : Float

        The short wavelength defining the cosine squared conforming filter.
    
    hi_lambda : Float

        The long wavelength defining the cosine squared conforming filter.
    
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
    
    - allow mask_polygon to be None, and default to smallest enclosing cardinal rectangle;
    - use convex hull instead of rectangle at 'pad out local with NaNs' stage;
    - investigate the use of windows in the FFTs, and optimise;
    
    """
    stoperror = False
    da_extremes = lambda da : f'; range = ({da.min().values:.0f}, {da.max().values:.0f})'

    # setup for the synthetic model
    if "synth" in str(loc_file) or "synth" in str(reg_file):
        local_in = harmonicSurvey(size=6000.)
        regional = harmonicSurvey(size=50000.)
        local_in.attrs['units'] = "gu"
        regional.attrs['units'] = "gu"
    
    # setup for data read from file
    else:
        local_in, _ = gridfile_to_xa(whizzFile=loc_file)
        regional, _ = gridfile_to_xa(whizzFile=reg_file)
        
        # check that valid units were provided
        if loc_units in ["mGal", "mgal"]:
            local_in = local_in * 10.0
            local_in.attrs['units'] = "gu"
        elif loc_units in ["gu", "ums2"]:
            local_in.attrs['units'] = "gu"
        else:
            print('ERROR - loc_units not one of "mGal", "ums2"')
            stoperror = True
            
        if reg_units in ["mGal", "mgal"]:
            regional = regional * 10.0
            regional.attrs['units'] = "gu"
        elif reg_units in ["gu", "ums2"]:
            regional.attrs['units'] = "gu"
        else:
            print('ERROR - reg_units not one of "mGal", "ums2"')
            stoperror = True
            
        if stoperror:
            return None
    
    # check for square grid cells in regional
    reg_cell_size = round((regional.x.data[1] - regional.x.data[0]) / 2.0, 1)
    tst_cell_size = round((regional.y.data[1] - regional.y.data[0]) / 2.0, 1)
    if reg_cell_size != tst_cell_size:
        print(f'ERROR regional grid does not have square cells: {reg_cell_size} != {(regional.y.data[1] - regional.y.data[0]) / 2.0}')
        return None

    # check for square grid cells in local
    loc_cell_size = round((local_in.x.data[1] - local_in.x.data[0]) / 2.0, 1)
    tst_cell_size = round((local_in.y.data[1] - local_in.y.data[0]) / 2.0, 1)
    if loc_cell_size != tst_cell_size:
        print(f'ERROR local grid does not have square cells: {loc_cell_size} != {(local_in.y.data[1] - local_in.y.data[0]) / 2.0}')
        return None

    if plot_flag:
        xdImage(local_in, f'starting local grid {da_extremes(local_in)}', hs=True)
        xdImage(regional,f'starting regional grid {da_extremes(regional)}', hs=True)
    
    # This mask is so we can mask the final data back to the original extents after conforming
    original_mask = np.ma.masked_array(local_in, ~np.isnan(local_in)).mask

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
    parr = local.pad(x=(padding, padding),
                     y=(padding, padding), 
                        keep_attrs=True,
                        mode='linear_ramp',
                        fill_value=grid_mean,
                        stat_length=None)
    dx = local.x.values[1] - local.x.values[0]
    dy = local.y.values[1] - local.y.values[0]
    for i in range(0,padding):
        parr.x.values[i] = parr.x.values[padding] - (padding - i) * dx
        parr.x.values[-1-i] = parr.x.values[-1-padding] + (padding - i) * dx
        parr.y.values[i] = parr.y.values[padding] - (padding - i) * dy
        parr.y.values[-1-i] = parr.y.values[-1-padding] + (padding - i) * dy
    if plot_flag:
        xdImage(parr,f'local, linear-ramp padded {da_extremes(parr)}', hs=True)
    
    # interpolate regional onto local (x,y) positions
    reg_match_pad = _grid_match(parr, regional)

    # (pad original_mask with `pad_cells` NaNs so that it matches the final dataarray shape)
    original_mask = np.pad(original_mask, pad_cells, mode='constant', constant_values=False)

    # smooth the padded grid in x, then y, directions
    test = parr.rolling(x=w, center=True).mean()#.dropna('x')
    test1 = test.rolling(y=w, center=True).mean()#.dropna('x')
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
    print(f'Equivalent wavelength = {2.0 * np.pi / (k[0,1] - k[0,0]):.3g} m.')
    k_zero = max(k_low, k_hi)
    k_unity = min(k_low, k_hi)
    if plot_flag:
        f2,a2 = plt.subplots(1,2,figsize=(6,3))
        a2[0].plot([k_zero, k_unity], [0., 1.])
        a2[0].set_title('low-pass')
        a2[0].set_xlabel('wavelength [m]')
    print(f'Low pass filter 2dB = {(k_zero + k_unity) / 2.:.3g} per m.')
    print(f'Low pass filter 2dB = {4. * np.pi / (k_zero + k_unity):.3g} m.')
    print(f'Low pass filtering at {k_unity}, {k_zero} per m.')
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
    print(f'Equivalent wavelength = {2.0 * np.pi / (k[0,0] - k[0,1]):.0f} m.')
    k_zero = min(k_low, k_hi)
    k_unity = max(k_low, k_hi)
    if plot_flag:
        a2[1].plot([k_zero, k_unity], [0., 1.])
        a2[1].set_title('high-pass')
        a2[0].set_xlabel('wavelength [m]')
    print(f'High pass filter 2dB = {(k_zero + k_unity) / 2.:.3g} per m.')
    print(f'High pass filter 2dB = {4. * np.pi / (k_zero + k_unity):.3g} m.')
    print(f'High pass filtering at {k_unity}, {k_zero} per m.')
    highfilter = cos_square_radial_filter(k, k_zero, k_unity)
    hi_fft = loc_fft * highfilter
    hi_grav = xrft.ifft(hi_fft, window=None, true_phase=True, true_amplitude=True)
    hi_grav.attrs = local_in.attrs
    if plot_flag:
        xdImage(hi_grav.real, f'high-pass filtered local grid {da_extremes(hi_grav.real)}', hs=True)
    
    # sum hi- and lo-pass filtered data
    sum_grav = low_grav.real# + reg_restore.real
    sum_grav.values = (sum_grav.values + hi_grav.real.values)
    
    # mask out data except where original data existed
    sum_masked = sum_grav.where(original_mask)

    # mask back to the given survey_mask, and trim to smallest enclosing rectangle
    geometries = [
        {
            'type': 'Polygon',
            'coordinates': [mask_polygon]
        }
    ]
    sum_masked.rio.write_crs(4283, inplace=True)
    sum_masked.rio.write_nodata(np.nan, encoded=True, inplace=True)
    sum_clipped = _trim_rectangle(sum_masked.real.rio.clip(geometries))
    sum_clipped.attrs = local_in.attrs

    # show before and after grids
    # with plt.ion():
    myfig, myax = plt.subplots(1,2, figsize=(12,6))
    vmin = min(local_in.min(), sum_clipped.min()).values
    vmax = max(local_in.max(), sum_clipped.max()).values
    xdImage(local_in, f'starting local grid {da_extremes(local_in)}', clipTo3Std=False, ax=myax[0], minClip=vmin, maxClip=vmax)
    xdImage(sum_clipped, f'final conformed local grid {da_extremes(sum_clipped)}', clipTo3Std=False, ax=myax[1], minClip=vmin, maxClip=vmax)
    
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

