from scipy.interpolate import RegularGridInterpolator
import numpy as np
import xarray as xr


def _grid_match(local, regional):
    """
    Interpolates the data in the regular, rectangular grid `regional`
    onto the (x,y) locations of the regular, rectangular grid `local`.
    expanding it by (2 x `pad_cells`) in length and width. The additional
    cells all have data value NaN. The grid is also mean-corrected.
   
    Parameters
    ----------
    local : xarray 2D DataArray
        Grid containing the reference locations (x,y).
    regional : xarray 2D DataArray
        Grid to be interpolated onto (x,y).

    Returns
    -------
    xarray 2D DataArray : Interpolated, mean-corrected grid.
        
    """
    interp = RegularGridInterpolator((regional.y.values,regional.x.values),
                                     regional.data,
                                     bounds_error=False,
                                     fill_value=None)
    X,Y = np.meshgrid(local.x.values, local.y.values, indexing='xy')
    reg_matched = xr.DataArray(
        data=(interp((Y,X))),
        dims=["y", "x"],
        coords={"y": local.y.values, "x": local.x.values})
    return reg_matched - reg_matched.mean()
    