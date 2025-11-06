import numpy as np
from pathlib import Path
import xarray as xr
import h5py

from pegasusQC.whizzFiles.updateLineSpacing import updateLineSpacing
import pegasusQC.config as config

groupName = config.groupName


def _calc_gridcell_size(whizzFile):
    """
    Returns a sensible grid cell size for gridding line-based data in xarray
    DataSet `ds`. Usefully provides a default for gridding functions.
    
    The algorithm is described here: "return 500.0"
    Later we will use the "TraverseLineSpacing", if available in the whizzfile.
    
    Parameters
    ----------
    whizzFile : Path or String
        The Path to, or String name of, the whizz file in HDF5 format.

    Returns
    -------
    cell_size : float
        the recommended cell size for gridding.
        
    """

    filename = str(whizzFile)
    with h5py.File(filename, 'r+') as f:
        if 'TraverseSpacing' in f[groupName]['CoordinateFrame'].attrs:
            trav_spacing = f[groupName]['CoordinateFrame'].attrs['TraverseSpacing']
            cell_size = _nice_number(trav_spacing / 4.0)

        else:
            updateLineSpacing(whizzFile)
            if 'TraverseSpacing' in f[groupName]['CoordinateFrame'].attrs:
                trav_spacing = f[groupName]['CoordinateFrame'].attrs['TraverseSpacing']
                cell_size = _nice_number(trav_spacing / 4.0)
            else:
                print('ERROR - could not calculate cell size for gridding.')
                print('    TraverseSpacing not set in whizzFile, and value not estimatable.')
                print('    Recommend running updateLineSpacing() with explicit values.')
                print('    And/or re-running craig_transform with explicit cell_size.')
                return None
    print(f'\nGridding data with cell size = {cell_size}.')
    return cell_size
    

def _nice_number(number):
    """
    TBD
    """
    return number