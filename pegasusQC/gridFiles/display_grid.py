import numpy as np
import pegasusQC.config as config
from pegasusQC.gridFiles.xdImage import xdImage
from pegasusQC.gridFiles.gridfiles import gridfile_to_xa

    
def display_grid(gridFile, mytitle, colormap=config.qc_colormap, cmap_norm='nonorm', 
    minClip=np.nan, maxClip=np.nan, cb_ticks='stats', nSigma=2,
    hs=True, azdeg=45, ax=None, clipTo3Std = True, whizzfile=None, e_chan='', n_chan=''):
    """
    Uses `xdImage()` to display the gridded data array in gridFile. All
    parameters after the name of the gridFile are just passed through
    to `xdImage()`.

    Parameters
    ----------
    gridFile : TYPE
        May be either an `ERS` or `NC` grid file. The default is ''.

    Returns
    -------
    None.

    """
    (xa, _) = gridfile_to_xa(gridFile, bandout=0)
    xdImage(xa, mytitle, colormap=colormap, cmap_norm=cmap_norm, 
        minClip=minClip, maxClip=maxClip, gridlines=True, cb_ticks=cb_ticks, nSigma=nSigma,
        hs=hs, azdeg=azdeg, ax=ax, clipTo3Std = clipTo3Std,
        whizzfile=whizzfile, e_chan=e_chan, n_chan=n_chan)
