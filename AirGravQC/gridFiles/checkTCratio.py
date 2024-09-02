from AirGravQC.gridFiles.xdImage import xdImage
from AirGravQC.gridFiles.gridfiles import gridfile_to_xa
import colorcet as cc


def checkTCratio(file000, filexxx, xxx, fileyyy, yyy, plotTitle):
    """
    Shows an image of the difference between data terrain-corrected
    at density y and the prediction based on un-corrected data and
    data terrain-corrected at x.

    Used to QC for consistency in terrain corrections.

    Parameters
    ----------
    file000 : String
        Name of the file containing the data before terrain correction. May be either an `ERS` or `NC` grid file.
    filexxx : String
        Name of the file containing the data after xpxx terrain correction. May be either an `ERS` or `NC` grid file.
    xxx : Float
        The density used in the filexxx terrain correction.
    fileyyy : String
        Name of the file containing the data after ypyy terrain correction. May be either an `ERS` or `NC` grid file.
    yyy : Float
        The density used in the fileyyy terrain correction.
    plotTitle : String
        A title for the plot of the difference between predicted and actual
        data in fileyyy.

    Returns
    -------
    None.

    """
    (g0, _) = gridfile_to_xa(file000, bandout=0)
    (gx, _) = gridfile_to_xa(filexxx, bandout=0)
    (gy, _) = gridfile_to_xa(fileyyy, bandout=0)
    
    tx = g0 - gx
    predy = g0 - yyy / xxx * tx
    predErr = predy - gy

    xdImage(predErr, plotTitle, colormap=cc.m_CET_L9, cmap_norm='no', hs=False, azdeg=90)
