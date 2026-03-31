#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check terrain correction density.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

from galileoQC.gridFiles.xdImage import xdImage
from galileoQC.gridFiles.gridfiles import gridfile_to_xa
import galileoQC.config as config


def checkTCratio(file000, filexxx, xxx, fileyyy, yyy, plotTitle, units=''):
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

    units : String, optional

        Plotted to the colorbar title. Default is empty string.

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
    predErr.attrs['Units'] = units

    xdImage(predErr, plotTitle, colormap=config.qc_colormap, cmap_norm='no', hs=False, azdeg=90)
