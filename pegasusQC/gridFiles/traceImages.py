#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Image the sum of three grids (trace of a tensor).

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

from pegasusQC.gridFiles.xdImage import xdImage
from pegasusQC.gridFiles.gridfiles import gridfile_to_xa
import pegasusQC.config as config


def traceImages(file1, file2, file3, plotTitle, units=''):
    """
    Plots an image of the trace of the three diagonal component grids of the tensor.

    Parameters
    ----------
    file1 : String

        A grid file containing the first diagonal tensor component. May be either an `ERS` or `NC` grid file.

    file2 : String

        A grid file containing the second diagonal tensor component. May be either an `ERS` or `NC` grid file.

    file3 : String

        A grid file containing the third diagonal tensor component. May be either an `ERS` or `NC` grid file.

    plotTitle : String

        A title for the plot of the trace image.

    units : String, optional

        Plotted to the colorbar title. Default is empty string.

    Returns
    -------
    None.

    """
    (x1, _) = gridfile_to_xa(file1, bandout=0)
    (x2, _) = gridfile_to_xa(file2, bandout=0)
    (x3, _) = gridfile_to_xa(file3, bandout=0)
    trace = x1 + x2 + x3
    trace.attrs['units'] = units
    xdImage(trace, plotTitle, colormap=config.qc_colormap, cmap_norm='no', hs=False, azdeg=90)

