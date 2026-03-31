#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Statistics from a grid file.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
from galileoQC.gridFiles.gridfiles import gridfile_to_xa


def imageStats(whizzFile=''):
    """
    Returns the basic statistics of the data in a grid file. The range,
    number of points, and spacings of northing and easting and the range, number
    of samples, mean, standard deviation, and number of real and missing values
    are reported.

    Parameters
    ----------
    whizzFile : Path, optional

        The Path to the grid file in either ERS or NC format. Default is '' in
        which case a file browser allows the user to select the grid file.

    Returns
    -------
    report (String) :

        the statistical summary report.

    """
    xa, fileUsed = gridfile_to_xa(whizzFile, bandout=0)
    if 'Units' in xa.attrs:
        data_units = f" {xa.attrs['Units']}"
    else:
        data_units = ""
        
    report = f'Statistics for {str(fileUsed.name)}'
    report += f'\n Datum  {xa.attrs["datum"]}; Projection {xa.attrs["projection"]}'
    for coord in xa.coords:
        report += f'\n {coord} from {np.nanmin(xa[coord])} to {np.nanmax(xa[coord])}; {xa[coord].shape[0]} samples at spacing {xa[coord].values[1] - xa[coord].values[0]}'
    report += f'\n Value from {np.nanmin(xa.values):.2f}{data_units} to {np.nanmax(xa.values):.2f}{data_units}; {xa.values.shape[0]} x {xa.values.shape[1]} samples'
    report += f'\n     Mean = {np.nanmean(xa.values):.2f}{data_units}, Stdev = {np.nanstd(xa.values):.2f}{data_units}'
    report += f'\n     Number of real values = {np.count_nonzero(~np.isnan(xa.values))}, Number of nans = {np.count_nonzero(np.isnan(xa.values))}'
    return report

