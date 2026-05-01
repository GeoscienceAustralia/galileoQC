#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Craig transform Gne/Guv to gD.

Author: Mark Helm Dransfield

Created: Oct 2025

License: CC BY-SA
"""

import numpy as np
import xarray as xr
import pathlib
import matplotlib.pyplot as plt

from galileoQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from galileoQC.transforms.spheretest import (sphereSurvey, _make_xr)
from galileoQC.transforms.gravity_from_curv import gravity_from_curv
from galileoQC.transforms._calc_gridcell_size import _calc_gridcell_size
from galileoQC.gridFiles.xarray_to_grid import xarray_to_grid
from galileoQC.gridFiles.xdImage import xdImage
from galileoQC.gridFiles.gridutility import report_gridStats
from galileoQC.gridFiles.sample_grid_to_line import sample_grid_to_line
from galileoQC.gridFiles.grid_to_xarray import gridfile_to_xa
from galileoQC.transforms.conform import conform
from galileoQC.gridFiles.write_ers import write_ers


def craig_transform(
    whizzFile=None, gne_chan=None, guv_chan=None, gd_chan=None,
    cell_size=None, result_units='um/s/s', survey_polygon=None,
    pad_cells=None, padding_mode="regional", regional_grid_file=None,
    regional_grav_units='mGal', numstns=None, firstorder=False,
    conforming=False, save_to_ers=False,
    plot_flag=False, verbose=False
):
    """
    Runs the Craig transform on gravity differential curvature data to
    calculate the vertical component of gravitational acceleration. Also
    reports statistics of results and plots relevant images.

    Parameters
    ----------
    whizzFile : Path or String, optional

        The Path to, or String name of, the whizz file in HDF5 format.
        Default None, in which case the transform is performed on internally
        generated synthetic data.

    gne_chan : String, optional

        The name of the channel in `whizzFile` containing the Gne survey data.
        Default None but must be provided if `whizzfile` is given.

    guv_chan : String, optional

        The name of the channel in `whizzFile` containing the Guv survey data.
        Default None but must be provided if `whizzfile` is given.

    gd_chan : String, optional

        The name of the channel in `whizzFile` to write the gD output. If the
        name already exists for a flight-line in the whizzFile, then the new
        data are NOT written.  Also used as the ERS filename if the grid is
        written to an ERS file and is the `long name` attribute in the output
        DataArray.) Default None in which case, the new data are NOT written.

    altitude_chan : String, optional

        The name of the channel in `whizzFile` containing the altitude data.
        Default None and ignored if provided. Later versions of s/w may use
        this input.

    cell_size : float, optional

        The grid cell size. Recommended as ~1/4 line spacing. Default None in
        which case an "optimum" (?!) number is calculated internally. Generally
        it is better to provide ~1/4 line spacing.

    result_units : String, optional

        The gravity units of the final resultant grid. Must be either "mGal" or
        "gu" or "um/s/s". Default "um/s/s".

    survey_polygon : numpy 2D array, optional

        The polygon vertices of the survey boundary, as an array or sequence of (x,y)
        pairs, in either clockwise or counter-clockwise order around the boundary.
        For example, survey_polygon = [(0, 0), (1, 0), (1,1), (0,1)]. Final output will be
        trimmed to this polygon if provided. Default None.

    pad_cells : int, optional

        The number of grid cells to pad the data. Default None in which case
        an optimum number is calculated internally.

    padding_mode : String, optional

        The method to be used in filling padding around the grid before DFT.
        Choices are "mean" and "regional". Default is "regional" which requires
        that `regional_grid_file` is provided.

    regional_grid_file : String, optional

        Name of the ERS file containing the regional grid. Required if
        `padding_mode` == "regional" or if `confroming` == True. Default None.

    regional_grav_units : String, optional

        The gravity units of the regional grid. Must be either "mGal" or
        "gu" or "um/s/s". Required if `padding_mode` == "regional".
        Default "mGal".

    firstorder : bool, optional

        If True, include first order Craig correction. Default False.
        
    conforming : bool, optional

        If True, then high-pass filter the gD data before returning it.
        Default False.

    save_to_ers : bool, optional

        If True, save the the transformed output grid as an ERS file. Default False.

    plot_flag : bool, optional

        If True, plot images of grids at intermediate stages of the processing.
        Meant for debugging the code but can be helpful for understanding the
        method. Default False.

    verbose : Bool, optional

        If True, then prints out details. Default False.

    RETURNS
    ----------
    
    gd_grid : xarray DataArray

        The resultant grid from the transform.
    """
    if not whizzFile is None and (padding_mode == "regional" or conforming == True):
        if regional_grid_file is None:
            print('\nERROR - regional padding mode and conforming require a regional grid file to be specified.')
            return
        elif regional_grav_units in ("mGal", "mgal"):
            regional_grav_units = "mGal"
        elif regional_grav_units in ("gu", "um/s/s"):
            regional_grav_units = "um/s/s"
        else:
            print('\nERROR - regional grid file gravity units not understood. Must be "mGal" or "um/s/s".')
            return
    else:
        regional_grid_file = None
        regional_grid = None

    if whizzFile is None:
        cell_size = 10.0
        pad_cells = 256
        result_units='um/s/s'
        padding_mode = "mean"
        
        n, e, d, Ann, Ane, And, Aee, Aed, Add, Auv, AD = sphereSurvey(numstns)
        # need data in xarray datasets
        Guv_syn = _make_xr(Auv, 'Guv', 'E', n, e)
        Gne_syn = _make_xr(Ane, 'Gne', 'E', n, e)
        gD_syn  = _make_xr(AD, 'gD', 'um/s/s', n, e)
        origD_grid_syn, _ = xarray_to_grid(gD_syn, cell_size, region=None, method='minc',
                  mask_polygon=[], mask_pixels=0, numneighbours=5, bdist=None, maxiters=100
                 )
        gD_grid_syn, gD_err_syn, Ane_grid = gravity_from_curv(
            Gne_syn, Guv_syn, cell_size, gd_chan=gd_chan, altitude=None,
            result_units=result_units, survey_polygon=survey_polygon,
            pad_cells=pad_cells, padding_mode=padding_mode, regional_grid=None,
            firstorder=False
        )

        if plot_flag:
            report_gridStats(gD_err_syn)
            xdImage(gD_err_syn, 'gD_err_syn (um/s/s)', hs=False)
            report_gridStats(gD_grid_syn)
            xdImage(gD_grid_syn, 'gD_grid_syn (um/s/s)', hs=False)
            report_gridStats(origD_grid_syn)
            xdImage(origD_grid_syn, 'origD_grid_syn (um/s/s)', hs=False)
            diff = origD_grid_syn.data[1:-1,1:-1] - gD_grid_syn.data
        diff_gD = gD_grid_syn.copy()
        diff_gD.data = diff
        if plot_flag:
            report_gridStats(diff_gD)
            xdImage(diff_gD, 'origD_grid_syn-gD_grid_syn (um/s/s)', hs=False)

        gD_grid = gD_grid_syn
    else:
        print('Reading located data from whizz file.')
        Ane = whizz_to_xarray(whizzFile, gne_chan)
        Auv = whizz_to_xarray(whizzFile, guv_chan)
        
        # Calculate the cell size while we know the whizz_file name
        if cell_size is None:
            cell_size = _calc_gridcell_size(whizzFile)
            if cell_size is None:
                print('ERROR - could not calculate cell size, please enter explicitly.')
                return None
        
        # If required, get the regional grid from file
        if conforming or padding_mode == 'regional':
            e_chan = Ane.attrs['x_channel']
            n_chan = Ane.attrs['y_channel']
            print('\nReading in the regional data.')
            reg_xa, _ = gridfile_to_xa(regional_grid_file)
            # Nicer if we use the same coordinate names for all!
            regional_grid = xr.DataArray(data=reg_xa.data,
                                         dims=['y', 'x'],
                                         coords={
                                             'y': reg_xa.y.values,
                                             'x': reg_xa.x.values
                                         },
                                         attrs={
                                         'author': 'Mark Dransfield',
                                         'x_channel': e_chan,
                                         'y_channel': n_chan,
                                         'z_channel': regional_grid_file,
                                         'Units': regional_grav_units
                                         })
            if verbose:
                report_gridStats(regional_grid, title='Regional Grid')
            if plot_flag:
                xdImage(regional_grid, 'Regional Grid')
            # regional_grid = _grid_clean(regional_grid)

        # main calculation
        gD_grid, gD_err, gD_raw = gravity_from_curv(
            Ane, Auv, cell_size, gd_chan=gd_chan, altitude=None, 
            result_units=result_units, survey_polygon=survey_polygon,
            pad_cells=pad_cells, padding_mode=padding_mode,
            regional_grid=regional_grid,
            firstorder=firstorder,
            plot_flag=plot_flag,
            verbose=verbose
        )
        if plot_flag:
            xdImage(gD_grid, 'Post-gfc: gD_grid (um/s/s)', hs=False)

        # conformed to regional if required.
        report_summary = '\n  Final transformed grids:'
        if conforming:
            report_summary = '\n  Final conformed and transformed grid stats:'
            if plot_flag:
                report_gridStats(gD_raw, title='Raw grid')
                report_gridStats(gD_grid, title='gD grid')
                report_gridStats(regional_grid, title='Regional grid')
                xdImage(gD_grid, ' Pre-conform: gD_grid (um/s/s)', hs=False)
                xdImage(gD_raw, ' Pre-conform: gD_raw (um/s/s)', hs=False)
            original_mask = np.ma.masked_array(gD_grid, ~np.isnan(gD_grid)).mask
            if plot_flag:
                plt.imshow(original_mask, origin='lower')
            if survey_polygon is None:
                survey_polygon = [
                    (float(gD_raw.x[0]), float(gD_raw.y[0])), 
                    (float(gD_raw.x[-1]), float(gD_raw.y[0])), 
                    (float(gD_raw.x[-1]), float(gD_raw.y[-1])), 
                    (float(gD_raw.x[0]), float(gD_raw.y[-1]))]

            gD_grid = conform(gD_raw, regional_grid, survey_polygon=survey_polygon, plot_flag=plot_flag, original_mask=original_mask)
            if gD_grid is None:
                return None
            if plot_flag:
                report_gridStats(gD_grid, title='conformed grid')

        # report statistics, and image grids (calculate clipping limits for images)
        print(report_summary)
        report_gridStats(gD_grid)
        im_min = np.nanmin(gD_grid.data)
        im_max = np.nanmax(gD_grid.data)
        print('  Final error grid stats:')
        report_gridStats(gD_err)
        im_mean = np.nanmean(gD_err.data)

        myfig, myax = plt.subplots(1,2, figsize=(12,6))
        gD_grid.plot(ax=myax[0], vmin=im_min, vmax=im_max)
        myax[0].set_xlabel(gD_grid.attrs['x_channel'])
        myax[0].set_ylabel(gD_grid.attrs['y_channel'])
        myax[0].set_title(gD_grid.attrs['title'])
        gD_err.plot(ax=myax[1], vmin=im_mean-(im_max-im_min)/2., vmax=im_mean+(im_max-im_min)/2.)
        myax[1].set_xlabel(gD_err.attrs['x_channel'])
        myax[1].set_ylabel(gD_err.attrs['y_channel'])
        myax[1].set_title(gD_err.attrs['title'])

        # Sample the result back into the whizzFile database.
        if not gd_chan is None:
            sample_grid_to_line(gD_grid, whizzFile)
    
        if save_to_ers:
            if isinstance(whizzFile, pathlib.Path):
                whizzFilePath = whizzFile
            elif isinstance(whizzFile, str):
                whizzFilePath = Path(whizzFile)
            else:
                print('Error - type of whizzFile not recognised. Must be Path or String')
                print('    Output grid not written to file.')
                return gD_grid

            my_long_name = "gD_craig"
            if 'long_name' in gD_grid.attrs:
                my_long_name = gD_grid.attrs['long_name']

            ersfilepath = whizzFilePath.with_name(my_long_name)
            write_ers(ersfilepath, gD_grid)

    return gD_grid
