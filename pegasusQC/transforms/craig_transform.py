import numpy as np

from pegasusQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from pegasusQC.transforms.spheretest import (sphereSurvey, _make_xr)
from pegasusQC.transforms.gravity_from_curv import gravity_from_curv
from pegasusQC.transforms._calc_gridcell_size import _calc_gridcell_size
from pegasusQC.gridFiles.xarray_to_grid import xarray_to_grid
from pegasusQC.gridFiles.xdImage import xdImage
from pegasusQC.gridFiles.gridutility import report_gridStats
from pegasusQC.gridFiles.sample_grid_to_line import sample_grid_to_line


def craig_transform(
    whizzFile=None, gne_chan=None, guv_chan=None, gd_chan=None,
    cell_size=None, result_units='um/s/s', mask_polygon=None,
    pad_cells=None, padding_mode="regional", regional_grid_file=None,
    regional_grav_units='mGal',
    numstns=None, firstorder=False
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
        The name of the channel in `whizzFile` to write the gD output to. If the
        name already exists in the whizzFile, then the new data are NOT written.
        Default None in which case, the code invents a name.
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
    mask_polygon : ?##?, optional
        The polygon of the survey boundary. Final output will be trimmed to
        this polygon if provided. Default None and the output is trimmed to
        the smallest rectangle containing all the input curvature data.
    pad_cells : int, optional
        The number of grid cells to pad the data. Default None in which case
        an optimum number is calculated internally.
    padding_mode : String, optional
        The method to be used in filling padding around the grid before DFT.
        Choices are "mean" and "regional". Default is "regional" which requires
        that `regional_grid_file` is provided.
    regional_grid_file : String, optional
        Name of the ERS file containing the regional grid. Required if
        `padding_mode` == "regional". Default None.
    regional_grav_units : String, optional
        The gravity units of the regional grid. Must be either "mGal" or
        "gu" or "um/s/s". Required if `padding_mode` == "regional".
        Default "mGal".
    firstorder : bool, optional
        If True, include first order Craig correction. Default False.
    """
    if not whizzFile is None and padding_mode == "regional":
        if regional_grid_file is None:
            print('\nERROR - regional padding mode requires a regional grid file to be specified.')
            return
        elif regional_grav_units in ("mGal", "mgal"):
            regional_grav_units = "mGal"
        elif regional_grav_units in ("gu", "um/s/s"):
            regional_grav_units = "um/s/s"
        else:
            print('\nERROR - regional grid file gravity units not understood.')
            return

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
            Gne_syn, Guv_syn, cell_size=cell_size, altitude=None,
            result_units=result_units, mask_polygon=None,
            pad_cells=pad_cells, padding_mode=padding_mode,
            firstorder=False
        )

        report_gridStats(gD_err_syn)
        xdImage(gD_err_syn, 'gD_err_syn (um/s/s)', hs=False)
        report_gridStats(gD_grid_syn)
        xdImage(gD_grid_syn, 'gD_grid_syn (um/s/s)', hs=False)
        report_gridStats(origD_grid_syn)
        xdImage(origD_grid_syn, 'origD_grid_syn (um/s/s)', hs=False)
        diff = origD_grid_syn.data[1:-1,1:-1] - gD_grid_syn.data
        diff_gD = gD_grid_syn.copy()
        diff_gD.data = diff
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
        
        # main calculation
        gD_grid, gD_err = gravity_from_curv(
            Ane, Auv, cell_size, gd_chan=gd_chan, altitude=None, 
            result_units=result_units, mask_polygon=mask_polygon,
            pad_cells=pad_cells, padding_mode=padding_mode,
            regional_grid_file=regional_grid_file,
            regional_grav_units=regional_grav_units,
            firstorder=firstorder
        )

        # report statistics, and image grids (calculate clipping limits for images)
        im_min = np.nanmin(gD_grid.data)
        im_max = np.nanmax(gD_grid.data)
        report_gridStats(gD_grid)
        xdImage(gD_grid, 'gD_grid (um/s/s)', minClip=im_min, maxClip=im_max, hs=False)
        report_gridStats(gD_err)
        xdImage(gD_err, 'gD_err (um/s/s)', minClip=im_min, maxClip=im_max, hs=False)

        #
        sample_grid_to_line(gD_grid, whizzFile)
    
    return gD_grid
