#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Write an xarray DataArray to an ERS file.

Author: Mark Helm Dransfield

Created: ca 2026

License: CC BY-SA
"""

import numpy as np
import pathlib
from pathlib import Path
from datetime import datetime, timezone


def write_ers(imagefile, data_array, nullcell=-999999):
    """
    Write a DataArray into an ERMapper image.

    Parameters
    ----------
    imagefile : String

        The name of the image file (includes path, but NOT the .ers suffix).

    data_array : xarray DataArray

        The array to be written to the image file.

    nullcell : Float, optional

        The value of any null cells (replaces any nans in the DataArray
        with this in the file). Default -999999.

    Returns
    -------
    NONE.

    """
    if isinstance(imagefile, pathlib.Path):
        imagefileStr = str(imagefile)
        imagefilePath = imagefile
    elif isinstance(imagefile, str):
        imagefileStr = imagefile
        imagefilePath = Path(imagefile)
    else:
        print('Error - type of imagefile not recognised. Must be Path or String')
        return

    if imagefilePath.suffix.upper() == ".ERS":
        imagefilePath = imagefilePath.with_suffix("")
        headerFilePath = imagefilePath.with_suffix(".ers")
    elif imagefilePath.suffix != "":
        print('Error - suffix of imagefile not recognised. Must be nothing.')
        return
    headerFilePath = imagefilePath.with_suffix(".ers")

    zz = data_array.values
    zz[np.isnan(zz)] = nullcell

    # check arguments here, check file doesn't exist
    # then modify zz as required for ERMapper

    fid = open(imagefilePath, 'wb')
    zz = np.flipud(data_array.values)
    zz.tofile(fid)
    fid.close()

    _write_ers_header(headerFilePath, data_array, nullcell=-999999)

    print(f'Grid written to ERS file {str(imagefilePath)}.')
    return


def _write_ers_header(headerFilePath, data_array, nullcell=-999999):
    """
    Write an ERMapper image header in crude form.

    Parameters
    ----------
    headerFilePath : Path

        The name of the header file (including .ers suffix).

    data_array : xarray DataArray

        The array written to the image file.

    nullcell : Float, optional

        The value of any null cells (replaces any nans in the DataArray
        with this in the file). Default -999999.

    Returns
    -------
    NONE.

    """
    have_units = False
    my_long_name = headerFilePath.with_suffix('').name
    my_datum = "NONE"
    my_projection = "NONE"
    my_CRS = "NONE"

    if 'long_name' in data_array.attrs:
        my_long_name = data_array.attrs['long_name']
    if 'Units' in data_array.attrs:
        my_units = data_array.attrs['Units']
        have_units = True
    if 'datum' in data_array.attrs:
        my_datum = data_array.attrs['datum']
    if 'projection' in data_array.attrs:
        my_projection = data_array.attrs['projection']
    if 'CRS' in data_array.attrs:
        my_CRS = data_array.attrs['CRS']

    my_fmt = '%a %d %b %Y %H:%M:%S'
    my_datime = datetime.now(timezone.utc).strftime(my_fmt) + " UTC"

    dx = data_array.x.values[1] - data_array.x.values[0]
    dy = data_array.y.values[1] - data_array.y.values[0]

    nrows = len(data_array.y)
    ncells = len(data_array.x)

    estr = f'{data_array.x.values[0]:.2f}'
    nstr = f'{data_array.y.values[-1]:.2f}'

    reg_cell_x = 0
    reg_cell_y = ncells - 1

    with open(headerFilePath, 'wt') as fid:
        fid.write(f'DatasetHeader Begin\n')
        fid.write(f'    Version        = "5.5"\n')
        fid.write(f'    Name           = "{my_long_name}"\n')
        fid.write(f'    LastUpdated    = {my_datime}\n')
        fid.write(f'    DataSetType    = ERStorage\n')
        fid.write(f'    DataType       = Raster\n')
        fid.write(f'    ByteOrder      = LSBFirst\n')
        if have_units:
            fid.write(f'    ZUnits          = "{my_units}"\n')
        fid.write(f'    Comments       = "written by galileoQC.write_ers"\n')
        fid.write(f'    CoordinateSpace Begin\n')
        fid.write(f'        Datum           = "{my_datum}"\n')
        fid.write(f'        Projection      = "{my_projection}"\n')
        fid.write(f'        CRS             = "{my_CRS}"\n')
        fid.write(f'        CoordinateType  = EN\n')
        fid.write(f'        Units           = "METERS"\n')
        fid.write(f'        Rotation        = 0:0:0.0\n')
        fid.write(f'    CoordinateSpace End\n')
        fid.write(f'    RasterInfo Begin\n')
        fid.write(f'        CellType      = IEEE8ByteReal\n')
        fid.write(f'        NullCellValue = {nullcell}\n')
        fid.write(f'        CellInfo Begin\n')
        fid.write(f'            Xdimension = {dx}\n')
        fid.write(f'            Ydimension = {dy}\n')
        fid.write(f'        CellInfo End\n')
        fid.write(f'        NrOfLines        = {nrows}\n')
        fid.write(f'        NrOfCellsPerLine = {ncells}\n')
        fid.write(f'        RegistrationCoord Begin\n')
        fid.write(f'            Eastings   = {estr}\n')
        fid.write(f'            Northings  = {nstr}\n')
        fid.write(f'        RegistrationCoord End\n')
        fid.write(f'        RegistrationCellX    = {reg_cell_x}\n')
        fid.write(f'        RegistrationCellY    = {reg_cell_y}\n')
        fid.write(f'        NrOfBands    = 1\n')
        fid.write(f'    RasterInfo End\n')
        fid.write(f'DatasetHeader End\n')
    return

