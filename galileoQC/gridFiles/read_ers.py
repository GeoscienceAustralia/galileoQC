#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read from an ermapper grid into numpy arrays.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
from pathlib import Path

def read_ers_image(ersFile, bandout=0):
    """
    Read eastings, northings, value from an ermapper grid into numpy arrays.

    Can only handle a restricted range of ERMapper images. In particular:

       DataSetType = ERStorage

       DataType = Raster

       CoordinateType = EN

       Rotation = 0:0:0.0

    Parameters
    ----------
    ersFile : String or Path

        The name of the ERMapper file, including the `.ers` extension.

    bandout : Int, optional

        The band. The default is 0, which gives the same result as 1.

    Returns
    -------
    1D numpy array

        Eastings - vector containing coordinates of MIDPOINT of each
        pixel in an image row.

    1D numpy array

        Northings - vector containing coordinates of MIDPOINT of each
        pixel in an image row.

    2D numpy array

        The grid values at (e,n).

    String

        The geographic datum as read from `ersFile`.

    String

        The geographic projection as read from `ersFile`.

    String

        The units string as read from `ersFile`, default "".

    """
    headerfile = str(ersFile)
    imagefile = str(ersFile.with_suffix(''))
    # should put a check arguments section here
    # imagefile = 'DrapeSurface_Fourier_FS'
    select_band = False            # the default, return all bands
    # want to be able to select band if bandout exists, select_band = True

    # get the header info we need to decode the image file
    [ncells, nrows, nbands, eastings, northings, nullcell, precision,
     headerbytes, originalnullcell, byteorder, datum,
     projection, units] = read_ers_header(headerfile)
    # [eastings, northings, zz, ncells, nrows, nbands, nullcell,
    # originalnullcell] = read_ers_header(imagefile)

    if (select_band):
        if (bandout > nbands):
            print('Asked for band # {0:2d}, but this imagefile only has \
                  {1:2d} bands!', bandout, nbands)
            print('Error in command usage for ermapper.read_ers_image')
            return  # check usage !!

        if (nbands == 1):
            select_band = False

    if nbands == 1:
        if precision == 'float64':
            zz = np.fromfile(imagefile, 'float64').reshape(nrows, ncells)
        elif precision == 'float32':
            zz = np.fromfile(imagefile, 'float32').reshape(nrows, ncells)
        else:
            print('unrecognised precision {}'.format(precision))
            return

        # want 'float32' a variable from byteorder
        # assume zero!! if headerbytes > 0: # skip this many bytes
        #     somehow ...

        # cif nrows x nbands is wrong, then report and break
        #   fprintf(1, 'Error in reading ERmapper imagefile: %s \n', imagefile)
        #   fprintf(1, 'Tried to read %d x %d = %d pixels \n', ncells, nrows,
        #            ncells x nrows)
        #   fprintf(1, 'Read = %d pixels \n', count)
        #   return # or break??

        # zz = flipud(zz')
    else:
        print('Error in reading ERmapper imagefile')
        print('Lazy Mark has not got around to multi-band reads.... Sorry.\n')
        return
    # now muck around with the nullcells  - convert these to NaN's regardless
    # of what they where on input. Should check for -inf as well!!

    if np.isinf(nullcell):
        zz[zz == np.inf] = np.nan
    else:
        if (nullcell != 0):
            zz[abs(zz - nullcell) < 0.01] = np.nan
        # elif (nullcell < 0):
        #   zz[abs(zz - nullcell) < 0.01] = np.nan
        else:
            zz[zz == nullcell] = np.nan

    nullcell = np.nan

    return eastings, northings, zz, datum, projection, units


def read_ers_header(imagefile):
    """
    Reads ERMapper ers header file information and returns same in a dict.

    Aims to extract key parameters from an ERMapper image header and return
    these to the caller. Attempts to parse the ERmapper header by extracting
    all the lines containing '=', and evaluating them. This approach may
    fail with unusual header types.

    Parameters
    ----------
    imagefile : string

        Name of ERMAPPER header file.

    Returns
    -------
    flag : Bool

        True if operation was successful.

    headerdict : dictionary

        Contains ncells (scalar number of image pixels in easting direction),
        nrows (scalar number of image pixels in northing direction), 
        nbands (number of bands in image),
        eastings (vector of coordinates of MIDPOINT of each pixel in an image row),
        northings (vector of coordinates of MIDPOINT of each pixel in an image column),
        nullcell (scalar value of NULL for the image),
        precision (string describing the image precision),
        headerbytes (number of bytes on header of image),
        byteorder (string equal to either 'ieee-be' (MSBFirst) or 'ieee-le' (LSBFirst)).

    """

    # should check input arguments

    headerdict = _parse_header(imagefile)

    # should put checks here that headerdict is ok

    ncells = headerdict['ncells']
    nrows = headerdict['nrows']
    nbands = headerdict['nbands']
    eastings = headerdict['eastings']
    northings = headerdict['northings']
    nullcell = headerdict['nullcell']
    precision = headerdict['precision']
    headerbytes = headerdict['headerbytes']
    originalnullcell = headerdict['originalnullcell']
    byteorder = headerdict['byteorder']
    datum = headerdict['datum']
    projection = headerdict['projection']
    units = headerdict['units']

    return [ncells, nrows, nbands, eastings, northings, nullcell,
            precision, headerbytes, originalnullcell, byteorder, datum,
            projection, units]


def _parse_header(imagefile, verbose=False):
    """
    Parses the information in an ERMapper header file and returns the information
    as a dictionary.

    Parameters
    ----------
    imagefile : String

        The name of the ERMapper header file to read.

    verbose : Bool, optional

        If True, then prints out details. Default = False.

    Returns
    -------
    headerdict : Dictionary

        Contains the metadata read from the header file.

    """
    import numpy as np
    if verbose:
        print(f'Parsing {imagefile}')

    # Prepare a dict containing all the fields we need to obtain from the
    # ermapper header file with defaults set.
    ersHeader = dict(
        NrOfCellsPerLine=0,
        NrOfLines=0,
        NrOfBands=0,
        Eastings=0.0,
        RegistrationCellX=0,
        Xdimension=30.0,
        Northings=0.0,
        RegistrationCellY=0,
        Ydimension=30.0,
        CellType='float32',
        NullCellValue=-9999.0,
        HeaderOffset=0,
        ByteOrder='ieee-le',
        Datum='RAW',
        Projection='RAW',
        Units='NONE')

    # open the filename
    headerfile = imagefile

    headerfilepath = Path(headerfile)
    if headerfilepath.with_suffix('.ers').exists():
        headerfilepath = headerfilepath.with_suffix('.ers')
    elif headerfilepath.with_suffix('.ERS').exists():
        headerfilepath = headerfilepath.with_suffix('.ERS')
    else:
        print(f'ERROR file {headerfile} does not exist.')
        return
    headerfile = str(headerfilepath)

    with open(headerfile, 'r') as fid:
        for myline in fid:
            if '=' in myline:
                [keyword, dummy, value] = myline.partition('=')

                # Remove leading and trailing blanks and translate
                keyword = keyword.strip()
                value = value.strip()
                value = _translate_value(value)
                rawvalue = value[0]

                if 'DATUM' in keyword.upper():
                    ersHeader['Datum'] = rawvalue
                elif 'PROJECTION' in keyword.upper():
                    ersHeader['Projection'] = rawvalue
                elif 'UNITS' in keyword.upper():
                    ersHeader['Units'] = rawvalue
                elif 'NROFCELLSPERLINE' in keyword.upper():
                    ersHeader['NrOfCellsPerLine'] = int(rawvalue)
                elif 'NROFLINES' in keyword.upper():
                    ersHeader['NrOfLines'] = int(rawvalue)
                elif 'NROFBANDS' in keyword.upper():
                    ersHeader['NrOfBands'] = int(rawvalue)
                elif 'XDIMENSION' in keyword.upper():
                    ersHeader['Xdimension'] = float(rawvalue)
                elif 'EASTINGS' in keyword.upper() or 'METERSX' in keyword.upper():
                    if len(rawvalue) == 0:
                        ersHeader['Eastings'] = float(0.0)
                    else:
                        ersHeader['Eastings'] = float(rawvalue)
                elif 'REGISTRATIONCELLX' in keyword.upper():
                    ersHeader['RegistrationCellX'] = float(rawvalue)
                elif 'NORTHINGS' in keyword.upper()  or 'METERSY' in keyword.upper():
                    if len(rawvalue) == 0:
                        ersHeader['Northings'] = float(0.0)
                    else:
                        ersHeader['Northings'] = float(rawvalue)
                elif 'REGISTRATIONCELLY' in keyword.upper():
                    ersHeader['RegistrationCellY'] = float(rawvalue)
                elif 'YDIMENSION' in keyword.upper():
                    ersHeader['Ydimension'] = float(rawvalue)
                elif 'CELLTYPE' in keyword.upper():
                    ersHeader['CellType'] = rawvalue
                elif 'NULLCELLVALUE' in keyword.upper():
                    ersHeader['NullCellValue'] = float(rawvalue)
                elif 'HEADEROFFSET' in keyword.upper():
                    ersHeader['HeaderOffset'] = int(rawvalue)
                elif 'BYTEORDER' in keyword.upper():
                    ersHeader['ByteOrder'] = rawvalue
    fid.close()
    # end with and while to read header file
    if verbose:
        print(ersHeader)

    ncells = ersHeader['NrOfCellsPerLine']
    nrows = ersHeader['NrOfLines']
    nbands = ersHeader['NrOfBands']
    xcell = ersHeader['Xdimension']
    ycell = ersHeader['Ydimension']
    xorigin = ersHeader['Eastings'] - ersHeader['RegistrationCellX'] \
        * ersHeader['Xdimension']
    yorigin = ersHeader['Northings'] + ersHeader['RegistrationCellY'] \
        * ersHeader['Ydimension'] - (ersHeader['NrOfLines']+1) \
        * ersHeader['Ydimension']
    precision = ersHeader['CellType']
    nullcell = ersHeader['NullCellValue']
    headerbytes = ersHeader['HeaderOffset']
    byteorder = ersHeader['ByteOrder']

    # row vector creation x = np.array([1,3,2])[None,:]
    # we use 1D array x = np.array([1,3,4])
    # want e,n of same Type as z, ie CellType
    eastings = np.arange(0, ncells)  # [None,:]
    eastings = xorigin + (eastings - 0.5) * xcell
    northings = np.arange(0, nrows)  # [None,:]
    northings = yorigin + (northings + 0.5) * ycell
    if verbose:
        print(precision, precision.__class__)
    if precision == 'float32':
        eastings = np.float32(eastings)
        northings = np.float32(northings)

    originalnullcell = nullcell
    # special test on Intrepid Nullcell's. Unless the image is stored as
    # float64, then -5e+75 will appear as -Inf
    # un-tested and not understood !!
    if ((nullcell < -4.0e+75) and (nullcell > -6.0e+75)):
        if not ('float64' in precision):
            nullcell = -np.inf

    if verbose:
        print('nEastings, nNorthings, nbands ', ncells, nrows, nbands)
        print('eastings[0], northings[0] ', eastings[0], northings[0])
        print('nullcell, precision, headerbytes ', nullcell, precision, headerbytes)
        print('originalnullcell, byteorder ', originalnullcell, byteorder)
        print('Projection, Datum', ersHeader['Projection'], ersHeader['Datum'])
        print('Coordinate Units ', ersHeader['Units'])

    headerdict = dict(ncells=ncells,
                      nrows=nrows,
                      nbands=nbands,
                      eastings=eastings,
                      northings=northings,
                      nullcell=nullcell,
                      precision=precision,
                      headerbytes=headerbytes,
                      originalnullcell=originalnullcell,
                      byteorder=byteorder,
                      datum=ersHeader['Datum'],
                      projection=ersHeader['Projection'],
                      units=ersHeader['Units'])

    return headerdict


def _translate_value(rawvalue):
    """
    Translates certain ERMapper header field names into something that
    Python understands.

    Parameters
    ----------
    rawvalue : TYPE

        DESCRIPTION.

    Returns
    -------
    list

        DESCRIPTION.

    """

    # remove quotes and replace parentheses from value string
    rawvalue = rawvalue.replace('\'', '')
    rawvalue = rawvalue.replace('"', '')
    # don't change without changing matching bracket
    rawvalue = rawvalue.replace('{', '[')
    rawvalue = rawvalue.replace('}', ']')

    # replacements that matlab understands
    rawvalue = rawvalue.replace('Unsigned8BitInteger', 'uint8')
    rawvalue = rawvalue.replace('Signed8BitInteger', 'int8')
    rawvalue = rawvalue.replace('Unsigned16BitInteger', 'uint16')
    rawvalue = rawvalue.replace('Signed16BitInteger', 'int16')
    rawvalue = rawvalue.replace('Signed32BitInteger', 'int32')
    rawvalue = rawvalue.replace('IEEE4ByteReal', 'float32')
    rawvalue = rawvalue.replace('IEEE8ByteReal', 'float64')
    rawvalue = rawvalue.replace('MSBFirst', 'ieee-be')
    rawvalue = rawvalue.replace('LSBFirst', 'ieee-le')
    # and their upper case equivalents
    rawvalue = rawvalue.replace('UNSIGNED8BITINTEGER', 'uint8')
    rawvalue = rawvalue.replace('SIGNED8BITINTEGER', 'int8')
    rawvalue = rawvalue.replace('UNSIGNED16BITINTEGER', 'uint16')
    rawvalue = rawvalue.replace('SIGNED16BITINTEGER', 'int16')
    rawvalue = rawvalue.replace('SIGNED32BITINTEGER', 'int32')
    rawvalue = rawvalue.replace('IEEE4BYTEREAL', 'float32')
    rawvalue = rawvalue.replace('IEEE8BYTEREAL', 'float64')
    rawvalue = rawvalue.replace('MSBFIRST', 'ieee-be')
    rawvalue = rawvalue.replace('LSBFIRST', 'ieee-le')

    return [rawvalue]


