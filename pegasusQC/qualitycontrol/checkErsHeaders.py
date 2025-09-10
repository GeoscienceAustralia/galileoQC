#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the ERS headers for grid files are consistent.
"""
import numpy as np

import pegasusQC.gridFiles.read_ers as ers


def checkErsHeaders(folderPath=r'\.'):
    """
    Compares all .ers files in the folder for certain key parameters and reports
    any that are different from the first file found in any parameter value.

    Parameters
    ----------
    folderPath : Path, optional
        The folder or directory containing the .ers files. The default is backslash dot.

    Returns
    -------
    None.

    """
    reportString = ''
    allOk = True
    fileOK = True
    
    # get a list of the ers file paths
    file_count = 0
    ersFiles = []
    folderFiles = folderPath.iterdir()
    for aFile in folderFiles:
        if aFile.is_file() and (aFile.suffix == '.ERS' or aFile.suffix == '.ers') and (aFile.name[0] != '.'):
            ersFiles.append(aFile)
            file_count += 1
    print(f'Found {file_count} .ers files ...')
    print(f'in: {str(folderPath)}')

    # get the header dict for each ers file
    ersDicts = []
    for aFile in ersFiles:
            [ncells, nrows, nbands, eastings, northings, nullcell,
                    precision, headerbytes, originalnullcell, byteorder, datum,
                    projection, units] = ers.read_ers_header(str(aFile))
            commonOK, reportStr = _commonErsHdrErrors(ncells, nrows, nbands, nullcell,
                    precision, headerbytes, originalnullcell, byteorder, datum,
                    projection, units)
            ersDicts.append([ncells, nrows, nbands, nullcell,
                    precision, headerbytes, originalnullcell, byteorder, datum,
                    projection, units])
            
    # compare the contents line-by-line
    print(f'Comparing ERS files against {ersFiles[0].name}.')
    firstDict = ersDicts[0]
    print(firstDict)
    for jj in range(1, len(ersDicts)):
        fileOK = True
        print(f'Checking file {ersFiles[jj].name}')
        for ii in range(0, len(firstDict)):
            if firstDict[ii] != ersDicts[jj][ii]:
                allOk = False
                fileOK = False
                print(f'  Different element {ii}. {firstDict[ii]} != {ersDicts[jj][ii]}')
        if fileOK:
            print('  Checked OK.')
    
    return


def _commonErsHdrErrors(ncells, nrows, nbands, nullcell, precision, headerbytes, originalnullcell, byteorder, datum, projection, units):
    """
    Checks a set of ERS header fields for various common errors.

    Parameters
    ----------
    ncells : Integer
        DESCRIPTION.
    nrows : Integer
        DESCRIPTION.
    nbands : Integer
        DESCRIPTION.
    nullcell : Integer
        DESCRIPTION.
    precision : Integer
        DESCRIPTION.
    headerbytes : Integer
        DESCRIPTION.
    originalnullcell : Integer
        DESCRIPTION.
    byteorder : Integer
        DESCRIPTION.
    datum : Integer
        DESCRIPTION.
    projection : Integer
        DESCRIPTION.

    Returns
    -------
    allOk : Bool
        True if no errors found, else False.
    reportStr : String
        A report listing all errors found.

    """
    allOk = True
    reportStr = ''
    
    if projection == 'WGS84':
        reportStr += 'ERROR: projection cannot be WGS84'
        allOk = False
    
    if units not in ['NONE', 'METERS', 'METER', 'METRES', 'METRE']:
        reportStr += 'ERROR: units cannot be {units}'
        allOk = False
    
    return allOk, reportStr


