#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check data gaps are within specification.
Author: Mark Helm Dransfield
Created: ca 2023
License: CC BY-SA
"""

import numpy as np
import h5py

import pegasusQC.config as config

groupName = config.groupName


def checkGaps(whizzFile, ignored_chans=[], maxGapSec=0.0, maxNumGaps=0, lines=[], verbose=True):
    """
    Checks every dataset for each channel and each survey line in filePath for
    gaps, and reports all gaps found.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    maxGapSec :  Float, optional
        The largest allowed gap measured in seconds. Default 0.0
    maxNumGaps : Integer, optional
        The maximum number of gaps allowed on any survey line. Default 0
    lines : Array{String}, optional
        Array of line numbers as strings. Default = [], meaning all lines are checked.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        lineGroups = list(g.values())
        channelNames = list(lineGroups[0].keys())
        num_lines_failed = 0
        total_num_lines = 0
        report = ''

        if lines == []:
            lines = list(g.keys())

        for line in lines:
            total_num_lines += 1
            gaps_on_line = 0
            num_channels = 0
            lineNo = line
            lineText = f'Line {lineNo}'
            for channel in channelNames:
                if channel in ignored_chans:
                    continue
                else:
                    num_channels += 1
                    numberMissing = np.count_nonzero(np.isnan(g[line][channel]))
                    if numberMissing > 0:
                        lineText += f'\n    {channel}, nans: {numberMissing}'
                        gaps_on_line += 1

            if gaps_on_line > 0:
                num_lines_failed += 1
                report += lineText + '\n'
        print(f'Checking for all gaps in {num_channels} channels on {total_num_lines} lines.')
        print(f'{num_lines_failed} lines failed.')
        if verbose:
            print(report)
