import numpy as np
import h5py

import AirGravQC.config as config

groupName = config.groupName


def checkGaps(whizzFile, maxGapSec=0.0, maxNumGaps=0):
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

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        _reportGaps(g, maxGapSec, maxNumGaps)
        
        
def _reportGaps(group, maxGapSec=0.0, maxNumGaps=0):
    """
    Checks every dataset for each channel and each survey line in the HDF5 group
    for gaps, and reports all gaps found.

    Parameters
    ----------
    group : HDF5 Whizz file 'Lines' group
        The group containing the survey line data.
    maxGapSec :  Float, optional
        The largest allowed gap measured in seconds. Default 0.0
    maxNumGaps : Integer, optional
        The maximum number of gaps allowed on any survey line. Default 0

    Returns
    -------
    None.

    """
    lineGroups = list(group.values())
    channelNames = list(lineGroups[0].keys())
    num_channels = len(channelNames)
    num_lines_failed = 0
    total_num_lines = 0
    message = ''

    for line in group.keys():
        total_num_lines += 1
        gaps_on_line = 0
        lineNo = line
        lineText = f'Line {lineNo}'
        for channel in channelNames:
            numberMissing = np.count_nonzero(np.isnan(group[line][channel]))
            if numberMissing > 0:
                lineText += f'\n    {channel}, nans: {numberMissing}'
                gaps_on_line += 1
        if gaps_on_line > 0:
            num_lines_failed += 1
            message += lineText + '\n'
    print(f'Checking for all gaps in all {num_channels} channels on all {total_num_lines} lines.')
    print(message)
    print(f'{num_lines_failed} lines failed.')
