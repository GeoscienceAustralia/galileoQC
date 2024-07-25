import numpy as np
from pathlib import Path
import h5py
import collections

import AirGravQC.utility.utility as util
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.config as config

groupName = config.groupName
projectName = config.projectName


def updateLineSpacing(whizzFile, lines=[], x='', y='', trav_spacing=0.0, ctrl_spacing=0.0, verbose=False):
    """
    Writes the mean traverse and mean control line-spacing as (float) attributes
    for the Lines group in `whizzFile`.

    If you know the spacings, set `trav_spacing` and `ctrl_spacing`. Otherwise,
    recommended usage is to set lines to an array of two adjacent lines that are straight
    and of the desired separation. This makes for the fastest analysis.

    Parameters
    ----------
    whizzFile : Path or String
        The Path to, or String name of, the whizz file in HDF5 format.
    lines : String Array, optional
        List of lines to use for the analysis. Default all lines.
    x : String, optional
        The name of the x (easting) channel in `whizz_file`. Default is to use
        the `XChannel` attribute.
    y : String, optional
        The name of the y (northing) channel in `whizz_file`. Default is to use
        the `YChannel` attribute.
    trav_spacing : Float, optional
        The traverse line spacing. Default (0.0) is to estimate from the `x` and `y` data.
    ctrl_spacing : Float, optional
        The control line spacing. Default (0.0) is to estimate from the `x` and `y` data.
    verbose : Bool, optional
        The verbosity of output. Default False.

    Returns
    -------
    Nothing.

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r+') as f:
        lines_group = f[groupName]['Lines']
        if lines == []:
            lines = list(g.keys())
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']

        no_calc_needed = False
        if trav_spacing > 0.0:
            f[groupName]['CoordinateFrame'].attrs['TravSpacing'] = trav_spacing
            no_calc_needed = True
        if ctrl_spacing > 0.0:
            f[groupName]['CoordinateFrame'].attrs['CtrlSpacing'] = ctrl_spacing
            no_calc_needed = True

        if no_calc_needed:
            return
        else:
            travOffset = []
            ctrlOffset = []
            meanTravTrack = 0.0

            lastslope = 0.0
            lastintercept = 0.0
            correlationlimit = 0.6

            for line in lines:
                linegroup = lines_group[line]
                dx = rd.getLineData(linegroup, x)
                dy = rd.getLineData(linegroup, y)

                equation = linregress(dx, dy)
                slope = equation.slope
                intercept = equation.equation
                rvalue = equation.rvalue

                # if slope inconsistent with track, report
                if not _slopeConsistent(linegroup, slope):
                    if verbose:
                        print(f'Line {line} skipped : slope inconsistent with MeanTrack attribute.')

                # skip the first line
                if lastslope == 0.0 and lastintercept == 0.0:
                    continue

                meanslope = np.mean([slope, lastslope])
                separation = (intercept - lastintercept) / np.sqrt(1 + meanslope * meanslope)

                try:
                    trackRad = linegroup.attrs['MeanTrack'] / 180.0 * np.pi
                except:
                    if verbose:
                        print(f'Line {line} skipped : MeanTrack attribute not found.')
                    continue
                try:
                    line_purpose = linegroup.attrs['LinePurpose']
                except:
                    if verbose:
                        print(f'Line {line} skipped : LinePurpose attribute not found.')
                    continue
                if line_purpose == 'Traverse':
                    travOffset.append(np.mean(dy - np.arctan(trackRad) * dx))
                    meanTravTrack += trackRad
                elif line_purpose == 'Control':
                    ctrlOffset.append(np.mean(dy - np.arctan(trackRad) * dx))
                else:
                    if verbose:
                        print(f'Line {line} skipped : LinePurpose not Traverse nor Control.')
                    continue

        meanTravTrack = meanTravTrack / (len(travOffset))
        travOffset = np.sort(travOffset)
        travSpacing = 0
        for idx, offset in enumerate(travOffset[1:]):
            travSpacing += np.cos(meanTravTrack) * (travOffset[idx] - travOffset[idx-1])
        travSpacing = travSpacing / (len(travOffset) - 1)

        ctrlOffset = np.sort(ctrlOffset)
        ctrlSpacing = 0
        for idx, offset in enumerate(ctrlOffset[1:]):
            ctrlSpacing += np.cos(trackRad) * (ctrlOffset[idx] - ctrlOffset[idx-1])
        ctrlSpacing = ctrlSpacing / (len(ctrlOffset) - 1)

        print(f'Spacing: Controls = {ctrlSpacing:.1f}, Traverses = {travSpacing:.1f}')

        return travOffset, travSpacing, meanTravTrack


def _slopeConsistent(linegroup, slope):
    return True


