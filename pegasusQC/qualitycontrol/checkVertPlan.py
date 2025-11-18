#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check that vertical exceedances of the aircraft from planned height are within specification.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import h5py
import pathlib
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from scipy import interpolate
import matplotlib.ticker as tkr

import pegasusQC.config as config
import pegasusQC.utility.utility as util

groupName = config.groupName


def checkVertPlan(planPaths, measPath, *, lines=[], planX='', planY='', planZ='', measX='', measY='', measZ='', allowance=30.0, maxCounter=13, maxDistance=0.0, known='', plot_flag=False, verbose=False):
    """
    Reports exceedances of actual vertical position from planned vertical positions
    (stored in a plan file) for an airborne survey Whizz database.

    The positions (`planX`, `planY`, `planY`) of each planned survey line
    are read from a `planPath`. The measured positions (`measX`, `measY`, `measZ`) are read
    from `measPath` and the vertical distance of each from the planned line is
    calculated. If this distance exceeds `allowance` for a distance greater than
    `maxDistance`, then an out-of-specification exceedance is reported for that line.
    If `maxDistance` is less than 1.0, then the test is instead against `maxCounter`
    consecutive positions. The default for `maxCounter` is 13.

    See also `checkClearance`, `checkSafeClearance`, `checkDrape`.

    Parameters
    ----------
    planPath : array of String or pathlib Path

        Names of HDF5 Whizz files, including path and extension, of survey plan.

    measPath : String or pathlib Path

        Name of a HDF5 Whizz file, including path and extension, of measured data.

    lines : Array{String}, optional

        Array of line numbers as strings. Default = [], meaning all lines are checked.

    planX : String, optional

        The name of the geoWhizz field or channel containing the planned x positions. The
        default is to read the xChannel field name from the Coordinate Frame.

    planY : String, optional

        The name of the geoWhizz field or channel containing the planned y positions. The
        default is to read the xChannel field name from the Coordinate Frame.

    measX : String, optional

        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.

    measY : String, optional

        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.

    allowance : Float, optional

        The maximum allowed deviation in height. The default is 30.0.

    maxCounter : Int, optional

        The maximum number of consecutive fids for which an exceedance
        greater than allowance is permitted. The default is 13.

    maxDistance : Float, optional

        The maximum number of consecutive metres for which an exceedance
        greater than allowance is permitted. If 0, then `maxCounter` is
        used instead of `maxDistance`. The default is 0.

    plot_flag : Bool, optional

        If True, plot exceedances for each failed line.

    verbose : Bool, optional

        If True, reports details for each failed line.

    Returns
    -------
    None.

    """
    if type(planPaths) is list:
        if not isinstance(planPaths[0], (str, pathlib.PurePath)):
            print("ERROR - the input parameter needs to be a list of planfile names but it is not.")
            return
    elif isinstance(planPaths, (pathlib.PurePath, str)):
        planPaths = [planPaths]
    else:
        print("ERROR - the first input parameter needs to be a list of planfile names but it is not.")
        return

    measFile = str(measPath)
    if maxDistance > 1.0:
        counting = False
    elif maxCounter > 0:
        counting = True
    else:
        print(f'ERROR. Require maxDistance > 1.0 m or maxCounter > 0')
        print(f'  maxDistance = {maxDistance}, maxCounter = {maxCounter}.')
        return

    exceedances_known = False
    this_exc_known = False
    number_exc_known = 0
    report = ''

    with h5py.File(measFile, 'r') as fm:
        gMeas = fm[groupName]['Lines']
        if measX == '':
            measX = fm[groupName]['CoordinateFrame'].attrs['XChannel']
        if measY == '':
            measY = fm[groupName]['CoordinateFrame'].attrs['YChannel']
        if measZ == '':
            measZ = fm[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        numErrors = int(0)
        numErrLines = 0
        num_lines_unplanned = 0

        if lines == []:
            lines = gMeas.keys()

        for line in lines:
            line_flagged = False
            gLineMeas = gMeas[line]
            planLineInPlan, gpline, planfile = util._get_line_planfiles(planPaths, gLineMeas, verbose=verbose)
            # planLineInPlan, gpline = util.getPlannedLine(gPlan, gLineMeas)
            with h5py.File(planfile, 'r') as fp:
                gPlan = fp[groupName]['Lines']
                # Need to deal with situation where these attributes are not set.
                # 1. Trap for possibility;
                # 2. Tell user to run updateCoordinateFrame
                # 3. End
                planX = fp[groupName]['CoordinateFrame'].attrs['XChannel']
                planY = fp[groupName]['CoordinateFrame'].attrs['YChannel']
                planZ = fp[groupName]['CoordinateFrame'].attrs['AltitudeChannel']

                lineName = util._get_lineName(gLineMeas)

                if planLineInPlan:
                    # print(gPlan, gpline, planX)
                    xP = np.array(gPlan[gpline][planX])
                    yP = np.array(gPlan[gpline][planY])
                    zP = np.array(gPlan[gpline][planZ])
                    xM = np.array(gMeas[line][measX])
                    yM = np.array(gMeas[line][measY])
                    zM = np.array(gMeas[line][measZ])

                    if known != '':
                        exceedances_known = True
                        exc_known = np.array(gMeas[line][known])
                        report_known = -1
                    
                    # Occasionally, some survey lines are flown at constant height,
                    # and the planned drape is all NaNs.
                    good = np.logical_not(np.isnan(zP))
                    if zP[good].size == 0:
                        report += f'Line {line} planned drape is all NaNs; skipping to next line.\n'
                        continue

                    # make life easier by transforming to a 2D problem
                    dirn = np.arctan2((yM[-1] - yM[0]), (xM[-1] - xM[0]))

                    xCentre = min(xP[0], xP[-1])
                    yCentre = min(yP[0], yP[-1])
                    (xm, ym) = util._rotateCoords(xM - xCentre, yM - yCentre, dirn)
                    (xp0, yp0) = util._rotateCoords(xP - xCentre, yP - yCentre, dirn)
                    xp = xp0
                    yp = yp0
                    zp = zP
                    
                    # interpolate (xm, zM) onto (xp, zmp)
                    if abs(xm[-1] - xm[0]) < abs(ym[-1] - ym[0]):
                        report += 'ERROR - expect xms > yms but this is not so.\n'
                    (zmp, zM_trim, xm, addreport) = _interpolateVert(xp, zp, xm, zM)
                    report += addreport
                    if zmp.size == 1 or zM_trim.size == 1 or xm.size == 1:
                        report += f'    interpolation failed on line {line}.\n'
                        continue

                    # calculate the deviation vector
                    z_dev = zM_trim - zmp
                
                    # check vertical deviations
                    # initialise fiducial counter and error report for line, assume no exceedances
                    fid = int(0)
                    exceeding = False
                
                    for one_x in z_dev:
                        fid += 1
                        in_spec = abs(one_x) < allowance
                        
                        if not (exceeding or in_spec):
                            start_fid = fid
                            exceeding = True
                            exc_fids = 0
                            this_exc_known = False
                        elif exceeding and not in_spec:
                            exc_fids += 1
                            if exceedances_known:
                                if exc_known[fid] > 0:
                                    report_known = exc_known[fid]
                                    this_exc_known = True
                        elif exceeding and in_spec:
                            num_fids = exc_fids
                            exceeding = False
                            ex0 = xM[start_fid]
                            ex1 = xM[start_fid + num_fids]
                            ey0 = yM[start_fid]
                            ey1 = yM[start_fid + num_fids]
                            exc_dist = util._displacement2(ex0, ex1, ey0, ey1)
                            if (counting and num_fids > maxCounter) or (not counting and exc_dist > maxDistance):
                                if not line_flagged:
                                    numErrLines += 1
                                    line_flagged = True
                                numErrors += 1
                                max_dev = np.nanmax(abs(z_dev[start_fid:start_fid + num_fids]))
                                report += f'\nL {lineName} deviates more than {allowance:.1f} m for'
                                report += f' {num_fids} fids ({exc_dist:.0f} m),'
                                report += f' max exceedance = {max_dev - allowance:.1f} m.'
                                report += f'\n  From ({ex0:.0f} E, {ey0:.0f} N) to ({ex1:.0f} E, {ey1:.0f} N).'
                                if this_exc_known:
                                    number_exc_known += 1
                                    report += f' Known exceedance {report_known}.'
                                    this_exc_known = False
                    if plot_flag and line_flagged:
                        if abs(np.cos(dirn)) > np.cos(np.pi/4.0):
                            _plot_vert_exceedance(xm, z_dev, xP, zP, xM, zM, measX, measZ, allowance, line, lineName, dirn)
                        else:
                            _plot_vert_exceedance(xm, z_dev, yP, zP, yM, zM, measY, measZ, allowance, line, lineName, dirn)
                else:
                    print(f'Line {lineName} not in plan.')
                    num_lines_unplanned += 1

        if num_lines_unplanned > 0:
            print(f'\n{num_lines_unplanned} lines not in plan and not checked.\n')
        else:
            print(f'\nAll lines in plan and checked.\n')
        print(f'Total number of exceedances = {numErrors} over {numErrLines} erroneous lines.')
        print(f'\n{number_exc_known} exceedances known in the database.\n')
        print(report)
        if plot_flag and numErrLines > 0:
            plt.show()


def _interpolateVert(xbase, ybase, xnew, ynew):
    """
    Interpolates ybase, sampled at xbase, onto the samples xnew.
    These three input arrays are pre-processed to ensure that xbase
    and xnew are monotonically increasing, whilst keeping ybase
    synchronised with xbase, and ynew synchronised with xnew.

    For use interpolating y data along a flown survey line onto positions
    on a planned survey line in preparation for comparison with planned y.

    Parameters
    ----------
    xbase : 1D numpy float array

        The independent variable of the inputs to be interpolated.

    ybase : 1D numpy float array

        The input dependent variable to be interpolated.

    xnew : 1D numpy float array

        The independent variable to interpolate onto.

    ynew : 1D numpy float array

        To be kept synchronised with xnew and returned.

    Returns
    -------
    yout : 1D numpy float array

        The values of ybase interpolated onto xnew.

    ynew : 1D numpy float array

        Synchronised with xnew (by trimming) and returned.

    """
    report = ''
    # clean out 'nan's'
    # good = ~np.isnan(xbase + ybase)
    good = np.logical_not(np.isnan(xbase + ybase))
    # if xbase[good].size < 10:
    #     print(good)
    #     print(xbase)
    #     print(ybase)
    #     print(xbase.size, ybase.size, len(good), xbase[good].size)
    xbase = xbase[good]
    ybase = ybase[good]
    if xbase.size < 10:
        report += f'ERROR - after trimming NaNs, the data vectors are too short for analysis.\n'
        return np.array([0]), np.array([0]), np.array([0]), report

    # ensure ordered in increasing x
    if xbase[1] < xbase[0]:
        xbase = xbase[::-1]
        ybase = ybase[::-1]
    if xnew[1] < xnew[0]:
        xnew = xnew[::-1]
        ynew = ynew[::-1]
    
    # trim base data and store
    keepsml = xbase < xnew[-1]
    keepbig = xbase > xnew[0]
    keep = keepsml & keepbig
    xbase = xbase[keep]
    ybase = ybase[keep]
    if xbase.size < 10:
        report += f'ERROR - trimmed planned data are too few. xnew [{xnew[0]:0.1f}:{xnew[-1]:0.1f}]\n'
        return np.array([0]), np.array([0]), np.array([0]), report

    # trim new data and store
    keepsml = xnew < xbase[-1]
    keepbig = xnew > xbase[0]
    keep = keepsml & keepbig
    xnew = xnew[keep]
    ynew = ynew[keep]
    if xnew.size < 10:
        report += f'ERROR - after trimming measured data, the vectors are too short for analysis.\n'
        return np.array([0]), np.array([0]), np.array([0]), report

    spl = interpolate.splrep(xbase, ybase, k=3, s=0)
    yout = interpolate.splev(xnew, spl)

    return yout, ynew, xnew, report


def _plot_vert_exceedance(xm, z_dev, xP, zP, xM, zM, measX, measZ, allowance, line, planLine, dirn):

    # if xm[1] < xm[0]:
    #     print(f'line {line}: xm decreasing.')
    # if xM[1] < xM[0]:
    #     print(f'line {line}: xM decreasing.')
    fig = plt.figure()
    plot_title = f'Line {line} (plan: {planLine}); Bearing {dirn * 180 / np.pi:.1f} deg E'
    fig.suptitle(plot_title, fontsize=10)
    
    ax = fig.add_subplot(2,1,1)
    ax.plot(xm, z_dev, 'b', lw=0.6)
    ax.plot(xm, -allowance * np.ones(z_dev.shape), 'r')
    ax.plot(xm, allowance * np.ones(z_dev.shape), 'r')
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax.set_xlim(xm[0], xm[-1])
    ax.set_ylabel('deviation from planned drape [m]', fontsize=8)
    ax.set_xlabel('distance along line [m]', fontsize=8)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    ax.grid()

    ax2 = fig.add_subplot(2,1,2)

    # Analysis was done without extrapolation so trim plot accordingly.
    xminlim = np.max([np.min(xM), np.min(xP)])
    xmaxlim = np.min([np.max(xM), np.max(xP)])
    ax2.plot(xP, zP, color='darkorange', label='Plan', lw=1.5, alpha=0.7)
    ax2.plot(xM, zM, color='blue', label='Measured', lw=1.5, alpha=0.7)
    ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax2.set_xlim(xminlim, xmaxlim)
    ax2.set_ylabel(f'{measZ} [m]', fontsize=8)
    ax2.set_xlabel(f'{measX} [m]', fontsize=8)
    for label in ax2.get_xticklabels(): label.set_fontsize(6)
    for label in ax2.get_yticklabels(): label.set_fontsize(6)
    ax2.legend(fontsize=8)
    ax2.grid()
    # The various rotations and swapping of x and y axes occasionally requires this.
    if xM[1] < xM[0]:
        ax2.invert_xaxis()
    return

