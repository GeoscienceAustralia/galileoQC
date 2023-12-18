import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import matplotlib.ticker as tkr

import AirGravQC.config as config
import AirGravQC.whizzFiles.pointfiles as gw
import AirGravQC.utility.utility as util

groupName = config.groupName


def checkVertPlan(planPath, measPath, lines=[], planX='', planY='', planZ='', measX='', measY='', measZ='', allowance=30.0, maxCounter=13, maxDistance=0.0, known='', plot_flag=False):
    """
    Reports exceedances of actual vertical position from planned vertical positions
    (stored in a plan file) for an airborne survey Whizz database.

    The positions (`planX`, `planY`, `planY`) of each planned survey line
    are read from `planPath`. The measured positions (`measX`, `measY`, `measZ`) are read
    from `measPath` and the vertical distance of each from the planned line is
    calculated. If this distance exceeds `allowance` for a distance greater than
    `maxDistance`, then an out-of-specification exceedance is reported for that line.
    If `maxDistance` is less than 1.0, then the test is instead against `maxCounter`
    consecutive positions. The default for `maxCounter` is 13.

    See also `checkClearance`, `checkSafeClearance`, `checkDrape`.

    Parameters
    ----------
    planPath : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, of survey plan.
    measPath : String or pathlib.PosixPath
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

    Returns
    -------
    None.

    TODO: If maxCounter unspecified, and maxDistance is specified, test against maxDistance.
          Use _exceedance_fail() to do this.

    """
    planfile = str(planPath)
    measFile = str(measPath)
    if maxDistance > 1.0:
        counting = False
    elif maxCounter > 0:
        counting = True
    else:
        print(f'ERROR. Require maxDistance > 1.0 m or maxCounter > 0')
        print(f'  maxDistance = {maxDistance}, maxCounter = {maxCounter}.')
    
    exceedances_known = False
    this_exc_known = False
    number_exc_known = 0

    with h5py.File(planfile, 'r') as fp:
        gPlan = fp[groupName]['Lines']
        if planX == '':
            planX = fp[groupName]['CoordinateFrame'].attrs['XChannel']
        if planY == '':
            planY = fp[groupName]['CoordinateFrame'].attrs['YChannel']
        if planZ == '':
            planZ = fp[groupName]['CoordinateFrame'].attrs['AltitudeChannel']

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
            report = ''

            if lines == []:
                lines = gMeas.keys()

            for line in lines:
                line_flagged = False
                lineName = util._get_lineName(gMeas[line])
                planLine = f"{gMeas[line].attrs['PlannedLine']:.1f}"
                if planLine in gPlan: # 5 DEC
                    xP = np.array(gPlan[planLine][planX])
                    yP = np.array(gPlan[planLine][planY])
                    zP = np.array(gPlan[planLine][planZ])
                    xM = np.array(gMeas[line][measX])
                    yM = np.array(gMeas[line][measY])
                    zM = np.array(gMeas[line][measZ])

                    if known != '':
                        exceedances_known = True
                        exc_known = np.array(gMeas[line][known])
                        report_known = -1
                    
                    # make life easier by transforming to a 2D problem
                    dirn = np.arctan2((yM[-1] - yM[0]), (xM[-1] - xM[0]))
                    (xm, ym) = util._rotateCoords(xM - xP[0], yM - yP[0], dirn)
                    (xp0, yp0) = util._rotateCoords(xP - xP[0], yP - yP[0], dirn)
                    xp = xp0
                    yp = yp0
                    zp = zP
                    
                    # interpolate (xm, zM) onto (xp, zmp)
                    if abs(xm[-1] - xm[0]) < abs(ym[-1] - ym[0]):
                        print('ERROR - expect xms > yms but this is not so.')
                    (zmp, zM_trim) = gw.interpolateLine(xp, zp, xm, zM, plot_flag=False)
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
                        if abs(np.cos(dirn)) > 0.5:
                            _plot_vert_exceedance(xm, z_dev, xP, zP, xM, zM, measX, measZ, allowance, line, planLine, dirn)
                        else:
                            _plot_vert_exceedance(xm, z_dev, yP, zP, yM, zM, measY, measZ, allowance, line, planLine, dirn)
                else:
                    print(f'Line {lineName} not in plan.')
                    num_lines_unplanned += 1

        print(f'\n{num_lines_unplanned} lines not in plan and not checked.\n')
        print(f'Total number of exceedances = {numErrors} over {numErrLines} erroneous lines.')
        print(f'\n{number_exc_known} exceedances known in the database.\n')
        print(report)
        if plot_flag and numErrLines > 0:
            plt.show()
 

def _plot_vert_exceedance(xm, z_dev, xP, zP, xM, zM, measX, measZ, allowance, line, planLine, dirn):
    fig = plt.figure()
    plot_title = f'Line {line} (plan: {planLine}); Bearing {dirn * 180 / np.pi:.1f} deg E'
    fig.suptitle(plot_title, fontsize=10)
    
    ax = fig.add_subplot(2,1,1)
    ax.plot(xm[1:], z_dev, 'b', lw=0.6)
    ax.plot(xm[1:], -allowance * np.ones(z_dev.shape), 'r')
    ax.plot(xm[1:], allowance * np.ones(z_dev.shape), 'r')
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax.set_xlim(xm[0], xm[-1])
    ax.set_ylabel('deviation from planned drape [m]', fontsize=8)
    ax.set_xlabel('distance along line [m]', fontsize=8)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    ax.grid()

    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(xP, zP, color='darkorange', label='Plan', lw=1.5, alpha=0.7)
    ax2.plot(xM, zM, color='blue', label='Measured', lw=1.5, alpha=0.7)
    ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax2.set_xlim(xM[0], xM[-1])
    ax2.set_ylabel(f'{measZ} [m]', fontsize=8)
    ax2.set_xlabel(f'{measX} [m]', fontsize=8)
    for label in ax2.get_xticklabels(): label.set_fontsize(6)
    for label in ax2.get_yticklabels(): label.set_fontsize(6)
    ax2.legend(fontsize=8)
    ax2.grid()
    return
