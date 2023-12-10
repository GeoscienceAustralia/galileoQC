import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr

import AirGravQC.config as config
import AirGravQC.utility.utility as util

groupName = config.groupName


def checkXYPlan(planPath, measPath, lines=[], planX='', planY='', measX='', measY='', allowance=200.0, maxCounter=14, maxDistance=0, known='', plot_flag=False, verbose=False):
    """
    Reports exceedances of actual horizontal position from planned horizontal
    positions for an airborne survey Whizz database.
    The positions (planX, planY) of the start and end of each planned survey line
    are read from planPath. The measured positions (measX, measY) are read from
    measPath and the perpendicular distance of each from the planned line is
    calculated. If this distance exceeds allowance for maxCounter consecutive
    fiducial, or maxDistance consecutive metres, then an out-of-specification
    exceedance is reported for that line.

    Parameters
    ----------
    planPath : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, with the survey
        positions plan.
    planX : String, optional
        The name of the geoWhizz field or channel containing the planned x position. The
        default is to read the xChannel field name from the Coordinate Frame.
    planY : String, optional
        The name of the geoWhizz field or channel containing the planned y position. The
        default is to read the yChannel field name from the Coordinate Frame.
    measPath : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, with the survey
        measured data.
    measX : String, optional
        The name of the geoWhizz field or channel containing the measured x position. The
        default is to read the xChannel field name from the Coordinate Frame.
    measY : String, optional
        The name of the geoWhizz field or channel containing the measured y position. The
        default is to read the yChannel field name from the Coordinate Frame.
    allowance : Float, optional
        The allowed horizontal distance for the measured line from the planned
        line. If any portion of a measured line is further than this from the
        planned position for more than the maximum allowed number of fids or the
        maximum allowed distance, then the line fails. The default is 200.0.
    maxCounter : Int, optional
        The maximum number of consecutive fids for which an exceedance
        greater than allowance is permitted. If 0, then the constraint is
        ignored. The default is 14.
    maxDistance : Float, optional
        The maximum number of consecutive metres for which an exceedance
        greater than allowance is permitted. If 0, then the constraint is
        ignored. The default is 0.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line.

    Returns
    -------
    None.

    """
    start_x = 0.0
    end_x = 0.0
    start_y = 0.0
    end_y = 0.0

    exceedances_known = False
    this_exc_known = False
    number_known = 0

    planfile = str(planPath)
    measFile = str(measPath)
    
    with h5py.File(planfile, 'r') as fp:
        gPlan = fp[groupName]['Lines']
        if planX == '':
            planX = fp[groupName]['CoordinateFrame'].attrs['XChannel']
        if planY == '':
            planY = fp[groupName]['CoordinateFrame'].attrs['YChannel']

        with h5py.File(measFile, 'r+') as fm:
            gMeas = fm[groupName]['Lines']
            if measX == '':
                measX = fm[groupName]['CoordinateFrame'].attrs['XChannel']
            if measY == '':
                measY = fm[groupName]['CoordinateFrame'].attrs['YChannel']
            numLines = len(gMeas.items())

            message = ''
            num_lines_exceeded = 0
            total_num_excs = 0
            num_lines_unplanned = 0


            if lines == []:
                lines = gMeas.keys()

            for line in lines:
                planLine = f"{gMeas[line].attrs['PlannedLine']:.1f}"
                lineName = util._get_lineName(gMeas[line])
                exceedance_in_line = False
                if planLine in gPlan:
                    xP = np.array(gPlan[planLine][planX])
                    yP = np.array(gPlan[planLine][planY])
                    xM = np.array(gMeas[line][measX])
                    yM = np.array(gMeas[line][measY])
                    max_deviation = 0.0

                    if known != '':
                        exceedances_known = True
                        exc_known = np.array(gMeas[line][known])
                        report_known = -1

                    # rotate to line of x ~ 0 using line direction in radians
                    # if abs(yP[-1] - yP[0]) > epsilon:
                    dirn = np.arctan2((yP[-1] - yP[0]), (xP[-1] - xP[0]))
                    # else:
                    #     dirn = np.pi / 2.0    

                    [x, y] = util._rotateCoords(xM - xP[0], yM - yP[0], dirn)

                    num_fids_in_exceedance = 0
                    exceedance_in_line = False

                    for fid in range(0, len(x)):
                        # There is an exceedance ...
                        if np.abs(y[fid]) > allowance: #x
                            # If a new exceedance, then initialise variables;
                            if num_fids_in_exceedance == 0:
                                start_x = xM[fid]
                                start_y = yM[fid]
                                num_fids_in_exceedance = 1
                                max_deviation = abs(y[fid]) - allowance #x
                                this_exc_known = False

                            # Else increment and update on the current exceedance.
                            else:
                                num_fids_in_exceedance += 1
                                max_deviation = max(np.abs(y[fid]) - allowance, max_deviation) #x
                                if exceedances_known:
                                    if exc_known[fid] > 0:
                                        report_known = exc_known[fid]
                                        this_exc_known = True

                        else:
                            if num_fids_in_exceedance > 0: # the current exceedance has ended
                                end_x = xM[fid]
                                end_y = yM[fid]
                                len_exceedance = util._length([start_x, end_x], [start_y, end_y])[1]
                                if util._exceedance_fail(num_fids_in_exceedance, len_exceedance, maxCounter, maxDistance):
                                    message += f'\nL {lineName} deviates more than {allowance} m for '
                                    message += f'{num_fids_in_exceedance} fids ({len_exceedance:.0f} m), max exceedance = {max_deviation:.0f} m.'
                                    message += f'\n  From ({start_x:.0f} E {start_y:.0f} N) to ({end_x:.0f} E {end_y:.0f} N).'
                                    total_num_excs += 1
                                    if this_exc_known:
                                        number_known += 1
                                        message += f' Known exceedance {report_known}.'
                                        this_exc_known = False
                                    exceedance_in_line = True
                                num_fids_in_exceedance = 0
                else:
                    print(f'Line {lineName} / {planLine} not in plan.')
                    num_lines_unplanned += 1
                if exceedance_in_line:
                    num_lines_exceeded += 1
                    if plot_flag:
                        if abs(np.cos(dirn)) > 0.5:
                            _plot_exceeding_line(x, y, xP, yP, xM, yM, measX, measY, allowance, line, planLine, dirn)
                        else:
                            _plot_exceeding_line(x, y, yP, xP, yM, xM, measY, measX, allowance, line, planLine, dirn)

            message = f'\n{num_lines_exceeded} lines with horizontal exceedances.\n' + message # 5 DEC
            message = f'\n{total_num_excs} horizontal exceedances.\n' + message # 5 DEC
            message = f'\n{num_lines_unplanned} lines not in plan and not checked.\n' + message # 5 DEC
            message = f'\n{number_known} exceedances known in the database.\n' + message # 5 DEC
            print(message)
            if plot_flag:
                plt.show()


def _plot_exceeding_line(x, y, xP, yP, xM, yM, measX, measY, allowance, line, planLine, dirn):
    """
    Plots a standard exceedance figure for checkXYPlan().

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    allowance : TYPE
        DESCRIPTION.
    line : TYPE
        DESCRIPTION.
    planLine : TYPE
        DESCRIPTION.
    dirn : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    fig = plt.figure()
    plot_title = f'Line {line} (plan: {planLine}); Bearing {dirn * 180 / np.pi:.1f} deg E'
    fig.suptitle(plot_title, fontsize=10)

    ax = fig.add_subplot(2,1,1)
    ax.plot(x, y, 'b')
    ax.plot(x, -allowance * np.ones(y.shape), 'r')
    ax.plot(x, allowance * np.ones(y.shape), 'r')
    ax.set_xlim(x[0], x[-1])
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax.set_xlabel('deviation from planned line [m]', fontsize = 8)
    ax.set_ylabel('distance along line', fontsize = 8)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    ax.grid()

    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(xP, yP, color='darkorange', label='Plan', lw=1.5, alpha=0.7)
    ax2.plot(xM, yM, color='blue', label='Measured', lw=1.5, alpha=0.7)
    ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax2.set_xlim(xM[0], xM[-1])
    ax2.set_ylabel(f'{measY} [m]', fontsize=8)
    ax2.set_xlabel(f'{measX} [m]', fontsize=8)
    for label in ax2.get_xticklabels(): label.set_fontsize(6)
    for label in ax2.get_yticklabels(): label.set_fontsize(6)
    ax2.legend(fontsize=8)
    ax2.grid()
    return
