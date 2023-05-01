import h5py


def checkGNSS(whizzFile, num_sats, pdop, vdop, hdop, nsats_min=4, max_pdop=6, max_vdop=4, max_hdop=4):
    """
    Checks that the data in a whizzFile meets the requirements for the minimum
    number of satellites, and maximum PDOP, VDOP and HDOP.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    num_sats : String
        Name of the channel containing the number of satellites visible for
        each measurement.
    pdop : String
        Name of the channel containing the PDOP for each measurement.
    vdop : String
        Name of the channel containing the VDOP for each measurement.
    hdop : String
        Name of the channel containing the HDOP for each measurement.
    nsats_min : Integer, optional
        The minimum number of satellites required, default 4.
    max_pdop : Integer, optional
        The maximum PDOP allowed, default 6
    max_vdop : Integer, optional
        The maximum VDOP allowed, default 6
    max_hdop : Integer, optional
        The maximum HDOP allowed, default 6

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']

        error_count = 0
        for line in g.keys():
            min_nsats_data = np.nanmin(np.array(g[line][num_sats]))
            max_pdop_data = np.nanmin(np.array(g[line][pdop]))
            max_vdop_data = np.nanmin(np.array(g[line][vdop]))
            max_hdop_data = np.nanmin(np.array(g[line][hdop]))
            if min_nsats_data < nsats_min:
                print(f'Line {line} fail: min sats = {min_nsats_data:.0f} < {nsats_min}')
                error_count += 1
            if max_pdop_data > max_pdop:
                print(f'Line {line} fail: max pdop = {max_pdop_data:.0f} > {max_pdop}')
                error_count += 1
            if max_vdop_data > max_vdop:
                print(f'Line {line} fail: max vdop = {max_vdop_data:.0f} > {max_vdop}')
                error_count += 1
            if max_hdop_data > max_hdop:
                print(f'Line {line} fail: max hdop = {max_hdop_data:.0f} > {max_hdop}')
                error_count += 1
        print(f'In {projName}, checked num sats, PDOP, VDOP and HDOP. Found {error_count} errors.')
        

def checkHeading(whizzFile, nominalHeading, x='', y='', tolerance=10.0):
    """
    Checks heading in degrees is within +/- tolerance (in degrees) of nominal (in degrees).

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    x : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    y : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    upper_limit = nominalHeading + tolerance
    lower_limit = nominalHeading - tolerance
    print(f'limits: {lower_limit}, {upper_limit}')

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']

        for line in g.keys():
            dx = np.diff(np.array(g[line][x]))
            dy = np.diff(np.array(g[line][y]))
            heading = np.arctan(dx / dy) * 180.0 / np.pi
            mean_heading = np.mean(heading)
            exceeds = any(h > upper_limit for h in heading)
            exceeds = exceeds or any(h < lower_limit for h in heading)
            
            if exceeds:
                print(f'Line {line}: heading limit exceeded. Mean {mean_heading}')
                fig = plt.figure()
                fig.suptitle(f'Heading Check Line {line}', fontsize=10)
                fig.subplots_adjust(top=0.85)
                
                ax = fig.add_subplot(1,1,1)
                ax.plot(np.array(g[line][x])[1:], heading, 'b', mfc='w')
                plt.ylabel('Estimated heading [deg]', fontsize = 6)
                plt.grid(True)
                for label in ax.get_xticklabels(): label.set_fontsize(6)
                for label in ax.get_yticklabels(): label.set_fontsize(6)
                plt.show()
    return
    

def checkLineLengths(whizzFile, min_len=50.0, measX='', measY=''):
    """
    Checks that all lines in whizzFile are at least min_len km long.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    min_len : TYPE, optional
        The minimum allowed line length in km. The default is 50.0.
    measX : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    measY : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.

    Returns
    -------
    None.

    """
    measFile = str(whizzFile)
    
    with h5py.File(measFile, 'r') as f:
        gMeas = f[groupName]['Lines']
        if measX == '':
            measX = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if measY == '':
            measY = f[groupName]['CoordinateFrame'].attrs['YChannel']
        
        num_failed_lines = 0
        for line in gMeas.keys():
            xM = np.array(gMeas[line][measX])
            yM = np.array(gMeas[line][measY])
            line_length = _displacement2(xM[0], xM[-1], yM[0], yM[-1])
            if line_length < min_len * 1000.0:
                num_failed_lines += 1
                print(f'Line {line} length = {line_length:.1f} less than allowed min {min_len*1000.0:.1f}')
        print(f'Number failed lines = {num_failed_lines}')
            

def checkOverlaps(whizzFile, min_overlap = 7.6, lines = [], verbose=False, plot_flag=False):
    """
    For every line in the file whizzFile, calculate the overlap with each other
    line that has the same prefix. Plot a map. Report overlaps.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    min_overlap : Float, optional
        The minimum overlap distance in km, default 7.6 km.
    lines : String list, optional.
        The line numbers to be checked. Default is all lines in the whizzFile.
    verbose : Bool, optional
        If True, report status of all overlaps, else only report errors. Default False.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line.

    Returns
    -------
    None.

    """
    from matplotlib.ticker import StrMethodFormatter
    measFile = str(whizzFile)
    if plot_flag:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    
    with h5py.File(measFile, 'r') as f:
        gMeas = f[groupName]['Lines']
        east = f[groupName]['CoordinateFrame'].attrs['XChannel']
        nrth = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if lines == []:
            lines = list(gMeas.keys())
        num_coinc_lines = 0
        report = ''

        for idx1 in range(0, len(lines)): #line1 in lines: #gMeas.keys():
            line1 = lines[idx1]
            line1_plan = gMeas[line1].attrs['PlannedLine']

            for idx2 in range (idx1 + 1, len(lines)):#line2 in gMeas.keys():
                line2 = lines[idx2]
                line2_plan = gMeas[line2].attrs['PlannedLine']
                # if the second line isn't the first line but has the same planned line no.
                if line1 != line2 and line1_plan == line2_plan:
                    num_coinc_lines += 1
                    # extract positions
                    n1 = np.array(gMeas[line1][nrth])
                    e1 = np.array(gMeas[line1][east])
                    n2 = np.array(gMeas[line2][nrth])
                    e2 = np.array(gMeas[line2][east])
                    # get line direction in radians
                    dirn = np.arctan((e1[-1] - e1[0])/(n1[-1] - n1[0]))
                    
                    # make life easier by transforming to a 2D problem
                    n0 = n1[0]
                    e0 = e1[0]
                    (y1, x1) = _rotateCoords(e1 - e0, n1 - n0, -dirn)
                    (y2, x2) = _rotateCoords(e2 - e0, n2 - n0, -dirn)

                    # find ends of lines
                    min1 = np.nanmin(x1)
                    max1 = np.nanmax(x1)
                    min2 = np.nanmin(x2)
                    max2 = np.nanmax(x2)
                    overlap = 0.0
                    whole_line = 0
                    
                    # find overlap length if any
                    if (min2 < min1 and min1 < max2) and not (min2 < max1 and max1 < max2):
                        overlap = abs(max2 - min1)
                    elif (min2 < max1 and max1 < max2) and not (min2 < min1 and min1 < max2):
                        overlap = abs(min2 - max1)
                    elif (min2 < max1 and max1 < max2) and (min2 < min1 and min1 < max2):
                        overlap = abs(max1 - min1)
                        whole_line = 1
                    elif (min1 < max2 and max2 < max1) and (min1 < min2 and min2 < max1):
                        overlap = abs(max2 - min2)
                        whole_line = 2

                    if plot_flag:
                        plotline1, = ax.plot(e1, n1, color='blue', lw=0.2)
                    
                    if whole_line == 1:
                        if verbose:
                            report += f'  OK: All of line {line1} in {line2}: length {overlap:.0f} m.\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='green', lw=0.2)
                    elif whole_line == 2:
                        if verbose:
                            report += f'  OK: All of line {line2} in {line1}: length {overlap:.0f} m.\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='green', lw=0.2)
                    elif overlap < 1:
                        if verbose:
                            report += f'  OK: Non-overlapping lines {line1}, {line2}: overlap {overlap:.0f} m < {min_overlap * 1000}.\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='green', lw=0.2)
                    elif overlap < min_overlap * 1000.0:
                        report += f'\n  ERROR: Repeat line {line1}, {line2}: overlap {overlap:.0f} m < {min_overlap * 1000}.\n\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='red', lw=0.2)
                    else:
                        if verbose:
                            report += f'  OK: Repeat line {line1}, {line2}: overlap by {overlap:.0f} m.\n'
                        if plot_flag:
                            plotline2, = ax.plot(e2, n2, color='green', lw=0.2)

    print(f'{num_coinc_lines} coincident lines found.')
    print(report)
    if num_coinc_lines > 0 and plot_flag:
        ax.set_aspect('equal')
        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        plt.xlabel('X [m]', fontsize = 10)
        plt.ylabel('Y [m]', fontsize = 10)
        plt.suptitle('Overlap Map', fontsize = 12)
        plt.title('[1st line blue, accepted line green, error line red]', fontsize = 10)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(10)
        for label in ax.get_yticklabels(): label.set_fontsize(10)
        plt.show()


def checkXYPlan(planPath, measPath, lines=[], planX='', planY='', measX='', measY='', allowance=200.0, maxCounter=14, maxDistance=0, known='', plot_flag=False):
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

                    # rotate to line of x ~ 0 using line direction in radians
                    # if abs(yP[-1] - yP[0]) > epsilon:
                    dirn = np.arctan2((yP[-1] - yP[0]), (xP[-1] - xP[0]))
                    # else:
                    #     dirn = np.pi / 2.0    

                    [x, y] = _rotateCoords(xM - xP[0], yM - yP[0], dirn)

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
                                        this_exc_known = True

                        else:
                            if num_fids_in_exceedance > 0: # the current exceedance has ended
                                end_x = xM[fid]
                                end_y = yM[fid]
                                len_exceedance = _dist([start_x, end_x], [start_y, end_y])[1]
                                if _exceedance_fail(num_fids_in_exceedance, len_exceedance, maxCounter, maxDistance):
                                    message += f'\nL {line} deviates more than {allowance} m for '
                                    message += f'{num_fids_in_exceedance} fids ({len_exceedance:.0f} m), max exceedance = {max_deviation:.0f} m.'
                                    message += f'\n  From ({start_x:.0f} E {start_y:.0f} N) to ({end_x:.0f} E {end_y:.0f} N).'
                                    total_num_excs += 1
                                    if this_exc_known:
                                        number_known += 1
                                        message += ' Known exceedance.'
                                        this_exc_known = False
                                    exceedance_in_line = True
                                num_fids_in_exceedance = 0
                else:
                    print(f'Line {line} / {planLine} not in plan.')
                    num_lines_unplanned += 1
                if exceedance_in_line:
                    num_lines_exceeded += 1
                    if plot_flag:
                        _plot_exceeding_line(x, y, xP, yP, xM, yM, measX, measY, allowance, line, planLine, dirn)

            message = f'\n{num_lines_exceeded} lines with horizontal exceedances.\n' + message # 5 DEC
            message = f'\n{total_num_excs} horizontal exceedances.\n' + message # 5 DEC
            message = f'\n{num_lines_unplanned} lines not in plan and not checked.\n' + message # 5 DEC
            message = f'\n{number_known} exceedances known in the database.\n' + message # 5 DEC
            print(message)
            if plot_flag:
                plt.show()


def _exceedance_fail(num_fids_in_exceedance, len_exceedance, maxCounter, maxDistance):
    """
    Given a specification, maxCounter, on number of fids and, maxDistance,
    on distance, checks to see if either  num_fids_in_exceedance > maxCounter
    or len_exceedance > maxDistance. If either is True, returns True.

    Parameters
    ----------
    num_fids_in_exceedance : Integer
        DESCRIPTION.
    len_exceedance : Float
        DESCRIPTION.
    maxCounter : Integer
        DESCRIPTION.
    maxDistance : Float
        DESCRIPTION.

    Returns
    -------
    Bool.

    """
    if maxCounter < 1 and maxDistance < 1:
        return False
    if maxCounter > 0:
        if num_fids_in_exceedance > maxCounter:
            return True
    if maxDistance > 0:
        if len_exceedance > maxDistance:
            return True
    return False


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
    ax.set_xlabel('deviation from planned line [m]', fontsize = 8)
    ax.set_ylabel('distance along line', fontsize = 8)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    ax.grid()

    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(xP, yP, color='magenta', label='Plan', lw=1.5, alpha=0.7)
    ax2.plot(xM, yM, color='cyan', label='Measured', lw=1.5, alpha=0.7)
    ax2.set_xlim(xM[0], xM[-1])
    ax2.set_ylabel(f'{measY} [m]', fontsize=8)
    ax2.set_xlabel(f'{measX} [m]', fontsize=8)
    for label in ax2.get_xticklabels(): label.set_fontsize(6)
    for label in ax2.get_yticklabels(): label.set_fontsize(6)
    ax2.legend(fontsize=8)
    ax2.grid()
    return

   
def checkVertPlan(planPath, measPath, lines=[], planX='', planY='', planZ='', measX='', measY='', measZ='', allowance=30.0, maxCounter=13, maxDistance=0.0, known='', plot_flag=False):
    """
    Reports exceedances of actual vertical position from planned vertical positions
    for an airborne survey Whizz database.
    The positions (`planX`, `planY`, `planY`) of each planned survey line
    are read from `planPath`. The measured positions (`measX`, `measY`, `measZ`) are read from
    `measPath` and the vertical distance of each from the planned line is
    calculated. If this distance exceeds `allowance` for a distance greater than
    `maxDistance`, then an out-of-specification exceedance is reported for that line.
    If `maxDistance` is less than 1.0, then the test is instead against `maxCounter`
    consecutive positions. The default for `maxCounter` is 13.

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
        DESCRIPTION. The default is 30.0.
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
                    
                    # make life easier by transforming to a 2D problem
                    dirn = np.arctan2((yM[-1] - yM[0]), (xM[-1] - xM[0]))
                    (xm, ym) = _rotateCoords(xM - xP[0], yM - yP[0], dirn)
                    (xp0, yp0) = _rotateCoords(xP - xP[0], yP - yP[0], dirn)
                    xp = xp0
                    yp = yp0
                    zp = zP
                    
                    # interpolate (xm, zM) onto (xp, zmp)
                    if abs(xm[-1] - xm[1]) < abs(ym[-1] - ym[0]):
                        print('ERROR - expect xms > yms but this is not so.')
                    (zmp, zM_trim) = mhd.interpolateLine(xp, zp, xm, zM, plot_flag=False)
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
                                    this_exc_known = True
                        elif exceeding and in_spec:
                            num_fids = exc_fids
                            exceeding = False
                            ex0 = xM[start_fid]
                            ex1 = xM[start_fid + num_fids]
                            ey0 = yM[start_fid]
                            ey1 = yM[start_fid + num_fids]
                            exc_dist = _displacement2(ex0, ex1, ey0, ey1)
                            if (counting and num_fids > maxCounter) or (not counting and exc_dist > maxDistance):
                                if not line_flagged:
                                    numErrLines += 1
                                    line_flagged = True
                                numErrors += 1
                                max_dev = np.nanmax(abs(z_dev[start_fid:start_fid + num_fids]))
                                report += f'\nL {line} deviates more than {allowance:.1f} m for'
                                report += f' {num_fids} fids ({exc_dist:.0f} m),'
                                report += f' max exceedance = {max_dev - allowance:.1f} m.'
                                report += f'\n  From ({ex0:.0f} E, {ey0:.0f} N) to ({ex1:.0f} E, {ey1:.0f} N).'
                                if this_exc_known:
                                    number_exc_known += 1
                                    report += ' Known exceedance.'
                                    this_exc_known = False
                    if plot_flag and line_flagged:
                        _plot_vert_exceedance(xm, z_dev, xP, zP, xM, zM, measX, measZ, allowance, line, planLine, dirn)
                else:
                    print(f'Line {line} not in plan.')
                    num_lines_unplanned += 1

        print(f'\n{num_lines_unplanned} lines not in plan and not checked.\n')
        print(f'Total number of exceedances = {numErrors} over {numErrLines} erroneous lines.')
        print(f'\n{number_exc_known} exceedances known in the database.\n')
        print(report)
        if plot_flag:
            plt.show()
                                       
   
def _rotateCoords(x, y, angle):
    """
    Rotates Cartesian vectors x and y, by `angle` to xr and yr.

    Parameters
    ----------
    x : numpy Float 1D array
        DESCRIPTION.
    y : numpy Float 1D array
        DESCRIPTION.
    angle : Float
        Angle in radians by which coordinates are to be rotated.

    Returns
    -------
    xr : numpy Float 1D array
        DESCRIPTION.
    yr : numpy Float 1D array
        DESCRIPTION.

    """
    xr = np.cos(angle) * x + np.sin(angle) * y
    yr = np.sin(angle) * x - np.cos(angle) * y
    return xr, yr
 

def _plot_vert_exceedance(xm, z_dev, xP, zP, xM, zM, measX, measZ, allowance, line, planLine, dirn):
    fig = plt.figure()
    plot_title = f'Line {line} (plan: {planLine}); Bearing {dirn * 180 / np.pi:.1f} deg E'
    fig.suptitle(plot_title, fontsize=10)
    
    ax = fig.add_subplot(2,1,1)
    ax.plot(xm[1:], z_dev, 'b', lw=0.6)
    ax.plot(xm[1:], -allowance * np.ones(z_dev.shape), 'r')
    ax.plot(xm[1:], allowance * np.ones(z_dev.shape), 'r')
    ax.set_xlim(xm[0], xm[-1])
    ax.set_ylabel('deviation from planned drape [m]', fontsize=8)
    ax.set_xlabel('distance along line [m]', fontsize=8)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    ax.grid()

    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(xP, zP, color='magenta', label='Plan', lw=1.5, alpha=0.7)
    ax2.plot(xM, zM, color='cyan', label='Measured', lw=1.5, alpha=0.7)
    ax2.set_xlim(xM[0], xM[-1])
    ax2.set_ylabel(f'{measZ} [m]', fontsize=8)
    ax2.set_xlabel(f'{measX} [m]', fontsize=8)
    for label in ax2.get_xticklabels(): label.set_fontsize(6)
    for label in ax2.get_yticklabels(): label.set_fontsize(6)
    ax2.legend(fontsize=8)
    ax2.grid()
    return


def checkClearance(whizzFile, nominalClearance, clearance_chan='', altitude_chan='', terrain_chan='', allowance=20.0, maxDistance=1000.0, xChannel='', yChannel='', only_low=False):
    """
    Checks the data from Whizz HDF5 file for height exceedances against a specification
    requiring heights to be within a relative range (allowance) about a nominal
    value (nominalClearance) over a particular distance (maxDistance).
    
    The clearance (calculated as altitude - terrain if clearance_chan=''), must be
    within +/- allowance of nominalClearance.
    For each line, the maximum absolute deviation of clearance from nominalClearance is
    reported. 
    If any samples along the line exceed the allowance, then profiles are plotted
    for visual analysis. No use is made of maxDistance (yet).
    The positions (xChannel and yChannel) are only used for plotting.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    clearance_chan : String, optional
        The name of the terrain clearance field in the Lines group
        of the Whizz HDF5 file.
    altitude_chan : String, optional
        The name of the absolute altitude or height field in the Lines group
        of the Whizz HDF5 file.
    terrain_chan : String, optional
        The name of the absolute terrain or DTM height field in the Lines group
        of the Whizz HDF5 file.
    nominalClearance : Float
        The planned clearance between above the terrain in metres.
    allowance : Float, optional
        The absolute maximum deviation allowed from the planned clearance in
        metres. The default is 20.0.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    only_low : String, optional
        If True, then only check for exceedances too close to the ground ("safety
        check"). The default is False.

    Returns
    -------
    None.

    """
    
    filename = str(whizzFile)
    if clearance_chan == '':
        if altitude_chan == '' or terrain_chan == '':
            print('ERROR - either clearance_chan, or both altitude_chan and terrain_chan must be specified.')

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        num_failed_lines = 0
        report = ''

        for line in g.keys():
            if clearance_chan == '':
                alt = np.array(g[line][altitude_chan])
                dtm = np.array(g[line][terrain_chan])
                clearance = alt - dtm
            else:
                clearance = np.array(g[line][clearance_chan])
            deviation = nominalClearance - clearance
            if only_low:
                maxDeviation = np.max(deviation)
            else:
                maxDeviation = np.max(abs(deviation))
            # print(f'Line {line}: Max deviation from {nominalClearance:.0f} m clearance = {maxDeviation:.0f} m.')
            if maxDeviation > allowance:
                num_failed_lines += 1
                report += f'\nClearance deviation of {maxDeviation:.0f} m on line {line}'
                x = np.array(g[line][xChannel])
                y = np.array(g[line][yChannel])
                distance = _dist(x, y)
                fig = plt.figure()

                ax = fig.add_subplot(2,1,1)
                if clearance_chan == '':
                    ax.plot(distance, alt)
                    ax.plot(distance, dtm)
                    plt.legend([altitude_chan, terrain_chan], fontsize=8)
                else:
                    ax.plot(distance, clearance)
                    plt.legend([clearance_chan], fontsize=8)
                plt.xlabel('distance along line [m]', fontsize = 6)
                plt.ylabel('height', fontsize = 6)
                if only_low:
                    plotTitle = projName + ': Minimum Clearance Check ' + line
                else:
                    plotTitle = projName + ': Absolute Clearance Check ' + line
                plt.title(plotTitle, fontsize = 8)
                plt.grid(True)
                for label in ax.get_xticklabels(): label.set_fontsize(6)
                for label in ax.get_yticklabels(): label.set_fontsize(6)
                ax2 = fig.add_subplot(2,1,2)
                ax2.plot(distance, deviation)
                ax2.plot(distance, allowance * np.ones(y.shape), 'r')
                ax2.plot(distance, -allowance * np.ones(y.shape), 'r')
                plt.ylabel('deviation', fontsize = 6)
                plt.xlabel('distance along line [m]', fontsize = 6)
                plt.grid(True)
                for label in ax2.get_xticklabels(): label.set_fontsize(6)
                for label in ax2.get_yticklabels(): label.set_fontsize(6)
        print(f'Number of failed lines = {num_failed_lines}.')
        print(report)
        if num_failed_lines > 0:
            plt.show()
        

def checkDrape(whizzFile, altitude, drape, warningClearance = 20.0, xChannel = '', yChannel = ''):
    """
    Checks actual altitude flown against planned drape in drape channel.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    altitude : TYPE
        The name of the geoWhizz field or channel containing the measured altitudes.
    drape : TYPE
        The name of the geoWhizz field or channel containing the measured drape heights.
    warningClearance : Float, optional
        DESCRIPTION. The default is 20.0.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        
        for line in g.keys():
            alt = g[line][altitude]
            drape = g[line][drape]
            clearance = np.abs(alt[()] - drape[()])
            maxDeviation = np.max(clearance)
            print(f'Drape {maxDeviation} on line {line}')
            if maxDeviation > warningClearance:
                print(f'Drape deviation of {maxDeviation} on line {line}')
                x = g[line][xChannel]
                y = g[line][yChannel]
                distance = _dist(x, y)
                fig = plt.figure()
                ax = fig.add_subplot(2,1,1)
                ax.plot(distance, alt)
                ax.plot(distance, drape)
                plt.xlabel('distance along line [m]', fontsize = 6)
                plt.ylabel('height', fontsize = 6)
                plotTitle = projName + ': Clearance Check ' + line
                plt.title(plotTitle, fontsize = 8)
                plt.legend([altitude, drape])
                plt.grid(True)
                for label in ax.get_xticklabels(): label.set_fontsize(6)
                for label in ax.get_yticklabels(): label.set_fontsize(6)
                ax2 = fig.add_subplot(2,1,2)
                ax2.plot(distance, clearance)
                plt.ylabel('deviation', fontsize = 6)
                plt.xlabel('distance along line [m]', fontsize = 6)
                plt.grid(True)
                for label in ax2.get_xticklabels(): label.set_fontsize(6)
                for label in ax2.get_yticklabels(): label.set_fontsize(6)
                plt.show()
        

def checkIntersectionHeights(whizzFile, controls, max_allowed_height=10.0):
    """
    Checks for .
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    max_allowed_height : Float, optional
        .

    Returns
    -------
    None

    """
    data_is_good = True
    report = ''
    filename = str(whizzFile)
    num_intersections_checked = 0
    num_failed_intersections = 0
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        zChannel = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        lines = list(g.keys())
        for linec in controls:
            x_ctrl = np.array(g[linec][xChannel])
            y_ctrl = np.array(g[linec][yChannel])
            z_ctrl = np.array(g[linec][zChannel])
            bear_ctrl = _calc_bearing(x_ctrl, y_ctrl)
            (y_ctrl, x_ctrl) = _rotateCoords(x_ctrl, y_ctrl, bear_ctrl)
            for linet in lines:
                if linet == linec:
                    continue
                x_trav = np.array(g[linet][xChannel])
                y_trav = np.array(g[linet][yChannel])
                z_trav = np.array(g[linet][zChannel])
                (y_trav, x_trav) = _rotateCoords(x_trav, y_trav, bear_ctrl)
                if _lines_cross(x_ctrl, y_ctrl, x_trav, y_trav):
                    # print(f'bearings: {bearingt}, {bearingt} -- {np.abs(np.cos(bearingc - bearingt))}')
                    dh = _intersection_height(x_trav, y_trav, z_trav, x_ctrl, y_ctrl, z_ctrl, bear_ctrl)
                    num_intersections_checked += 1
                    if dh > max_allowed_height:
                        num_failed_intersections += 1
                        # print(dh)
                        report += f'\n  {linet} : {linec} intersection height difference = {dh:.1f} > {max_allowed_height:.1f}.'
                        data_is_good = False
                # else:
                #     report += f'\n  {linet} : {linec} un-tested since not perpendicular.'
    if data_is_good:
        report += f'All {num_intersections_checked} intersection heights were less than {max_allowed_height:.1f}.'
    else:
        tmpstr = f'Of {num_intersections_checked} intersections checked'
        tmpstr += f', {num_failed_intersections} exceeded the '
        tmpstr += f'{max_allowed_height} m allowed height difference.\n'
        report = tmpstr + report
    print(report)
    return


def _intersection_height(x_trav, y_trav, z_trav, x_ctrl, y_ctrl, z_ctrl, bearingc):
    """
    """

    y = np.abs(y_trav - _mean_1std(y_ctrl))
    it_arr = np.where(y == y.min())
    it = it_arr[0]

    x = np.abs(x_ctrl - x_trav[it])
    ic_arr = np.where(x == x.min())
    ic = ic_arr[0]

    return np.abs(z_trav[it] - z_ctrl[ic])[0]


def _mean_1std(x):
    """
    Calculate the mean of the values in x that fall in the range of +/- stdev(x).
    This is a simplistic "mean excluding outliers".
    """
    mean1 = np.mean(x)
    std1 = np.std(x)
    idx = np.argwhere(np.logical_and(x > (mean1 - std1), x < (mean1 + std1)))
    return np.mean(x[idx])


def _lines_cross(x_ctrl, y_ctrl, x_trav, y_trav):
    """
    Assumes x,y coordinates rotated so that y_ctrl is approximately constant (y-axis
    approximately parallel to the control line).
    """

    # The traverse is 'north' or 'south' of the control line
    if y_trav.min() > y_ctrl.max() or y_trav.max() < y_ctrl.min():
        return False

    # We don't allow shallow angle crossovers (and won't count lines that were supposed to be parallel!).
    min_cosine = 0.1
    bear_ctrl = _calc_bearing(x_ctrl, y_ctrl)
    bear_trav = _calc_bearing(x_trav, y_trav)
    if np.abs(np.cos(bear_ctrl - bear_trav)) < min_cosine:
        return False

    # Find the index to the traverse sample whose y coordinate is closest to the mean control y coordinate
    y = np.abs(y_trav - _mean_1std(y_ctrl))
    it_arr = np.where(y == y.min())
    if np.size(it_arr):
        itc = it_arr[0]
    else:
        return False

    # Now can check if the traverse is 'east' or 'west' of the control line.
    if x_ctrl.min() > x_trav[itc] or x_ctrl.max() < x_trav[itc]:
        return False

    return True


def _calc_bearing(x, y):
    """
    arctan(mean(diff(x) / mean(diff(y))))
    """
    return np.arctan(_mean_1std(np.diff(x)) / _mean_1std(np.diff(y)))


def checkSpeeds(whizzFile, xChannel='', yChannel='', tChannel='', vel_north='', vel_east='', nominalSpeed=60.0, 
    maxDuration=0.0, maxDistance=0.0, allowance=0.1, minSafeSpeed=42.0, plot_flag=False):
    """
    Checks the data from Whizz HDF5 file for speed exceedances against a specification
    requiring ground speeds to be within a relative range (allowance) about a nominal
    value (nominalSpeed) over a particular distance (maxDistance).
    
    The positions (xChannel and yChannel) are assumed to be sampled uniformly in time
    and the first two time (tChannel) values for the first line in the file are 
    differenced to obtain the sampling interval. From this, the number of samples, N, 
    for the aircraft to travel maxDistance is calculated. The algorithm compares the
    actual distance flown in N samples with maxDistance. If the actual distance is
    greater than

    *maxDistance + allowance * maxDistance*
    
    or less than 

    *maxDistance - allowance * maxDistance*
    
    an error is printed.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    tChannel : String, optional
        The name of the geoWhizz field or channel containing the measured times. The
        default is to read the tChannel field name from the Coordinate Frame.
    nominalSpeed : Float, optional
        The specified ground speed in m/s. The default is 60.0.
    maxDuration : Float, optional
        The time in seconds over which the speed estimate is determined.
        The default is 13.3.
    maxDistance : Float, optional
        The distance in metres over which the speed estimate is determined.
        The default is 1000.0.
    allowance : Float, optional
        The magnitude, relative to the nominalSpeed, of the range of allowed
        speeds. The default is 0.1 (i.e. +/- 10% of nominal).
    minSafeSpeed : Float, optional
        The minimum allowed instantaneous safe speed in m/s. The default is 42.0.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line. Default False.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    if maxDistance > 1.0:
        check_by_time = False
    elif maxDuration > 1.0:
        check_by_time = True
    else:
        print(f'ERROR. Require maxDistance > 1.0 m or maxDuration > 1.0')
        print(f'  maxDistance = {maxDistance}, maxDuration = {maxDuration}.')
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if tChannel == '':
            tChannel = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
        title_str = ''
        project = f[groupName].attrs['ProjectName']
        block = f[groupName].attrs['BlockID']
        if project != '':
            title_str += f'{project}'
        if block != '':
            title_str += f' {block}'
        _reportSpeeds(g, maxDuration=maxDuration, maxDistance=maxDistance, xChannel=xChannel, yChannel=yChannel,
                tChannel=tChannel, vel_north=vel_north, vel_east=vel_east, nominalSpeed=nominalSpeed, 
                     allowance=allowance, minSafeSpeed=minSafeSpeed, title_str=title_str,
                     plot_flag=plot_flag)
        

def _reportSpeeds(group, maxDuration=0.0, maxDistance=0.0, xChannel='X', yChannel='Y',
                tChannel='time', vel_north='', vel_east='', nominalSpeed=60.0,
                allowance=0.1, minSafeSpeed=42.0, title_str='', plot_flag=False):
    """
    Checks the data from Whizz Line group for speed exceedances against a specification
    requiring ground speeds to be within a relative range (`allowance`) about a nominal
    value (`nominalSpeed`) over a specified distance (`maxDistance`). If `maxDistance` is
    not provided, and `maxDuration` is provided, then the check is instead over the
    specified duration.
    
    If `vel_north` and `vel_east` are both provided, then the data in those channels
    is used to calculate speed along the line. If not, then speeds are calculated
    by differencing the positions (`xChannel` and `yChannel`) and dividing by the
    difference of the first two time (`tChannel`) values for the first line in the
    group. This assumes that data are sampled uniformly in time.
    
    From this, the number of samples, N, 
    for the aircraft to travel maxDistance is calculated. The algorithm compares the
    actual distance flown in N samples with maxDistance. If the actual distance is
    greater than
        maxDistance + allowance * maxDistance
    or less than 
        maxDistance - allowance * maxDistance
    an error is printed.

    Parameters
    ----------
    group : HDF5 Group
        The Whizz line group containing the survey line data.
    maxDuration : Float, optional
        The time in seconds over which the speed estimate is determined.
        The default is 0.0.
    maxDistance : Float, optional
        The distance in metres over which the speed estimate is determined.
        The default is 0.0.
    xChannel : String, optional
        The field in the line containing X positions. The default is 'X'.
    yChannel : String, optional
        The field in the line containing Y positions. The default is 'Y'.
    tChannel : String, optional
        The field in the line containing sample times. The default is 'time'.
    vel_north : String, optional
        The field in the line containing velocity north. The default is ''.
    vel_east : String, optional
        The field in the line containing velocity east. The default is ''.
    nominalSpeed : Float, optional
        The specified ground speed in m/s. The default is 60.0.
    allowance : Float, optional
        The magnitude, relative to the nominalSpeed, of the range of allowed
        speeds. The default is 0.1 (i.e. +/- 10% of nominal).
    minSafeSpeed : Float, optional
        The minimum allowed instantaneous safe speed in m/s. The default is 42.0.
    title_str : String, optional
        A title string for the plots. Default ''.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line. Default False.

    Returns
    -------
    None.

    """
    calc_from_pos = False
    if vel_north == '' or vel_east == '':
        calc_from_pos = True
        print('Velocities not known - will calculate from positions')

    check_against_dist = False
    max_allowed_str = f'{maxDuration:.1f} s'
    if maxDuration < 1.0:
        check_against_dist = True
        max_allowed_str = f'{maxDistance:.1f} m'

    settings = f'Nominal ground speed {nominalSpeed:.1f} m/s; '
    settings += f'allowed {nominalSpeed-allowance*nominalSpeed:.1f} : '
    settings += f'{nominalSpeed+allowance*nominalSpeed:.1f} for < {max_allowed_str}.\n'
    print(settings)
    
    num_failed_lines = 0
    num_failures = 0
    num_exceed_lines = 0
    total_num_lines = 0
    report = ''

    lines = list(group.keys())
    total_num_lines = len(lines)

    for line in lines:
        if title_str == '':
            plot_title = f'Line {line}'
        else:
            plot_title = f'{title_str}: Line {line}'

        x, y, dist, t, speed = _get_data(group[line], xChannel, yChannel, tChannel, vel_north, vel_east)

        if speed[speed < minSafeSpeed].size > 0:
            lineUnsafeSlow = True
            print(f' For at least one reading in L{line}, the ground speed was < {minSafeSpeed} (might be unsafe).')
        
        rel_speed = (speed - nominalSpeed) / nominalSpeed

        speed_extreme = 0.0
        num_fids_in_exceedance = 0
        exceedance_in_line = False
        too_slow = False
        line_fails = False

        for fid in range(0, len(x)):
            # There is a speed exceedance ...
            if rel_speed[fid] > allowance or rel_speed[fid] < -allowance:
                # If a new exceedance, then initialise variables;
                if num_fids_in_exceedance == 0:
                    num_exceed_lines += 1
                    if rel_speed[fid] < -allowance:
                        too_slow = True
                    else:
                        too_slow = False
                    start_x = x[fid]
                    start_y = y[fid]
                    start_t = t[fid]
                    num_fids_in_exceedance = 1
                    speed_extreme = speed[fid]
                # Else increment and update on the current exceedance.
                else:
                    # check we haven't swapped from too slow to too fast or vice versa
                    if (too_slow and rel_speed[fid] > allowance) or ((not too_slow) and rel_speed[fid] < allowance):
                        print(f'WARNING: Exceedance reversed speed in one fid. NOT BELIEVABLE. At time={t[fid]:.3f}')
                    num_fids_in_exceedance += 1
                    if too_slow:
                        speed_extreme = min(speed[fid], speed_extreme)
                    else:
                        speed_extreme = max(speed[fid], speed_extreme)
            else:
                if num_fids_in_exceedance > 0: # the current exceedance has ended
                    end_x = x[fid]
                    end_y = y[fid]
                    end_t = t[fid]
                    dist_exceedance = _dist([start_x, end_x], [start_y, end_y])[1]
                    durn_exceedance = end_t - start_t
                    if too_slow:
                        speed_msg = "too slow"
                    else:
                        speed_msg = "too fast"
                    if _exceedance_fail(durn_exceedance, dist_exceedance, maxDuration, maxDistance):
                        report += f'\nL {line} {speed_msg} for {durn_exceedance} sec '
                        report += f'({dist_exceedance:.0f} m), peak exceedance = {speed_extreme:.0f} m/s.'
                        report += f'\n  From ({start_x:.0f} E {start_y:.0f} N) to ({end_x:.0f} E {end_y:.0f} N).'
                        exceedance_in_line = True
                        num_failures += 1
                    num_fids_in_exceedance = 0
                    too_slow = False
        if exceedance_in_line:
            num_failed_lines += 1
            if plot_flag:
                _plot_speed(t, speed, nominalSpeed * (1.0 - allowance), 
                   nominalSpeed * (1.0 + allowance), plot_title=plot_title)

    print(f' Checked {total_num_lines} lines and {num_exceed_lines} had some short exceedance(s).')
    print(f' {num_failed_lines} lines failed for exceedance > allowed.')
    print(f' Total number of full exceedances = {num_failures}.')
    print(report)
    if plot_flag:
        plt.show()


def _get_data(line_group, xChannel, yChannel, tChannel, vel_north, vel_east):
    """
    """
    x = np.array(line_group[xChannel])
    y = np.array(line_group[yChannel])
    t = np.array(line_group[tChannel])
    distance = np.sqrt(x * x + y * y)
    if vel_north == '' or vel_east == '':
        sampleTime = t[1] - t[0]
        xVel = np.diff(x) / sampleTime
        yVel = np.diff(y) / sampleTime
        temp = np.sqrt(xVel * xVel + yVel * yVel)
        speed = np.append(temp, np.mean(temp))
    else:
        xVel = np.array(line_group[vel_east])
        yVel = np.array(line_group[vel_north])
        speed = np.sqrt(xVel * xVel + yVel * yVel)

    return x, y, distance, t, speed


def _plot_speed(t, speed, min_speed=54, max_speed=66, plot_title=''):
    """
    Plots the speed as a function of t and adds limit lines showing the allowed
    tolerance

    Parameters
    ----------
    t : TYPE
        Time (sec) data vector.
    speed : TYPE
        The speed data vector.
    min_speed : TYPE, optional
        DESCRIPTION. The default is 54.
    max_speed : TYPE, optional
        DESCRIPTION. The default is 66.
    plot_title : TYPE, optional
        DESCRIPTION. The default is ''.

    Returns
    -------
    None.

    """
    fig = plt.figure()
    fig.suptitle(plot_title, fontsize=10)
    fig.subplots_adjust(top=0.85)
    
    ax = fig.add_subplot(1,1,1)
    ax.plot(t, speed, 'b', t, np.ones(t.size) * min_speed, 'r', t, np.ones(t.size) * max_speed, 'r', mfc='w')
    plt.xlabel('Time [s]', fontsize = 6)
    plt.ylabel('Speed [m/s]', fontsize = 6)
    plotTitle = 'Speed' + ' Stats'
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    
    
def _displacement2(x0, x1, y0, y1):
    """
    The displacement distance from (xo, y0) to (x1, y1).

    Parameters
    ----------
    x0 : TYPE
        DESCRIPTION.
    x1 : TYPE
        DESCRIPTION.
    y0 : TYPE
        DESCRIPTION.
    y1 : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return np.sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1))
            

