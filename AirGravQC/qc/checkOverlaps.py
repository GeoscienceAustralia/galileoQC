import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr

import AirGravQC.config as config
import AirGravQC.whizzFiles.pointfiles as gw
import AirGravQC.utility.utility as util
import AirGravQC.whizzPlots.whizzPlot as wpl

groupName = config.groupName
            

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
                    dirn = np.arctan2((e1[-1] - e1[0]), (n1[-1] - n1[0]))
                    
                    # make life easier by transforming to a 2D problem
                    n0 = n1[0]
                    e0 = e1[0]
                    (y1, x1) = util._rotateCoords(e1 - e0, n1 - n0, -dirn)
                    (y2, x2) = util._rotateCoords(e2 - e0, n2 - n0, -dirn)

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
    if report == '':
        report = f'All overlaps meet requirement (>{min_overlap:.1f} km).'

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
