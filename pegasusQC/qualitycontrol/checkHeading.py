#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check aircraft heading is within specification.
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.utility.utility as util
import pegasusQC.whizzPlots.whizzPlot as wpl

groupName = config.groupName
        

def checkHeading(whizzFile, nominalHeadings=[], headingchan='', x='', y='', tolerance=10.0, known='', lines=[], plot_flag=False):
    """
    Checks heading in degrees is within +/- tolerance (in degrees) of nominal (in degrees). Actually
    checks against `sin(nominalHeading +/- tolerance)`.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    nominalHeadings : [Float], optional
        The desired headings in degrees from north. The default is to check against the mean heading.
    headingchan : String, optional
        The name of the geoWhizz channel containing the headings. The
        default is to calculate the heading from the x and y channels.
    x : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    y : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    tolerance : Float, optional
        Headings within +/- tolerance degrees of nominalHeading are ok.
    known : String, optional
        If present, the name of the channel containing the "known error" flag.
        This is reported against any error so that known errors can be distinguished
        from unknown errors.
    lines : Array{String}, optional
        Array of line numbers as strings. Default = [], meaning all lines are checked.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    report = ''
    num_failed_lines = 0

    exceedances_known = False
    this_exc_known = False
    number_known = 0

    checkmean = True
    testmode = 'mean'
    if len(nominalHeadings) > 0:
        checkmean = False
        testmode = 'nominal'

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if lines == []:
            lines = g.keys()
        numLines = len(lines)

        for line in lines:
            x_data = rd.getLineData(g[line], x)
            y_data = rd.getLineData(g[line], y)
            distance = util._length(x_data, y_data)

            if known != '':
                exceedances_known = True
                exc_known = rd.getLineData(g[line], known)
                report_known = -1

            if headingchan == '':
                dx = np.diff(x_data)
                dy = np.diff(y_data)
                heading = np.arctan2(dx, dy) * 180.0 / np.pi
                plot_x = distance[1:]
            else:
                heading = rd.getLineData(g[line], headingchan)
                plot_x = distance
            allok = True
            min_heading = np.nanmin(heading % 360)
            max_heading = np.nanmax(heading % 360)
            mean_heading = np.nanmean(heading % 360)

            if checkmean:
                allok = all(angle_in_range(heading, mean_heading, tolerance))
                if allok:
                    break
            else:
                for nomhead in nominalHeadings:
                    allok = all(angle_in_range(heading, nomhead, tolerance))
                    if allok:
                        break
            
            if not allok:
                num_failed_lines += 1
                report += f'Line {line}: at least one sample failed. '
                report += f'Min {min_heading:.2f}, Max {max_heading:.2f} deg, Mean {mean_heading:.2f} deg.'
                if exceedances_known:
                    if np.max(exc_known) > 0:
                        report += f'Exceedance known on line: {exc_known[fid]:.0f}'
                report += '\n'
                if plot_flag:
                    fig = plt.figure()
                    fig.suptitle(f'Heading Check Line {line}', fontsize=10)
                    fig.subplots_adjust(top=0.85)
                    
                    ax = fig.add_subplot(1,1,1)
                    thou_format = tkr.FuncFormatter(util._space_thou)
                    ax.plot(plot_x, heading % 360, 'b', mfc='w')
                    ax.xaxis.set_major_formatter(thou_format)
                    plt.ylabel('Estimated heading [deg]', fontsize = 6)
                    plt.xlabel(f'{x} [m]', fontsize = 6)
                    plt.grid(True)
                    for label in ax.get_xticklabels(): label.set_fontsize(6)
                    for label in ax.get_yticklabels(): label.set_fontsize(6)

    # print(f'Heading limits: [{nominalHeading}, +/-{tolerance}] deg or equivalent.')
    print(f'  Checked {numLines} lines for heading - {testmode} > tolerance {tolerance}; {num_failed_lines} failed.\n')
    print(report)
    if plot_flag and num_failed_lines > 0:
        plt.show()
    return


def angle_in_range(alpha, nominal, tolerance):
    lower = nominal - abs(tolerance)
    upper = nominal + abs(tolerance)
    return (alpha - lower) % 360 <= (upper - lower) % 360
