import numpy as np
import matplotlib.pyplot as plt
import h5py
# import xarray as xr
# import verde as vd
# import pooch

import AirGravQC.config as config
import AirGravQC.whizzFiles.pointfiles as gw
import AirGravQC.whizzPlots.whizzPlot as wpl
import AirGravQC.utility.utility as util
import matplotlib.ticker as tkr
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica Neue']})

groupName = config.groupName


def linesMap(whizzFiles=[], easting='', northing='', whizzPlanFile='', planLines=[], planEast='', planNorth=''):
    """
    Plots a line map of the survey contained in the HDF5 Whizz file.

    Parameters
    ----------
    whizzFiles : Array of String or pathlib.PosixPath
        Each element is the name of a HDF5 Whizz file, including path and extension.
    easting : String, optional
        The name of the field containing eastings. The default is the name
        stored in the Coordinates attribute XChannel.
    northing : String, optional
        The name of the field containing eastings. The default is the name
        stored in the Coordinates attribute YChannel.
    whizzPlanFile : String or pathlib.PosixPath
        The name of a HDF5 Whizz file, including path and extension.
    planLines : Array of String
        Each element is the name of a survey line in whizzPlanFile, flown lines will
        only be plotted if their planned line number is in this list. Default [] ignores
        this limit.
    planEast : String, optional
        The name of the field containing eastings. The default is the name
        stored in the Coordinates attribute XChannel.
    planNorth : String, optional
        The name of the field containing eastings. The default is the name
        stored in the Coordinates attribute YChannel.

    Returns
    -------
    None.

    """
    if whizzPlanFile == '' and whizzFiles == []:
        print("No files provided so no Line Map can be made.")
        return

    plot_subtitle = ''
    if whizzPlanFile != '' and whizzFiles != []:
        plot_subtitle = '[planned (red); flown (blue)]'

    fig = plt.figure(figsize=(8, 6))
    thou_format = tkr.FuncFormatter(util._space_thou)
    ax = fig.add_subplot(1,1,1)
    plotTitle = ''

    if whizzPlanFile != '':
        planname = str(whizzPlanFile)

        with h5py.File(planname, 'r') as f:
            if planEast == '':
                planEast = f[groupName]['CoordinateFrame'].attrs['XChannel']
            if planNorth == '':
                planNorth = f[groupName]['CoordinateFrame'].attrs['YChannel']
            g = f[groupName]['Lines']
            plotTitle = wpl.make_plot_title(f[groupName]) + ': Line Map'
                    
            if planLines == []:
                planLines = list(g.keys())
            for line in planLines:
                lX = g[line][planEast][0:]
                lY = g[line][planNorth][0:]
                planline, = ax.plot(lX, lY, color='red', lw=0.6, alpha=0.7)

    if whizzFiles != []:
        for file in whizzFiles:
            filename = str(file)

            with h5py.File(filename, 'r') as f:
                if easting == '':
                    easting = f[groupName]['CoordinateFrame'].attrs['XChannel']
                if northing == '':
                    northing = f[groupName]['CoordinateFrame'].attrs['YChannel']
                g = f[groupName]['Lines']
                plotTitle = wpl.make_plot_title(f[groupName]) + ': Line Map'
                        
                for line in list(g.keys()):
                    lineAttrs = list(g[line].attrs)
                    plot_line_flag = True
                    if 'PlannedLine' in lineAttrs:
                        planned_line = f"{g[line].attrs['PlannedLine']:.3f}"
                        if not (planned_line in planLines or whizzPlanFile == ''):
                            continue

                    lX = gw.getLineData(g[line], easting)[0:]
                    lY = gw.getLineData(g[line], northing)[0:]
                    flownline, = ax.plot(lX, lY, color='blue', lw=0.6, alpha=0.7)
    # If we get one survey line flown due north or east, we need this:        
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    x_range = abs(xmax - xmin)
    y_range = abs(ymax - ymin)
    if x_range < 0.1 * y_range:
        xmax += 0.05 * y_range
        xmin += -0.05 * y_range
    if y_range < 0.1 * x_range:
        ymax += 0.05 * x_range
        ymin += -0.05 * x_range
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    ax.set_aspect('equal')
    ax.xaxis.set_major_formatter(thou_format)
    ax.yaxis.set_major_formatter(thou_format)
    plt.xlabel(f'{easting} [m]', fontsize = 10)
    plt.ylabel(f'{northing} [m]', fontsize = 10)
    plt.suptitle(plotTitle, fontsize = 12)
    plt.title(plot_subtitle, fontsize = 10)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(8)
    for label in ax.get_yticklabels(): label.set_fontsize(8)
    plt.show()
