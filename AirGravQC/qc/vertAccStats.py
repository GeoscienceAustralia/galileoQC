import numpy as np
import h5py

import AirGravQC.config as config
import AirGravQC.whizzFiles.pointfiles as gw
import AirGravQC.whizzPlots.whizzPlot as wpl

groupName = config.groupName

        
def vertAccStats(whizzFile, vertaccel='', vertvel='', vertpos='', samplePeriod=0.0):
    """
    Uses numpy.diff to estimate the vertical acceleration. Plots the
    min, max, mean and stdev for each survey line as a function of
    line.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    vertaccel : String, optional
        The name of the channel containing the vertical accelerations. Default is
        to assume it does not exist.
    vertvel : String, optional
        The name of the channel containing the vertical velocities. Default is
        to assume it does not exist.
    vertpos : String, optional
        The name of the channel containing the vertical positions. Default is
        to use the channel in ['CoordinateFrame'].attrs['AltitudeChannel'].
    samplePeriod : Float, optional
        The time between samples in the `whizzFile`. The default is to calculate
        this using the data in ['CoordinateFrame'].attrs['TimeChannel'].

    Returns
    -------
    None

    """
    outlierRatio = 6.0
    specAcc = 1.0
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        lines_group = f[groupName]['Lines']
        gProject = f[groupName]

        numLines = len(g.items())
        chMin = np.zeros((numLines,))
        chMax = np.zeros((numLines,))
        chMean = np.zeros((numLines,))
        chStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            accel_units = 'm/s'
            if samplePeriod == '':
                time_chan = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
                time = gw.getLineData(lines_group[line], time_chan)
                samplePeriod = abs(time[1] - time[0])
            if not vertaccel == '':
                accel = gw.getLineData(lines_group[line], vertaccel)
                titlestr = 'Acceleration from {vertaccel}'
                accel_units = getLineDataUnits(linegroup, channel)
            elif not vertvel == '':
                accel = np.diff(gw.getLineData(lines_group[line], vertvel), n = 1)
                titlestr = 'Acceleration from {vertvel}'
            elif not vertpos == '':
                accel = np.diff(gw.getLineData(lines_group[line], vertpos), n = 2)
                titlestr = 'Acceleration from {vertpos}'
            else:
                vertpos = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
                accel = np.diff(gw.getLineData(lines_group[line], vertpos), n = 2)
                titlestr = 'Acceleration from {vertpos}'

            lineNo[count] = line
            chMin[count] = np.min(accel)
            chMax[count] = np.max(accel)
            chMean[count] = np.mean(accel)
            chStd[count] = np.std(accel)
            maxDevRatio = (chMax[count] - chMean[count]) / chStd[count]
            minDevRatio = (chMean[count] - chMin[count]) / chStd[count]
            if chStd[count] > specAcc:
                print(lineNo[count], ': 100 mill-g exceedance ', chStd[count])
            if maxDevRatio > outlierRatio:
                print(lineNo[count], ': max outlier ratio ', maxDevRatio)
            if minDevRatio > outlierRatio:
                print(lineNo[count], ': min outlier ratio', minDevRatio)
            count += 1

    figtitle = wpl.make_plot_title(gProject)
    xlabelstr = 'Line number'
    ylabelstr = 'Vertical acceleration [{accel_units}]'
    wpl.plotBoxWhisker(chMin, chMax, chMean, chStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)


