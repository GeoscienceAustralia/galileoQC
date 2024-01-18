import numpy as np
import h5py

import AirGravQC.config as config
import AirGravQC.whizzFiles.pointfiles as gw
import AirGravQC.whizzPlots.whizzPlot as wpl

groupName = config.groupName


def checkDiffStats(whizzFile, channel1, channel2, flightLines=[]):
    """
    Plot the statistics of channel1 - channel2 in each flightLine. 

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    flightLines : String List, optional
        A list of flightline, e.g. ['1000110.0']. Default is all lines in whizzFile.
    channel1 : String
        The name of a channel.
    channel2 : String
        The name of a channel.

    Returns
    -------
    None.

    """
    
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if flightLines == []:
            flightLines = list(g.keys())
        corr_units = g[flightLines[0]][channel1].attrs['Units']
        if not (g[flightLines[0]][channel1].attrs['Units'] == corr_units):
            print('Error: {channel1} and {channel2} do not have the same units.')
            return

        # initialise variables
        chMin = np.zeros((numLines,))
        chMax = np.zeros((numLines,))
        chMean = np.zeros((numLines,))
        chStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in flightLines:
            if line != 'CoordinateFrame':
                lineNo[count] = line
                dd = gw.getLineData(g[line], channel1) - gw.getLineData(g[line], channel2)

                if np.sum(~np.isnan(dd)) > 3:
                    chMin[count] = np.nanmin(dd)
                    chMax[count] = np.nanmax(dd)
                    chMean[count] = np.nanmean(dd)
                    chStd[count] = np.nanstd(dd)
                else:
                    chMin[count] = 0.0
                    chMax[count] = 0.0
                    chMean[count] = 0.0
                    chStd[count] = 0.0
                    print(f'Less than three real values in {lineNo[count]:.2f} for {channel1}, {channel2}, no statistics.')
                count += 1

        figtitle = wpl.make_plot_title(gProject)
        titlestr = channel1 + ' - ' + channel2 + ' Stats'
        xlabelstr = 'Line number'
        ylabelstr = f'{channel1} - {channel2} [{corr_units}]'
        wpl.plotBoxWhisker(chMin, chMax, chMean, chStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)
    return
