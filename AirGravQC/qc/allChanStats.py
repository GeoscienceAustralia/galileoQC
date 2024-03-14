import numpy as np
import h5py

import AirGravQC.config as config
import AirGravQC.whizzPlots.whizzPlot as wpl
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.utility.utility as util

groupName = config.groupName


def allChanStats(whizzFile, allChannels=[], lines=[], d1_chans=[], mr_chans=[], sin_chans=[]):
    """
    Generate statistical plots for the channels across all lines. The plots show
    the min, mean, max and stdev for each channel as a function of line number.

    The statistics can be optionally calculated on the data after first differencing,
    or mean removal, or both. In the latter case, first differencing is done first.

    Parameters
    ----------
    whizzFile : String or pathlib Path.
        Name of a HDF5 Whizz file, including path and extension.
    allChannels : [String], optional.
        A list of the channels or fields to plot. Default is all in whizzFile.
    lines : String list, optional.
        The line numbers to be checked. Default is all lines in the whizzFile.
    d1_chans : [String], optional.
        An array of names of channels from `allChannels` whose first difference
        along each survey line should be calculated before the statistics.
    mr_chans : [String], optional.
        An array of names of channels from `allChannels` whose mean along each
        survey line should be subtracted before calculating statistics.
    sin_chans : [String], optional.
        An array of names of channels from `allChannels` whose sine
        along each survey line should be calculated before the statistics.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        gProject = f[groupName]
        g = gProject['Lines']
        if lines == []:
            lines = list(g.keys())
        if allChannels == []:
            lineGroups = list(g.values())
            allChannels = list(lineGroups[0].keys())
        numLines = len(lines)
        
        for channel in allChannels:
            remove_mean = False
            diff_one = False
            remove_sine = False
            if channel in mr_chans:
                remove_mean = True
            if channel in d1_chans:
                diff_one = True
            if channel in sin_chans:
                remove_sine = True

            # initialise variables
            chMin = np.zeros((numLines,))
            chMax = np.zeros((numLines,))
            chMean = np.zeros((numLines,))
            chStd = np.zeros((numLines,))
            lineNo = np.zeros((numLines,))
            count = 0

            # get the units for the y axis label
            xlabelstr = 'Line number'
            my_units = rd.getChannelAttrs(g[lines[0]], channel)
            if my_units == '':
                ylabelstr = f'{channel}'
            else:
                ylabelstr = f'{channel} [{my_units}]'
        
            for line in lines:
                if line != 'CoordinateFrame':
                    lineNo[count] = line
                    dd = rd.getLineData(g[line], channel)
                    if remove_sine:
                        dd = np.sin(dd * (np.pi / 180.0))
                    if diff_one:
                        dd = np.append(np.diff(dd), dd[-1]-dd[-2])
                    if remove_mean:
                        dd = dd - np.mean(dd)

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
                        print(f'Less than three real values in {lineNo[count]:.2f} for {channel}, no statistics.')
                    count += 1

            figtitle = wpl.make_plot_title(gProject)
            titlestr = ''
            if remove_mean:
                titlestr += 'mr('
            if diff_one:
                titlestr += ' d1('
            if remove_sine:
                titlestr += 'sin('

            titlestr += channel

            if remove_mean:
                titlestr += ')'
            if diff_one:
                titlestr += ')'
            if remove_sine:
                titlestr += ')'
            titlestr += ' Stats'
            wpl.plotBoxWhisker(chMin, chMax, chMean, chStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)
    return
