import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import AirGravQC.config as config
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.utility.utility as util

groupName = config.groupName

    
def checkHighFreq(whizzFile, lines=[], noiseLimit=50, channels=[], cutoffs=[0.15, 3.6], tChannel='', vertaccel='', vertvelocity='', vertdispl='', verbose=False, plot_flag=False):
    """
    Reports the high frequency noise for each line in lines from filename (a whizz file)
    which exceeds noiseLimit. See Mark Dransfield's documentation for details of method.

    For AGG (Falcon), the recommended channels are ANE, AUV, BNE, BUV. For FTG, the 6 raw cross and inline components.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    lines : String list, optional
        The line numbers to be checked. Default is all lines in the whizzFile.
    noiselimit : Float
        The maximum allowable high frequency noise on a line.
    cutoffs : Array of Float, optional
        The low pass, and high pass, cutoff frequencies that define "HighFreq".
        Recommended values: AGG - [0.15, 3.6], FTG - [0.1, 0.48]. Default AGG.
    channels : Array[String]
        An array of channel names containing the gradient component data.
    vertaccel : String, optional
        The name of the channel containing the vertical velocity field. Default ''.
    vertvelocity : String, optional
        The name of the channel containing the vertical velocity field. Default ''.
    vertdispl : String, optional
        The name of the channel containing the vertical velocity field. Default ''.
    verbose : Bool, optional
        If True, report status of all overlaps, else only report errors. Default False.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line.

    Returns
    -------
    None

    """
    filename = str(whizzFile)
    if channels == []:
        return

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if tChannel == '':
            tChannel = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
        if lines == []:
            lines = g.keys()

        num_subplots = 4
        num_lines = len(list(lines))
        num_failed_lines = 0
        summary = f'Checked {num_lines} lines; no line had high frequency signal above {noiseLimit}.'
        reportStr = ''
        for line in lines:
            time = rd.getLineData(g[line], tChannel)
            time = time - time[0]
            fs = 1.0 / abs((time[1] - time[0]))
            if vertaccel != '':
                turb = rd.getLineData(g[line], vertaccel)
                time1 = time
            elif vertvelocity != '':
                data = rd.getLineData(g[line], vertvelocity)
                turb = np.diff(data, n = 1)
                time1 = time[1:]
            elif vertdispl != '':
                data = rd.getLineData(g[line], vertdispl)
                turb = np.diff(data, n = 2)
                time1 = time[1:-1]
            else:
                num_subplots = 3

            for channel in channels:
                data = rd.getLineData(g[line], channel)
                noSlope = np.zeros((len(data),))
                filtered = np.zeros((len(data),))
                myStd = np.zeros((len(data)-50,))
        
                stData = np.mean(data[:10])
                enData = np.mean(data[-10:])
                slope = (enData - stData) / len(data)
                for ii in range(0, len(data)):
                    noSlope[ii] = data[ii] - stData - ii * slope
                
                filtered = util._butter_bandpass_filter(noSlope, cutoffs[0], cutoffs[1], fs, order = 6)
                # filtered = util._butter_bandpass_filter(noSlope, 0.1, 0.48, fs, order = 6)
                
                for ii in range(0, len(data)-50):
                    myStd[ii] = np.std(filtered[ii:ii+50])
                    
                if np.max(myStd) > noiseLimit:
                    num_failed_lines += 1
                    reportStr += f'Line {line}: peak HF noise in {channel} = {np.max(myStd):.1f}.\n'
                    if plot_flag:
                        fig = plt.figure()
                        y_maxscale = 5.0 * noiseLimit
                        plotTitle = f'Line {line}, Channel {channel}: orig and slope removed'
                        _subplot_hiF_analysis(fig, num_subplots, 1, plotTitle, time, data, x2=time, y2=noSlope)

                        plotTitle = 'filtered'
                        _subplot_hiF_analysis(fig, num_subplots, 2, plotTitle, time, filtered)

                        sTime = time[25:25+len(myStd)]
                        plotTitle = 'rolling stdev'
                        _subplot_hiF_analysis(fig, num_subplots, 3, plotTitle, sTime, myStd, bounds=[time[0], time[-1], 0.0, noiseLimit])

                        if num_subplots == 4:
                            plotTitle = 'Turbulence'
                            _subplot_hiF_analysis(fig, num_subplots, 4, plotTitle, time1, turb, bounds=[time[0], time[-1], 0.0, 2.0])

                        fig.tight_layout()
        if num_failed_lines == 1:
            summary = f'Checked {num_lines} lines; 1 line had high frequency signal above {noiseLimit}.'
        elif num_failed_lines > 1:
            summary = f'Checked {num_lines} lines; {num_failed_lines} lines had high frequency signal above {noiseLimit}.'
        print(summary)
        if verbose:
            print(reportStr)
        if plot_flag and num_failed_lines > 0:
            plt.show()


def _subplot_hiF_analysis(fig, num_subplots, plotIdx, plotTitle, x1, y1, x2=np.array([]), y2=np.array([]), bounds=[]):
    if bounds != []:
        xmin = bounds[0]
        xmax = bounds[1]
        ymin = bounds[2]
        ymax = bounds[3]
    else:
        xmin = x1[0]
        xmax = x1[-1]
        ymin = 0.0
        ymax = ymin
    ax = fig.add_subplot(num_subplots, 1, plotIdx)
    if np.array(x2).size == 0 or np.array(y2).size == 0:
        ax.plot(x1, y1, lw=0.5)
    else:
        ax.plot(x1, y1, x2, y2, lw=0.5)
    ax.set_xlim(xmin, xmax)
    if ymin != ymax:
        ax.set_ylim(ymin, ymax)
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)

