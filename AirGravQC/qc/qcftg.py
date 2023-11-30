#=========================
#
# FTG
#
#=========================

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import AirGravQC.config as config
import AirGravQC.whizzFiles.pointfiles as gw
import AirGravQC.utility.utility as util

groupName = config.groupName

def checkInlineSum(whizzFile, inline1='', inline2='', inline3='', dontfilter=False):
    """
    Estimates the inline sum for each sample on each line in an FTG whizzFile.
    Plots the min, max, mean and stdev (units of eotvos) for each survey line
    as a function of line.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.

    Returns
    -------
    None

    """
    if inline1 == '':
        inline1 = 'Inline1_raw'
    if inline2 == '':
        inline2 = 'Inline2_raw'
    if inline3 == '':
        inline3 = 'Inline3_raw'

    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']

        numLines = len(g.items())
        chMin = np.zeros((numLines,))
        chMax = np.zeros((numLines,))
        chMean = np.zeros((numLines,))
        chStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            data = g[line][inline1]
            data1 = data[()]
            data = g[line][inline2]
            data2 = data[()]
            data = g[line][inline3]
            data3 = data[()]
            ils_BP = _inLineSum(data1, data2, data3, dontfilter=dontfilter)
            print(f'Line {line}, standard deviation of band-pass filtered in-line sum = {np.std(ils_BP):.2g}')
            lineNo[count] = line
            chMin[count] = np.min(ils_BP)
            chMax[count] = np.max(ils_BP)
            chMean[count] = np.mean(ils_BP)
            chStd[count] = np.std(ils_BP)
            count += 1

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(lineNo, chMin, 'bo', mfc='w')
        ax.plot(lineNo, chMax, 'bo', mfc='w')
        ax.errorbar(lineNo, chMean, chStd, capsize=3, marker='s', c='blue')
        plt.ylabel('Inline Sum', fontsize = 6)
        plotTitle = groupName + ': ' + ' Inline Sum Stats'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        plt.show()


def ilsNoiseVturb(whizzFile, diagComponent1, diagComponent2, diagComponent3, noiseSpec=17.0, vertaccel='', vertvelocity='', vertdispl='', labelLines=False):
    """
    For a Bell Air-FTG. For each line, reports the standard deviation of the in-line sums,
    and plots these as a scatter plot against the standard deviation of the vertical
    acceleration (if the acceleration is not supplied, it is estimated as the time difference
    of the vertical velocity, or if that is not supplied, the second time difference of the
    vertical displacement).

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    diagComponent1 : String
        The name of the channel containing the first tensor diagonal (in-line) component.
    diagComponent2 : String
        The name of the channel containing the second tensor diagonal (in-line) component.
    diagComponent3 : String
        The name of the channel containing the third tensor diagonal (in-line) component.
    noiseSpec : Float, optional
        The noise specification (largest allowed in-line sum for any flight line). Default 17.0 E.
    vertaccel : String, optional
        The name of the channel containing the vertical velocity field. Default ''.
    vertvelocity : String, optional
        The name of the channel containing the vertical velocity field. Default ''.
    vertdispl : String, optional
        The name of the channel containing the vertical velocity field. Default ''.
    labelLines : Bool, optional
        if True, label (with the line number) all points on the plot where the
        line failed the specification. Defaults to False

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    report = ''

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        numLines = len(g.items())
        accStd = np.zeros((numLines,))
        ilsStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        labelx = []
        labely = []
        labelt = []
        failed_lines = 0
        count = 0

        for line in g.keys():
            linegroup = g[line]
            if vertaccel != '':
                accel = gw.getLineData(linegroup, vertaccel)
            elif vertvelocity != '':
                data = gw.getLineData(linegroup, vertvelocity)
                accel = np.diff(data, n = 1)
            elif vertdispl != '':
                data = gw.getLineData(linegroup, vertdispl)
                accel = np.diff(data, n = 2)
            else:
                print("ERROR - need one of vertical acceleration, velocity or displacement (height/altitude).")
                return

            accStd[count] = np.std(accel)

            data1 = gw.getLineData(linegroup, diagComponent1)
            data2 = gw.getLineData(linegroup, diagComponent2)
            data3 = gw.getLineData(linegroup, diagComponent3)
            ilsStd[count] = np.std(_inLineSum(data1, data2, data3))
            if ilsStd[count] > noiseSpec:
                if labelLines:
                    labelx.append(accStd[count])
                    labely.append(ilsStd[count])
                    labelt.append(line)
                failed_lines += 1
                report += f'Line {line}: in-line sum = {ilsStd[count]:.1f} exceeds specification of {noiseSpec}.\n'
                # print(f'Line {line}: in-line sum = {ilsStd[count]:.1f} exceeds specification of {noiseSpec}.')
            lineNo[count] = line
            count += 1
        
        # x = {'label': 'Vertical acceleration [m/s/s]', 'data': accStd}
        # y = {'label': 'Inline Sum [E]', 'data': ilsStd}
        # wpl.plotxy(y, x, plotTitle = 'In-line Sum Noise versus Turbulence', xOffset=False, plot_symbol='+')

        fig = plt.figure()
        fig.suptitle(f'In-line Sum Noise vs Turbulence - {projName}', fontsize=12)
        fig.subplots_adjust(top=0.85)
        ax = fig.add_subplot(1,1,1)
        ax.plot(accStd, ilsStd, 'go')
        if labelLines:
            for ii in range(failed_lines):
                plt.text(labelx[ii], labely[ii], labelt[ii], va='top', ha='right', size=6.0)
        plt.ylabel(f'Inline Sum [E]', fontsize = 8)
        plt.xlabel(f'Turbulence [m/s/s]', fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(8)
        for label in ax.get_yticklabels(): label.set_fontsize(8)
        fig.tight_layout()
        plt.show()
    print(report)



def checkRawFTG(whizzFile, lines=[], noiseLimit=50, gradients=[], timechan='', vertaccel='', vertvelocity='', vertdispl=''):
    """
    Reports the high frequency noise for each line in lines from filename (a whizz file)
    which exceeds noiseLimit. See Mark Dransfield's documentation for details of method.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    lines : String list, optional
        The line numbers to be checked. Default is all lines in the whizzFile.
    noiselimit : Float
        The maximum allowable high frequency noise on a line.
    gradients : Array[String]
        An array of channel names containing the gradient component data.
    vertaccel : String, optional
        The name of the channel containing the vertical velocity field. Default ''.
    vertvelocity : String, optional
        The name of the channel containing the vertical velocity field. Default ''.
    vertdispl : String, optional
        The name of the channel containing the vertical velocity field. Default ''.

    Returns
    -------
    None

    """
    filename = str(whizzFile)
    if gradients == []:
        gradients = ['Cross1_raw', 'Cross2_raw', 'Cross3_raw', 'Inline1_raw', 'Inline2_raw', 'Inline3_raw']

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if timechan == '':
            timechan = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
        if lines == []:
            lines = g.keys()
        
        for line in lines:
            reportStr = f'Line {line} Noise: '
            time = g[line][timechan]
            time = time - time[0]
            if vertaccel != '':
                turb = g[line][vertaccel]
                time1 = time
            elif vertvelocity != '':
                data = g[line][vertvelocity]
                turb = np.diff(data, n = 1)
                time1 = time[1:]
            elif vertdispl != '':
                data = g[line][vertdispl]
                turb = np.diff(data, n = 2)
                time1 = time[1:-1]
            else:
                print("ERROR - need one of vertical acceleration, velocity or displacement (height/altitude).")
                return

            for channel in gradients:
                data = gw.getLineData(g[line], channel)
                noSlope = np.zeros((len(data),))
                filtered = np.zeros((len(data),))
                myStd = np.zeros((len(data)-50,))
        
                stData = np.mean(data[:10])
                enData = np.mean(data[-10:])
                slope = (enData - stData) / len(data)
                for ii in range(0, len(data)):
                    noSlope[ii] = data[ii] - stData - ii * slope
                        
                filtered = util._butter_bandpass_filter(noSlope, 0.1, 0.48, 1, order = 6)
                
                for ii in range(0, len(data)-50):
                    myStd[ii] = np.std(filtered[ii:ii+50])
                    
                reportStr += f'{channel} = {np.max(myStd):.1f}; '
                
                if np.max(myStd) > noiseLimit:
                    fig = plt.figure()
                    ax1 = fig.add_subplot(4,1,1)
                    ax1.plot(time, data, time, noSlope, lw=0.5)
                    ax1.set_xlim(time[0], time[-1])
                    plotTitle = f'Line {line}, Channel {channel}: orig and slope removed'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax1.get_xticklabels(): label.set_fontsize(6)
                    for label in ax1.get_yticklabels(): label.set_fontsize(6)
                    ax2 = fig.add_subplot(4,1,2)
                    ax2.plot(time, filtered, lw=0.5)
                    ax2.set_xlim(time[0], time[-1])
                    ax2.set_ylim(-100.0, 100.0)
                    plotTitle = 'filtered'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax2.get_xticklabels(): label.set_fontsize(6)
                    for label in ax2.get_yticklabels(): label.set_fontsize(6)
                    ax3 = fig.add_subplot(4,1,3)
                    sTime = time[25:25+len(myStd)]
                    ax3.plot(sTime, myStd, lw=0.5)
                    ax3.set_xlim(time[0], time[-1])
                    ax3.set_ylim(0.0, 50.0)
                    plotTitle = 'rolling stdev'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax3.get_xticklabels(): label.set_fontsize(6)
                    for label in ax3.get_yticklabels(): label.set_fontsize(6)
                    ax4 = fig.add_subplot(4,1,4)
                    ax4.plot(time1, turb, lw=0.5)
                    ax4.set_xlim(time[0], time[-1])
                    ax4.set_ylim(0.0, 2.0)
                    plotTitle = 'turb'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax4.get_xticklabels(): label.set_fontsize(6)
                    for label in ax4.get_yticklabels(): label.set_fontsize(6)
                    fig.tight_layout()
                    plt.show()
            # print(reportStr)

    
def checkHighFreq(whizzFile, lines=[], noiseLimit=50, channels=[], tChannel='', verbose=False, plot_flag=False):
    """
    Reports the high frequency noise for each line in lines from filename (a whizz file)
    which exceeds noiseLimit. See Mark Dransfield's documentation for details of method.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    lines : String list, optional
        The line numbers to be checked. Default is all lines in the whizzFile.
    noiselimit : Float
        The maximum allowable high frequency noise on a line.
    channels : Array[String]
        An array of channel names containing the gradient component data.

    Returns
    -------
    None

    """
    filename = str(whizzFile)
    if channels == []:
        return

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if lines == []:
            lines = g.keys()
        
        num_lines = len(list(lines))
        num_failed_lines = 0
        summary = f'Checked {num_lines} lines; no line had high frequency signal above {noiseLimit}.'
        reportStr = ''
        for line in lines:
            if tChannel == '':
                tChannel = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
            time = gw.getLineData(g[line], tChannel)
            time = time - time[0]

            for channel in channels:
                data = gw.getLineData(g[line], channel)
                noSlope = np.zeros((len(data),))
                filtered = np.zeros((len(data),))
                myStd = np.zeros((len(data)-50,))
        
                stData = np.mean(data[:10])
                enData = np.mean(data[-10:])
                slope = (enData - stData) / len(data)
                for ii in range(0, len(data)):
                    noSlope[ii] = data[ii] - stData - ii * slope
                        
                filtered = util._butter_bandpass_filter(noSlope, 0.1, 0.48, 1, order = 6)
                
                for ii in range(0, len(data)-50):
                    myStd[ii] = np.std(filtered[ii:ii+50])
                    
                if np.max(myStd) > noiseLimit:
                    num_failed_lines += 1
                    reportStr += f'Line {line}: peak HF noise in {channel} = {np.max(myStd):.1f}.\n'
                    if plot_flag:
                        fig = plt.figure()
                        y_maxscale = 5.0 * noiseLimit
                        ax1 = fig.add_subplot(3,1,1)
                        ax1.plot(time, data, time, noSlope, lw=0.5)
                        ax1.set_xlim(time[0], time[-1])
                        plotTitle = f'Line {line}, Channel {channel}: orig and slope removed'
                        plt.title(plotTitle, fontsize = 8)
                        plt.grid(True)
                        for label in ax1.get_xticklabels(): label.set_fontsize(6)
                        for label in ax1.get_yticklabels(): label.set_fontsize(6)
                        ax2 = fig.add_subplot(3,1,2)
                        ax2.plot(time, filtered, lw=0.5)
                        ax2.set_xlim(time[0], time[-1])
                        ax2.set_ylim(-y_maxscale, y_maxscale)
                        plotTitle = 'filtered'
                        plt.title(plotTitle, fontsize = 8)
                        plt.grid(True)
                        for label in ax2.get_xticklabels(): label.set_fontsize(6)
                        for label in ax2.get_yticklabels(): label.set_fontsize(6)
                        ax3 = fig.add_subplot(3,1,3)
                        sTime = time[25:25+len(myStd)]
                        ax3.plot(sTime, myStd, lw=0.5)
                        ax3.set_xlim(time[0], time[-1])
                        ax3.set_ylim(0.0, y_maxscale)
                        plotTitle = 'rolling stdev'
                        plt.title(plotTitle, fontsize = 8)
                        plt.grid(True)
                        for label in ax3.get_xticklabels(): label.set_fontsize(6)
                        for label in ax3.get_yticklabels(): label.set_fontsize(6)
                        fig.tight_layout()
                        plt.show()
        if num_failed_lines == 1:
            summary = f'Checked {num_lines} lines; 1 line had high frequency signal above {noiseLimit}.'
        elif num_failed_lines > 1:
            summary = f'Checked {num_lines} lines; {num_failed_lines} lines had high frequency signal above {noiseLimit}.'
        print(summary)
        if verbose:
            print(reportStr)

    
def _FTGeigen(Txx, Txy, Txz, Tyy, Tyz, Tzz, line = "", noiselimit=30.0):
    """
    Returns the standard deviation of the Frobenius norm of the gravity gradient.
    Plots an analysis if this is larger than `noiseLimit`.

    Parameters
    ----------
    Txx : numpy 1D array
        The xx component of the gravity gradient in ENU coordinates.
    Txy : numpy 1D array
        The xy component of the gravity gradient in ENU coordinates.
    Txz : numpy 1D array
        The xz component of the gravity gradient in ENU coordinates.
    Tyy : numpy 1D array
        The yy component of the gravity gradient in ENU coordinates.
    Tyz : numpy 1D array
        The yz component of the gravity gradient in ENU coordinates.
    Tzz : numpy 1D array
        The zz component of the gravity gradient in ENU coordinates.
    line : String, optional.
        The line number containing the above data, just used for
        reporting. Default ''.
    noiselimit : Float, optional.
        The maximum allowed standard deviation of the Frobenius norm.
        Default 30.0 E.

    Returns
    -------
    Float
        The standard deviation of the Frobenius norm of the gravity gradient.

    """
    numSamples = len(Txx)
    trace = np.zeros((numSamples,))
    det = np.zeros((numSamples,))
    I2 = np.zeros((numSamples,))
    frob = np.zeros((numSamples,))
    
    for ii in range(0, numSamples):
        a = np.array([[Txx[ii], Txy[ii], Txz[ii]], [Txy[ii], Tyy[ii], Tyz[ii]], [Txz[ii], Tyz[ii], Tzz[ii]]])
        w, v = np.linalg.eig(a)
        trace[ii] = w.sum()
        det[ii] = np.cbrt(w[0] * w[1] * w[2])
        I2[ii] = w[0] * w[1] - (w[0] + w[1]) * (w[0] + w[1])
        frob[ii] = np.linalg.norm(a, 'fro')

    if np.std(frob) > noiselimit:
        myTitle = 'Trace for Line ' + line
        fig = plt.figure(figsize=(5,8))
        ax1 = fig.add_subplot(4,1,1)
        ax1.plot(trace, '.', ms=2)
        plt.ylabel('Trace', fontsize = 6)
        plotTitle = myTitle
        plt.title(plotTitle, fontsize = 6)
        plt.grid(True)
        for label in ax1.get_xticklabels(): label.set_fontsize(6)
        for label in ax1.get_yticklabels(): label.set_fontsize(6)

        ax2 = fig.add_subplot(4,1,2)
        ax2.plot(det, '.', ms=2)
        plt.ylabel('Det', fontsize = 6)
        plotTitle = 'Det'
        plt.title(plotTitle, fontsize = 6)
        #ax.set_ylim(0.0, 0.0)
        plt.grid(True)
        for label in ax2.get_xticklabels(): label.set_fontsize(6)
        for label in ax2.get_yticklabels(): label.set_fontsize(6)
        
        ax3 = fig.add_subplot(4,1,3)
        ax3.plot(-np.sqrt(-I2), '.', ms=2)
        plt.ylabel('I2', fontsize = 6)
        plotTitle = 'I2'
        plt.title(plotTitle, fontsize = 6)
        #ax.set_ylim(0.0, 0.0)
        plt.grid(True)
        for label in ax3.get_xticklabels(): label.set_fontsize(6)
        for label in ax3.get_yticklabels(): label.set_fontsize(6)
        
        ax4 = fig.add_subplot(4,1,4)
        ax4.plot(frob, '.', ms=2)
        plt.ylabel('Frob Norm', fontsize = 6)
        plotTitle = f'Frob Norm: std = {np.std(frob):.1f}'
        plt.title(plotTitle, fontsize = 6)
        ax4.set_ylim(0.0, 100.0)
        plt.grid(True)
        for label in ax4.get_xticklabels(): label.set_fontsize(6)
        for label in ax4.get_yticklabels(): label.set_fontsize(6)
        
        fig.tight_layout()
        plt.show()
    return np.std(frob)
        
        
def checkFrobenius(whizzFile, lines = [], il1='Inline1_raw', il2='Inline2_raw', il3='Inline3_raw', cr1='Cross1_raw', cr2='Cross2_raw', cr3='Cross3_raw', noiselimit=30.0):
    """
    Reports the noise for each line in lines from filename (a whizz file)
    which exceeds noiseLimit. Here the noise is calculated by `_FTGeigen` as
    the Frobenius norm.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    lines : String list, optional.
        The line numbers to be checked. Default is all lines in the whizzFile.
    noiselimit : Float
        The maximum allowable in-line noise on a line.

    Returns
    -------
    None

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if lines == []:
            lines = g.keys()
        
        for line in lines:
            i1 = gw.getLineData(g[line], il1)
            i2 = gw.getLineData(g[line], il2)
            i3 = gw.getLineData(g[line], il3)
            c1 = gw.getLineData(g[line], cr1)
            c2 = gw.getLineData(g[line], cr2)
            c3 = gw.getLineData(g[line], cr3)
            (Gxx, Gxy, Gxz, Gyy, Gyz, Gzz) = _FTGTransform(i1, i2, i3, c1, c2, c3)
            Txx = util._butter_bandpass_filter(Gxx, 0.1, 0.49, 1.0, order = 6)
            Txy = util._butter_bandpass_filter(Gxy, 0.1, 0.49, 1.0, order = 6)
            Txz = util._butter_bandpass_filter(Gxz, 0.1, 0.49, 1.0, order = 6)
            Tyy = util._butter_bandpass_filter(Gxy, 0.1, 0.49, 1.0, order = 6)
            Tyz = util._butter_bandpass_filter(Gyz, 0.1, 0.49, 1.0, order = 6)
            Tzz = util._butter_bandpass_filter(Gzz, 0.1, 0.49, 1.0, order = 6)
            noise = _FTGeigen(Txx, Txy, Txz, Tyy, Tyz, Tzz, line, noiselimit)
            print(f'Check line {line}. Noise = {noise:.1f}')
    

def _inLineSum(il1, il2, il3, fs=1.0, lowcut=0.03, highcut=0.1, dontfilter=False):
    """
    Calculates the filtered in-line sum of the in-line components.

    Parameters
    ----------
    il1 : numpy 1D array
        The first in-line component data.
    il2 : numpy 1D array
        The second in-line component data.
    il3 : numpy 1D array
        The third in-line component data.
    fs : Float
        The sample frequency.
    lowcut : Float
        The low-pass frequency in Hz of the filter.
    highcut : Float
        The high-pass frequency in Hz of the filter.

    Returns
    -------
    numpy 1D array
        The filtered in-line sum.

    """

    order = 6
    ils = (il1 + il2 + il3) / np.sqrt(3.0)
    ils = ils - np.mean(ils)
    if dontfilter:
        return ils

    return util._butter_bandpass_filter(ils, lowcut, highcut, fs, order = order)

    
def _FTGTransform(il1, il2, il3, cr1, cr2, cr3):
    """
    Returns Txx, Txy, Txz, Tyy, Tyz, Tzz from an algebraic transform of
    the three inline and three cross components of the FTG gradiometer.
    Assumes, I believe, that the FTG is oriented so that the horizontal
    projection of the spin axis of GGI3 is north (the y-axis). From equation (5)
    in J. Brewster, Comparison of gravity gradiometer designs using the 3D
    sensitivity function. In SEG International Exposition and 86th Annual
    Meeting, 2016.

    Parameters
    ----------
    il1 : numpy 1D array
        The first in-line component data.
    il2 : numpy 1D array
        The second in-line component data.
    il3 : numpy 1D array
        The third in-line component data.
    cr1 : numpy 1D array
        The first cross component data.
    cr2 : numpy 1D array
        The second cross component data.
    cr3 : numpy 1D array
        The third cross component data.

    Returns
    -------
    Txx : numpy 1D array
        The xx component of the gravity gradient in ENU coordinates.
    Txy : numpy 1D array
        The xy component of the gravity gradient in ENU coordinates.
    Txz : numpy 1D array
        The xz component of the gravity gradient in ENU coordinates.
    Tyy : numpy 1D array
        The yy component of the gravity gradient in ENU coordinates.
    Tyz : numpy 1D array
        The yz component of the gravity gradient in ENU coordinates.
    Tzz : numpy 1D array
        The zz component of the gravity gradient in ENU coordinates.

    """
    sq2 = np.sqrt(2.0)
    sq3 = np.sqrt(3.0)
    Txx = 2.0 / 3.0 * ( il3 - sq2 * cr1)
    Txy = 0.75/np.sqrt(2) * (np.sqrt(2/3) * il3 - np.sqrt(1.5) * (cr1 - cr2))
    Txz = np.sqrt(6)/(2+3*np.sqrt(2)) * (cr1 - cr2 + 2 * il3)
    Tyy = sq2 / 3.0 * cr1 - (sq2 * cr2 + cr3) / sq3  - il3 / 3.0
    Tyz = 3/np.sqrt(2) * (2/3 * (cr3 - il1 + il2) - (cr1 + cr2)/3)
    Tzz = sq2 / 3.0 * cr1 + (sq2 * cr2 + cr3) / sq3  - il3 / 3.0
    
    return Txx, Txy, Txz, Tyy, Tyz, Tzz


def checkRawAGG(whizzFile, ane, auv, bne, buv, turb, time='', lines=[], noiseLimit=10.0):
    """
    Reports the high frequency noise for each line in lines (read from a whizz file)
    which exceeds noiseLimit.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    aNE : String
        The name of the channel containing the A_NE field.
    aUV : String
        The name of the channel containing the A_UV field.
    bNE : String
        The name of the channel containing the B_NE field.
    bUV : String
        The name of the channel containing the B_UV field.
    turbulence : String
        The name of the channel containing the turbulence field.
    lines : String list, optional
        The line numbers to be checked. Default is all lines in the whizzFile.
    noiselimit : Float, optional
        The maximum allowable high frequency noise on a line. Default=10.0.

    Returns
    -------
    None

    """
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if lines == []:
            lines = g.keys()
        
        for line in lines:
            reportStr = f'Line {line} Noise: '
            if time == '':
                time = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
            time_data = gw.getLineData(g[line], time)
            time_data = time_data - time_data[0]
            fs = 1.0 / abs((time_data[1] - time_data[0]))
            turb_data = gw.getLineData(g[line], turb)
            for channel in [ane, auv, bne, buv]:
                data = gw.getLineData(g[line], channel)
                noSlope = np.zeros((len(data),))
                filtered = np.zeros((len(data),))
                myStd = np.zeros((len(data)-50,))
        
                stData = np.mean(data[:10])
                enData = np.mean(data[-10:])
                slope = (enData - stData) / len(data)
                for ii in range(0, len(data)):
                    noSlope[ii] = data[ii] - stData - ii * slope
                        
                filtered = util._butter_bandpass_filter(noSlope, 0.15,
                                                  3.6, fs, order = 6)
                
                for ii in range(0, len(data)-50):
                    myStd[ii] = np.std(filtered[ii:ii+50])
                    
                reportStr += f'{channel} = {np.max(myStd):.1f}; '
                
                if np.max(myStd) > noiseLimit:
                    fig = plt.figure()
                    ax1 = fig.add_subplot(4,1,1)
                    ax1.plot(time_data, data, time_data, noSlope, lw=0.5)
                    ax1.set_xlim(time_data[0], time_data[-1])
                    plotTitle = f'Line {line}, Channel {channel}: orig and slope removed'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax1.get_xticklabels(): label.set_fontsize(6)
                    for label in ax1.get_yticklabels(): label.set_fontsize(6)
                    ax2 = fig.add_subplot(4,1,2)
                    ax2.plot(time_data, filtered, lw=0.5)
                    ax2.set_xlim(time_data[0], time_data[-1])
                    ax2.set_ylim(-10.0, 10.0)
                    plotTitle = f'filtered [0.15,3.6]; fs={fs:.1f}'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax2.get_xticklabels(): label.set_fontsize(6)
                    for label in ax2.get_yticklabels(): label.set_fontsize(6)
                    ax3 = fig.add_subplot(4,1,3)
                    sTime = time_data[25:25+len(myStd)]
                    ax3.plot(sTime, myStd, lw=0.5)
                    ax3.set_xlim(time_data[0], time_data[-1])
                    ax3.set_ylim(0.0, 15.0)
                    plotTitle = 'rolling stdev'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax3.get_xticklabels(): label.set_fontsize(6)
                    for label in ax3.get_yticklabels(): label.set_fontsize(6)
                    ax4 = fig.add_subplot(4,1,4)
                    ax4.plot(time_data, turb_data, lw=0.5)
                    ax4.set_xlim(time_data[0], time_data[-1])
                    ax4.set_ylim(0.0, 2.0)
                    plotTitle = 'turb'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax4.get_xticklabels(): label.set_fontsize(6)
                    for label in ax4.get_yticklabels(): label.set_fontsize(6)
                    fig.tight_layout()
                    plt.show()
            print(reportStr)


