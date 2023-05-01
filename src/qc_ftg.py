def checkInlineSum(whizzFile):
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
            data = g[line]['Inline1_raw']
            data1 = data[()]
            data = g[line]['Inline2_raw']
            data2 = data[()]
            data = g[line]['Inline3_raw']
            data3 = data[()]
            ils_BP = _inLineSum(data1, data2, data3)
            print(line, ' STD BPF = ', np.std(ils_BP))
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


def ilsNoiseVturb(whizzFile, diagComponent1, diagComponent2, diagComponent3, vertvelocity):
    """
    For a Bell Air-FTG. For each line, reports the standard deviation of the in-line sums,
    and plots these as a scatter plot against the standard deviation of the vertical
    acceleration (estimated as the time difference of the vertical velocity).

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
    vertvelocity : String
        The name of the channel containing the vertical velocity field.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        numLines = len(g.items())
        accStd = np.zeros((numLines,))
        ilsStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            data = g[line][vertvelocity]
            accel = np.diff(data, n = 1)
            accStd[count] = np.std(accel)

            data1 = np.array(g[line][diagComponent1])
            data2 = np.array(g[line][diagComponent2])
            data3 = np.array(g[line][diagComponent3])
            ilsStd[count] = np.std(_inLineSum(data1, data2, data3))
            print(line, ' ', ilsStd[count])
            lineNo[count] = line
            count += 1
        
        x = {'label': 'Vertical acceleration [m/s/s]', 'data': accStd}
        y = {'label': 'Inline Sum [E]', 'data': ilsStd}
        wpl.plotxy(y, x, plotTitle = 'In-line Sum Noise versus Turbulence', xOffset=False, plot_symbol='+')


def checkRawFTG(whizzFile, lines=[], noiseLimit=50):
    """
    Reports the high frequency noise for each line in lines from filename (a whizz file)
    which exceeds noiseLimit. See QC Report for Lawin for details of method.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    lines : String list, optional
        The line numbers to be checked. Default is all lines in the whizzFile.
    noiselimit : Float
        The maximum allowable high frequency noise on a line.

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
            time = g[line]['time']
            time = time - time[0]
            vv = g[line]['vertvelocity']
            turb = np.diff(vv, n = 1)
            for channel in ['Cross1_raw', 'Cross2_raw', 'Cross3_raw', 'Inline1_raw', 'Inline2_raw', 'Inline3_raw']:
                data = np.array(g[line][channel])
                noSlope = np.zeros((len(data),))
                filtered = np.zeros((len(data),))
                myStd = np.zeros((len(data)-50,))
        
                stData = np.mean(data[:10])
                enData = np.mean(data[-10:])
                slope = (enData - stData) / len(data)
                for ii in range(0, len(data)):
                    noSlope[ii] = data[ii] - stData - ii * slope
                        
                filtered = butter_bandpass_filter(noSlope, 0.1, 0.48, 1, order = 6)
                
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
                    ax4.plot(time[1:], turb, lw=0.5)
                    ax4.set_xlim(time[0], time[-1])
                    ax4.set_ylim(0.0, 2.0)
                    plotTitle = 'turb'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax4.get_xticklabels(): label.set_fontsize(6)
                    for label in ax4.get_yticklabels(): label.set_fontsize(6)
                    fig.tight_layout()
                    plt.show()
            print(reportStr)

    
def _FTGeigen(Txx, Txy, Txz, Tyy, Tyz, Tzz, line = "", noiselimit=30.0):
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
        
        
def eigenPlot(whizzFile, lines = [], noiselimit=30.0):
    """
    Reports the noise for each line in lines from filename (a whizz file)
    which exceeds noiseLimit. Here the noise is calculated by `_FTGeigen` as
    the Frobenius norm.
    TODO : Why is the function name misleading?

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
            i1 = np.array(g[line]['Inline1_raw'])
            i2 = np.array(g[line]['Inline2_raw'])
            i3 = np.array(g[line]['Inline3_raw'])
            c1 = np.array(g[line]['Cross1_raw'])
            c2 = np.array(g[line]['Cross2_raw'])
            c3 = np.array(g[line]['Cross3_raw'])
            (Gxx, Gxy, Gxz, Gyy, Gyz, Gzz) = _FTGTransform(i1, i2, i3, c1, c2, c3)
            Txx = butter_bandpass_filter(Gxx, 0.1, 0.49, 1.0, order = 6)
            Txy = butter_bandpass_filter(Gxy, 0.1, 0.49, 1.0, order = 6)
            Txz = butter_bandpass_filter(Gxz, 0.1, 0.49, 1.0, order = 6)
            Tyy = butter_bandpass_filter(Gxy, 0.1, 0.49, 1.0, order = 6)
            Tyz = butter_bandpass_filter(Gyz, 0.1, 0.49, 1.0, order = 6)
            Tzz = butter_bandpass_filter(Gzz, 0.1, 0.49, 1.0, order = 6)
            noise = _FTGeigen(Txx, Txy, Txz, Tyy, Tyz, Tzz, line, noiselimit)
            print(f'Check line {line}. Noise = {noise:.1f}')
    

def _inLineSum(il1, il2, il3, fs=1.0, lowcut=0.03, highcut=0.1):
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
    fs = 1.0
    lowcut = 0.03
    highcut = 0.1
    order = 6
    ils = (il1 + il2 + il3) / np.sqrt(3.0)
    ils = ils - np.mean(ils)
    return butter_bandpass_filter(ils, lowcut, highcut, fs, order = order)

    
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


