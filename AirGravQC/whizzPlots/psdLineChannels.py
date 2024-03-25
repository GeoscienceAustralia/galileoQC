import numpy as np
import matplotlib.pyplot as plt
import h5py
    import scipy.signal as sig

import AirGravQC.config as config
import AirGravQC.qc.qualityAnalysis as qc
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.utility.utility as util
import matplotlib.ticker as tkr

groupName = config.groupName


def psdLineChannels(whizzFile, flightLine, channels, time='', plotTitle = ''):
    """
    Plot the PSD (log-log Sqrt(Power) from welch method) of channels in flightLine. 

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension.
    flightLine : String
        A flightline, e.g. '1000110.0'.
    channels : String Array
        The names of the channels to plot.
    time : String, optional
        The name of the channel containing the time data. Defaults to 'TimeChannel'.
    plotTitle : String, optional
        A title for the plot. The default is '' in which case the title will be Project Line Channel.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if time == '':
            time = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
        t = rd.getLineData(g[flightLine], time)
        f_sample = 1.0 / (t[1] - t[0])
        for channel in channels:
            data = rd.getLineData(g[flightLine], channel)
            freq, Pxx = sig.welch(data, nfft=2048, fs = f_sample)
            period = 1.0 / freq[1:]
            plt.plot(period, np.sqrt(Pxx[1:]), 'g', lw=0.6)

    plt.xlim([0, 200])
    
    plt.xlabel('Period [s]', fontsize = 6)
    plt.ylabel(f'sqrt(PSD)', fontsize = 6)
    if plotTitle == '':
        plotTitle = f'{projName} L{flightLine}: sqrt(Pwr))'
    plt.title(plotTitle, fontsize = 8)
    plt.grid(True)
    for label in ax.get_xticklabels(): label.set_fontsize(6)
    for label in ax.get_yticklabels(): label.set_fontsize(6)
    plt.show()

    return
    
