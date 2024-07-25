import numpy as np
import h5py
from scipy.signal import correlate
import matplotlib.pyplot as plt
# from scipy.signal import butter, lfilter
# from matplotlib.ticker import StrMethodFormatter
# import matplotlib.ticker as tkr

import AirGravQC.config as config
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.whizzFiles.pointfiles as gw
# import AirGravQC.gridFiles.read_ers as ers
# import AirGravQC.gridFiles.gridfiles as grd
# import AirGravQC.utility.utility as util
# import AirGravQC.whizzPlots.whizzPlot as wpl

groupName = config.groupName

def checkPhase(whizzFile, channel1, channel2, tChannel='', tolerance=1.0, lines=[], verbose=False, plot_flag=False):
    """
    For every survey line in a geoWhizz HDF5 file, given two channels, calculate
    the phase shift (in seconds if time data is available, or number of samples otherwise)
    required to maximise the correlation between the two channels.

    If the result is greater than the tolerance and plot_flag=True, the correlation is plotted.

    Parameters
    ----------
    whizzFile : String
        The name of a geoWhizz HDF5 file.
    channel1 : String
        The name of the first channel.
    channel2 : String
        The name of the second channel.
    tChannel : String, optional
        The name of the time channel. The default, '', uses the 'TimeChannel' in
        the whizzFile attributes.
    tolerance : Float, optional
        The maximum allowed phase shift (in seconds if time data are available, else in number of samples).
    lines : String list, optional.
        The line numbers to be checked. Default is all lines in the whizzFile.
    verbose : Bool, optional
        If True, report status of all lines, else only report errors. Default False.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line.
    
    Returns
    -------
    None.

    """
    with h5py.File(whizzFile, 'r') as f:
        g = f[groupName]['Lines']
        numLines = len(g.items())
        offsets = np.zeros((numLines,))
        count = 0
        # TODO: trap case where tChannel does not exist and report "number of samples"
        if lines == []:
            lines = g.keys()
        if tChannel == '':
            tChannel = f[groupName]['CoordinateFrame'].attrs['TimeChannel']

        for line in lines:
            linegroup = g[line]
            A = rd.getLineData(linegroup, channel1)
            B = rd.getLineData(linegroup, channel2)
            time = rd.getLineData(linegroup, tChannel)
            fs = 1.0 / abs(time[1] - time[0])
            
            nsamples = len(A)

            # regularize datasets by subtracting mean and dividing by s.d.
            A -= A.mean()
            A /= A.std()
            B -= B.mean()
            B /= B.std()

            # Find cross-correlation
            xcorr = correlate(A, B)
            
            # delta time array to match xcorr
            dt = np.arange(1-nsamples, nsamples)
            recovered_time_shift = dt[xcorr.argmax()]
            offsets[count] = recovered_time_shift
            count += 1

            time = np.arange(dt[0], dt[-1], 0.1)
            # Now interpolate through gaps by cubic spline
            (xcorrInt, _) = gw.interpolateLine(dt, xcorr, time)
            recovered_time_shift2 = time[xcorrInt.argmax()]
            if verbose:
                print(f'Line {line}: Recovered time shift = {recovered_time_shift2 / fs:.1f} sec.')
            
            if abs(recovered_time_shift2 / fs) > tolerance and plot_flag:
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                ax.plot(time[0:len(xcorrInt)] / fs, xcorrInt)
                plt.ylabel('Correlation [arbitrary units]', fontsize = 6)
                plt.xlabel('FID difference [s]', fontsize = 6)
                plotTitle = f'Line {line}: Correlation of {channel1} v {channel2}'
                plt.title(plotTitle, fontsize = 8)
                plt.grid(True)
                for label in ax.get_xticklabels(): label.set_fontsize(6)
                for label in ax.get_yticklabels(): label.set_fontsize(6)
                plt.show()
                print(f'Offset MEAN = {np.mean(offsets):.2f}; STD = {np.std(offsets):.3f}')
        
        
