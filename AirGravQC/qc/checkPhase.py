import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
from matplotlib.ticker import StrMethodFormatter
# from scipy.signal import butter, lfilter
from scipy.signal import correlate

import AirGravQC.config as config
import AirGravQC.whizzFiles.pointfiles as gw
# import AirGravQC.gridFiles.read_ers as ers
# import AirGravQC.gridFiles.gridfiles as grd
import AirGravQC.utility.utility as util
import AirGravQC.whizzPlots.whizzPlot as wpl

groupName = config.groupName


def checkPhase(filename, channel1, channel2):
    """
    For every survey line in a geoWhizz HDF5 file, given two fields, calculate
    the phase shift (in number of samples) required to maximise the correlation
    between the two fields

    Parameters
    ----------
    filename : String
        The name of a geoWhizz HDF5 file.
    channel1 : String
        The name of the first field.
    channel2 : String
        The name of the second field.
    
    Returns
    -------
    None._time_frequency(group)

    """
    with h5py.File(filename, 'r') as f:
        gProject = f[groupName]
        f_units = 's'
        f_sample = util._time_frequency(gProject)
        if f_sample < 0.0001:
            f_sample = 1.0
            f_units = 'samples'
        g = gProject["Lines"]
        numLines = len(g.items())
        offsets = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            linegroup = g[line]
            A = gw.getLineData(linegroup, channel1)
            B = gw.getLineData(linegroup, channel2)
            nsamples = len(A)

            # regularize datasets by subtracting mean and dividing by s.d.
            A -= np.mean(A)
            A /= np.std(A)
            B -= np.mean(B)
            B /= np.std(B)

            # Put in an artificial time shift between the two datasets
            #time_shift = 20
            #A = numpy.roll(A, time_shift)

            # Find cross-correlation
            xcorr = correlate(A, B)
            
            # delta time array to match xcorr
            dt = np.arange(1-nsamples, nsamples)
            recovered_time_shift = dt[xcorr.argmax()]
            offsets[count] = recovered_time_shift
            count += 1

            time = np.arange(dt[0], dt[-1], 0.1) 
            # Now interpolate through gaps by cubic spline
            (xcorrInt, _) = util._interpolateLine(dt, xcorr, time)
            recovered_time_shift2 = time[xcorrInt.argmax()] / f_sample
            print(f'Line {line}: Recovered time shift = {recovered_time_shift2:.1f} {f_units}.')
            
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(time[1:] / f_sample, xcorrInt)
        plt.ylabel('Correlation [arbitrary units]', fontsize = 6)
        plt.xlabel(f'FID difference [{f_units}]', fontsize = 6)

        figtitle = wpl.make_plot_title(gProject)
        plotTitle = f'{figtitle} Line {line}: Correlation of {channel1} v {channel2}'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        plt.show()
        print(f'Offset MEAN = {np.mean(offsets):.2f}; STD = {np.std(offsets):.3f}')
        