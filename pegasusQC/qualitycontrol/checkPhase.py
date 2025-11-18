#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Report the phase shift required to match two channels.

Author: Mark Helm Dransfield

Created: ca 2023

License: CC BY-SA
"""

import numpy as np
import h5py
from scipy.signal import correlate
from scipy import interpolate
import matplotlib.pyplot as plt

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd

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

        if lines == []:
            lines = list(g.keys())
        numLines = len(lines)
        offsets = np.zeros((numLines,))
        count = 0
        # TODO: trap case where tChannel does not exist and report "number of samples"
        if tChannel == '':
            tChannel = f[groupName]['CoordinateFrame'].attrs['TimeChannel']

        for line in lines:
            print(f'analysing {line} of {numLines}.')
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
            xcorrInt = _interpolateCorr(dt, xcorr, time)
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
        

def _interpolateCorr(xbase, ybase, xnew):
    """
    Interpolates `ybase`, sampled at `xbase`, onto the new vector `xnew`.
    These three input arrays are pre-processed to ensure that `xbase`
    and `xnew` are monotonically increasing, whilst keeping `ybase`
    synchronised with `xbase`.

    For use interpolating phase correlations in `checkPhase`.

    Parameters
    ----------
    xbase : 1D numpy float array

        The independent variable of the inputs to be interpolated.

    ybase : 1D numpy float array

        The dependent variable of the inputs to be interpolated.

    xnew : 1D numpy float array

        The independent variable to interpolate onto.

    Returns
    -------
    out : 1D numpy float array

        The values of ybase interpolated onto xnew.

    """
    # clean out 'nan's'
    good = ~np.isnan(xbase + ybase)
    xbase = xbase[good]
    ybase = ybase[good]
    if xbase.size < 10:
        print(f'ERROR - after trimming NaNs, the data vectors are too short for analysis.')
        return

    # ensure ordered in increasing x
    if xbase[1] < xbase[0]:
        xbase = xbase[::-1]
        ybase = ybase[::-1]
    if xnew[1] < xnew[0]:
        xnew = xnew[::-1]
    
    # trim varying data and store
    keepsml = xbase < xnew[-1]
    keepbig = xbase > xnew[0]
    keep = keepsml & keepbig
    xbase = xbase[keep]
    ybase = ybase[keep]
    if xbase.size < 10:
        print(f'ERROR - after trimming measured data, the vectors are too short for analysis.')
        return

    # trim base data and store
    keepsml = xnew < xbase[-1]
    keepbig = xnew > xbase[0]
    keep = keepsml & keepbig
    xnew = xnew[keep]
    if xnew.size < 10:
        print(f'ERROR - after trimming plan data, the vectors are too short for analysis.')
        return

    spl = interpolate.splrep(xbase, ybase, k=3, s=0)
    out = interpolate.splev(xnew, spl)

    return out

        
