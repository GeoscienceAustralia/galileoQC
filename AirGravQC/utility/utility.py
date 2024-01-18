#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Aug 14 20:28:20 2021

@author: markdransfield

"""

# import necessary modules
import numpy as np
from scipy.signal import butter, lfilter
from scipy.interpolate import CloughTocher2DInterpolator
from scipy import interpolate


def _distance(x, y):
    return np.sqrt(x * x + y * y)


def _length(x, y):
    """
    Returns a 1D numpy array where the i-th element is the length of the vector from
    (x[0], y[0]) to (x[i], y[i]).

    Parameters
    ----------
    x : numpy Float array
        DESCRIPTION.
    y : numpy Float array
        DESCRIPTION.

    Returns
    -------
    numpy Float array
        DESCRIPTION.

    """
    return  np.sqrt((x - x[0]) * (x - x[0]) + (y - y[0]) * (y - y[0]))


def _displacement2(x0, x1, y0, y1):
    """
    The displacement distance from (xo, y0) to (x1, y1).

    Parameters
    ----------
    x0 : TYPE
        DESCRIPTION.
    x1 : TYPE
        DESCRIPTION.
    y0 : TYPE
        DESCRIPTION.
    y1 : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return np.sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1))
            

def _butter_bandpass(lowcut, highcut, fs, order=5):
    """
    Calculates the parameters of a Butterworth bandpass filter.

    Parameters
    ----------
    lowcut : Float
        The lowpass cutoff frequency or wavenumber.
    highcut : Float
        The highpass cutoff frequency or wavenumber.
    fs : Float
        The sample frequency or wavenumber.
    order : Integer, optional
        The filter order. Default 5.

    Returns
    -------
    b : Float
        DESCRIPTION.
    a : Float
        DESCRIPTION.

    """
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def _butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    """
    Uses a Butterworth bandpass filter [`lowcut`, `highcut`] to filter data
    (sampled at `fs`).

    Parameters
    ----------
    data : Numpy 1D array
        The input data.
    lowcut : Float
        The lowpass cutoff frequency or wavenumber.
    highcut : Float
        The highpass cutoff frequency or wavenumber.
    fs : Float
        The sample frequency or wavenumber.
    order : Integer, optional
        The filter order. Default 5.

    Returns
    -------
    Numpy 1D array
        The output data.

    """
    b, a = _butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def _space_thou(x, pos):  # formatter function takes tick label and tick position
    s = '%d' % x
    groups = []
    while s and s[-1].isdigit():
        groups.append(s[-3:])
        s = s[:-3]
    return s + ' '.join(reversed(groups)) #u'\2009'


def _get_lineName(linegroup):
    lineNo = linegroup.attrs['LineNumber']
    lineName = f'{lineNo:.2f}'
    if 'Flight' in linegroup.attrs:
        flightNo = linegroup.attrs['Flight']
        lineName += ":" + f'{flightNo:.0f}'
    return lineName


def _rotateCoords(x, y, angle):
    """
    Rotates Cartesian vectors x and y, by `angle` to xr and yr.

    Parameters
    ----------
    x : numpy Float 1D array
        DESCRIPTION.
    y : numpy Float 1D array
        DESCRIPTION.
    angle : Float
        Angle in radians by which coordinates are to be rotated.

    Returns
    -------
    xr : numpy Float 1D array
        DESCRIPTION.
    yr : numpy Float 1D array
        DESCRIPTION.

    """
    xr = np.cos(angle) * x + np.sin(angle) * y
    yr = np.sin(angle) * x - np.cos(angle) * y
    return xr, yr


def _exceedance_fail(num_fids_in_exceedance, len_exceedance, maxCounter, maxDistance):
    """
    Given a specification, maxCounter, on number of fids and, maxDistance,
    on distance, checks to see if either  num_fids_in_exceedance > maxCounter
    or len_exceedance > maxDistance. If either is True, returns True.

    Parameters
    ----------
    num_fids_in_exceedance : Integer
        DESCRIPTION.
    len_exceedance : Float
        DESCRIPTION.
    maxCounter : Integer
        DESCRIPTION.
    maxDistance : Float
        DESCRIPTION.

    Returns
    -------
    Bool.

    """
    if maxCounter < 1 and maxDistance < 1:
        return False
    if maxCounter > 0:
        if num_fids_in_exceedance > maxCounter:
            return True
    if maxDistance > 0:
        if len_exceedance > maxDistance:
            return True
    return False


def _failsDeviation(data, limit, nSamples):
    """
    Checks data, testing for more than `nSamples` consecutive instances where
    the data values are outside the range [-limit, limit].
    
    Parameters
    ----------
    data : Numpy 1D Float array
        The data to be checked.
    limit : Float
        The limit, assumed > 0.
    nSamples : Float
        Maximum allowed number of deviations allowed.

    Returns
    -------
    Bool
        True if the data failed.
    Float
        The number of exceedances (regardless of whether consecutive).
    
    """
    
    if np.max(data) < limit and np.min(data) > -limit:
        # print('all ok')
        return False, 0
    
    indxSpikes = np.argwhere(np.logical_or(data > limit, data < -limit))
    numExceedances = indxSpikes.shape[0]
    if numExceedances < nSamples:
        # print(f'only {numExceedances} exceedances')
        return False, numExceedances
    
    else:
        count = 1
        for jj in range(1, len(indxSpikes)):
            if indxSpikes[jj][0] == indxSpikes[jj-1][0] + 1:
                count += 1
            else:
                if count >= nSamples:
                    return True, count
                else:
                    count = 1
        if count >= nSamples:
            return True, count
        # print(f'{numExceedances} but not consecutive, most {count}')
        return False, numExceedances
        

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

    return _butter_bandpass_filter(ils, lowcut, highcut, fs, order = order)


def _interpolateLine(timeIn, dataIn, timeOut, spare=[], plot_flag=False):
    """
    Interpolates dataIn, sampled at timeIn, onto the samples timeOut.
    These three input arrays are pre-processed to ensure that timeIn
    and timeOut are monotonically increasing, whilst keeping dataIn
    synchronised with timeIn, and spareOut synchronised with timeOut.

    Parameters
    ----------
    timeIn : 1D numpy float array
        The independent variable of the inputs to be interpolated.
    dataIn : 1D numpy float array
        The input dependent variable to be interpolated.
    timeOut : 1D numpy float array
        The independent variable to interpolate onto.
    spare : 1D numpy float array
        To be kept synchronised with timeOut and returned.
    plot_flag : Bool, optional
        If True, plot the np.diff() of dataIn and timeIn.
        Default False.

    Returns
    -------
    out : 1D numpy float array
        The values of dataIn interpolated onto timeOut.
    spareOut : 1D numpy float array
        To be kept synchronised with timeOut and returned.

    """
    
    timeIn = timeIn[np.logical_not(np.isnan(dataIn))]
    dataIn = dataIn[np.logical_not(np.isnan(dataIn))]
    spareOut = spare
    # print(np.min(np.diff(timeIn)), np.max(np.diff(timeIn)))
    # model = interpolate.InterpolatedUnivariateSpline(timeIn, dataIn)
    # return model(timeOut)
    
    min_length = 100
    if timeIn.size < min_length or dataIn.size < min_length or timeOut.size < min_length:
        out = np.zeros(timeOut.shape)
        out[:] = np.nan
        return out, spare
    if len(spare) < min_length:
        spare = []

    (timeIn_trim, dataIn_trim) = _trim_monotonic(timeIn, sync=dataIn)
    (timeOut_trim, spare_trim) = _trim_monotonic(timeOut, sync=spare)
    if len(spare) > min_length:
        spareOut = spare_trim

    if plot_flag:
        print(f'Len t: {len(timeIn)}; Len d: {len(dataIn)}')
        print(f'')
        plt.plot(timeIn[1:], np.diff(dataIn), 'b', timeIn[1:], np.diff(timeIn), 'g')
        plt.show()
    spl = interpolate.splrep(timeIn_trim, dataIn_trim, k=3, s=0)
    # out = interpolate.splev(timeOut, spl)
    out = interpolate.splev(timeOut_trim, spl)
    return out, spareOut


def _trim_monotonic(data_in, sync=[]):
    """
    Force input to be monotonic by removing samples; keep
    sync aligned by removing corresponding samples there 
    as well.

    """
    b = np.diff(data_in)
    # print(len(b[b > 0]), len(b[b < 0]))
    if len(b[b > 0]) < len(b[b < 0]):
        # if it is mostly decreasing, reverse, ...
        data_temp = data_in[::-1]
        if sync.size != 0:
            sync_temp = sync[::-1]
    else:
        data_temp = data_in
        if len(sync) != 0:
            sync_temp = sync
    # Now mostly increasing, keep only the increases
    b = np.concatenate(([0],np.diff(data_temp)))
    data_trim = data_temp[b > 0]
    if len(sync) != 0:
        sync_out = sync_temp[b > 0]
    else:
        sync_out = []
    return data_trim, sync_out






def _time_frequency(group):
    """
    Returns the sample frequency of the data in the project group.
    Simply calculated as the inverse of the difference of the first
    two sample times in the first line. Relies on the TimeChannel
    attribute being set. 

    Parameters
    ----------
    group : HDF5 group
        Must be the project group (top level of hierarchy in whizz File).

    Returns
    -------
    Float. The sample frequency in Hz.

    """
    flightLines = list(group['Lines'].keys())
    if 'TimeChannel' in group['CoordinateFrame'].attrs:
        time = group['CoordinateFrame'].attrs['TimeChannel']
        t = np.array(group['Lines'][flightLines[0]][time])
        return 1.0 / np.abs(t[1] - t[0])
    else:
        return 0.0


