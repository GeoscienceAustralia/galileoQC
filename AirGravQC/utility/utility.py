#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Aug 14 20:28:20 2021

@author: markdransfield

"""

# import necessary modules
import numpy as np
from scipy.signal import butter, lfilter


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


def _space_thou(x, pos):
    """
      Formatter function takes tick label and tick position

    Parameters
    ----------
    x : 
        tick label.
    y : 
        tick position.

    Returns
    -------

    Integer number formatted as string with spaces as a thousands separator.

    """
    s = '%d' % x
    groups = []
    while s and s[-1].isdigit():
        groups.append(s[-3:])
        s = s[:-3]
    return s + ' '.join(reversed(groups)) #u'\2009'


def _get_lineName(linegroup):
    """
    Returns the line number and, if extant, the flight number as a
    formatted string, suitable for reporting.

    Parameters
    ----------
    linegroup : HDF5 group
        The subgroup containing the line.

    Returns
    -------
    lineName : String
        The line number and flight number formatted like LLLL.LL:FFF.

    """
    smallnumber = 0.0001
    lineNo = linegroup.attrs['LineNumber']
    if lineNo - int(lineNo) < smallnumber:
        lineName = f'{lineNo:.0f}'
    elif lineNo * 10 - int(lineNo * 10) < smallnumber:
        lineName = f'{lineNo:.1f}'
    else:
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
    if maxCounter < 1 and maxDistance < 1.0E-6:
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


def getPlannedLine(gPlan, gLineMeas):
    """
    Given the single line group, `gLineMeas`, for a line in a whizzfile of measured
    data, returns the corresponding line group, 'gpline', from the lines group, 
    `gPlan`, of a whizzfile of planned data.

    If successful, `planLineInPlan` is True, else it is False.

    Parameters
    ----------
    gPlan : HDF5 group
        The ['Lines'] group for all the survey lines in a survey plan file.
    gLineMeas : HDF5 group
        The ['Lines'][line] group for a single survey line in a measured survey file.

    Returns
    -------
    planLineInPlan : Bool
        True if `gpline` found.
    gpline : HDF5 group
        The desired line group in the plan lines group.
    """
    planLineInPlan = False
    gpline = ''
    if not 'PlannedLine' in gLineMeas.attrs:
        print(f'ERROR. No PlannedLine attribute set for line.')
        return None, None
    plannedLineNo = gLineMeas.attrs['PlannedLine'] # a Float double
    for pline in gPlan.keys():
        if gPlan[pline].attrs['LineNumber'] == plannedLineNo:
            planLineInPlan = True
            gpline = pline
            break
    return planLineInPlan, gpline
 

def getLineNumberInGroup(linesgroup, lineNumber):
    """
    Given the single line group, `gLineMeas`, for a line in a whizzfile of measured
    data, returns the corresponding line group, 'gpline', from the lines group, 
    `gPlan`, of a whizzfile of planned data.

    If successful, `lineInGroup` is True, else it is False.

    Parameters
    ----------
    linesgroup : HDF5 group
        The ['Lines'] group for all the survey lines in a whizzfile.
    lineNumber : HDF5 group
        The ['Lines'][line] group for a single survey line in a whizzfile.

    Returns
    -------
    lineInGroup : Bool
        True if `gpline` found.
    gpline : HDF5 group
        The desired line group in the plan lines group.
    """
    lineInGroup = False
    gpline = ''
    for line in linesgroup.keys():
        if linesgroup[line].attrs['LineNumber'] == lineNumber:
            lineInGroup = True
            gpline = line
            break
    return lineInGroup, gpline
 

def e2norm(x, y):
    """
    Returns the Euclidean norm of `x` and `y`.
    """
    return np.sqrt(x * x + y * y)


def controls_lessthan_1000(all_lines):
    """
    Separates line numbers on the basis that those less than 999.99
    are control lines, and the rest are traverse lines.

    Parameters
    ----------
    all_lines : list of strings
        The line numbers to be separated.

    Returns
    -------
    (trav_strs, ctrl_strs) : two lists of strings

    """
    ctrl_strs = []
    trav_strs = []
    for line in all_lines:
        if float(line) < 999.99:
            ctrl_strs.append(line)
        else:
            trav_strs.append(line)
    return trav_strs, ctrl_strs


def controls_nineteen(all_lines):
    """
    Separates line numbers on the basis that those that start with
    "19" are control lines, and the rest are traverse lines.

    Parameters
    ----------
    all_lines : list of strings
        The line numbers to be separated.

    Returns
    -------
    (trav_strs, ctrl_strs) : two lists of strings

    """
    ctrl_strs = []
    trav_strs = []
    for line in all_lines:
        if line[:2] == "19":
            ctrl_strs.append(line)
        else:
            trav_strs.append(line)
    return trav_strs, ctrl_strs


def convertUTMtoGeo(crs_epsg, x_utm, y_utm):
    """
    Return the latitude and longitude for projected coordinates with a common datum.

    Parameters
    ----------
    crs_epsg : Integer
        The common datum CRS EPSG code.
    x_utm : Float
        The projected x or easting.
    y_utm : Float
        The projected y or northing.

    Returns
    -------
    (latitude, longitude) : two Floats

    """
    from pyproj import CRS
    from pyproj import Transformer

    mycrs = CRS.from_epsg(crs_epsg)
    proj = Transformer.from_crs(mycrs, mycrs.geodetic_crs)
    return proj.transform(x_utm, y_utm)


def _mean_1std(x):
    """
    Calculate the mean of the values in x that fall in the range of +/- stdev(x).
    This is a simplistic "mean excluding outliers".
    """
    mean1 = np.mean(x)
    std1 = np.std(x)
    idx = np.argwhere(np.logical_and(x > (mean1 - std1), x < (mean1 + std1)))
    return np.mean(x[idx])


def _calc_bearing(x, y):
    """
    arctan(mean(diff(x) / mean(diff(y))))
    """
    return np.arctan2(_mean_1std(np.diff(x)), _mean_1std(np.diff(y)))

