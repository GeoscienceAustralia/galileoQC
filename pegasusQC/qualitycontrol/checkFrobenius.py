#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Checks the FTG Frobenius norm against a specification.
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import pegasusQC.config as config
import pegasusQC.whizzFiles.retrieveData as rd
import pegasusQC.utility.utility as util

groupName = config.groupName
        
        
def checkFrobenius(whizzFile, 
    il1='Inline1_raw', il2='Inline2_raw', il3='Inline3_raw', 
    cr1='Cross1_raw', cr2='Cross2_raw', cr3='Cross3_raw', 
    noiselimit=30.0, lines=[], verbose=True, plot_flag=False):
    """
    Reports the noise for each line in lines from filename (a whizz file)
    which exceeds noiseLimit. Here the noise is calculated by `_FTGeigen` as
    the Frobenius norm.

    Parameters
    ----------
    whizzFile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension, to be checked.
    noiselimit : Float
        The maximum allowable in-line noise on a line.
    lines : String list, optional.
        The line numbers to be checked. Default is all lines in the whizzFile.

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
            i1 = rd.getLineData(g[line], il1)
            i2 = rd.getLineData(g[line], il2)
            i3 = rd.getLineData(g[line], il3)
            c1 = rd.getLineData(g[line], cr1)
            c2 = rd.getLineData(g[line], cr2)
            c3 = rd.getLineData(g[line], cr3)
            (Gxx, Gxy, Gxz, Gyy, Gyz, Gzz) = _FTGTransform(i1, i2, i3, c1, c2, c3)
            Txx = util._butter_bandpass_filter(Gxx, 0.1, 0.49, 1.0, order = 6)
            Txy = util._butter_bandpass_filter(Gxy, 0.1, 0.49, 1.0, order = 6)
            Txz = util._butter_bandpass_filter(Gxz, 0.1, 0.49, 1.0, order = 6)
            Tyy = util._butter_bandpass_filter(Gxy, 0.1, 0.49, 1.0, order = 6)
            Tyz = util._butter_bandpass_filter(Gyz, 0.1, 0.49, 1.0, order = 6)
            Tzz = util._butter_bandpass_filter(Gzz, 0.1, 0.49, 1.0, order = 6)
            noise = _FTGeigen(Txx, Txy, Txz, Tyy, Tyz, Tzz, line, noiselimit, plot_flag=plot_flag)
            if verbose:
                print(f'Check line {line}. Noise = {noise:.1f}')
    
    
def _FTGeigen(Txx, Txy, Txz, Tyy, Tyz, Tzz, line = "", noiselimit=30.0, plot_flag=False):
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

    if plot_flag and np.std(frob) > noiselimit:
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


