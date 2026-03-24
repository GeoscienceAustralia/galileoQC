#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot flight-line location over scatter plot of Australian ground gravity data.

Author: Mark Helm Dransfield

Created: 2023

License: CC BY-SA
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import xarray as xr
import pooch

import galileoQC.config as config
import galileoQC.whizzFiles.retrieveData as rd

groupName = config.groupName


def plotLinesOnGroundStns(whizzFile, line, minlon=-360, maxlon=360, minlat=-90, maxlat=90, min_reliability=0, fig_title=''):
    """
    Plots the location of the line (in latitude, longitude) overlain on
    scatter plots of the ground gravity station information.

    Parameters
    ----------
    whizzFile : Path

        Each element is the name of a HDF5 Whizz file, including path and extension.

    line : String

        Line identifier for line whose location is to be plotted.

    minlon : Float, optional

        Minimum longitude of data to be extracted. Defaults to -360 deg.

    maxlon : Float, optional

        Maximum longitude of data to be extracted.

    minlat : Float, optional

        Minimum latitude of data to be extracted.

    maxlat : Float, optional

        Maximum latitude of data to be extracted.

    min_reliability : Float, optional

        Only stations with reliability > `min_reliability` will be plotted. Default 0.

    fig_title : String, optional

        The figure title.


    Returns
    -------
    None.

    """

    # Get the coordinates of the points on the flight-line.
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        latitude = f[groupName]['CoordinateFrame'].attrs['LatitudeChannel']
        longitude = f[groupName]['CoordinateFrame'].attrs['LongitudeChannel']
        g = f[groupName]['Lines']
        if fig_title == '':
            fig_title = f[groupName].attrs['ProjectName'] + ': Line(s) on Ground Gravity'
        lon = rd.getLineData(g[line], longitude)[0:]
        lat = rd.getLineData(g[line], latitude)[0:]

    # Extract geographic range of data
    margin = 0.15
    if minlat == -90:
        minlat = np.min(lat) - margin
    if maxlat == 90:
        maxlat = np.max(lat) + margin
    if minlon == -360:
        minlon = np.min(lon) - margin
    if maxlon == 360:
        maxlon = np.max(lon) + margin

    # ensure coords are in Australia
    if (minlat > -11.0) | (minlat < -44.0):
        print(f'Min lat is {minlat} but must be in range [-11, -44]')
        return
    if (maxlat > -11.0) | (maxlat < -44.0):
        print(f'Max lat is {maxlat} but must be in range [-11, -44]')
        return
    if (minlon < 113.0) | (minlon > 154.0):
        print(f'Min lon is {minlon} but must be in range [113, 154]')
        return
    if (maxlon < 113.0) | (maxlon > 154.0):
        print(f'Min lon is {maxlon} but must be in range [113, 154]')
        return
    
    # ensure coords are ordered
    if minlon >= maxlon:
        print(f'{minlon} is not < {maxlon} as required')
        return
    if minlat >= maxlat:
        print(f'{minlat} is not < {maxlat} as required')
        return
    
    # get ground station data
    doi = "10.6084/m9.figshare.13643837"
    # Known MD5 checksum (from figshare)
    checksum = "md5:16c94a792003714efee2bdb4f3181d3a"
    fname = pooch.retrieve(
        url=f"doi:{doi}/australia-ground-gravity.nc",
        known_hash=checksum,
    )
    data = xr.load_dataset(fname)
    ds_5371 = data.where((data.longitude > minlon) 
        & (data.longitude < maxlon) 
        & (data.latitude > minlat) 
        & (data.latitude < maxlat) 
        & (data.reliability_index > min_reliability))
    new = ds_5371.dropna('point')
    
    SMALL_SIZE = 12
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 16

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

    fig = plt.figure(figsize=(14,7))
    fig.suptitle(fig_title, fontsize=BIGGER_SIZE)
    fig.subplots_adjust(top=0.95)
    
    ax = fig.add_subplot(2,2,1)
    new.plot.scatter(x='longitude',y='latitude',s=8,hue='gravity_accuracy',ax=ax)
    ax.plot(lon, lat, 'b', lw=1.0)
    ax.axis('equal')
    plt.grid()
    ax = fig.add_subplot(2,2,2)
    new.plot.scatter(x='longitude',y='latitude',s=8,hue='reliability_index',ax=ax)
    ax.plot(lon, lat, 'b', lw=1.0)
    ax.axis('equal')
    plt.grid()
    ax = fig.add_subplot(2,2,3)
    new.plot.scatter(x='longitude',y='latitude',s=8,hue='height_error',ax=ax)
    ax.plot(lon, lat, 'b', lw=1.0)
    ax.axis('equal')
    plt.grid()
    ax = fig.add_subplot(2,2,4)
    new.plot.scatter(x='longitude',y='latitude',s=8,hue='gravity',ax=ax)
    ax.plot(lon, lat, 'b', lw=1.0)
    ax.axis('equal')
    plt.grid()
    plt.tight_layout()
                                       
   
