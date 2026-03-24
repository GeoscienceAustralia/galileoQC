#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a synthetic data set for testing.

Author: Mark Helm Dransfield

Created: Oct 2025

License: CC BY-SA
"""

import numpy as np
import xarray as xr

def sphereSurvey(numstns=20):
    """
    Creates an airborne survey over 3 point mass sources, calculates all
    gravity acceleration vector and curvature tensor outputs.

    All information is referenced to an NED coordinates frame. Accelerations
    are in um/s/s and curvatures are in eotvos.

    One-off function for testing Craig transform code on synthetic data.
    
    Parameters
    ----------
    numstns : int, optional

        The number of stations (samples) in the grid's x and y directions.
        Default 20.

    Returns
    -------
    n, e, d, Ann, Ane, And, Aee, Aed, Add, Auv, AD : 

        a list of numpy 1D arrays

    """
    # survey parameters - note that the loop order below is important.
    nn = 201
    ne = 201
    de = numstns
    dn = numstns
    oe = -(divmod(ne,2)[0]) * de
    on = -(divmod(nn,2)[0]) * dn
    d = -50

    xm = np.zeros((nn*ne, 3))

    for n in range(0, nn):
        for m in range(0, ne):
            xm[n*ne+m, 0] = on+dn*n
            xm[n*ne+m, 1] = oe+de*m
            xm[n*ne+m, 2] = d

    # sphere 1 parameters and calculate

    ps = 8
    rs = 50
    xs = np.zeros((3))

    Ann, Ane, And, Aee, Aed, Add, Auv, AD = sphere(xm, xs, rs, ps)

    # sphere 2 parameters and calculate

    ps = 4
    rs = 80
    xs = np.array((0.0, 230.0, 0.0))

    A2nn, A2ne, A2nd, A2ee, A2ed, A2dd, A2uv, AD2 = sphere(xm, xs, rs, ps)

    return xm[:,0], xm[:,1], xm[:,2], Ann+A2nn, Ane+A2ne, And+A2nd, Aee+A2ee, Aed+A2ed, Add+A2dd, Auv+A2uv, AD+AD2


def sphere(xm, xs, rs, ps):
    """
    Calculates the vertical gravity, gD, and all gravity gradient
    components (in NED coordinates), on sample points simulating an
    airborne survey above a massive sphere at the centre of the survey area.
    
    One-off function for testing Craig transform code on synthetic data.
    
    Parameters
    ----------
    xm : numpy 2D array

        The NED coordinates of the synthetic survey measurement
        points in metres. Array size is [n,3] for `n` measurements.

    xs : numpy 1D array

        The NED coordinates of the centre of the sphere in metres.

    rs : str

        The radius of the sphere in metres.

    ps : float

        The density of the sphere in gm/cc (== tonnes per cubic metre).

    Returns
    -------
    Gnn, Gne, Gnd, Gee, Ged, Gdd, Guv, gD : a list of numpy 1D arrays

        gravity gradients and vertical acceleration.

    """

    import numpy as np

    G = 6.673E-11							# SI
    M = ps * 4 / 3 * np.pi * rs**3 * 1E3	# kg, *1E3 to convert ps: gm/cc to kg/m^3
    GM3 = G * M * 3 * 1E9					# Eo
    GM = G * M * 1E9						# Eo
    GMu = G * M * 1E6						# um/s/s

    nO = len(xm)                            # the number of observations
    Gnn = np.zeros(nO)
    Gne = np.zeros(nO)
    Gnd = np.zeros(nO)
    Gee = np.zeros(nO)
    Ged = np.zeros(nO)
    Gdd = np.zeros(nO)
    Guv = np.zeros(nO)
    gD = np.zeros(nO)

    for n in range(0, nO):
        dx = xm[n, :] - xs
        if n < 3:
            print(dx)
        r = np.sqrt(dx[0]**2 + dx[1]**2 + dx[2]**2)
        Gne[n] = GM3 * dx[0] * dx[1] / r**5
        Gnd[n] = GM3 * dx[0] * dx[2] / r**5
        Ged[n] = GM3 * dx[1] * dx[2] / r**5
        Guv[n] = GM3 * (dx[0]**2 - dx[1]**2) / r**5 / 2
        Gdd[n] = GM * (3 * dx[2]**2 - r**2) / r**5
        Gee[n] = GM * (3 * dx[1]**2 - r**2) / r**5
        Gnn[n] = GM * (3 * dx[0]**2 - r**2) / r**5
        gD[n] = - GMu * dx[2] / r**3

    return Gnn, Gne, Gnd, Gee, Ged, Gdd, Guv, gD


def _make_xr(data, name, units, n, e, n_chan='northing', e_chan='easting'):
    """
    Puts the data in 1D numpy arrays (n, e, data) into a fiducial-
    dimensioned 1D xarray DataSet. The n and e channels are assumed
    to have units 'm'.
    
    One-off function for testing Craig transform code on synthetic data.
    
    Parameters
    ----------
    data : numpy 1D array

        Airborne survey line-based data.

    name : str

        Name of the `data`.

    units : str

        Name of the units of the `data`.

    n : numpy 1D array

        The northings (or 'y' coordinates) for the `data`. Must be in metres.

    e : numpy 1D array

        The eastings (or 'x' coordinates) for the `data`. Must be in metres.

    n_chan : str, optional

        The name to be given to the `n` data, usually 'northing' or 'y'.
        Default 'northing'.

    e_chan : str, optional

        The name to be given to the `e` data, usually 'easting' or 'x'.
        Default 'easting'.

    Returns
    -------
    ds : xarray 1D DataSet

        DataSet containing the provided data, dimensioned by fiducial.

    """
    z_chan = name
    fiducials = np.arange(0, len(e))
    ds = xr.Dataset({
                e_chan: xr.DataArray(
                    data = e,
                    coords={'fiducials': fiducials}, 
                    dims = ['fiducials'],
                    attrs = {
                        'units': 'm'
                    }
                ),
                n_chan: xr.DataArray(
                    data = n,
                    coords={'fiducials': fiducials}, 
                    dims = ['fiducials'],
                    attrs = {
                        'units': 'm'
                    }
                ),
                z_chan: xr.DataArray(
                    data = data,
                    coords={'fiducials': fiducials}, 
                    dims = ['fiducials'],
                    attrs = {
                        'units': units
                    }
                )},
                attrs = {
                'author': 'Mark Dransfield',
                'x_channel': e_chan,
                'y_channel': n_chan,
                'z_channel': z_chan,
                'title': z_chan
                }
            )
    return ds
