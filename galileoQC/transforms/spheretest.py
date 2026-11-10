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

def sphereSurvey(numstns=None):
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
    if numstns is None:
        numstns = 20
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

    a_nn, a_ne, a_nd, a_ee, a_ed, a_dd, a_uv, a_D = sphere(xm, xs, rs, ps)

    # sphere 2 parameters and calculate

    ps = 4
    rs = 80
    xs = np.array((0.0, 230.0, 0.0))

    a_2nn, a_2ne, a_2nd, a_2ee, a_2ed, a_2dd, a_2uv, a_D2 = sphere(xm, xs, rs, ps)

    return xm[:,0], xm[:,1], xm[:,2], a_nn+a_2nn, a_ne+a_2ne, a_nd+a_2nd, a_ee+a_2ee, a_ed+a_2ed, a_dd+a_2dd, a_uv+a_2uv, a_D+a_D2


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
    GM_u = G * M * 1E6						# um/s/s

    nO = len(xm)                            # the number of observations
    G_nn = np.zeros(nO)
    G_ne = np.zeros(nO)
    G_nd = np.zeros(nO)
    G_ee = np.zeros(nO)
    G_ed = np.zeros(nO)
    G_dd = np.zeros(nO)
    G_uv = np.zeros(nO)
    gD = np.zeros(nO)

    for n in range(0, nO):
        dx = xm[n, :] - xs
        if n < 3:
            print(dx)
        r = np.sqrt(dx[0]**2 + dx[1]**2 + dx[2]**2)
        G_ne[n] = GM3 * dx[0] * dx[1] / r**5
        G_nd[n] = GM3 * dx[0] * dx[2] / r**5
        G_ed[n] = GM3 * dx[1] * dx[2] / r**5
        G_uv[n] = GM3 * (dx[0]**2 - dx[1]**2) / r**5 / 2
        G_dd[n] = GM * (3 * dx[2]**2 - r**2) / r**5
        G_ee[n] = GM * (3 * dx[1]**2 - r**2) / r**5
        G_nn[n] = GM * (3 * dx[0]**2 - r**2) / r**5
        gD[n] = - GM_u * dx[2] / r**3

    return G_nn, G_ne, G_nd, G_ee, G_ed, G_dd, G_uv, gD


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
                        'Units': 'm'
                    }
                ),
                n_chan: xr.DataArray(
                    data = n,
                    coords={'fiducials': fiducials}, 
                    dims = ['fiducials'],
                    attrs = {
                        'Units': 'm'
                    }
                ),
                z_chan: xr.DataArray(
                    data = data,
                    coords={'fiducials': fiducials}, 
                    dims = ['fiducials'],
                    attrs = {
                        'Units': units
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
