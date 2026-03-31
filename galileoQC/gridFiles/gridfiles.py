#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rarely used but possibly useful grid functions.

Author: Mark Helm Dransfield

Created: Sat Jul 18 16:43:31 2020

License: CC BY-SA
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from pathlib import Path
import colorcet as cc
import xarray as xr
import netCDF4 as nc4
import filebrowser as fb
# import rioxarray
import h5py
import matplotlib.ticker as tkr

import galileoQC.gridFiles.graphics as graphics
import galileoQC.utility.utility as util
import galileoQC.gridFiles.read_ers as ers
import galileoQC.whizzFiles.retrieveData as rd
import galileoQC.config as config
from galileoQC.gridFiles.graphicsShaded import graphicsShaded
from galileoQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from galileoQC.gridFiles.xarray_to_grid import xarray_to_grid
from galileoQC.gridFiles.xdImage import xdImage
import galileoQC.gridFiles.gridutility as gut

groupName = config.groupName
projectName = config.projectName

# GRID REPORTING

def subtractImages(imagefile1, imagefile2, scale=1.0, band1=0, band2=0):
    """
    Takes in the paths to two grid files Each may be a single band ERMapper image
    file or a netCDF4 grid file) and returns the difference as an xarray:

    `imagefile1 - scale x imagefile2`.

    If the images do not cover the same space, then the return flag is set False
    and only the first image is returned.

    Parameters
    ----------
    imagefile1 : Path

        The Path to the first grid file, must have extension `ers` or `nc`.

    imagefile2 : Path

        The Path to the second grid file, must have extension `ers` or `nc`.

    scale : (Float, optional)

        imagefile2 is multiplied by scale before subtraction. Default 1.0.

    Returns
    -------
    flag : Bool

        True if operation was successful.

    x : xarray

        The difference if both span the same coordinates, otherwise the first.

    """
    
    x1, _ = gridfile_to_xa(imagefile1, bandout=band1)
    x2, _ = gridfile_to_xa(imagefile2, bandout=band2)
    
    try:
        xr.testing.assert_equal(x1.coords, x2.coords) #x1.coords.all() == x2.coords.all():
        return True, x1 - scale * x2
    except:
        return False, x1


def gridfile_to_xr(whizzFile='', bandout=0):
    """
    Returns an xarray Dataset containing the geographically located data from a
    gridfile in either ERS or NC format.

    Parameters
    ----------
    whizzFile : Path

        The Path to the grid file, must have extension `ers` or `nc`.

    bandout : Int, optional

        The band to be read if the grid file is `ERS`. The default is 0.

    Returns
    -------
    xd : xarray Dataset

        The data from `whizzFile`.

    """
    if whizzFile == '':
        whizzFile = Path(fb.get_grid_filename(filetypes=(('NetCdf4 grid', '*.nc'),
                                                         ('ERMapper grid', '*.ers'))))
    if whizzFile.suffix.upper() == '.ERS':
        e, n, z, datum, projection, units = ers.read_ers_image(whizzFile, bandout=bandout)
        xa = xr.DataArray(data=z,#np.flip(z, 0), # DANGER!!!
                          dims=["N", "E"],
                          coords={"N": n, "E": e})
        xa.dropna(dim='N',how='all')
        xa.dropna(dim='E',how='all')
        fname = whizzFile.with_suffix('').name
        xd = xr.Dataset(data_vars={fname: xa})
        xd.attrs["long_name"] = fname
        xd.attrs["datum"] = datum
        xd.attrs["projection"] = projection
        if datum == 'WGS84' and projection == 'SUTM55':
            xd.rio.write_crs("epsg:32755", inplace=True)
        xd.attrs["Units"] = units
    elif whizzFile.suffix.upper() == '.NC':
        nc = nc4.Dataset(str(whizzFile), mode='r')
        xd = xr.open_dataset(xr.backends.NetCDF4DataStore(nc))
    else:
        print(f'ERROR: whizzFile has suffix {Path(whizzFile).suffix.upper()} but must be `.nc` or `.ers`.')
    return xd, whizzFile


def gridfile_to_xa(whizzFile='', bandout=0):
    """
    Returns an xarray DataArray containing the geographically located data from a
    gridfile in either ERS or NC format.

    Parameters
    ----------
    whizzFile : Path

        The Path to the grid file, must have extension `ers` or `nc`.

    bandout : Int, optional

        The band to be read if the grid file is `ERS`. The default is 0.

    Returns
    -------
    xd : xarray Dataset

        The data from `whizzFile`.

    """
    if whizzFile == '':
        whizzFile = Path(fb.get_grid_filename(filetypes=(('NetCdf4 grid', '*.nc'),
                                                         ('ERMapper grid', '*.ers'))))
    if whizzFile.suffix.upper() == '.ERS':
        e, n, z, datum, projection, units = ers.read_ers_image(whizzFile, bandout=bandout)
        xa = xr.DataArray(data=np.flip(z, 0), # DANGER!!!
                          dims=["N", "E"],
                          coords={"N": n, "E": e})
        xa.dropna(dim='N',how='all')
        xa.dropna(dim='E',how='all')
        fname = whizzFile.with_suffix('').name
        xa.attrs["long_name"] = fname
        xa.attrs["datum"] = datum
        xa.attrs["projection"] = projection
        if datum == 'WGS84' and projection == 'SUTM55':
            xa.rio.write_crs("epsg:32755", inplace=True)
        xa.attrs["Units"] = units
    elif whizzFile.suffix.upper() == '.NC':
        # nc = nc4.Dataset(str(whizzFile), mode='r')
        xa = xr.load_dataarray(str(whizzFile))#xr.backends.NetCDF4DataStore(nc))
    else:
        print(f'ERROR: whizzFile has suffix {Path(whizzFile).suffix.upper()} but must be `.nc` or `.ers`.')
    return xa, whizzFile


def ers_to_netcdf4(ersFile='', ncFile='', datum='', projection='', long_name='', units=''):
    """
    Writes an `ERS` grid to a netCDF4 file.

    Parameters
    ----------
    ersFile : pathlib Path, optional

        The ERS grid file to be converted. The default is ''. If '', then file browser allows
        the user to select the file.

    ncFile : pathlib Path, optional

        The NC file to be written. The default is ''. if '', then the name of the ERS file is
        used, with the 'NC' extension.

    datum : Str, optional

        The geographic datum. The default is ''. If '', it is read from the ERS file.

    projection : Str, optional

        The geographic projection. The default is ''. If '', it is read from the ERS file.

    long_name : Str, optional

        The NC Dataset long_name attr. The default is ''. If '', then the ERS filename is used.

    units : Str, optional

        The units of the Dataarray values. The default is ''.

    Returns
    -------
    None.

    """
    if ersFile == '':
        ersFile = Path(fb.get_grid_filename(filetypes=[('ERMapper grid', '*.ers')]))
    if ncFile == '':
        ncFile = ersFile.with_suffix('.nc')
    xd, fname = gridfile_to_xr(ersFile)
    if datum == '':
        datum = xd.attrs['datum']
    if long_name == '':
        long_name = xd.attrs['long_name']
    if projection == '':
        projection = xd.attrs['projection']
    xd.attrs = {'datum' : datum,
                'projection' : projection,
                'long_name' : long_name,
                'Units' : units}
    print(xd)
    xd.to_netcdf(ncFile)
     
    
def update_grid(whizzFile='', datum='', projection='', long_name='', units=''):
    """
    NOT FUNCTIONAL. Adds any specified (ignores the defaults) attributes to a NetCdf4 grid files DataSet.

    Parameters
    ----------
    whizzFile : String or Path, optional

        The name of the netCDF4 grid file, including the `.nc` extension. Defaults to
        putting up a filebrowser dialog so the user can select a file.

    datum : String, optional

        The geographic datum. The default is ''.

    projection : String, optional

        The geographic projection. The default is ''.

    long_name : String, optional

        A descriptive name for the grid. The default is ''.

    units : String, optional

        The units of the grid `z` field. The default is ''.

    Returns
    -------
    None.

    """
    if whizzFile == '':
        whizzFile = Path(fb.get_grid_filename([('NetCdf4 grid', '*.nc')]))
    nc = nc4.Dataset(str(whizzFile), mode='r')
    xa = xr.open_dataset(xr.backends.NetCDF4DataStore(nc))
    nc.close()
    xs = xa
    # xs.attrs = {'datum' : datum,
    #             'projection' : projection,
    #             'long_name' : long_name,
    #             'Units' : units}
    suffix = whizzFile.suffix
    fname = whizzFile.with_suffix('').name + '+'
    newWhizz = whizzFile.with_name(fname).with_suffix(suffix)
    print(xs)
    xs.to_netcdf(newWhizz)
    

def cmap_exists(name):
    """
    Returns True if the Matplotlib Colormap `name` exists.

    Parameters
    ----------
    name : Colormap instance

        The name of the Colormap.

    Returns
    -------
    bool

    """
    try:
         cm.get_cmap(name)
         return True
    except ValueError:
         pass
    return False


def graphicsTernary(e, n, red, green, blue, mytitle):
    """
    Uses Matplotlib's `imshow()` to create a ternary image from three channels
    of located information. Each channel is histogram-equalised before imaging.

    Parameters
    ----------
    e : np.array(Float, 1D)

        The easting vector.

    n : np.array(Float, 1D)

        The northing vector.

    z : np.array(Float, 2D)

        The data to be imaged (referenced to the easting and northing positions).

    red : np.array(Float, 2D)

        The data for the red channel.

    green : np.array(Float, 2D)

        The data for the green channel.

    blue : np.array(Float, 2D)

        The data for the blue channel.

    mytitle : String

        The title for the image plot.

    Returns
    -------
    None.

    """
    # histogram equalise each channel
    red2 = _histEqual(red)
    green2 = _histEqual(green)
    blue2 = _histEqual(blue)
    # stack into one 3D array
    z = np.dstack((red2, green2, blue2))
    
    fig, ax = plt.subplots(figsize=(12,6))
    thou_format = tkr.FuncFormatter(util._space_thou)
    ax.imshow(e, n, z)#, cmap = 'Reds')
    ax.set_title(mytitle, fontsize=12)
    ax.xaxis.set_major_formatter(thou_format)
    ax.yaxis.set_major_formatter(thou_format)
    ax.set_xlabel('Eastings [m]', fontsize=9)
    ax.set_ylabel('Northings [m]', fontsize=9)
    ax.grid(True)
    plt.axis('equal')

    plt.show()

    
def _histEqual(img):
    """
    Performs histogram equalisation on an image.

    Parameters
    ----------
    img : 2D image array

        DESCRIPTION.

    Returns
    -------
    2D image array
    
        DESCRIPTION.

    """
    cleanImg = img[np.logical_not(np.isnan(img))]
    hist,bins = np.histogram(cleanImg.flatten(),256)#,[0,256])
    cdf = hist.cumsum()
    cdf_m = np.ma.masked_equal(cdf,0)
    cdf_m = (cdf_m - cdf_m.min()) / (cdf_m.max() - cdf_m.min()) * 255
    cdf = np.ma.filled(cdf_m,0).astype('uint8')
    img255 = (img - np.min(cleanImg)) / (np.max(cleanImg) - np.min(cleanImg)) * 255.0
    return cdf[img255.astype('uint8')]

    


