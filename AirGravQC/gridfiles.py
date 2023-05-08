"""
Created on Sat Jul 18 16:43:31 2020

@author: markdransfield

A library of geophysical QA/QC functions for gridded data.

"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from pathlib import Path
import colorcet as cc
import xarray as xr
import netCDF4 as nc4
import filebrowser as fb
import rioxarray
import h5py
import pygmt

import AirGravQC.graphics as graphics
import AirGravQC.config as config

groupName = config.groupName
projectName = config.projectName


def whizz_to_xarray(whizz_file, z_chan, n_chan='', e_chan='', remove_mean=False, diff_one=False):
    """
    Return a point-located xArray Dataset of (northing, easting, z), over the `fiducials` dimension,
    from a whizz_file.

    Parameters
    ----------
    whizzFile : Path or String
        The Path to, or String name of, the whizz file in HDF5 format.
    z_chan : String
        The name of the channel in `whizz_file` to be imaged.
    n_chan : String, optional
        The name of the channel in `whizz_file` containing the northings (y).
        Default "" causes the name of the "YChannel" in `whizz_file` to be used.
    e_chan : String, optional
        The name of the channel in `whizz_file` containing the eastings (x).
        Default "" causes the name of the "XChannel" in `whizz_file` to be used.
    remove_mean : Bool, optional
        If true, the mean is subtracted from each survey line of data before
        writing to `my_dataset`. Default False.
    diff_one : Bool, optional
        If true, the first difference along each survey line of data is
        written to `my_dataset`. Default False.

    Returns
    -------
    my_dataset : (xArray Dataset)
        Contains the data. If `whizz_to_xarray()` is unable to return data,
        it returns an empty xArray Dataset (test with `len(aa.attrs) == 0`).

    """
    filename = str(whizz_file)
    
    with h5py.File(filename, 'r') as f:
        lines_group = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if e_chan == '':
            e_chan = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if n_chan == '':
            n_chan = f[groupName]['CoordinateFrame'].attrs['YChannel']

        if z_chan == e_chan or z_chan == n_chan:
            print(f'Cannot process {z_chan}, same as {e_chan} or {n_chan}.')
            return xr.Dataset()
        totalNumFids = 0
        
        for line in lines_group.keys():
            xData = np.array(lines_group[line][e_chan])
            totalNumFids += len(xData)
        print(f'Total number of fids in whizz file = {totalNumFids}.')

        # initialise an xarray to take the data, with name and units
        fiducials = np.arange(0, totalNumFids)
        data = np.zeros((totalNumFids,))
        eastings = np.zeros((totalNumFids,))
        northings = np.zeros((totalNumFids,))

        my_dataset = xr.Dataset({
            e_chan: xr.DataArray(
                data = eastings,
                coords={'fiducials': fiducials}, 
                dims = ['fiducials'],
                attrs = {
                    'units': 'm'
                }
            ),
            n_chan: xr.DataArray(
                data = northings,
                coords={'fiducials': fiducials}, 
                dims = ['fiducials'],
                attrs = {
                    'units': 'm'
                }
            ),
            z_chan: xr.DataArray(
                data = np.zeros((totalNumFids,)),
                coords={'fiducials': fiducials}, 
                dims = ['fiducials'],
                attrs = {
                    'units': '-'
                }
            )},
            attrs = {
            'author': 'Mark Dransfield',
            'x_channel': e_chan,
            'y_channel': n_chan,
            'z_channel': z_chan
            }
        )

        sfid = 0
        efid = 0
        for line in lines_group.keys():
            sfid = efid
            xData = np.array(lines_group[line][e_chan])
            yData = np.array(lines_group[line][n_chan])
            zData = np.array(lines_group[line][z_chan])
            efid += len(yData)
            if remove_mean:
                zData = zData - np.mean(zData)
            if diff_one:
                zData = np.append(np.diff(zData), zData[-1]-zData[-2])

            my_dataset[e_chan][sfid:efid] = xData
            my_dataset[n_chan][sfid:efid] = yData
            my_dataset[z_chan][sfid:efid] = zData

            if 'Units' in lines_group[line][e_chan].attrs:
                my_dataset[e_chan].attrs['units'] = lines_group[line][e_chan].attrs["Units"]
            if 'Units' in lines_group[line][n_chan].attrs:
                my_dataset[n_chan].attrs['units'] = lines_group[line][n_chan].attrs["Units"]
            if 'Units' in lines_group[line][z_chan].attrs:
                my_dataset[z_chan].attrs['units'] = lines_group[line][z_chan].attrs["Units"]
            if remove_mean and diff_one:
                my_dataset.attrs['title'] = f'{z_chan} (mr) (d1)'
            elif remove_mean:
                my_dataset.attrs['title'] = f'{z_chan} (mr)'
            elif diff_one:
                my_dataset.attrs['title'] = f'{z_chan} (d1)'
            else:
                my_dataset.attrs['title'] = z_chan
            
    return my_dataset


def xarray_to_grid(my_data, grid_space):
    """
    Uses `PyGMT` to interpolate the `my_data` onto a regular grid. Method
    uses data to 5 x the grid spacing to focus on QC issues.
    TODO:
    1. allow change in size of grid spacing.

    Parameters
    ----------
    my_data : xArray Dataset
        Contains x, y, and z data dimensioned by fiducial.

    Returns
    -------
    grid : 
        Contains the gridded data.
    region : tuple
        Four floating point values for xmin, xmax, ymin, ymax.
    grid_space : Float
        The distance between grid cell centres in grid distance units.

    """

    # grid_space = 500.0

    x_chan = my_data.attrs['x_channel']
    y_chan = my_data.attrs['y_channel']
    z_chan = my_data.attrs['z_channel']
    myunits = my_data[z_chan].attrs['units']
    print(f'Processing (x, y, z) = ({x_chan}, {y_chan}, {z_chan}). {z_chan} in {myunits}.')

    # grid spacing and search radius
    inspc = f'{grid_space:.0f}+e'
    maxradius = 5.0 * grid_space

    region = pygmt.info(data=my_data[[x_chan, y_chan]], spacing=1)  # West, East, South, North
    print(f"Data points cover region: {region}")

    # Preprocess z_chan data using blockmedian
    data_trm = pygmt.blockmedian(
        data=my_data[[x_chan, y_chan, z_chan]],
        spacing=inspc,
        region=region,
    )

    # make grid
    grid = pygmt.surface(
        x=data_trm[0],
        y=data_trm[1],
        z=data_trm[2],
        spacing=inspc,
        M=maxradius,
        region=region,  # xmin, xmax, ymin, ymax
        T=0.35,  # tension factor
    )

    grid.attrs['units'] = myunits
    grid.attrs['long_name'] = z_chan
    grid.attrs['title'] = my_data.attrs['title']
    grid['x'].attrs['orig_name'] = x_chan
    grid['y'].attrs['orig_name'] = y_chan

    return grid, region


def image_pygmt(grid, region):
    """
    Shows a figure of an image of `grid` over a `region, using
    `PyGMT`'s 'grdimage()`'. The inputs are generated by `xarray_to_grid()`.

    Parameters
    ----------
    grid : 
        Contains the gridded data.
    region : tuple
        Four floating point values for xmin, xmax, ymin, ymax.

    Returns
    -------
    Nothing.

    """

    plot_width = 0.12 # ie 12 cm

    # colour bar labeling
    yc_label = f'af+l{grid.attrs["long_name"]}'
    uc_label = f'y+l{grid.attrs["units"]}'

    # axis labeling and control
    xa_label = f"xag+l{grid['x'].attrs['orig_name']}"
    ya_label = f"yag+l{grid['y'].attrs['orig_name']}"
    my_title = f'+t{grid.attrs["title"]}'

    # Colour scaling
    cmin = np.nanpercentile(grid.data, 1)
    cmax = np.nanpercentile(grid.data, 99)
    print(f'z range (1st to 99th percentile) ({cmin:.4g}, {cmax:.4g})')
    if abs(cmax - cmin) < 1.0E-10:
        cmin = -1.0E-10
        cmax = 1.0E-10

    # fit the plot to a width of 12 cm
    hscale = (grid['x'].attrs['actual_range'][1] - grid['x'].attrs['actual_range'][0]) / plot_width
    myproj = f'x1:{hscale:.0f}'

    fig = pygmt.Figure()
    fig.basemap(
        frame=[xa_label, ya_label, my_title],
        region=region,
        projection=myproj,
    )
    pygmt.makecpt(cmap="viridis", series=(cmin, cmax))
    fig.grdimage(grid=grid, nan_transparent=True, cmap=True)
    fig.colorbar(position="JMR", frame=[yc_label, uc_label])
    fig.show()


def grid_n_image(whizz_file, z_chans, mr_chans, d1_chans, grid_space):
    """
    Every channel in `z_chans` from `whizz_file` is interpolated onto a grid and imaged.
    Channels listed in `mr_chans` have the mean value of each survey line subtracted first.
    Channels listed in `d1_chans` are first differenced along each survey line first.

    Parameters
    ----------
    whizzFile : Path or String
        The Path to, or String name of, the whizz file in HDF5 format.
    z_chans : [String]
        An array of names of channels in `whizz_file` to be imaged.
    mr_chans : [String]
        An array of names of channels in `whizz_file` to be imaged.
    d1_chans : [String]
        An array of names of channels in `whizz_file` to be imaged.
    grid_space : Float
        The distance between grid cell centres in grid distance units.

    Returns
    -------
    Nothing.

    """

    for z_chan in z_chans:
        remove_mean = False
        diff_one = False
        if z_chan in mr_chans:
            remove_mean = True
        if z_chan in d1_chans:
            diff_one = True

        print(f'Gridding and imaging {z_chan}')
        my_data = whizz_to_xarray(whizz_file, z_chan, remove_mean=remove_mean, diff_one=diff_one)
        if len(my_data.attrs) == 0:
            continue
        my_grid, my_region = xarray_to_grid(my_data, grid_space)
        image_pygmt(my_grid, my_region)


def imageStats(whizzFile=''):
    """
    Returns the basic statistics of the data in a grid file. The range,
    number of points, and spacings of northing and easting and the range, number
    of samples, mean, standard deviation, and number of real and missing values
    are reported.

    Parameters
    ----------
    whizzFile : Path, optional
        The Path to the grid file in either ERS or NC format. Default is '' in
        which case a file browser allows the user to select the grid file.

    Returns
    -------
    report (String) : the statistical summary report.

    """
    xs, fileUsed = gridfile_to_xr(whizzFile, bandout=0)
    
    # Need to extract a DataArray from xs
    xa = xs
    
    report = f'Statistics for {str(fileUsed.name)}'
    report += f'\n Datum  {xa.attrs["datum"]}; Projection {xa.attrs["projection"]}'
    for coord in xa.coords:
        report += f'\n {coord} from {np.nanmin(xa[coord])} to {np.nanmax(xa[coord])}; {xa[coord].shape[0]} samples at spacing {xa[coord].values[1] - xa[coord].values[0]}'
    report += f'\n Value from {np.nanmin(xa.values):.2f} to {np.nanmax(xa.values):.2f}; {xa.values.shape[0]} x {xa.values.shape[1]} samples'
    report += f'\n     Mean = {np.nanmean(xa.values):.2f}, Stdev = {np.nanstd(xa.values):.2f}'
    report += f'\n     Number of real values = {np.count_nonzero(~np.isnan(xa.values))}, Number of nans = {np.count_nonzero(np.isnan(xa.values))}'
    return report


def checkTCratio(file000, filexxx, xxx, fileyyy, yyy, plotTitle):
    """
    Shows an image of the difference between data terrain-corrected
    at density y and the prediction based on un-corrected data and
    data terrain-corrected at x.
    Used to QC for consistency in terrain corrections.

    Parameters
    ----------
    file000 : String
        Name of the file containing the data before terrain correction.
    filexxx : String
        Name of the file containing the data after xpxx terrain correction.
    xxx : Float
        The density used in the filexxx terrain correction.
    fileyyy : String
        Name of the file containing the data after ypyy terrain correction.
    yyy : Float
        The density used in the fileyyy terrain correction.
    plotTitle : String
        A title for the plot of the difference between predicted and actual
        data in fileyyy

    Returns
    -------
    None.

    """
    n, e, g0 = read_ers_image(file000)
    n, e, gx = read_ers_image(filexxx)
    n, e, gy = read_ers_image(fileyyy)
    
    tx = g0 - gx
    predy = g0 - yyy / xxx * tx
    predErr = predy - gy
    graphicsShaded(n, e, predErr, plotTitle,  hs=False, colormap=cc.m_CET_L9, cmap_norm='no', azdeg=90)


def traceImages(file1, file2, file3, plotTitle):
    """
    Plots the trace of the three diagonal components of the tensor.

    Parameters
    ----------
    file1 : String
        A grid file containing the first diagonal tensor component.
    file2 : TYPE
        A grid file containing the second diagonal tensor component.
    file3 : TYPE
        A grid file containing the third diagonal tensor component.
    plotTitle : String
        A title for the plot of the trace image

    Returns
    -------
    None.

    """
    n1, e1, z1 = read_ers_image(file1)
    n2, e2, z2 = read_ers_image(file2)
    n3, e3, z3 = read_ers_image(file3)
    trace = z1 + z2 + z3
    graphicsShaded(n1, e1, trace, plotTitle,  hs=False, colormap=cc.m_CET_L9, cmap_norm='no', azdeg=90)


def subtractImages(imagefile1, imagefile2, scale=1.0, band1=0, band2=0):
    """
    Takes in the paths to two grid files Each may be a single band ERMapper image
    file or a netCDF4 grid file) and returns the difference as an xarray:
        (imagefile1 - scale *  imagefile2).
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
    flag (Bool): True if operation was successful.
    
    x : xarray
        The difference if both span the same coordinates, otherwise the first.

    """
    
    x1, _ = gridfile_to_xr(imagefile1, bandout=band1)
    x2, _ = gridfile_to_xr(imagefile2, bandout=band2)
    
    if x1.coords == x2.coords:
        return True, x1 - scale * x2
    else:
        return False, x1


def gridfile_to_xr(whizzFile='', bandout=0):
    """
    Returns an xarray Dataset containing the geographically located data from a
    gridfile in either ERS or NC format.

    Parameters
    ----------
    whizzFile : Path
        The Path to the grid file, must have extension `ers` or `nc
    bandout : TYPE, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    xd : xarray Dataset
        DESCRIPTION.

    """
    if whizzFile == '':
        whizzFile = Path(fb.get_grid_filename(filetypes=(('NetCdf4 grid', '*.nc'),
                                                         ('ERMapper grid', '*.ers'))))
    if whizzFile.suffix.upper() == '.ERS':
        e, n, z, datum, projection = read_ers_image(whizzFile, bandout=bandout)
        xa = xr.DataArray(data=z,#np.flip(z, 0), # DANGER!!!
                          dims=["N", "E"],
                          coords={"N": n, "E": e})
        # x2 = xa.dropna('N',how='all')
        # xa = x2.dropna('E',how='all')
        fname = whizzFile.with_suffix('').name
        xd = xr.Dataset(data_vars={fname: xa})
        xd.attrs["long_name"] = fname
        xd.attrs["datum"] = datum
        xd.attrs["projection"] = projection
        if datum == 'WGS84' and projection == 'SUTM55':
            xd.rio.write_crs("epsg:32755", inplace=True)
    elif whizzFile.suffix.upper() == '.NC':
        nc = nc4.Dataset(str(whizzFile), mode='r')
        xd = xr.open_dataset(xr.backends.NetCDF4DataStore(nc))
    else:
        print(f'ERROR: whizzFile has suffix {Path(whizzFile).suffix.upper()} but must be `.nc` or `.ers`.')
    return xd, whizzFile


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
        The units of the Dataarray values. The default is ''

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
                'units' : units}
    print(xd)
    xd.to_netcdf(ncFile)
     
    
def update_grid(whizzFile='', datum='', projection='', long_name='', units=''):
    """
    Adds attributes to a NetCdf4 grid files DataSet.

    Parameters
    ----------
    whizzFile : TYPE, optional
        DESCRIPTION. The default is ''.
    datum : TYPE, optional
        DESCRIPTION. The default is ''.
    projection : TYPE, optional
        DESCRIPTION. The default is ''.
    long_name : TYPE, optional
        DESCRIPTION. The default is ''.
    units : TYPE, optional
        DESCRIPTION. The default is ''.

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
    #             'units' : units}
    suffix = whizzFile.suffix
    fname = whizzFile.with_suffix('').name + '+'
    newWhizz = whizzFile.with_name(fname).with_suffix(suffix)
    print(xs)
    xs.to_netcdf(newWhizz)

    
def read_ers_image(whizzFile, bandout=0):
    """
    Read eastings, northings, value from an ermapper grid into numpy arrays.

    Can only handle a restricted range of ERMapper images. In particular:
                       DataSetType = ERStorage
                       DataType = Raster
                       CoordinateType = EN
                       Rotation = 0:0:0.0
    Parameters
    ----------
    whizzFile : TYPE, optional
        DESCRIPTION. The default is ''.
    datum : TYPE, optional
        DESCRIPTION. The default is ''.
    projection : TYPE, optional
        DESCRIPTION. The default is ''.
    long_name : TYPE, optional
        DESCRIPTION. The default is ''.
    units : TYPE, optional
        DESCRIPTION. The default is ''.

    Returns
    -------
    1D numpy array
        Eastings - vector containing coordinates of MIDPOINT of each
        pixel in an image row
    1D numpy array
        Northings - vector containing coordinates of MIDPOINT of each
        pixel in an image row
    2D numpy array
        The grid values at (e,n).

    """
    headerfile = str(whizzFile)
    imagefile = str(whizzFile.with_suffix(''))
    # should put a check arguments section here
    # imagefile = 'DrapeSurface_Fourier_FS'
    select_band = False            # the default, return all bands
    # want to be able to select band if bandout exists, select_band = True

    # get the header info we need to decode the image file
    [ncells, nrows, nbands, eastings, northings, nullcell, precision,
     headerbytes, originalnullcell, byteorder, datum,
     projection] = read_ers_header(headerfile)
    # [eastings, northings, zz, ncells, nrows, nbands, nullcell,
    # originalnullcell] = read_ers_header(imagefile)

    if (select_band):
        if (bandout > nbands):
            print('Asked for band # {0:2d}, but this imagefile only has \
                  {1:2d} bands!', bandout, nbands)
            print('Error in command usage for ermapper.read_ers_image')
            return  # check usage !!

        if (nbands == 1):
            select_band = False

    if nbands == 1:
        if precision == 'float64':
            print('64 bit precision')
            zz = np.fromfile(imagefile, 'float64').reshape(nrows, ncells)
        elif precision == 'float32':
            print('32 bit precision')
            zz = np.fromfile(imagefile, 'float32').reshape(nrows, ncells)
        else:
            print('unrecognised precision {}'.format(precision))
            return

        # want 'float32' a variable from byteorder
        # assume zero!! if headerbytes > 0: # skip this many bytes
        #     somehow ...

        # cif nrows * nbands is wrong, then report and break
        #   fprintf(1, 'Error in reading ERmapper imagefile: %s \n', imagefile)
        #   fprintf(1, 'Tried to read %d x %d = %d pixels \n', ncells, nrows,
        #            ncells * nrows)
        #   fprintf(1, 'Read = %d pixels \n', count)
        #   return # or break??

        # zz = flipud(zz')
    else:
        print('Error in reading ERmapper imagefile')
        print('Lazy Mark has not got around to multi-band reads.... Sorry.\n')
        return
    # now muck around with the nullcells  - convert these to NaN's regardless
    # of what they where on input. Should check for -inf as well!!

    if np.isinf(nullcell):
        zz[zz == np.inf] = np.nan
    else:
        if (nullcell != 0):
            zz[abs(zz - nullcell) < 0.01] = np.nan
        # elif (nullcell < 0):
        #   zz[abs(zz - nullcell) < 0.01] = np.nan
        else:
            zz[zz == nullcell] = np.nan

    nullcell = np.nan

    return eastings, northings, zz, datum, projection


def display_grid(whizzFile, mytitle, colormap=cc.m_CET_L9, cmap_norm='nonorm', 
                   minClip=np.nan, maxClip=np.nan, cb_ticks='stats', nSigma=2,
                   hs=True, azdeg=45, ax=None, clipTo3Std = True):
    """
    Uses `xrImage()` to display the gridded data array in whizzFile. All
    parameters after the name of the whizzFile are just passed through
    to `xrImage()`.

    Parameters
    ----------
    whizzFile : TYPE, optional
        May be either an `ERS` or `NC` grid file. The default is ''.

    Returns
    -------
    None.

    """
    (xa, _) = gridfile_to_xr(whizzFile, bandout=0)
    xrImage(xa, mytitle, colormap=colormap, cmap_norm=cmap_norm, 
                       minClip=minClip, maxClip=maxClip, cb_ticks=cb_ticks, nSigma=nSigma,
                       hs=hs, azdeg=azdeg, ax=ax, clipTo3Std = clipTo3Std)
    
    
def displayImage(e, n, z, mytitle):
    """
    Given a 2D array z located by 1D arrays e and n, plot an image of z

    Any invalid values (NaN, Inf) of z are masked out. The plot has grid lines,
    colorbar, x and y labels (Eastings and Northings respectively) and a title.

    Parameters
    ----------
    e : numpy array
        The 1D array of eastings.
    n : numpy array
        The 1D array of northings.
    z : numpy array
        The 2D array of z values.
    mytitle : String
        The image plot title.

    Returns
    -------
    None.

    """
    zm = np.ma.masked_invalid(z)
    # might change to image histogram-equalised zm or limit the range to 95%
    # z indices ought to be fixed in ermread
    plt.figure(dpi=100)
    plt.pcolormesh(e, n, zm[::-1, :], cmap='gray')
    plt.axis('equal')
    plt.xlabel('Eastings [m]')
    plt.ylabel('Northings [m]')
    plt.title(mytitle)
    plt.grid(True)
    plt.colorbar()
    plt.show()

    return


def cmap_exists(name):
    """
    Returns True if the Matplotlib Colormap `name` exists.

    Parameters
    ----------
    name : Colormap instance
        DESCRIPTION.

    Returns
    -------
    bool
        DESCRIPTION.

    """
    try:
         cm.get_cmap(name)
         return True
    except ValueError:
         pass
    return False


def imageAllInDir(path_name):
    """
    Quick and dirty method to image all the grids in a 
    directory for QC.

    Parameters
    ----------
    path_name : Path
        The Path where the ERS grid files are located.

    Returns
    -------
    None.

    """
    # get a list of the ers file paths
    file_count = 0
    ersFiles = []
    folderFiles = path_name.iterdir()
    for aFile in folderFiles:
        if aFile.is_file() and (aFile.suffix == '.ERS' or aFile.suffix == '.ers') and (aFile.name[0] != '.'):
            ersFiles.append(aFile)
            file_count += 1
    print(f'Found {file_count} .ers files ...')
    print(f'in: {str(path_name)}')

    for f in ersFiles:
        (dxt, _) = gridfile_to_xr(f)
        xrImage(dxt, str(f.name), colormap=cc.m_CET_R1)


def xrImage(data_array, mytitle, colormap=cc.m_CET_L9, cmap_norm='nonorm', 
                   minClip=np.nan, maxClip=np.nan, cb_ticks='stats', nSigma=2,
                   hs=True, azdeg=45, ax=None, clipTo3Std = True):
    """
    Uses `graphicsShaded()` to display the gridded data in data_array. All
    parameters after the name of the whizzFile are just passed through
    to `graphicsShaded()`.
    Assumes only one DataArray in the xarray dataset

    Parameters
    ----------
    data_array : 2D xarray dataset
        DESCRIPTION.

    Returns
    -------
    None.

    """
    vmin = np.nan
    vmax = np.nan
    first_data = list(data_array.keys())[0]
    fd_mean = data_array[first_data].mean()
    fd_std = data_array[first_data].std()
    if clipTo3Std:
        a = fd_mean - 3.0 * fd_std
        vmin = a.data * 1
        a = fd_mean + 3.0 * fd_std
        vmax = 1 * a.data
    elif ~np.isnan(minClip + maxClip):
        vmin = minClip
        vmax = maxClip
    
    graphicsShaded(data_array.E, data_array.N, data_array[first_data], mytitle, colormap, cmap_norm, minClip=vmin, maxClip=vmax, 
                   cb_ticks=cb_ticks, nSigma=nSigma, hs=hs, azdeg=azdeg, ax=ax)


def graphicsShaded(e, n, z, mytitle, colormap=cc.m_CET_L9, cmap_norm='nonorm', 
                   minClip=np.nan, maxClip=np.nan, cb_ticks='stats', nSigma=2,
                   hs=True, azdeg=45, ax=None):
    """
    Creates a colour image of a data array, with colour bar and grid-lines. The
    shape of z must be shape(e) X shape(n). The (e, n, z) typically are the 
    output of erm.read_ers_image().

    Parameters
    ----------
    e : np.array(Float, 1D
        The easting vector.
    n : np.array(Float, 1D)
        The northing vector.
    z : np.array(Float, 2D)
        The data to be imaged (referenced to the easting and northing positions).
    mytitle : String
        The figure title.
    colormap : Colormap, optional
        A colour map, eg cc.m_CET_L9. The default is cc.m_CET_L9.
    cmap_norm : String, optional
        Must be one of 'nonorm' (no normalisation, ie linear stretch); 'equalize'
        (equlaization stretch); 'auto'. The default is 'nonorm'.
    minClip : Float, optional
        z -> z < minClip : minClip: z. The default is np.nan - no clipping
    maxClip : Float, optional
        z -> z > maxClip : maxClip: z. The default is np.nan - no clipping
    cb_ticks : TYPE, optional
        DESCRIPTION. The default is 'stats'.
    nSigma : TYPE, optional
        DESCRIPTION. The default is 2.
    hs : Bool, optional
        hill-shading. The default is True.

    Returns
    -------
    None.

    """
    from matplotlib.ticker import StrMethodFormatter
    
    if not np.isnan(minClip) and not np.isnan(maxClip):
        z = np.clip(z, minClip, maxClip)
    elif np.isnan(minClip) and (not np.isnan(maxClip)):
        z = np.clip(z, np.min(z), maxClip)
    elif (not np.isnan(minClip)) and np.isnan(maxClip):
        z = np.clip(z, minClip, np.max(z))
    
    # register the supplied colormap for access via its name
    if not 'myCmap' in plt.colormaps():
        plt.register_cmap('myCmap', colormap)
    
    if ax == None:
        fig, ax = plt.subplots(figsize=(12,6))
    ax.set_title(mytitle, fontsize=14)
    ax.set_xlabel('Eastings [m]', fontsize=12)
    ax.set_ylabel('Northings [m]', fontsize=12)
    ax.grid(True)
    ax.axes.set_aspect('equal')
    plt.tight_layout()
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    graphics.imshow_hs(z, ax, cmap='myCmap',  cmap_norm=cmap_norm, hs=hs, colorbar=True,
                   azdeg=45, altdeg = 45, blend_mode = 'alpha', alpha = 0.7,
                   extent=(e[0], e[-1], n[0], n[-1]),origin='upper')


def graphicsTernary(e, n, red, green, blue, mytitle):
    """
    Creates a ternary image from three channels of located information.

    Parameters
    ----------
    e : TYPE
        DESCRIPTION.
    n : TYPE
        DESCRIPTION.
    red : TYPE
        DESCRIPTION.
    green : TYPE
        DESCRIPTION.
    blue : TYPE
        DESCRIPTION.
    mytitle : TYPE
        DESCRIPTION.

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
    ax.imshow(e,n,z)#, cmap = 'Reds')
    ax.set_title(mytitle, fontsize=12)
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

    
def displayShadedImage(e, n, z, mytitle):
    """
    Given a 2D array z located by 1D arrays e and n, plot an image of z.

    Any invalid values (NaN, Inf) of z are masked out. The plot has grid lines,
    colorbar, x and y labels (Eastings and Northings respectively) and a title.

    Parameters
    ----------
    e : 1D numpy array
        Eastings.
    n : 1D numpy array
        Northings.
    z : 2D numpy array
        Data to be imaged.
    mytitle : String
        Plot title.

    Returns
    -------
    None

    """
    from matplotlib.colors import LightSource

    my_ls = LightSource(azdeg=315, altdeg=50)
    cmap = plt.cm.nipy_spectral  # gist_earth

    zm = np.ma.masked_invalid(z)
    rgb = my_ls.shade(zm, cmap, fraction=0.2)
    # might change to image histogram-equalised zm or limit the range to 95%
    # z indices ought to be fixed in ermread
    fig, ax = plt.subplots()
    fig.dpi = 200

    ax.imshow(rgb, cmap=cmap)
    ax.set_title(mytitle, fontsize=8)
    ax.set_xlabel('Eastings [m]', fontsize=6)
    ax.set_ylabel('Northings [m]', fontsize=6)
    ax.grid(True)

    #   Much better is to run pcolormesh, steal its labels, remove it and set the
    #   labels
    #    ax.xaxis.set_ticklabels(('1', '692000', '694000', '696000', '698000',
    #                             '700000', '702000', '704000', '706000'),
    #                            fontsize=6)
    #    labels = ax.yaxis.get_ticklabels(which='major')
    #    ax.yaxis.set_ticklabels(('1', '3676000', '3674000', '3672000', '3670000',
    #                            '3668000', '3666000', '3664000', '3662000'),
    #                            fontsize=6)

    plt.axis('equal')
    plt.show()
    im = ax.imshow(zm, cmap=cmap)
    im.remove()
    fig.colorbar(im)

    return


def read_ers_header(imagefile):
    """
    reads ERMapper ers header file information and returns same in a dict

    NAME: READ_ERS_HEADER

    CALL AS: headerdict = read_ers_header(imagefile)

    FUNCTION: Python function to extract key parameters from an ERMapper image
    header and return these to the caller. Python attempts to parse the
    ERmapper header by extracting all the lines containing '=', and evaluating
    them. This is approach may fail with unusual header types.

    INPUT:  'imagefile' = string name of ERMAPPER image file 

    OUTPUT:  None

    RETURNS: headerdict = a dictionary with the following fields:
             ncells = scalar containing number of image pixels in easting
                           direction
             nrows = scalar containing number of image pixels in northing
                           direction
             nbands = number of bands in image
             eastings = vector containing coordinates of MIDPOINT of each
                           pixel in an image row
             northings = vector containing coordinates of MIDPOINT of each
                           pixel in an image column
             nullcell = scalar containing the value of NULL for the image
             precision = string for MATLAB's fread containing description of
                           image precision
             headerbytes = Number of bytes on header of image
             byteorder = string for fopen containing 'ieee-be' (MSBFirst) or
                           'ieee-le' (LSBFirst)

    BLAME: Mark Dransfield (based on Tim Monks Matlab code, version 2.1)
    VERSION: 1.0
    Copyright (c) 1998 by CGG Aviation (Australia) Pty Ltd.
    """

    # should check input arguments

    headerdict = Parse_header(imagefile)

    # should put checks here that headerdict is ok

    ncells = headerdict['ncells']
    nrows = headerdict['nrows']
    nbands = headerdict['nbands']
    eastings = headerdict['eastings']
    northings = headerdict['northings']
    nullcell = headerdict['nullcell']
    precision = headerdict['precision']
    headerbytes = headerdict['headerbytes']
    originalnullcell = headerdict['originalnullcell']
    byteorder = headerdict['byteorder']
    datum = headerdict['datum']
    projection = headerdict['projection']

    return [ncells, nrows, nbands, eastings, northings, nullcell,
            precision, headerbytes, originalnullcell, byteorder, datum,
            projection]


def Parse_header(imagefile, verbose=False):
    """
    parses the information in an ERMapper header file and returns the information
    as a dictionary.

    Parameters
    ----------
    imagefile : TYPE
        DESCRIPTION.
    verbose : Bool, optional
        If True, then prints out details. Default = False.

    Returns
    -------
    headerdict : TYPE
        DESCRIPTION.

    """
    # do ERMapper parsing internally with Python

    import numpy as np
    if verbose:
        print(f'Parsing {imagefile}')

    # Prepare a dict containing all the fields we need to obtain from the
    # ermapper header file with defaults set.
    ersHeader = dict(
        NrOfCellsPerLine=0,
        NrOfLines=0,
        NrOfBands=0,
        Eastings=0.0,
        RegistrationCellX=0,
        Xdimension=30.0,
        Northings=0.0,
        RegistrationCellY=0,
        Ydimension=30.0,
        CellType='float32',
        NullCellValue=-9999.0,
        HeaderOffset=0,
        ByteOrder='ieee-le',
        Datum='RAW',
        Projection='RAW')

    # open the filename
    headerfile = imagefile
    with open(headerfile, 'r') as fid:
        for myline in fid:
            if '=' in myline:
                [keyword, dummy, value] = myline.partition('=')

                # Remove leading and trailing blanks and translate
                keyword = keyword.strip()
                value = value.strip()
                value = _translate_value(value)
                rawvalue = value[0]

                if 'DATUM' in keyword.upper():
                    ersHeader['Datum'] = rawvalue
                elif 'PROJECTION' in keyword.upper():
                    ersHeader['Projection'] = rawvalue
                elif 'NROFCELLSPERLINE' in keyword.upper():
                    ersHeader['NrOfCellsPerLine'] = int(rawvalue)
                elif 'NROFLINES' in keyword.upper():
                    ersHeader['NrOfLines'] = int(rawvalue)
                elif 'NROFBANDS' in keyword.upper():
                    ersHeader['NrOfBands'] = int(rawvalue)
                elif 'XDIMENSION' in keyword.upper():
                    ersHeader['Xdimension'] = float(rawvalue)
                elif 'EASTINGS' in keyword.upper() or 'METERSX' in keyword.upper():
                    if len(rawvalue) == 0:
                        ersHeader['Eastings'] = float(0.0)
                    else:
                        ersHeader['Eastings'] = float(rawvalue)
                elif 'REGISTRATIONCELLX' in keyword.upper():
                    ersHeader['RegistrationCellX'] = float(rawvalue)
                elif 'NORTHINGS' in keyword.upper()  or 'METERSY' in keyword.upper():
                    if len(rawvalue) == 0:
                        ersHeader['Northings'] = float(0.0)
                    else:
                        ersHeader['Northings'] = float(rawvalue)
                elif 'REGISTRATIONCELLY' in keyword.upper():
                    ersHeader['RegistrationCellY'] = float(rawvalue)
                elif 'YDIMENSION' in keyword.upper():
                    ersHeader['Ydimension'] = float(rawvalue)
                elif 'CELLTYPE' in keyword.upper():
                    ersHeader['CellType'] = rawvalue
                elif 'NULLCELLVALUE' in keyword.upper():
                    ersHeader['NullCellValue'] = float(rawvalue)
                elif 'HEADEROFFSET' in keyword.upper():
                    ersHeader['HeaderOffset'] = int(rawvalue)
                elif 'BYTEORDER' in keyword.upper():
                    ersHeader['ByteOrder'] = rawvalue
    fid.close()
    # end with and while to read header file
    if verbose:
        print(ersHeader)

    ncells = ersHeader['NrOfCellsPerLine']
    nrows = ersHeader['NrOfLines']
    nbands = ersHeader['NrOfBands']
    xcell = ersHeader['Xdimension']
    ycell = ersHeader['Ydimension']
    xorigin = ersHeader['Eastings'] - ersHeader['RegistrationCellX'] \
        * ersHeader['Xdimension']
    yorigin = ersHeader['Northings'] + ersHeader['RegistrationCellY'] \
        * ersHeader['Ydimension'] - (ersHeader['NrOfLines']+1) \
        * ersHeader['Ydimension']
    precision = ersHeader['CellType']
    nullcell = ersHeader['NullCellValue']
    headerbytes = ersHeader['HeaderOffset']
    byteorder = ersHeader['ByteOrder']

    # row vector creation x = np.array([1,3,2])[None,:]
    # we use 1D array x = np.array([1,3,4])
    # want e,n of same Type as z, ie CellType
    eastings = np.arange(0, ncells)  # [None,:]
    eastings = xorigin + (eastings - 0.5) * xcell
    northings = np.arange(0, nrows)  # [None,:]
    northings = yorigin + (northings + 0.5) * ycell
    if verbose:
        print(precision, precision.__class__)
    if precision == 'float32':
        eastings = np.float32(eastings)
        northings = np.float32(northings)

    originalnullcell = nullcell
    # special test on Intrepid Nullcell's. Unless the image is stored as
    # float64, then -5e+75 will appear as -Inf
    # un-tested and not understood !!
    if ((nullcell < -4.0e+75) and (nullcell > -6.0e+75)):
        if not ('float64' in precision):
            nullcell = -np.inf

    if verbose:
        print('nEastings, nNorthings, nbands ', ncells, nrows, nbands)
        print('eastings[0], northings[0] ', eastings[0], northings[0])
        print('nullcell, precision, headerbytes ', nullcell, precision, headerbytes)
        print('originalnullcell, byteorder ', originalnullcell, byteorder)
        print('Projection, Datum', ersHeader['Projection'], ersHeader['Datum'])

    headerdict = dict(ncells=ncells,
                      nrows=nrows,
                      nbands=nbands,
                      eastings=eastings,
                      northings=northings,
                      nullcell=nullcell,
                      precision=precision,
                      headerbytes=headerbytes,
                      originalnullcell=originalnullcell,
                      byteorder=byteorder,
                      datum=ersHeader['Datum'],
                      projection=ersHeader['Projection'])

    return headerdict


def _translate_value(rawvalue):
    """
    Returns a string that can be eval'ed after protecting strings and other
    things

    Parameters
    ----------
    rawvalue : TYPE
        DESCRIPTION.

    Returns
    -------
    list
        DESCRIPTION.

    """
    """"""""

    # remove quotes and replace parentheses from value string
    rawvalue = rawvalue.replace('\'', '')
    rawvalue = rawvalue.replace('"', '')
    # don't change without changing matching bracket
    rawvalue = rawvalue.replace('{', '[')
    rawvalue = rawvalue.replace('}', ']')

    # replacements that matlab understands
    rawvalue = rawvalue.replace('Unsigned8BitInteger', 'uint8')
    rawvalue = rawvalue.replace('Signed8BitInteger', 'int8')
    rawvalue = rawvalue.replace('Unsigned16BitInteger', 'uint16')
    rawvalue = rawvalue.replace('Signed16BitInteger', 'int16')
    rawvalue = rawvalue.replace('Signed32BitInteger', 'int32')
    rawvalue = rawvalue.replace('IEEE4ByteReal', 'float32')
    rawvalue = rawvalue.replace('IEEE8ByteReal', 'float64')
    rawvalue = rawvalue.replace('MSBFirst', 'ieee-be')
    rawvalue = rawvalue.replace('LSBFirst', 'ieee-le')
    # and their upper case equivalents
    rawvalue = rawvalue.replace('UNSIGNED8BITINTEGER', 'uint8')
    rawvalue = rawvalue.replace('SIGNED8BITINTEGER', 'int8')
    rawvalue = rawvalue.replace('UNSIGNED16BITINTEGER', 'uint16')
    rawvalue = rawvalue.replace('SIGNED16BITINTEGER', 'int16')
    rawvalue = rawvalue.replace('SIGNED32BITINTEGER', 'int32')
    rawvalue = rawvalue.replace('IEEE4BYTEREAL', 'float32')
    rawvalue = rawvalue.replace('IEEE8BYTEREAL', 'float64')
    rawvalue = rawvalue.replace('MSBFIRST', 'ieee-be')
    rawvalue = rawvalue.replace('LSBFIRST', 'ieee-le')

    return [rawvalue]


def write_ers_image(imagefile, zz):
    """
    CALL AS: write_ers_image(imagefile, zz)

    FUNCTION: Python function to write an array into an ERMapper image.
              Does a simple numpy tofile(), without writing a header (*.ers)
              file. Up to the user to create the header!!!

    INPUT: imagefile = string name of ERMAPPER image file (includes path,
               but NOT THE .ers SUFFIX)
           zz = Numpy array with image values of size (nrows, ncells)
               (2D array) for single band images, and a 3D array of size
               (nrows, ncells, nbands) for multiband images. Ordering of array
               is flipped NS compared with ERMAPPER storage system, IE
               southernmost row 1st in matrix.

    RETURNS: None

    Based loosely on Tim Monks' WRITE_ERMAPPER_IMAGE (BHP, 8/12/1998).

    Mark Dransfield Oct 2017
    CGG Aviation (Australia) Pty Ltd

    """
    import numpy as np

    zz[np.isnan(zz)] = -1.0E32  # should pass nullcell!!??

    # check arguments here, check file doesn't exist
    # then modify zz as required for ERMapper

    fid = open(imagefile, 'wb')
    zz.tofile(fid)
    fid.close()

    return


def write_ers_header(imagefile, headerdict):
    """
    CALL AS: write_ers_header(imagefile, headerdict)

    FUNCTION: Python function to write an ERMapper image header in crude form

    INPUT:  imagefile = string name of ERMAPPER image file (includes path,
                but NOT THE .ers SUFFIX)
            headerdict = dictionary as returned by erm.Parse_header(imagefile)

    RETURNS: None

    Suggested usage is to call Parse_header() to return a headerdict from an
    ERMapper file with the same or similar properties to those desired. Modify
    this headerdict to suit requirements and then pass it to write_ers_header()
    to create the new header.

    Based on Tim Monks' WRITE_ERMAPPER_IMAGE (BHP, 2/12/1998).

    Mark Dransfield Oct 2017
    CGG Aviation (Australia) Pty Ltd

    """

    # would be nice to upgrade to multi-band version as Tim had.
    # check arguments
    # if (nargin < 2 || nargin > 3)
    #    fprintf(1, 'Error in call to write_ermapper_header\n')
    #    fprintf(1, 'Usage:    write_ermapper_header(imagefile, \
    #           headerdict [,bandnames])\n')
    #    error('Error in write_ermapper_header usage')
    # end

    # open the filename
    headerfile = imagefile + '.ers'
    with open(headerfile, 'wt') as fid:
        #  if fid == -1
        #      fprintf('Error opening ERmapper header file (%s)\n', headerfile)
        #      error(message)
        #  end

        # cell size
        dx = headerdict['eastings'][1] - headerdict['eastings'][0]
        dy = headerdict['northings'][1] - headerdict['northings'][0]
        estr = '{0:.2f}\n'.format(headerdict['eastings'][0])
        nstr = '{0:.2f}\n'.format(headerdict['northings'][-1])

        # now print out a standard ERmapper header
        fid.write('DatasetHeader Begin\n')
        fid.write('  Version        = "5.5"\n')
        fid.write('  DataSetType    = ERStorage\n')
        fid.write('  DataType   = Raster\n')
        fid.write('  ByteOrder  = LSBFirst\n')
        fid.write('  Comments   = "MHDs Python module"\n')
        fid.write('  CoordinateSpace Begin\n')
        fid.write('    Datum        = "' + headerdict['datum'] + '"\n')
        fid.write('    Projection   = "' + headerdict['projection'] + '"\n')
        fid.write('    CoordinateType   = EN\n')
        fid.write('    Rotation = 0:0:0.0\n')
        fid.write('  CoordinateSpace End\n')
        fid.write('  RasterInfo Begin\n')
        # always IEEE4ByteReal if using write_ers_image.py
        fid.write('    CellType = IEEE4ByteReal\n')
        fid.write('    NullCellValue = ' + str(headerdict['nullcell']) + '\n')
        fid.write('    CellInfo Begin\n')
        fid.write('      Xdimension = ' + str(dx) + '\n')
        fid.write('      Ydimension = ' + str(dy) + '\n')
        fid.write('    CellInfo End\n')
        fid.write('    NrOfLines    = ' + str(headerdict['nrows']) + '\n')
        fid.write('    NrOfCellsPerLine = ' + str(headerdict['ncells']) + '\n')
        fid.write('    RegistrationCoord Begin\n')
        fid.write('      Eastings   = ' + estr)
        fid.write('      Northings  = ' + nstr)
        fid.write('    RegistrationCoord End\n')
        fid.write('    RegistrationCellX    = 0\n')
        fid.write('    RegistrationCellY    = ' + str(headerdict['nrows']) + '\n')
        fid.write('    NrOfBands    = ' + str(headerdict['nbands']) + '\n')

        #    if (nargin == 3)
        #        if (length(bandnames) == headerdict.nbands)
        #       for i=1:headerdict.nbands
        #           fid.write('     BandId Begin\n')
        #           fid.write('         Value = "%s"\n', char(bandnames{i}))
        #           fid.write('     BandId End\n')
        #       end
        #        else
        #       disp(sprintf('BAND NAMES NOT WRITTEN OUT: %d names supplied, \
        #                  %d bands!', length(bandnames), headerdict.nbands))
        #        end
        #    end
        fid.write('  RasterInfo End\n')
        fid.write('DatasetHeader End\n')

    fid.close()

    return


