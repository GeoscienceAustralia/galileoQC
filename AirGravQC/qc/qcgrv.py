#=========================
#
# Gravimetry
#
#=========================

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import AirGravQC.config as config
import AirGravQC.whizzFiles.pointfiles as gw
import AirGravQC.gridFiles.read_ers as grd
import AirGravQC.whizzPlots.whizzPlot as wpl
import AirGravQC.utility.utility as util

groupName = config.groupName

def checkAtmosEffect(whizzFile, atmosCorr, GRS80_height=''):
    """
    Subtracts the atomspheric correction in the data file from one calculated using
    Hinze et al (2005) and the height data in the data file.
    The min, max, mean and standard deviation of the difference is calculated
    for each line and presented in an extended whisker plot.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    atmosCorr : String
        The name of the geoWhizz field or channel containing the atmospheric correction.
    GRS80_height : String, optional
        The name of the geoWhizz field or channel containing the GRS80_height. The
        default is to read the altitudeChannel field name from the Coordinate Frame.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if GRS80_height == '':
            GRS80_height = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']

        flightLine = list(g.keys())[0]
        corr_units = g[flightLine][atmosCorr].attrs['Units']
        if corr_units == 'mGal':
            unit_scale = 10.0
        elif corr_units == 'gu' or corr_units == 'µm/s/s':
            unit_scale = 1.0
        else:
            print('ERROR - correction units not recognised, expected mGal or µm/s/s (gu)')
            return

        numLines = len(g.items())
        diffMin = np.zeros((numLines,))
        diffMax = np.zeros((numLines,))
        diffMean = np.zeros((numLines,))
        diffStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            ht_data = gw.getLineData(g[line], GRS80_height)
            cor_data = gw.getLineData(g[line], atmosCorr)
           
            cal_data = _atmosEffect(ht_data)
            err_data = cor_data * unit_scale - cal_data 
            diffMin[count] = np.min(err_data)
            diffMax[count] = np.max(err_data)
            diffMean[count] = np.mean(err_data)
            diffStd[count] = np.std(err_data)
            lineNo[count] = line
            count += 1
            
        figtitle = wpl.make_plot_title(f[groupName]) + ' Atmospheric Correction Check'
        titlestr = atmosCorr + ' Stats'
        xlabelstr = 'Line number'
        ylabelstr = 'Atmospheric Correction Error [um/s/s]'
        wpl.plotBoxWhisker(diffMin, diffMax, diffMean, diffStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)

    return

    
def checkLatCorr(whizzFile, latCorr, latitude=''):
    """
    Subtracts the latitude correction in the data file from one calculated using
    Hinze et al (2005) and the latitude data in the data file. The min, max, mean
    and standard deviation of the difference is calculated for each line and
    presented in an extended whisker plot.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    latCorr : String
        The name of the geoWhizz field or channel containing the latitude correction
        (sometimes called "normal gravity").
    latitude : String, optional
        The name of the geoWhizz field or channel containing the latitude. The
        default is to read the latitude field name from the Coordinate Frame.

    Returns
    -------
    None.

    """         
    
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if latitude == '':
            latitude = f[groupName]['CoordinateFrame'].attrs['LatitudeChannel']

        flightLine = list(g.keys())[0]
        corr_units = g[flightLine][latCorr].attrs['Units']
        if corr_units == 'mGal':
            unit_scale = 10.0
        elif corr_units == 'gu' or corr_units == 'µm/s/s':
            unit_scale = 1.0
        else:
            print('ERROR - correction units not recognised, expected mGal or µm/s/s (gu)')
            return

        numLines = len(g.items())
        diffMin = np.zeros((numLines,))
        diffMax = np.zeros((numLines,))
        diffMean = np.zeros((numLines,))
        diffStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            lat_data = gw.getLineData(g[line], latitude)
            cor_data = gw.getLineData(g[line], latCorr)
            if line == '8474.0':
                fig = plt.figure()
                ax = fig.add_subplot(3,1,1)
                ax.plot(cor_data * 10)
                plt.title('lat corr in data')
                ax = fig.add_subplot(3,1,2)
                ax.plot(lat_data)
                plt.title('Latitude')
                ax = fig.add_subplot(3,1,3)
                ax.plot(_normalGravity(lat_data))
                plt.title('My Estimate')
                
            err_data = cor_data * unit_scale + _normalGravity(lat_data)
            diffMin[count] = np.min(err_data)
            diffMax[count] = np.max(err_data)
            diffMean[count] = np.mean(err_data)
            diffStd[count] = np.std(err_data)
            lineNo[count] = line
            count += 1
        
        figtitle = wpl.make_plot_title(f[groupName]) + ' Latitude Correction Check'
        titlestr = latCorr + ' Stats'
        xlabelstr = 'Line number'
        ylabelstr = 'Latitude Correction Error [um/s/s]'
        wpl.plotBoxWhisker(diffMin, diffMax, diffMean, diffStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)
    return
    

def checkEotvosCorr(whizzFile, eotCorr, latitude='', x='', y='', GRS80_height='', time='', east_vel='', north_vel=''):
    """
    Subtracts the eotvos correction in the data file from one calculated using
    Hinze et al (2005) and the latitude and position data in the data file.
    The min, max, mean and standard deviation of the difference is calculated
    for each line and presented in an extended whisker plot.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    eotCorr : String
        The name of the geoWhizz field or channel containing the eotvos correction.
    latitude : String, optional
        The name of the geoWhizz field or channel containing the latitude. The
        default is to read the latitude field name from the Coordinate Frame.
    x : String, optional
        The name of the geoWhizz field or channel containing the x position. The
        default is to read the xChannel field name from the Coordinate Frame.
    y : String, optional
        The name of the geoWhizz field or channel containing the y position. The
        default is to read the yChannel field name from the Coordinate Frame.
    GRS80_height : String, optional
        The name of the geoWhizz field or channel containing the GRS80_height. The
        default is to read the altitudeChannel field name from the Coordinate Frame.
    time : String, optional
        The name of the geoWhizz field or channel containing the time. The
        default is to read the timeChannel field name from the Coordinate Frame.
    east_vel : String, optional
        The name of the geoWhizz field or channel containing the velocity in the
        east direction. The default ('') is to calculate this from the x and time
        channels.
    north_vel : String, optional
        The name of the geoWhizz field or channel containing the velocity in the
        north direction. The default ('') is to calculate this from the y and time
        channels.

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if latitude == '':
            latitude = f[groupName]['CoordinateFrame'].attrs['LatitudeChannel']
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if y == '':
            y = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if GRS80_height == '':
            GRS80_height = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        if time == '':
            time = f[groupName]['CoordinateFrame'].attrs['TimeChannel']            

        flightLine = list(g.keys())[0]
        corr_units = g[flightLine][eotCorr].attrs['Units']
        if corr_units == 'mGal':
            unit_scale = 10.0
        elif corr_units == 'gu' or corr_units == 'µm/s/s':
            unit_scale = 1.0
        else:
            print('ERROR - correction units not recognised, expected mGal or µm/s/s (gu)')
            return

        numLines = len(g.items())
        diffMin = np.zeros((numLines,))
        diffMax = np.zeros((numLines,))
        diffMean = np.zeros((numLines,))
        diffStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            lat_data = gw.getLineData(g[line], latitude)
            x_data = gw.getLineData(g[line], x)
            y_data = gw.getLineData(g[line], y)
            ht_data = gw.getLineData(g[line], GRS80_height)
            time_data = gw.getLineData(g[line], time)
            cor_data = gw.getLineData(g[line], eotCorr)
            if (east_vel == '')  | (north_vel == ''):
                (n_speed, e_speed) = _calc_speed(x_data, y_data, time_data)
            else:
                n_speed = gw.getLineData(g[line], north_vel)
                e_speed = gw.getLineData(g[line], east_vel)
            cal_data = _eotvosCorrection(e_speed, n_speed, lat_data, ht_data)
            err_data = cor_data * unit_scale + cal_data
            diffMin[count] = np.min(err_data)
            diffMax[count] = np.max(err_data)
            diffMean[count] = np.mean(err_data)
            diffStd[count] = np.std(err_data)
            lineNo[count] = line
            
            if np.abs(diffMin[count] - diffMax[count]) > 10.0:
                fig = plt.figure()
                fig.suptitle(f'{projName} L{lineNo[count]}', fontsize=10)
                fig.subplots_adjust(top=0.85)
                ax = fig.add_subplot(2,1,1)
                ax.plot(x_data, cor_data*unit_scale, x_data, -cal_data)
                ax = fig.add_subplot(2,1,2)
                ax.plot(x_data, err_data)
                
            count += 1
        
        figtitle = wpl.make_plot_title(f[groupName]) + ' Eotvos Correction Check'
        titlestr = eotCorr + ' Stats'
        xlabelstr = 'Line number'
        ylabelstr = 'Eotvos Correction Error [um/s/s]'
        wpl.plotBoxWhisker(diffMin, diffMax, diffMean, diffStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)
    return


def _calc_speed(e, n, t):
    """
    Returns the (e_speed, n_speed) velocity given the (e, n) positions as
    a function of t (time). Uses a simple numpy.diff approach.
 
    Parameters
    ----------
    e : Numpy 1D array
        The easting data.
    n : Numpy 1D array
        The northing data.
    t : Numpy 1D array
        The time data.

    Returns
    -------
    e_speed : Numpy 1D array
        The velocity in the east direction.
    n_speed : Numpy 1D array
        The velocity in the north direction.

    """
    n_speed = np.diff(n) / np.diff(t)
    last = n_speed[-1]
    n_speed = np.append(n_speed, last)
    e_speed = np.diff(e) / np.diff(t)
    last = e_speed[-1]
    e_speed = np.append(e_speed, last)
    return (e_speed, n_speed)
    

def checkFreeAirCorr(whizzFile, faCorr, latitude='', GRS80_height=''):
    """
    Subtracts the free-air correction in the data file from one calculated using
    Hinze et al (2005) and the latitude and height data in the data file.
    The min, max, mean and standard deviation of the difference is calculated
    for each line and presented in an extended whisker plot.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    faCorr : String
        The name of the geoWhizz field or channel containing the free-air correction.
    latitude : String, optional
        The name of the geoWhizz field or channel containing the latitude. The
        default is to read the latitude field name from the Coordinate Frame.
    GRS80_height : String, optional
        The name of the geoWhizz field or channel containing the GRS80_height. The
        default is to read the altitudeChannel field name from the Coordinate Frame.

    Returns
    -------
    None.

    """
    
    filename = str(whizzFile)

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if latitude == '':
            latitude = f[groupName]['CoordinateFrame'].attrs['LatitudeChannel']
        if GRS80_height == '':
            GRS80_height = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']

        flightLine = list(g.keys())[0]
        corr_units = g[flightLine][faCorr].attrs['Units']
        if corr_units == 'mGal':
            unit_scale = 10.0
        elif corr_units == 'gu' or corr_units == 'µm/s/s':
            unit_scale = 1.0
        else:
            print('ERROR - correction units not recognised, expected mGal or µm/s/s (gu)')
            return

        numLines = len(g.items())
        diffMin = np.zeros((numLines,))
        diffMax = np.zeros((numLines,))
        diffMean = np.zeros((numLines,))
        diffStd = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        for line in g.keys():
            lat_data = gw.getLineData(g[line], latitude)
            ht_data = gw.getLineData(g[line], GRS80_height)
            cor_data = gw.getLineData(g[line], faCorr)
           
            cal_data = _freeAirCorrection(ht_data, lat_data)
            err_data = cor_data * unit_scale + cal_data 
            diffMin[count] = np.min(err_data)
            diffMax[count] = np.max(err_data)
            diffMean[count] = np.mean(err_data)
            diffStd[count] = np.std(err_data)
            lineNo[count] = line
            count += 1
            
            if line == '8334.0':
                fig = plt.figure()
                fig.suptitle(f'{projName}', fontsize=10)
                fig.subplots_adjust(top=0.85)
                ax = fig.add_subplot(3,1,1)
                ax.plot(cor_data*1e6)
                ax = fig.add_subplot(3,1,2)
                ax.plot(cal_data)
                ax = fig.add_subplot(3,1,3)
                ax.plot(err_data)
        
        figtitle = wpl.make_plot_title(f[groupName]) + ' Free Air Correction Check'
        titlestr = faCorr + ' Stats'
        xlabelstr = 'Line number'
        ylabelstr = 'Free Air Correction Error [um/s/s]'
        wpl.plotBoxWhisker(diffMin, diffMax, diffMean, diffStd, lineNo, figtitle, titlestr, xlabelstr, ylabelstr)
    return
    

def _eotvosCorrection(eSpeed, nSpeed, latitude, height=0):
    """
    Calculates the Eotvos correction for moving-base gravimetry. Uses the exact
    equation (see for example Jekeli (2016), slide 51)

    Parameters
    ----------
    eSpeed (Float) : the aircraft speed in the East direction in m/s/s.
    
    nSpeed (Float) : the aircraft speed in the North direction in m/s/s.
    
    latitude (Float) : latitude in degrees (N pos, S neg).
    
    height (Float) : the height in metres above the GRS80 ellipsoid; optional, default 0.0.

    Returns
    -------
    eotvos (Float) : the eotvos correction for gravity in um/s/s.

    """
    
    radius = 6378137.0
    ellipticity = 0.0818191908426
    ellSquared = ellipticity * ellipticity
    angularVelocity = 7.2921158553E-5
    cosLat = np.cos(latitude * np.pi / 180)
    sinLat = np.sin(latitude * np.pi / 180)
    ellSinLatSquared = ellSquared * sinLat * sinLat
    
    angTerm = -2.0 * angularVelocity * eSpeed * cosLat
    eastTerm = - (eSpeed * eSpeed) / (height + radius / (np.sqrt(1 - ellSinLatSquared)))
    northTerm = - (nSpeed * nSpeed) / (height + (radius * (1 - ellSquared)) / ((1 - ellSinLatSquared) * np.sqrt(1 - ellSinLatSquared)))
    
    eotvos = 1.0E6 * (angTerm + eastTerm + northTerm)
    return eotvos
                    

def _normalGravity(latitude):
    """
    Returns the ellipsoid theoretical gravity in µm/s/s of Hinze et al (2005), eqn 2,
    based on the GRS80 ellipsoid.    

    Parameters
    ----------
    latitude (Float) : latitude in degrees (N pos, S neg).
    

    Returns
    -------
    latitude correction (Float) : ellipsoid theoretical gravity in um/s/s.

    """
    
    gNormal = 9780326.7715
    k = 0.001931851353
    eSquared = 0.0066943800229
    sinLatSquared = np.sin(latitude * np.pi / 180.0)
    sinLatSquared = sinLatSquared * sinLatSquared
    
    
    ratio = (1 + k * sinLatSquared) / np.sqrt(1 - eSquared * sinLatSquared)
    return gNormal * (ratio)# - 1.0)


def _atmosEffect(height):
    """
    Returns the atmospheric gravity correction in um/s/s of Aiken et al (2005), eqn 3.

    Parameters
    ----------
    height (Float) : the height in metres above the GRS80 ellipsoid; optional, default 0.0.


    Returns
    -------
    atmospheric correction (Float) : in m/s/s.

    """
    return 10.0 * (0.874 - 9.9E-5 * height + 3.56E-9 * height * height)


def _freeAirCorrection(height, latitude):
    """
    Returns the free-air gravity correction in um/s/s of Hinze et al (2005), eqn 3.
    
    Parameters
    ----------
    height (Float) : the height in metres above the GRS80 ellipsoid; optional, default 0.0.

    latitude (Float) : latitude in degrees (N pos, S neg).

    Returns
    -------
    free-air gravity correction in um/s/s.

    """
    
    sinLatSquared = np.sin(latitude * np.pi / 180.0)
    sinLatSquared = sinLatSquared * sinLatSquared
    freeAir = -(0.3087691 - 0.0004398 * sinLatSquared) * height
    freeAir += 7.2125E-8 * height * height
    # Hinze et al work in mGal but we want um/s/s so x 10
    return 10 * freeAir    
    
                   
def diffGroundGrid(whizzFile, whizzChannel, whizzLine, gridPath, plot_title='Ground & Airborne Comparison'):
    """
    Samples data from the grid file onto the line from whizzFile
    and displays those data and the whizzChannel on a plot.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    whizzChannel : String
        Name of the channel of line data to be compared.
    whizzLine : String
        The identifier (eg '10030.0') for the line from which the
        whizzChannel data are taken for comparison.
    gridPath : Path
        The path (PathLib) to the ERS gridfile from which data is
        to be extracted for comparison.
    plot_title : String, optional
        A title for the plot. Default is 'Ground & Airborne Comparison'.

    Returns
    -------
    None.

    """

    # retrieve the grid information
    eg, ng, zg, datum, projection = grd.read_ers_image(gridPath)
    ngmin = np.min(ng)
    ngmax = np.max(ng)
    ngd = np.abs(ng[1] - ng[0])
    egmin = np.min(eg)
    egmax = np.max(eg)
    egd = np.abs(eg[1] - eg[0])
    zg = zg[::-1, :]
    newChannelName = gridPath.stem
    
    print('\nGrid file read for channel ', newChannelName)
    print('  SW Corner = ',ngmin, egmin, '  NE Corner = ',ngmax, egmax, '. Spacings = ', ngd, egd)

    # retrieve the line information
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        x = f[groupName]['CoordinateFrame'].attrs['XChannel']
        y = f[groupName]['CoordinateFrame'].attrs['YChannel']
        em = np.array(g[whizzLine][x])
        nm = np.array(g[whizzLine][y])
        zm = np.array(g[whizzLine][whizzChannel])

    # Some useful summary info, particularly if the line is not inside
    # the grid area.
    lineNo = whizzLine
    lineText = 'Line ' + lineNo
    print(lineText)
    print('  North min, max ', np.min(nm), np.max(nm), '  East min, max ', np.min(em), np.max(em))

    # Check that the endpoints of the line are both within the grid area.
    i_float, i_disp = divmod(np.nanmin(nm) - ngmin, ngd)
    j_float, j_disp = divmod(np.nanmin(em) - egmin, egd)
    if i_float < 0 or j_float < 0 or i_float > len(ng) or j_float > len(eg):
        print('ERROR - the flight line is not wholly within the area of the grid')
        return
    i_float, i_disp = divmod(np.nanmax(nm) - ngmin, ngd)
    j_float, j_disp = divmod(np.nanmax(em) - egmin, egd)
    if i_float < 0 or j_float < 0 or i_float > len(ng) or j_float > len(eg):
        print('ERROR - the flight line is not wholly within the area of the grid')
        return

    # A point on the line is close to the projection of a grid cell centre if
    # it is closer than eps (half the sample spacing along the line).
    sample_spacing = util._distance(em[1] - em[0], nm[1] - nm[0])
    eps = sample_spacing / 2.0

    # initialise the arrays to store the sampled data
    num_grid_samples = int(float(len(em)) * np.sqrt(2) * sample_spacing / max(ngd, egd))
    ns = np.zeros((num_grid_samples,))
    es = np.zeros((num_grid_samples,))
    zs = np.zeros((num_grid_samples,))
    diffs = np.zeros((num_grid_samples,))
    m = 0

    # traverse along the line extracting grid data as it occurs
    for k in range(0, len(em)):
        # find the closest grid cell centres (i,j), (i+1,j), (i,j+1), (i+1,j+1)
        i_float, i_disp = divmod(nm[k] - ngmin, ngd)
        j_float, j_disp = divmod(em[k] - egmin, egd)

        if i_float < 0 or j_float < 0 or i_float > len(ng) or j_float > len(eg):
            print('ERROR - the flight line is not wholly within the area of the grid')
            return

        # if the nearest grid cell centre is on the line, use the grid data
        if i_disp < eps and j_disp < eps:
            i0 = int(i_float)
            j0 = int(j_float)
            ns[m] = nm[k]
            es[m] = em[k]
            zs[m] = zg[i0,j0]
            diffs[m] = zs[m] - zm[k]
            m += 1
        # if the nearest grid cell centre is on the same northing, interpolate the grid data over easting
        elif i_disp < eps:
            i0 = int(i_float)
            j0 = int(j_float)
            j1 = j0 + 1
            ns[m] = nm[k]
            es[m] = em[k]
            zs[m] = zg[i0,j0] + (zg[i0,j1] - zg[i0,j0]) * (em[k] - eg[j0]) / (eg[j1] - eg[j0])
            diffs[m] = zs[m] - zm[k]
            m += 1
        # similarly ... for same easting, interpolate over northing
        elif j_disp < eps:
            i0 = int(i_float)
            i1 = i0 + 1
            j0 = int(j_float)
            ns[m] = nm[k]
            es[m] = em[k]
            zs[m] = zg[i0,j0] + (zg[i1,j0] - zg[i0,j0]) * (nm[k] - ng[i0]) / (ng[i1] - ng[i0])
            diffs[m] = zs[m] - zm[k]
            m += 1

    # clean up a couple of things before plotting
    if (np.ptp(em) > np.ptp(nm)):
        plotx_s = es[:m-1]
        plotx_m = em
    else:
        plotx_s = ns[:m-1]
        plotx_m = nm
    plotz_s = zs[:m-1]
    plotz_m = zm * 10.0
    # results ...
    print(f'Stdev(diff) = {np.nanstd(diffs):.1f}')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(plotx_s, plotz_s, 'b', label=newChannelName)
    ax.plot(plotx_m, plotz_m, 'g', label=lineText + ' ' + whizzChannel)
    ax.set_title(plot_title)
    ax.set_xlabel('Easting [m]')
    ax.set_ylabel('Bouguer Gravity [um/s/s]')
    ax.legend()
    ax.grid()


def checkRepeatLines(whizzFiles, channel, repeatLines, x='', z='', xOffset=True):
    """
    For all repeatLines, plot (x, channel) and report stats of differences to mean.
    This will require trimming to [minX, maxX] and interpolating to common x.
    Repeat the analysis for the `z` channel (height).

    Parameters
    ----------
    whizzFiles : array of HDF5 Whizz file pathlib Paths
        The pathlib Paths to the Whizz HDF5 files containing the survey repeat line data.
    channel : String
        The name of the channel or field to analyse and plot. Usually a gravity channel.
    repeatLines : [String], optional
        A list of flightlines, e.g. ['1000110.0', '1000210.0', '1000310.0']. 
    x : String, optional
        The name of the independent variable for the plot. Defaults to the `XChannel`.
    z : String, optional
        The name of the height variable for the analysis and plot. Defaults to the `XChannel`.
    xOffset : Bool, optional
        If True, map x to x - x[0] before plotting. The default is True.

    Returns
    -------
    None

    """

    # build the arrays to store the data
    temp_repeats = repeatLines.copy()
    xBase, xData, yData, zData, minBigX, maxSmallX, deltaX  = _xBaseInterpolant(whizzFiles, channel, temp_repeats, x, z)
    temp_repeats = repeatLines.copy()

    # Interpolate the data to common x and store in arrays
    lineCount = 0
    for whizzFile in whizzFiles:

        filename = str(whizzFile)
        with h5py.File(filename, 'r') as f:
            g = f[groupName]['Lines']
            if x == '':
                x = f[groupName]['CoordinateFrame'].attrs['XChannel']
                north = f[groupName]['CoordinateFrame'].attrs['YChannel']
            if z == '':
                z = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
            all_flightLines = list(g.keys())

            # if the channel has an attribute 'Units'
            dd = g[all_flightLines[0]][channel]
            chan_y_label = channel
            if 'Units' in dd.attrs.keys():
                chan_y_units = dd.attrs['Units']
                chan_y_label += ' ' + chan_y_units
            ddz = g[all_flightLines[0]][z]
            if 'Units' in ddz.attrs.keys():
                chan_z_units = ddz.attrs['Units']

            # read the data into the arrays
            for line in all_flightLines:
                if line in temp_repeats:
                    baseLine = g[line].attrs['PlannedLine']
                    xd = gw.getLineData(g[line], x)
                    yd = gw.getLineData(g[line], channel)
                    zd = gw.getLineData(g[line], z)

                    # Get the heading TODO: use this to check RMS(mean difference vs heading direction)
                    dx = np.diff(xd)
                    dy = np.diff(gw.getLineData(g[line], north))
                    heading = np.arctan2(dx, dy) * 180.0 / np.pi
                    mean_heading = np.mean(heading)
                    print(f'Line {line} heading = {mean_heading:.1f} deg.')

                    # ensure ordered in increasing x
                    if xd[1] < xd[0]:
                        xd = xd[::-1]
                        yd = yd[::-1]
                        zd = zd[::-1]

                    xStart = 0
                    xEnd = xd.size - 1
                    
                    # trim data and store
                    for xSample in range(0, xd.size):
                        if xd[xSample] < (maxSmallX - deltaX / 2.0):
                            xStart = max(xSample, xStart)
                        else:
                            break
                    for xSample in range(xd.size-1, 0, -1):
                        if xd[xSample] > (minBigX + deltaX / 2.0):
                            xEnd = min(xSample, xEnd)
                        else:
                            break
                            
                    # interpolate data
                    (yOut, _) = gw.interpolateLine(xd-xBase[0], yd, xBase-xBase[0])
                    (zOut, _) = gw.interpolateLine(xd-xBase[0], zd, xBase-xBase[0])

                    vec_len = len(xBase)-1 # interpolateLine has lost a datapoint in outputs
                    # print(f'line {line}, shapes: xBase {xBase.shape}, xData {xData.shape}')
                    xData[lineCount, 0:vec_len] = xBase[1:]
                    yData[lineCount, 0:vec_len] = yOut
                    zData[lineCount, 0:vec_len] = zOut
                    lineCount += 1
                    # In case the line is in more than one geoWhizz file
                    temp_repeats.remove(line)
        
    # analyse statistics and report with plots
    wpl._plotRepeatAnalysis(xBase, xOffset, lineCount, xData, yData, zData, channel, repeatLines, baseLine, z, chan_z_units, chan_y_label, chan_y_units)
            
    return


def _xBaseInterpolant(whizzFiles, channel, repeatLines, x='', z=''):

    nSamples = 0
    minBigX = 1.0E12
    maxSmallX = -1.0E12
    nLines = len(repeatLines)
    linecount = 0
    
    for whizzFile in whizzFiles:
        filename = str(whizzFile)
        with h5py.File(filename, 'r') as f:
            g = f[groupName]['Lines']
            if x == '':
                x = f[groupName]['CoordinateFrame'].attrs['XChannel']
                north = f[groupName]['CoordinateFrame'].attrs['YChannel']
            if z == '':
                z = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
            all_flightLines = list(g.keys())

            # nSamples is the array width for data storage
            for line in all_flightLines:
                if line in repeatLines:
                    linecount += 1
                    xs = gw.getLineData(g[line], x)
                    nSamples = max(nSamples, xs.size)
                    minBigX = min(max(xs), minBigX)
                    maxSmallX = max(min(xs), maxSmallX)
                    deltaX = np.abs(xs[1] - xs[0])
                    repeatLines.remove(line)
                
    if minBigX < maxSmallX:
        return 0.0
    xBase = np.linspace(maxSmallX, minBigX, num=nSamples, endpoint=True)
    print(f'{linecount} of {nLines} lines analysed, each with {nSamples} samples.')
    xData = np.empty((nLines, nSamples))
    xData[:] = np.nan
    yData = np.empty((nLines, nSamples))
    yData[:] = np.nan
    zData = np.empty((nLines, nSamples))
    zData[:] = np.nan

    return xBase, xData, yData, zData, minBigX, maxSmallX, deltaX

