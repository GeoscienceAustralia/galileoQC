# Methods

This section provides a brief description of airborne gravity QC methods using the `AirGravQC` functions.

The examples show the necessary and common parameters passed to each function. There are other parameters available for many of the functions providing additional options. Further information can be seen in the [API](modules.rst) and examples are available in the [Cookbooks](#cookbook-target).

There are a number of types of different gravity sensor flying airborne gravity surveys on several different types of aircraft. Occasionally, gravity survey aircraft also collect magnetic data and some functions are also provided for QC of aeromagnetic survey data, although these are less thoroughly tested. The methods described here have been tested on the following aircraft and gravity sensors:

- Aircraft - BT-67, Cessna C208, Twin Otter
- Gravimeters - Sander AirGrav
- Gravity Gradiometers - Falcon AGG (both analog and digital), Bell Geo FTG

The acquired data referred to in the following Python function calls are stored in a geoWhizz file, `dh`. The planned survey line positions are stored in a different geoWhizz file, `ph`.

## geoWhizz Files

The `AirGravQC` package includes a non-proprietary binary data file format, called `geoWhizz`. The `geoWhizz` format uses HDF5 [^HDF5] which is a fast hierarchical format that is well established and supported by good open-source packages in a wide variety of languages.

There is no intent to widely disburse the `geoWhizz` format or files in this format. It is provided only to allow very fast data access of large data files during QC.

More information on the `geoWhizz` format is TODO!

A user of `AirGravQC` will typically receive located line data in one of a variety of industry format files, either proprietary or open. If proprietary, then it is up to the user to transform the data into an open format. `AirGravQC` supports `XYZ` line format and, very poorly at this time, `ASEG-GDF2` [^GDF2] format.

It is strongly recommended that the user have their data delivered in `XYZ` format at this time. They can then follow the work-flow demonstrated in the tutorial or any of the notebooks to write the data, and associated meta-data, to `geoWhizz` format.

A typical work-flow to prepare the data for QC would use the following main steps (using the Canobie project data as an example). More information is provided in the tutorial and notebooks. TODO LINK to these!!

### Convert `XYZ` to `geoWhizz`

Create a `geoWhizz` file in HDF5 format containing the data held in the XYZ file, `dx`, and give it a project name. In the following sections, the output `geoWhizz` file is `dh`.

```python
mhd.xyzToHDF(Path(dx), projectName='Canobie')
```

### Project Meta-data

These meta-data describe all the contained data. There are other options possible but these are the most useful. All of these meta-data must be delivered by the acquirer together with the data.

```python
block_name = 'Prelim Canobie Data'
mhd.updateProject(dh, acquirer='Xcalibur', blockID=block_name)
mhd.updateCoordFrame(dh, lat='LATITUDE', lon='LONGITUDE', 
	x='EASTING', y='NORTHING', time='Time_1980', alt='HEIGHT')
mhd.updateCoordFrame(dh, geoDatum='WGS84', htDatum='WGS84', 
	projection='UTM', utmz='54')
```

### Line numbering system

The line number typically encodes some or all of the following information: whether the line is a traverse or control line, the number of the corresponding planned line, the segment number for when a line is deliberately to be flown in separate flights, the re-flight number for when a part of a planned is re-flown for some reason, whether the line is a survey, test, or repeat line. Different acquisition companies use different numbering schemes and these may also vary for different projects flown by the same company.

A better system would be good.

This example has lines numbered according to a system, `Xcal_can`, used on the Canobie data.

```python
mhd.updateLineAttributes(dh, line_type='Xcal_can')
```

### Channel meta-data

As with the project meta-data, the channel meta-data must be delivered with the data. Particularly important are the units, some of which are used, and many assumed, by `AirGravQC`. The following command should be run on all data channels which is somewhat tedious (but necessary).

```python
mhd.updateChannelAttributes(dh, 'ANE_TC_2p67', units='eotvos', 
	description='')
```

### Report Contents

A summary of the contents of the created `geoWhizz` data file is useful.

```python
mhd.reportWhizz(dh)
```

### Report Flights

Summarises the flight numbers and, optionally, the lines flown on each.

```python
mhd.reportFlights(dh, detailed=True)
```

### Report Sampling

This is done in time (where the sampling rate should be constant) and space (we expect variations due to varying aircraft ground speed).

```python
mhd.reportSampling(dh)
```

### The Plan

There should be an `XYZ` file, `px`, that provides the easting, northing, and altitude (or drape) for each planned survey line. We perform a similar process to that used on the acquired data to make a `geoWhizz` plan file, `ph`. This file is used in navigation checking. TODO LINK TO NAV !!

```python
mhd.xyzToHDF(Path(px), projectName='Canobie')
block_name = 'Survey Plan'
mhd.updateProject(ph, acquirer='Xcalibur', blockID=block_name)
mhd.updateCoordFrame(ph, x='EASTING', y='NORTHING')
mhd.updateCoordFrame(ph, geoDatum='WGS84', projection='UTM', 
	utmz='54')
mhd.updateLineAttributes(ph, line_type='Xcal_can')
mhd.updateChannelAttributes(ph, 'EASTING', units='m')
mhd.updateChannelAttributes(ph, 'NORTHING', units='m')
mhd.reportWhizz(ph)
```

### Lines Map


```python
wp.linesMap([dh], whizzPlanFile=ph)
```

```{figure} linesMap.png

A map showing the flown lines (blue) against the planned lines (red). This provides a visual check that the lines are in about the right location and shows the amount of the survey flown so far. The map title, and the x and y axes are labelled using meta-data stored in the `geoWhizz` file.
```

## Navigation and Positioning

There are a number of requirements on the navigation and positioning of the aircraft that must be met.

### GNSS Satellites

There is a minimum number of GNSS satellites visible to the GNSS receiver, and there are maximum allowable values of horizontal, vertical, and position dilution of precision (HDOP, VDOP, and PDOP) that must be met for all acceptable data. The values of these variables (number of satellites, HDOP, VDOP and PDOP) must be recored as channels in the data, for each data point. The values in these channels are checked:

```python
qc.checkGNSS(dh, 'NumSats', 'PDOP', 'VDOP', 'HDOP',
	nsats_min=5, max_pdop=4, max_hdop=4, max_vdop=4)
```

### Horizontal Position

The aircraft position is allowed to deviate horizontally from the planned survey line by more than `allowance` metres for a distance no greater than `maxDistance` metres.

```python
qc.checkXYPlan(ph, dh, allowance=40.0, maxDistance=1000.0)
```

### Drape Following

The aircraft height is allowed to deviate vertically from the planned survey line by more than `allowance` metres for a distance no greater than `maxDistance` metres.

```python
qc.checkVertPlan(ph, dh, allowance=20.0, maxDistance=1000.0)
```

### Segment Length

All flown line segments must be no shorter than the required minimum. This varies with sensor but is typically 11.2 km for AGG.

```python
qc.checkLineLengths(dh, min_len=11.2)
```

### Minimum Overlap

Segments flown on the same line must overlap by some minimum amount. In this example, the minimum is 0.6 km.

```python
qc.checkOverlaps(dh, min_overlap=0.6)
```

### Minimum Clearance

There is no technical specification on minimum height but the acquisition providers should be informed if the aircraft flies very close to the ground surface on any occasion. This example reports any occasions closer than 60.0 m (more than 20 m below 80.0 m) from the ground.

```python
qc.checkClearance(dh, 80.0, altitude_chan='HEIGHT',
	terrain_chan='DTM', allowance=20.0, only_low=True)
```

### Ground Speed

If the aircraft ground speed is different from the `nominalSpeed` by more than `allowance` $\times$ `nominalSpeed` for a distance greater than `maxDistance`, then that portion of the line must be re-flown.

```python
qc.checkSpeeds(dh, nominalSpeed=60., allowance=0.15,
	maxDistance=1000.0)
```

### 

```python
	
```

### 

```python
	
```

### 

```python
	
```



## Gravimetry

Four functions provide explicit checks of the corrections applied to the gravimeter data which are all required to meet the standards in the referenced publications.

### Atmospheric Correction

The atmospheric effect is checked against [^HinzeEtAl], equation (3), using the `checkAtmosEffect()` function.

> $$ \delta g_{atm} = 0.874 - 9.9\times10^{-5} h + 3.56 \times 10^{-9} h^{2}$$

```python
qc.checkAtmosEffect(dh, 'ATMCOR', GRS80_height='MGA_Z')
```

```{figure} atmos_check.png

A statistical analysis plot of the atmospheric correction check. The mean (solid squares), range (open circles) and one standard deviation (horizontal line) summarise the statistics. In this example, all errors are less than $0.005\mu m s^{-2}$ which can be attributed to the limited precision of the data in a text-based `XYZ` format.
```

### Latitude Correction

The latitude correction (normal gravity) in the data is checked against the 1980 Geodetic Reference System ellipsoid using the Somigliana formula for theoretical gravity, $g_{T}$, on this ellipsoid [^HinzeEtAl], equation (2).

> $$ g_{T} = g_{e} \frac{\left( 1 + k \sin^{2} \phi\right)}{\left(1 - e^{2}\sin^{2} \phi\right)^{\frac{1}{2}}}$$

```python
qc.checkLatCorr(dh, 'LATCOR')
```

```{figure} latcor_check.png

A statistical analysis plot of the latitude correction check. The mean (solid squares), range (open circles) and one standard deviation (horizontal line) summarise the statistics. In this example, all errors are less than $0.005\mu m s^{-2}$ which can be attributed to the limited precision of the data in a text-based `XYZ` format.
```

### Free-air Correction

The free-air (or height) correction is checked against the second-order approximation formula for the height above the GRS80 ellipsoid ([^HinzeEtAl], equation (5)):

> $$ \delta g_{h} = -\left(0.3087691 - 0.0004389 \sin^{2} \phi \right) h + 7.2125 \times 10^{-8} h^{2} $$

```python
qc.checkFreeAirCorr(dh, 'FACOR_GRS80')
```

```{figure} freeAirCor_check.png

A statistical analysis plot of the free-air correction check. The mean (solid squares), range (open circles) and one standard deviation (horizontal line) summarise the statistics. In this example, all errors are less than $0.02\mu m s^{-2}$ which can be attributed to the limited precision of the data in a text-based `XYZ` format.
```

### Eotvos Correction

The Eotvos effect is checked against the equation in slide 30 of [^Jekeli] which is more accurate than the older approximation by Harlan.

> $$ \delta g_{Eotvos} = 2 \omega_{E} \nu_{E} \cos\phi + \frac{\nu_{N}^{2}}{M+h} + \frac{\nu_{E}^{2}}{N+h} $$

```python
qc.checkEotvosCorr(dh, 'EOTCOR', east_vel='V_EAST',
	north_vel='V_NORTH')
```

```{figure} eotcor_check.png

A statistical analysis plot of the Eotvos correction check. The mean (solid squares), range (open circles) and one standard deviation (horizontal line) summarise the statistics. In this example, all errors are less than $0.8\mu m s^{-2}$ which can be attributed to the limited precision of the speed data.
```

## AGG Gradiometry

The primary method for checking AGG data is `diffNoiseVturb()` which plots the error in each complement for each survey line against the mean turbulence experienced by the AGG along the line.	


### AGG Difference Noise

The AGG is a dual-complement, single-axis gradiometer more well known under its "Falcon" brand. The delivered gradients ($\Gamma_{NE}$ and $\Gamma_{UV}$) are the average of the gradients from each complement. Half the difference between the measured gradient components provides an error estimate of the noise in each delivered component everywhere on the survey.

The Deed requires that the standard deviation of this difference over an entire flown line segment be less than some agreed value (5.0 eotvos is typical). This value is called the *difference noise*. [^Dransfield2013]

The difference noise depends on turbulence so it is instructive to plot each complement's noise against turbulence. This also shows any outliers from the general trend.

```python
qc.diffNoiseVturb(dh, eNE='Noise_NE', eUV='Noise_UV',turbulence='TURBULENCE')
```

```{figure} AGG_noiseVturb.png
Difference Noise versus turbulence for each flown survey line. There is an outlier at a turbulence just less than 0.8 $m s^{-2}$ with noise above 3.5 E, although this is not sufficiently different from the trend to be a concern.
```

## FTG Gradiometry

### In Line Sum (ILS)

The FTG is a single-complement, three-axis gradiometer. On each axis, it measures ($\Gamma_{xy}$ and $\Gamma_{uv}^{i}$) relative to that ($i$-th) axis. Theoretically, the sum of the three $\Gamma_{uv}^{i}$ should equal zero and variations from zero reflect error in the measurement. The actual quantity, called the *in-line sum*, must be scaled for the number of axes and is expressed as:

> $$ \eta = \frac{\Gamma_{uv}^{1} + \Gamma_{uv}^{2} + \Gamma_{uv}^{3}}{\sqrt{3}} $$

The in-line sum tends generally to increase with turbulence so it is useful to plot it against turbulence. The FTG data used in testing did not include a turbulence channel but it did have a vertical velocity channel. This is supplied to `ilsNoiseVturb()` and differenced to form an acceleration channel.

The Deed specificationallows for the in-line sum to be filtered and, currently, this happens behind the scenes in the code.

```python
qc.ilsNoiseVturb(dh, diagComponent1, diagComponent2, 
	diagComponent3, vertvelocity)
```

## Aeromagnetics

:::{warning}
The aeromagnetic QC functions are based on old code, are poorly tested and are not recommended for production use. They are provided here primarily as a starting point for developers.
:::

### Fourth Difference

This check has not been used for surveys under a Deed but it is based on a common specification in the aeromagnetic survey business. The fourth difference of the total magnetic field may not exceed the agreed peak-to-peak value more than an agreed number of consecutive times along any given survey line.

Usually a fourth difference channel (here `TCDiff4`) is supplied with the data but, if not, then the raw magnetic data (`rawMag`) channel may be used.

```python
qc.checkTCDiff4(dh, TCDiff4='', rawMag='', limit = 0.02,
	nSamples = 3000, plotAll = False)
```

## General

There are a number of statistical and plotting checks that are general in nature and collected together here.

### Constant Slope

Some channels in the data are expected to vary uniformly with sampling along the survey line. Obvious examples are channels containing the date, line number, flight number, project number and so forth which should be constant; and channels containing the time or fiducial.

Errors in these channels are usually trivial and minor but nevertheless require fixing when they occur. `checkConstantSlope()` checks all the requested `fields`.

A common reported fault is when a channel containing the local time of day in seconds past midnight fails to be of constant slope at or very near the value of 86400.0. That value corresponds to midnight and, since the time of day resets to 0.0 at midnight, the slope will dramatically change. The cause is usually that the clock is not set to local time (since the survey flights do not take place at night).

```python
qc.checkConstantSlope(dh, fields=[])
```

### Gaps in Data

The Deed sets bounds on the largest permissible gap in certain data channels. With modern data acquisition systems, there ought to be no gaps at all. Consequently, the `checkGaps()` function checks all channels for all survey lines and reports any gaps in data found.

```python
qc.checkGaps(dh)
```

### Channel Statistics

It is useful to be able to quickly review the statistics of most, or all, channels in the data on a line-by-line basis in order to find any unusual behaviour in the data.

The function `allChanStats()` reports these statistics for all requested channels in the form of an expanded whisker plot.

The mean, standard deviation and range of every channel are plotted for every line. This summarises a lot of information and one can just run an eye quickly over the plots looking for outliers, and check the vertical scales to ensure that the values are in about the right range.

```python
qc.allChanStats(dh, allChannels=[])
```

```{figure} allChanStats_height.png
An example plot from `allChanStats()` for the HEIGHT channel. The mean HEIGHT for each survey line is shown by a solid blue square; the vertical blue line indicates +/- 1 standard deviation; and the open blue circles indicate the full range of the data.
```

### ERS Header Files

All delivered grids are required to be in ERS format and should be to the grid cell locations. Assuming that all the grid files are in one directory, it is simple to compare their contents against any one (here the first appearing in the directory listing) and report any substantial differences.

Some common errors in `.ers` header files are also checked for.

```python
qc.checkErsHeaders(folderPath='\.')
```

## Grid and Image

Most of the QC work is done on a line-by-line basis across the survey data. It is also useful to review images of the survey data interpolated to a regular grid. There are no specifications in the Deed but problems or artefacts in the data are often easier to see in an image.

### Grid and Image

Errors in the line number, flight number, fiducial, latitude, or any channel at all are possible and can often be easily seen in an image. Accordingly, `grid_n_image()` is provided to interpolate every named channel to a regular grid and image it to the Jupyter-lab notebook.

The work is all done by calls to `GMT` via the Python `pyGMT` package after some pre-work.

Channels (for example, line number, flight number, day of year) might vary between lines but be constant along a line. Others might vary at a constant rate along a line but change dramatically between lines flown on different days or in different directions. And some might change sign depending on line direction (bearing, velocity, Eotvos correction).

Two pre-filters are provided to deal with these situations: a mean-removal filter subtracts the mean of the channel for each line from the data on that line; and a first difference filter returns the difference between successive samples for the gridding algorithm.

```python
z_chans = ['ANE_TC_2p67', 'AUV_TC_2p67', 'BNE_TC_2p67', 'BUV_TC_2p67', 'Bearing', 'CLEARANCE', 'DTM',
           'EASTING', 'FIDUCIAL', 'FLIGHT', 'GDD_Fourier_2p67', 'GNE_Fourier_2p67', 'GUV_Fourier_2p67', 'HDOP',
           'HEIGHT', 'JOB_ID', 'LATITUDE', 'LINE', 'LONGITUDE', 'NORTHING', 'Noise_NE', 'Noise_UV', 'NumSats',
           'PDOP', 'TURBULENCE', 'T_DD', 'T_NE', 'T_UV', 'Time_1980', 'Time_Day', 'VDOP', 'gD_Fourier_2p67']
mr_chans = ['Bearing']
d1_chans = ['FIDUCIAL', 'FLIGHT', 'JOB_ID', 'LINE', 'Time_1980', 'Time_Day']
erm.grid_n_image(dh, z_chans, mr_chans, d1_chans, 500.0)
```

Only a few of the output images are shown here just to illustrate the typical output.

```{figure} img_GDD.png
All images are rendered to a linear colour scale from the 1st to 99th percentile. The white and black regions are where the data are in the extreme percentiles.
```


```{figure} img_numsats.png

caption
```


```{figure} img_flight.png
The first difference of the flight number should be, and is, exactly zero, as expected. Without the first differencing, this would be a very complex image.
```


```{figure} img_bearing.png
Removing the mean bearing for each survey line leaves just the variations in aircraft heading along each survey line. The bearing only varies by about 2.5 degrees. Noticeably, some lines have much less variation than others. A comparison with a turbulence image (not shown here) might explain this.
```

### Display Grid

When grids are provided by the acquirer, there is no need to use `grid_n_image()`. Instead, the grids can be simplay imaged to the Jupyter-lab notebook by `display_grid()`.

```python
erm.display_grid(boug_grid_path, 'Bouguer Gravity')
```

```{figure} img_Altitude.png
A shaded image of a grid. To be replaced by one of better quality!!!
```


[^Dransfield2013]: M. H. Dransfield and A. N. Christensen. Performance of airborne gravity gradiometers. The Leading Edge, 32(8):908–922, Aug. 2013

[^HinzeEtAl]: W. J. Hinze, C. Aiken, J. M. Brozena, B. Coakley, D. Dater, G. Flanagan, R. Forsberg, T. G. Hildenbrand, G. R. Keller, J. Kellogg, R. Kucks, X. Li, A. Mainville, R. Morin, M. Pilkington, D. Plouff, D. Ravat, D. Roman, J. Urrutia-Fucugauchi, M. Veronneau, M. Webring, and D. Winester. New standards for reducing gravity data: The North American gravity database. Geophysics, 70(4):J25, 2005

[^Jekeli]: C. Jekeli, Theoretical fundamentals of airborne gradiometry. In Airborne Gravity for Geodesy Summer School, 23-27 May, 2016.

[^HDF5]: See <https://www.h5py.org> or Andrew Collette, 2013, Python and HDF5, O'Reilly Media, Inc., 
ISBN: 9781449367831.

[^GDF2]: `AirGravQC` is using the code from <https://github.com/kinverarity1/aseg_gdf2>. The standard can be found at <https://www.aseg.org.au/sites/default/files/pdf/ASEG-GDF2-REV4.pdf>.

