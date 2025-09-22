(geowhizz-target)=
# geoWhizz

## Introduction

The __pe*ga*susQC__ package includes a non-proprietary binary data file format, called `geoWhizz`. The `geoWhizz` format uses HDF5 [^HDF5] which is a fast hierarchical format that is well established and supported by good open-source packages in a wide variety of languages.

All point-located data must be converted to `geoWhizz` format at the beginning of the QC work-flow because all the QC functions expect it.

There is no intent to widely distribute the `geoWhizz` format or any files in this format. It is provided only to allow very fast data access of large data files during QC. There is at least one effort underway to develop a standard for airborne geophysical survey data files [^GSpy] and it is possible that __pe*ga*susQC__ might, in the future, replace `geoWhizz` with this or some similar format.

A user of __pe*ga*susQC__ will typically receive located line data in one of a variety of industry format files, either proprietary or open. If proprietary, then it is up to the user to transform the data into an open format. __pe*ga*susQC__ supports `XYZ` line format [^XYZ] and `ASEG-GDF2` [^GDF2] format, the latter is with some use of the `aseg-gdf2` [^Kent] package.

Both are demonstrated in the tutorial notebooks.

HDF5 is structured into `groups`, `data` and `attributes`. Groups can be thought of as containers for subsidiary groups (the entire file is the root group), data and attributes. Attributes provide meta-data to groups and data.

## The root group

Starting from the top of the hierarchy, the root group (the data-file) has a group name which is its version number. This provides for the possibility of future changes in file structure. It also has the following attributes:

1. ProjectName. The Project Name and Block ID are used in titles for plots and reports.
2. BlockID. Projects might be broken into several survey blocks and this allows for them to be in separate geoWhizz files. More typically, field data is delivered in data blocks and the Block ID allows for easy identification of the delivery, for example "Delivery008".
3. Acquirer. The name of the company who are flying the survey.
4. AcquirerProjectID. The project or job number or identifier used by the Acquirer. This is usually different to the clients project ID. The geoWhizz file does not store the latter since QC is of data delivered by the Acquirer.
5. ReportName. The citation details for the acquisition and processing report. Used for final data only.

The root group has two sub-groups:
`CoordinateFrame` and `Lines`.

## The Coordinate Frame

The Coordinate Frame is a list of attributes that describe the coordinates. For many purposes, QC is of some data as a function of position, fiducial, or time and it is convenient to store the names of the channels containing these reference data. Thus LatitudeChannel is a String attribute containing the name of the channel which stores the latitudes. The following attributes store the obvious channel names:

`LatitudeChannel`, `LongitudeChannel`, `AltitudeChannel`, `XChannel`, `YChannel`, `TimeChannel`, `FidChannel`.

Additionally, the Coordinate Frame contains the following datum attributes:

`GeoDatum`, `HeightDatum`, `Projection`, `UTMZone`, `TimeDatum`.

It would be good to allow for the use of EPSG codes as well but this has not yet happened. Currently, the datum attributes are not used by any function in __pe*ga*susQC__ but they are reported so that they can be checked against contractual requirements.

## The Lines Group

The Lines Group contains an array of subsidiary line groups, one for each survey line segment.

Each line group has the following attributes (see [here](#linenumbers-target) for background information):

`LineNumber`, `HasBeenFlown`, `PlannedLine`, `Segment`, `ReflightNumber`.

The PlannedLine attribute is the most important since the aircraft positions are required to follow the planned positions to within certain specified parameters. Consequently __pe*ga*susQC__ needs to know which planned line was being followed for any flown line. The HasBeenFlown parameter is not used.

Each line group also has an array of data-set groups, one data-set group for each channel of data provided by the acquirer.

## The Dataset Group

A dataset group contains one channel of data for one flown survey line segment. It has the following attributes:

`Name`, `Units`, `Alias`, `Description`, `chan_precision`.

The Name and Units are obvious and are used often. Alias and Description are provided by analogy with other data formats (netCDF4, XArray for example). Finally, the data are typically imported from an ASCII (or possibly UTF-8) file where each number has a generally small precision (the number of digits after the decimal point in the ASCII or UTF representation). This number is recorded in chan_precision so that errors in the data caused by the limited precision can be traced by checking the attribute.

The data are stored in a 1D array, accessed as a Numpy array by __pe*ga*susQC__.


[^HDF5]: See <https://www.h5py.org> or Andrew Collette, 2013, Python and HDF5, O'Reilly Media, Inc., ISBN: 9781449367831.

[^Kent]: See <https://github.com/kinverarity1/aseg_gdf2/tree/main> or https://pypi.org/project/aseg-gdf2/.

[^GSpy]: James, S. R., Foks, N.L., and Minsely, B. J. 2022. GSPy: A new toolbox and data standard for Geophysical Datasets. Frontiers in Earth Science. 10. doi:10.3389/feart.2022.907614

[^GDF2]: Pratt, D.A., 2003. The ASEG-GDF2 standard for point located data. ASEG Standards Committee, 33. Viewed 17 November 2022, https://www.aseg.org.au/sites/default/files/pdf/ASEG-GDF2-REV4.pdf 

[^XYZ]: A reference that defines this very well known and widely used format has not been found.


