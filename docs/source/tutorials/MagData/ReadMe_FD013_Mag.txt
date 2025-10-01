==========================================================================================================
Project Code: GeoAus21.AUS
Project:      Airborne Gravity Survey - Victoria and South Australia
Client:       Geoscience Australia

Date:                          2024-10-17
Field Gravity Data Delivery:   FD013
Data Description:              Magnetic Flight Line Data
Flight Range:                  2102-2217
==========================================================================================================

LINE DATA
==============

File Name:  FD013_Mag.xyz
Data Rate:  10Hz
File format as follows:

   Name  Size  Units   Null  Description

   LINE     8    -     *.**  Line Number LLLL.AB where L is line number, A is segment number B is reflight number
 FLIGHT     7    -      *    Flight Number
   YEAR     6    -      *    Year
    DOY     5    -      *    Day of year
  FTIME    10    s      *    Fiducial Seconds Past Midnight UTC
  MGA-X    12    m      *    X coordinate, GDA2020 / MGA zone 55
  MGA-Y    12    m      *    Y coordinate, GDA2020 / MGA zone 55
  MGA-Z    11    m      *    GPS Elevation (above GRS80 Ellipsoid)
  MSL-Z    11    m      *    GPS Elevation (above AHD Geoid)
    LAT    14  degree   *    Latitude (GDA2020 datum)
   LONG    14  degree   *    Longitude (GDA2020 datum)
    DEM    12    m      *    Topography from COP30 (height above GDA2020)
DIURNAL    12    nT     *    Diurnal / ground magnetics base station
COMPMAG    12    nT     *    Raw compensated, lag corrected magnetic total field
  DCMAG    12    nT     *    Magnetic total field diurnal corrected
IGRFMAG    12    nT     *    Magnetic total field diurnal and IGRF (International Geomagnetic Reference Field) corrected
LVLDMAG    12    nT     *    Magnetic total field, levelled to survey

Line Number convention:
LLLL.AB

'LLLL' is the planned line number.
Line numbers with four digits are traverse lines.
Line numbers with three digits are tie/control lines.
'A' represents the preplanned segment number of the line.
A planned line that is not continunous will have separate segemnts with different segment numbers.
'B' represents the particular occasion that a line/line segment was flown, starting from zero.
For example 1245.20 would be the planned segment '2' of the planned line 1245 flown for the first time.
If this line segment were flown as two partial lines, their identifying line numbers would be 1245.20 and 1245.21.

For the regular survey lines, the planned segment numbers are associated with the particular region and
aircraft with which they are to be flown, specifically:


AREA:                     Planned Segment Numbers:
Caravan Area A            3
Caravan Area B            4
Caravan Area C            5 and 6
Twin Ottawa Area          0 and 1


GRIDS
==============
Datum:           WGS84
Projection:      UTM55S
Grid cell size:  200m
Formats:         ERDAS ERMapper

Name               Units         Description
--------------     ------        ----------------------------------------------
MAG                nT            Preliminary magnetic total field, levelled to survey
FVM                nT            First vertical derivative of preliminary magnetic total field, levelled to survey

==========================================================================================================


Flown and compiled by:

Sander Geophysics Ltd.
260 Hunt Club Road
Ottawa ON  K1V 1C1
CANADA
Phone : +1 613 521 9626
Fax   : +1 613 521 0215
Email : info@sgl.com
Web   : www.sgl.com
