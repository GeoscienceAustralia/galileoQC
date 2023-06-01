==========================================================================================================
Project Code: GeoAus21.AUS
Project:      Airborne Gravity Survey – Victoria and South Australia
Client:       Geoscience Australia

Date:                          05 March 2023
Field Gravity Data Delivery:   DLV013
Data Description:              Gravity Repeat Line Data - Latrobe Repeat Line (Twin-Otter survey)
Number of Test Line Passes:    10 passes of Latrobe Repeat Line (7001)
==========================================================================================================

LINE DATA
==============

File Name:  DLV013_Grv.xyz
Data Rate:  2Hz
File format as follows:

  Size     Name   Units     Null    Description

  08             LINE    -        -       Line Number XXXX.YY where XXXX is line number and YY is segment or reflight number
  06           FLIGHT    -        -       Flight Number
  05             YEAR    -        -       Year
  05              DOY    -        -       Day of year
  10            FTIME    s        *       Fiducial Seconds Past Midnight UTC
  11            MGA-X    m        *       X coordinate, GDA2020 / MGA zone 55
  11            MGA-Y    m        *       Y coordinate, GDA2020 / MGA zone 55
  10            MGA-Z    m        *       GPS Elevation (above GRS80 Ellipsoid)
  10            MSL-Z    m        *       GPS Elevation (above AHD Geoid)
  13              LAT    degree   *       Latitude (GDA2020 datum)
  13             LONG    degree   *       Longitude (GDA2020 datum)
   8             PDOP    -        *       Position dilution of precision (PDOP)
   8             VDOP    -        *       Vertical dilution of precision (VDOP)
   8             HDOP    -        *       Horizontal dilution of precision (HDOP)
   8            NSATS    -        *       Number of satellites (DGPS solution)
  11              DEM    m        *       Topography from COP30 (height above GDA2020)
  11             RALT    m        *       TRT Radar altitude (Terrain Clearance)
  11             LALT    m        *       Laser altitude (Terrain Clearance)
  12               FX    mGal     *       X acceleration
  12               FY    mGal     *       Y acceleration
  12               FZ    mGal     *       Z acceleration
  12           V_EAST    m/s      *       East Velocity Component
  12          V_NORTH    m/s      *       North Velocity Component
  12           EOTCOR    mGal     *       Eotvos correction
  12           LATCOR    mGal     *       Latitude correction
  12          STATCOR    mGal     *       Static correction
  12           ATMCOR    mGal     *       Atmospheric Correction
  12      FACOR_GRS80    mGal     *       Free-Air correction (wrt GRS80)
  12      FACOR_GEOID    mGal     *       Free-Air correction (wrt GEOID)
  12            TACOR    mGal     *       Topographic Attraction Correction
  12      FA56s_GRS80    mGal     *       Preliminary Filtered Free-Air Gravity, 56s full-wavelength filter (wrt GRS80)
  12      FA56s_GEOID    mGal     *       Preliminary Filtered Free-Air Gravity, 56s full-wavelength filter (wrt GEOID)
  12     FA100s_GRS80    mGal     *       Preliminary Filtered Free-Air Gravity, 100s full-wavelength filter (wrt GRS80)
  12     FA100s_GEOID    mGal     *       Preliminary Filtered Free-Air Gravity, 100s full-wavelength filter (wrt GEOID)
  12   B56s-267_GRS80    mGal     *       Preliminary Filtered Bouguer Gravity, 56s full-wavelength filter, rock density 2.67 g/cm^3 (wrt GRS80)
  12   B56s-267_GEOID    mGal     *       Preliminary Filtered Bouguer Gravity, 56s full-wavelength filter, rock density 2.67 g/cm^3 (wrt GEOID)
  12  B100s-267_GRS80    mGal     *       Preliminary Filtered Bouguer Gravity, 100s full-wavelength filter, rock density 2.67 g/cm^3 (wrt GRS80)
  12  B100s-267_GEOID    mGal     *       Preliminary Filtered Bouguer Gravity, 100s full-wavelength filter, rock density 2.67 g/cm^3 (wrt GEOID)

==========================================================================================================

Note on Line Number convention "LLLL.AB" :

- 'LLLL' is the planned line number (4 digits: traverse lines, 3 digits: control/tie and border lines).
- 'A' represents the pre-planned segment number of the line. A planned line that is not continunous will have
      separate segments with different segment numbers.
- 'B' represents the particular occasion that a line/line segment was flown, starting from zero.
      For example, 1245.20 would be the planned segment '2' of the planned line 1245 flown for the first time.
      If this line segment were flown as two partial lines, their identifying line numbers would be 1245.20 and 1245.21.

The Latrobe Test/Repeat Line has been given the designated line number 7001.
The Lethbridge Test/Repeat Line has been given the designated line number 7002.
The Devenish Test Line has been given the designated line numner 8001.
The Otway Test Line has been given the designated line numner 8002.
For the regular survey lines, the planned segment numbers are associated with the particular region and
aircraft with which they are to be flown, specifically:

AREA:                     Planned Segment Numbers:
Caravan Area A            3
Caravan Area B            4
Caravan Area C            5 and 6
Twin Ottawa Area          0 and 1

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
