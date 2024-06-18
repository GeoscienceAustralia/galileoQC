==========================================================================================================
Project Code: GeoAus21.AUS
Project:      Airborne Gravity Survey - Victoria and South Australia
Client:       Geoscience Australia

Date:                          18 April 2022
Field Gravity Data Delivery:   DLV003
Data Description:              Gravity Flight Line Data
Flight Range:                  1001-1014
==========================================================================================================

LINE DATA
==============

File Name:  DLV003_Grv.xyz
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


GRIDS
==============
Datum:           WGS84
Projection:      UTM55S
Grid cell size:  200m
Formats:         ERDAS ERMapper

Name               Units         Description
--------------     ------        ----------------------------------------------
F056a5000           mGal         Preliminary Free-Air Gravity Anomaly, 5000m full-wavelength spatial filter
B056a5000           mGal         Preliminary Bouguer Gravity, 5000m full-wavelength spatial filter, rock density 2.67 g/cm^3
F056a5000-fvd     Eotvos         First Vertical Derivative of Preliminary Free-Air Gravity Anomaly, 5000m full-wavelength spatial filter
B056a5000-fvd     Eotvos         First Vertical Derivative of Preliminary Bouguer Gravity, 5000m full-wavelength spatial filter, rock density 2.67 g/cm^3

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
