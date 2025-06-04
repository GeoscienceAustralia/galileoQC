(tutorials-target)=
# Tutorials

This section provides links to the tutorial notebooks (in Jupyter-lab format).

You can download any, or all, of the Jupyterlab notebooks from github and run them yourself once you have AirGravQC installed and running.

## Data Acknowledgements

The tutorials use field data from airborne gravity surveys to demonstrate the use of AirGravQC on real data. Permission to use field data has been kindly given by:

	Stephan Sander (Co-President, Sander Geophysics),
	Chris van Galder (Chief Geophysicist, Xcalibur Smart Mapping), and
	Colm Murphy (Chief Geoscientist, Bell Geospace)

The data sets used are described in:

The Kauring AirGrav data were collected by Sander Geophysics in 2012 and are available, with report, at [^Kauring2012].

The Canobie Falcon data were collected by Xcalibur Multiphysics in 2021. The data and report are available at [^CanobieRefs].

The Vinton Dome FTG data were supplied by Bell Geospace. The data are described in [^Murphy2004].

The Eastern Victoria airborne gravimeter survey was flown by Sander Geophysics in 2022-25. Details can be found at [^EastVicRefs].

## Prepare XYZ data

The Geosoft XYZ data format file is commonly used in geophysics to transmit geophysical survey data in a text-readable form.

We need all the data in HDF5 geoWhizz format because all the QC functions expect that format. (More on the geoWhizz format elsewhere in the `AirGravQC` documentation.)

For three of the example data sets (Eastern Victoria, Canobie, and Vinton Dome), the field data were supplied in Geosoft XYZ format. This tutorial works through the preparation of these XYZ data files as geowhizz HDF5 files for quality control checks.

The appropriate `Prepare_*` tutorial should be run before running any of the other tutorials that utilise the relevant data. 

```{toctree}
Prepare_EastVicData.ipynb
Prepare_VintonDomeData.ipynb
Prepare_CanobieData.ipynb
```

## Prepare ASEG-GDF2 data

One of the example data sets is supplied in ASEG-GDF2 format, a geophysical standard format used to transmit geophysical survey data in a text-readable form[^Pratt2003]. This tutorial demonstrates the preparation of these data as a geowhizz HDF5 files for quality control checks.

For this example we are using the data from the Kauring airborne gravity survey flown by SGL in 2012[^Kauring2012].

This tutorial should be run before running any of the other tutorials that utilise these data.

```{toctree}
Prepare_2012KauringData.ipynb
```

## Navigation and Position

A variety of requirements are set on navigation and position as demonstrated here.

```{toctree}
Navigation_and_Position.ipynb
```

## Intersection Analysis

An important analysis tool in airborne geophysics generally is intersection analysis which analyses the differences in some channel of data at the intersections of traverse and control flight-lines.

```{toctree}
Intersection_Analysis.ipynb
```

## Odd - Even Analysis

For surveys where the flight-lines are more closely spaced than the low-pass filter wavelength, odd-even grid analysis ([^Sander_2002]) is useful for assessing the noise in gravimeter data.

```{toctree}
Odd_Even_Analysis.ipynb
```

## Gravimetry Corrections

There are four standard corrections that must be made to airborne gravity data: atmospheric correction, free-air correction, Eotvos correction and latitude correction. In this notebook we check the calculation of these corrections.

```{toctree}
Gravimeter_Corrections.ipynb
```

## Falcon Noise Analysis

The primary measure of noise in a Falcon survey is difference noise[^Christensen2015], and its relationship to turbulence. High frequency noise[^Sunderland2022], and demodulation phase analysis are also useful checks for the QC of Falcon airborne survey data.

```{toctree}
Falcon_Noise_Analysis.ipynb
```

## FTG Noise Analysis

In-line noise[^Murphy2010], and its relationship to turbulence, is the standard noise measure for FTG surveys. High frequency noise[^Sunderland2022] and possibly the Frobenius norm analysis are also useful checks.

```{toctree}
FTG_Noise_Analysis.ipynb
```

## Statistical Analyses

Airborne survey data typically contain a large number of channels of data collected along a large number of flight-lines. It is useful, when reviewing any airborne geophysical survey, to be able to rapidly analyse and plot simple statistics so that any discrepancies or unusual behaviour can be rapidly found.

```{toctree}
Statistical_Analysis.ipynb
```

## Grid Imaging

The following tutorial demonstrates some of the gridding and imaging functions.

```{toctree}
Grid_Imaging.ipynb
```

[^Pratt2003]: D. A. Pratt. Geophysicists, 2003.The ASEG-GDF2 standard for point located data. In
https://www.aseg.org.au/public/200/files/ASEG-GDF2-REV4.pdf. Australian Society of Exploration

[^Murphy2004]: C. A. Murphy and G. R. Mumaw. 3D Full Tensor Gradiometry: a high resolution gravity measuring instrument resolving ambiguous geological interpretations. In ASEG-PESA 17th Geophysical Conference and Exhibition, pages 1–4, Sydney, July 2004.

[^Kauring2012]: https://researchdata.edu.au/kauring-airborne-gravity-wa-2012/3431742

[^CanobieRefs]: https://geoscience.data.qld.gov.au/data/gravity-gradiometry/gg100099, or https://researchdata.edu.au/canobie-airborne-gravity-survey-2021/1959965, or https://ecat.ga.gov.au/geonetwork/srv/api/records/744d3146-af5a-4354-8edd-9633073d51c5

[^EastVicRefs]: https://www.land.vic.gov.au/surveying/projects-and-initiatives/airborne-gravity-survey, or http://earthresources.efirst.com.au/product.asp?pID=1337&cID=37

[^Sander_2002]: S. Sander, S. Ferguson, L. Sander, V. Lavoie, and R. B. Charters. Measurement of noise in airborne gravity data using even and odd grids. First Break, 20(8):525–528, 2002.

[^Christensen2015]: A. N. Christensen, M. H. Dransfield, and C. van Galder. Noise and repeatability of airborne gravity gradiometry. First Break, 33:55–63, April 2015.

[^Murphy2010]: C. A. Murphy. Recent developments with Air-FTG®. In R. Lane, editor, Airborne Gravity 2010 - Abstracts from theASEG-PESA Airborne Gravity 2010 Workshop: Published jointly by Geoscience Australia and the Geological Survey of New South Wales, Geoscience Australia Record 2010/23 and GSNSW File GS2010/0457., pages 142–151, June 2010.

[^Sunderland2022]: A. Sunderland, Y. Naveh, L. Ju, D. G. Blair, B. Anderson, and M. Dransfield. Acoustic and vibration isolator for a gravity gradiometer. Review of Scientific Instruments, 93(6), 2022.
