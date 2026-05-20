(tutorials-target)=
# Tutorials

This section provides links to the tutorial notebooks (in Jupyter-lab format).

The examples show the necessary and common parameters passed to each function. There are other parameters available for many of the functions providing additional options. Further information can be seen in the [API](modules-target).

You can download any, or all, of the Jupyterlab notebooks from github and run them yourself once you have __galileoQC__ installed and running. Once you have copied a notebook, you can modify it and use it on your own data.

The notebooks are at <https://github.com/GeoscienceAustralia/galileoQC/tree/main/docs/source/tutorials>.

## Prepare XYZ data

The Geosoft XYZ data format file is commonly used in geophysics to transmit geophysical survey data in a text-readable form.

We need all the data in HDF5 geoWhizz format because all the QC functions expect that format. (More on the geoWhizz format [here](#geowhizz-target).)

For three of the example data sets (Eastern Victoria gravimetry, Melbourne magnetometry, and Vinton Dome), the field data are supplied in Geosoft XYZ format. This tutorial demonstrates the preparation of these XYZ data files as geowhizz HDF5 files for later use in quality control checks.

The appropriate `Prepare_*` tutorial should be run before running any of the other tutorials that utilise the relevant data.

```{toctree}
Prepare_EastVicData.ipynb
Prepare_VintonDomeData.ipynb
Prepare_AeromagData.ipynb
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

## Repeat Lines

Another important analysis tool in airborne gravimetry in particular is repeat line analysis which analyses the repeatability in some channel of data over a series of flight-lines flown over the same planned line.

```{toctree}
CheckRepeats_Analysis.ipynb
```

## Odd - Even Grids

For surveys where the flight-lines are more closely spaced than the low-pass filter wavelength, odd-even grid analysis ([^Sander_2002]) is useful for assessing the noise in gravimeter data.

```{toctree}
Odd_Even_Analysis.ipynb
```

## Gravimetry Corrections

There are four standard corrections that must be made to airborne gravity data: atmospheric correction, free-air correction, Eotvos correction and latitude correction.  These corrections must be made using the modern standard formulae appropriate to airborne gravimetry. The formulae for the atmospheric, free-air and latitude corrections used here were taken from [^Hinze2005] where they appear in units of $mGal$. Here they are modified to be in $\mu ms^{-2}$ but __galileoQC__ uses the units of the measured data, either $mGal$ or $\mu ms^{-2}$. The Eotvos correction formula was taken from [^Jekeli2016] and also appears in [^Zhao2015].


```{toctree}
Gravimeter_Corrections.ipynb
```

## Falcon Noise Analysis

The primary measure of noise in a Falcon survey is difference noise[^Dransfield2013], and its relationship to turbulence. High frequency noise[^Sunderland2022] is also a useful check for the QC of Falcon airborne survey data.

```{toctree}
Falcon_Noise_Analysis.ipynb
```

## FTG Noise Analysis

The FTG is a single-complement, three-axis gradiometer. On each axis, it measures ($\Gamma_{xy}^{i}$ and $\Gamma_{uv}^{i}$) relative to that ($i$-th) axis. Theoretically, the sum of the three $\Gamma_{uv}^{i}$ should equal zero and variations from zero reflect error in the measurement. The actual quantity, called the *in-line sum*[^Murphy2010], must be scaled for the number of axes and is expressed as:

> $$ \eta = \frac{\Gamma_{uv}^{1} + \Gamma_{uv}^{2} + \Gamma_{uv}^{3}}{\sqrt{3}} $$

The in-line sum tends generally to increase with turbulence so it is useful to plot it against turbulence. 

In-line noise, and its relationship to turbulence, is the standard noise measure for FTG surveys.

There are mechanisms [^Sunderland2022] by which gravity gradiometer noise at a frequency higher than the data frequency band can be down-converted to the data frequency band, resulting in errors in the data. This can occur in both AGG and FTG systems.

The checkHighFreq function checks for periods of high amplitude, high frequency signal in the raw gradiometer channels. It highlights sections of a survey line with excess high-frequency signal which might result in high gradient error. It is important to know that these erroneous gradients are true acceleration gradients and can be difficult to see by other methods.

```{toctree}
FTG_Noise_Analysis.ipynb
```

## Ground Gravity Comparisons

Airborne gravity test lines are flown to allow comparison with ground gravity. When making such comparisons, it is useful to know the sampling and quality of the ground gravity and, in Australia, this can be assessed by `plotLinesOnGroundStns`. This code relies heavily on the data collection[^Wynne] in its compiled form[^Uieda].

```{toctree}
Compare_Gravity.ipynb
```

## Statistical Analyses

Airborne survey data typically contain a large number of channels of data collected along a large number of flight-lines. It is useful, when reviewing any airborne geophysical survey, to be able to rapidly analyse and plot simple statistics so that any discrepancies or unusual behaviour can be rapidly found.

```{toctree}
Statistical_Analysis.ipynb
```

## Grid Imaging

Most of the QC work is done on a line-by-line basis across the survey data. It is also useful to review images of the survey data interpolated to a regular grid because problems or artefacts in the data can be easier to see in an image. The following tutorial demonstrates some of the gridding and imaging functions.

```{toctree}
Grid_Imaging.ipynb
```

## Craig transform

The Craig transform method transforms differential curvature components $G_{UV}$ and $G_{NE}$ into vertical gravity[^Dransfield2019]. This is useful because gravity interpreters are used to working with the gravity data, and because there are often regional or other gravity data that cover part or all of the survey area with which the transformed gravity data can be compared.

If regional gravity data are available, it is also possible to conform [^Dransfield2010] the survey gravity, derived from the Craig transform , to the regional gravity.

The following tutorial demonstrates the use of the Craig transform and the conforming process.

```{toctree}
Craig_transform.ipynb
```

## Aeromagnetics

The following tutorial demonstrates some of the aeromag QC functions.

```{toctree}
Aeromag_Analysis.ipynb
```

## Data Acknowledgements

The tutorials use field data from airborne gravity surveys to demonstrate the use of __galileoQC__ on real data. Permission to use field data has been kindly given by:

Stephan Sander (Co-President, Sander Geophysics),

Chris van Galder (Chief Geophysicist, Xcalibur Smart Mapping), and

Colm Murphy (Chief Geoscientist, Bell Geospace)


The Kauring AirGrav data were collected by Sander Geophysics in 2012 and are available, with report, at [^Kauring2012].

The Canobie Falcon field data were collected by Xcalibur Multiphysics in 2021. The final data and report are available at [^CanobieRefs].

The Vinton Dome FTG final data were supplied by Bell Geospace. The data are described in [^Murphy2004].

The Eastern Victoria airborne gravimeter survey field data were supplied by Sander Geophysics in 2022-25. Details can be found at [^EastVicRefs].

## References

[^Pratt2003]: D. A. Pratt. Geophysicists, 2003.The ASEG-GDF2 standard for point located data. In
https://www.aseg.org.au/public/200/files/ASEG-GDF2-REV4.pdf. Australian Society of Exploration

[^Murphy2004]: C. A. Murphy and G. R. Mumaw. 3D Full Tensor Gradiometry: a high resolution gravity measuring instrument resolving ambiguous geological interpretations. In ASEG-PESA 17th Geophysical Conference and Exhibition, pages 1–4, Sydney, July 2004.

[^Kauring2012]: https://researchdata.edu.au/kauring-airborne-gravity-wa-2012/3431742

[^CanobieRefs]: https://geoscience.data.qld.gov.au/data/gravity-gradiometry/gg100099, or https://researchdata.edu.au/canobie-airborne-gravity-survey-2021/1959965, or https://ecat.ga.gov.au/geonetwork/srv/api/records/744d3146-af5a-4354-8edd-9633073d51c5

[^EastVicRefs]: https://www.land.vic.gov.au/surveying/projects-and-initiatives/airborne-gravity-survey, or http://earthresources.efirst.com.au/product.asp?pID=1337&cID=37

[^Sander_2002]: S. Sander, S. Ferguson, L. Sander, V. Lavoie, and R. B. Charters. Measurement of noise in airborne gravity data using even and odd grids. First Break, 20(8):525–528, 2002.

[^Murphy2010]: C. A. Murphy. Recent developments with Air-FTG®. In R. Lane, editor, Airborne Gravity 2010 - Abstracts from theASEG-PESA Airborne Gravity 2010 Workshop: Published jointly by Geoscience Australia and the Geological Survey of New South Wales, Geoscience Australia Record 2010/23 and GSNSW File GS2010/0457., pages 142–151, June 2010.

[^Sunderland2022]: A. Sunderland, Y. Naveh, L. Ju, D. G. Blair, B. Anderson, and M. Dransfield. Acoustic and vibration isolator for a gravity gradiometer. Review of Scientific Instruments, 93(6), 2022.

[^Hinze2005]: W. J. Hinze, C. Aiken, J. M. Brozena, B. Coakley, D. Dater, G. Flanagan, R. Forsberg, T. G. Hildenbrand, G. R. Keller, J. Kellogg, R. Kucks, X. Li, A. Mainville, R. Morin, M. Pilkington, D. Plouff, D. Ravat, D. Roman, J. Urrutia-Fucugauchi, M. Véronneau}, M. Webring, and D. Winester. New standards for reducing gravity data: The North American gravity database. Geophysics, 70(4):J25, 2005.

[^Jekeli2016]: C. Jekeli. Theoretical Fundamentals of Airborne Gravimetry. In Airborne Gravity for Geodesy Summer School, 23-27 May, May 2016 (accessed 23 May 2025 at https://www.ngs.noaa.gov/grav-d/2016SummerSchool/history-fundamentals.shtml).

[^Zhao2015]: L. Zhao, R. Forsberg, M. Wu, A. V. Olesen, K. Zhang, and J. Cao. A flight test of the strapdown airborne gravimeter sga-wz in greenland. Sensors, 15:13258–13269, 2015.

[^Dransfield2013]: M. H. Dransfield and A. N. Christensen. Performance of airborne gravity gradiometers. The Leading Edge, 32(8):908–922, Aug. 2013.

[^Dransfield2019]: M. H. Dransfield and T. Chen. Heli-borne gravity gradiometry in rugged terrain. Geophysical Prospecting, 67(6):1626–1636, July 2019.

[^Dransfield2010]: M. H. Dransfield. Conforming Falcon gravity and the global gravity anomaly. Geophysical Prospecting, 58(3):469--483, May 2010.

[^Uieda]: Uieda, L. (2021). Ground gravity data compilation for Australia version 2.0. figshare. https://doi.org/10.6084/m9.figshare.13643837

[^Wynne]: Wynne, P. (2018). NetCDF Ground Gravity Point Surveys Collection (Version 1.0). Commonwealth of Australia (Geoscience Australia). https://doi.org/10.26186/5C1987FA17078

