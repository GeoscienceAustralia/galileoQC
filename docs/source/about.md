# About

__pe*ga*susQC__ is a collection of functions for use in checking airborne gravity data under a quality control procedure. A number of assumptions underlie its development and use.

The airborne data in question might be from a gravimeter, or gravity gradiometer, or both, carried on an aircraft together with other instruments as part of a survey to map the changes in gravity over some area. This is done by flying a series of nominally parallel and equally spaced survey lines covering the area. The flying height is either constant with respect to some height datum, or else it approximately and reasonably smoothly drapes over the terrain.

__pe*ga*susQC__ provides functions to:

- check the navigation and positioning of the aircraft on the survey lines;
- calculate and display the statistics of the collected data;
- verify the accuracy of the calculated corrections to the gravimetry data;
- report on the noise characteristics of gradiometry data;
- interpolate data to regular grids, and image and display those grids;
- check consistency of supplied data in various forms.

These various checks are designed largely in line with the Airborne Gravity Deed from Geoscience Australia. The Deed provides a standard technical specification for airborne gravity surveys and is intended for use in contracts for airborne gravity surveys. pe*ga*susQC provides a number of additional checks beyond those explicitly laid out in the Deed; these check the data against best practice and/or provide additional information for the person performing the QC.

The tested, and preferred, method of using __pe*ga*susQC__ is via a Jupyter-lab notebook. Some examples are available in [Tutorials](tutorials-target). The notebooks can be used as a template for new QC projects if that is useful. They provide a method of keeping the reporting of the QC together with its analysis.


*Mark Dransfield*
> {sub-ref}`today`

## Acknowledgements

Lots of help from various people at GA.
The minimum curvature code was gratefully copied from PyGMI (https://patrick-cole.github.io/pygmi/index.html) in July 2025.

## Citing

One day there might be a publication to cite.

## Contributing

Once this is properly set up, people are encouraged to contribute.

## Licence

Whatever GA want.

## The Airborne Gravity Deed

Once published, there should be a reference here.
