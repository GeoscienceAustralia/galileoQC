# About

__galileoQC__ is a collection of functions for use in checking airborne gravity data under a quality control procedure. A number of assumptions underlie its development and use.

The airborne data in question might be from a gravimeter, or gravity gradiometer, or both, carried on an aircraft together with other instruments as part of a survey to map the changes in gravity over some area. This is done by flying a series of nominally parallel and equally spaced survey lines covering the area. The flying height is either constant with respect to some height datum, or else it approximately and reasonably smoothly drapes over the terrain.

__galileoQC__ provides functions to:

- check the navigation and positioning of the aircraft on the survey lines;
- calculate and display the statistics of the collected data;
- verify the accuracy of the calculated corrections to the gravimetry data;
- report on the noise characteristics of gradiometry data;
- report on the quality of aeromagnetic data;
- interpolate data to regular grids, and image and display those grids;
- check consistency of supplied data in various forms.

These various checks are designed largely in line with the Airborne Gravity Deed from Geoscience Australia. The Deed provides a standard technical specification for airborne gravity surveys and is intended for use in contracts for airborne gravity surveys. galileoQC provides a number of additional checks beyond those explicitly laid out in the Deed; these check the data against best practice and/or provide additional information for the person performing the QC.

The tested, and preferred, method of using __galileoQC__ is via a JupyterLab notebook. Some examples are available as JupyterLab notebooks in [Tutorials](tutorials-target). The notebooks can be used as templates for new QC projects if that is useful. They provide a method of keeping the reporting of the QC together with its analysis.

*Mark Dransfield* {sub-ref}`today`

:::{note}
__galileoQC__ is at early stage of development and will keep improving in the future. The commonly used functions in the API should be quite stable, but minor utilities could change in the next version. If you find any bugs or would like to request any enhancements, please raise an issue on GitHub.
:::

## Acknowledgements

A massive amount of help, encouragement and support from many people at Geoscience Australia, particularly Roger Miller, Yvette Poudjom Djomani, Negin Moghaddam, Mike Barlow, Anandaroop Ray, and Jack McCubbine. The minimum curvature code was copied from PyGMI (https://patrick-cole.github.io/pygmi/index.html) in July 2025 and the graphical imaging code from graphics (https://github.com/jobar8/graphics) sometime in 2020. Grateful thanks to the authors of those codes as well as all the contributors to all the packages used. The idea of using tutorials as the key technical documentation is partly based on Des Fitzgerald's "cookbooks" for Intrepid.

## Testing

The __galileoQC__ package has been tested on the following operating systems: Windows 10, Windows 11, macos 15, Debian 12, Fedora 38, and Ubuntu 22. All testing has been done using Jupyter Notebook (the JupyterLab Desktop app was used often in development but its maintenance ceased in Aug2025).

Data from a wide variety of gravity instruments, mounted on a wide variety of aircraft, have been reviewed from many surveys flown in several countries.

## Citing

M. H. Dransfield and R. Miller. Python code for quality control of airborne gravimetry and gravity
gradiometry field data. Preview, 237:78, August 2025.

R. Miller and M. H. Dransfield. Metrics for quality control of airborne gravity gradiometry field data. Preview, 237:91, August 2025.

## Contributing

Once this is properly set up, people are encouraged to contribute.

## Licence

Creative Commons Attribution-ShareAlike 4.0 International

This license requires that reusers give credit to the creator. It allows reusers to distribute, remix, adapt, and build upon the material in any medium or format, even for commercial purposes. If others remix, adapt, or build upon the material, they must license the modified material under identical terms. The full text is in [License](license-target).

```{toctree}
:maxdepth: 2
:hidden:
license.md
```

## The Airborne Gravity Deed

Once published, there should be a reference here.
