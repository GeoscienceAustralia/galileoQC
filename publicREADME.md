# __pe*ga*susQC__
__pe*ga*susQC__ is a python package to check the quality of the airborne gravity dataset delivered by suppliers.

## Dependencies

__pe*ga*susQC__ uses the following packages and has been successfully tested against the version(s) listed.

	python 3.10
	colorcet 3.0
	filebrowser 1.1
	gspy 0.1
	gmt >= 6.3
	h5py 3.7
	ipympl 0.9
	matplotlib 3.5
	netCDF4 1.6
	numpy 1.22
	pathlib 1.0
	pooch 1.6
	pygmt 0.8
	rioxarray 0.11
	scikit-image 0.19
	scipy 1.8
	verde 1.7
	xarray 2022.12

	jupyterlab

Running __pe*ga*susQC__ is usually done in jupyterlab and this approach is highly recommended. All testing has been done using the app [JupyterLab Desktop](https://github.com/jupyterlab/jupyterlab-desktop) but the on-line [JupyterLab](https://jupyter.org) should work equally well.

You should ensure you have all of the above installed before installing and running __pe*ga*susQC__. Further [below](#more-installation-details) are some additional tips that might help in installing the dependencies.

## Install

These installation instructions are for the /dev version and will need changing for the release version.

First, download the /dev version of __pe*ga*susQC__ from GitHub onto a useful location on your computer. Somewhere in your Documents folder or its equivalent is typical.

Ensure you have already installed the [dependencies](#dependencies). It is usually best to have first created an environment and have all the dependencies installed there.

Then install __pe*ga*susQC__. The example has the path to __pe*ga*susQC__ on Bob's computer. Replace this with the equivalent path on your computer. If you want to be able to edit the code, use:

	python3 -m pip install -e /Users/bob/Documents/GitHub/pegasusQC

If not, you can use:

	python3 -m pip install /Users/bob/Documents/GitHub/pegasusQC

Now you can run python from the command line:
	
	python3

... and import __pe*ga*susQC__. If it imports successfully, we are done and can now go to the [JupyterLab workbook](#test).
	
	>>> import pegasusQC as qc

## More Installation Details

This sections provides some additional hints and help for the installation of the dependencies of __pe*ga*susQC__. Create and activate your virtual environment first, then install python, then GMT, and then the rest of the dependencies.

### Virtual Environment

If you don't already have a [virtual environment](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment) in which you want to do your QC work, set one up. It can be done with venv or [conda](https://docs.conda.io/en/latest/). Here is a conda example that creates a working environment called 'whizz' and installs the latest versions of python and pip into that environment. You don’t have to call it whizz.

	conda create --name whizz python pip

Make the environment the current one by activating it (again, I use conda):

	conda activate whizz

When you have finished your work, you ought to deactivate the virtual environment:

	conda deactivate whizz

JupyterLab Desktop shows you the current virtual environment in the top right of the menu bar and, once you have all the installation done and all is working, you can use that to set the environment as you wish.

### Tips

Some of the dependencies install best with pip:

	pip install mydependencypackage

and others best with conda:

	conda install yourdependencypackage

Install GMT:

	conda install gmt

GMT and pygmt may have some installation issues. Help can be had from the [pygmt](https://www.pygmt.org/latest/install.html) installation URL.

Then, since the environment is new, we need to install the following packages (some seem to prefer conda, others prefer pip; this is currently a mystery):

	conda install matplotlib h5py pathlib colorcet xarray netCDF4 rioxarray pygmt scipy pooch verde jupyterlab ipympl

	pip install numpy netCDF4 filebrowser gspy scikit-image

... and, for now, revert to version 0.1.1 of gspy.

	pip install --upgrade gspy==0.1.1

Install [JupyterLab Desktop](https://github.com/jupyterlab/jupyterlab-desktop) and now we are ready to go!

## Test

Functionality can be tested by running the Jupyter notebook `Test_pegasusQC_overall.ipynb` in the top-level of the __pe*ga*susQC__ installation folder.
