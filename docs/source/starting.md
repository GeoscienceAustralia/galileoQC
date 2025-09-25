# Getting Started

The __pe*ga*susQC__ package is run from the Jupyter-lab notebook. You should create a working environment (called `whizz` in this documentation) with python3, pip, conda or miniconda

## Installation

To use __pe*ga*susQC__, first install it into your `whizz` environment using pip (conda should work as well but conda installations of __pe*ga*susQC__ are untested):

```bash
(whizz) $ pip install pegasusQC
```

```bash
(whizz) $ conda install pegasusQC
```

Then start up a Jupyter-lab notebook. For more information about this, see [^Jupyter] . The Jupyter-lab app [^JupyterApp] is recommended.

:::{tip}
While __pe*ga*susQC__ remains on the GA private repository, it can be installed with the following steps:

Clone the repository using Git (or download the zip file from github, and extract the contents):

```bash
git clone https://github.com/GeoscienceAustralia/pegasusQC.git
```
Navigate to the cloned directory ("pegasusQC" for this example):

```bash
cd pegasusQC
```
Install the package via setup.py:

```bash
pip install .\name_of_directory_containing_setup_py
```
:::

## Session Setup

In a Jupyter-lab notebook, set up your session with the following commands.

Firstly import the necessary Python packages ...

```python
from pathlib import Path # useful for file names
%matplotlib widget # needed to, for example, zoom into plots
```

... then import __pe*ga*susQC__ ... 

```python
import pegasusQC as qc
```

... and finally set the path to your field and plan data:

```python
# Setup the path to the acquired data, ...
data_root = "/Users/username/Documents/GitHub/pegasusQC/examples/"
dx = Path(data_root + r'AGG/Canobie/20211130.xyz')
dh = dx.with_suffix('.hdf5')

# ..., and plan data.
plan_root = data_root
px = Path(plan_root + r'AGG/Canobie/902212_1.xyz')
ph = px.with_suffix('.hdf5')
```

The variables `dx`, `dh`, `px` and `ph` are often used to represent the four key data filenames. Of course, you can use any variable names for your projects but it is worth noting these because you will see them occasionally in the documentation.

If you want to see how to use __pe*ga*susQC__, go to [Tutorials](#tutorials-target) and follow along with the examples. You can also download the python tutorial notebooks and example data from github and experiment with the functions yourself.

Now you are ready to QC your data!

[^Jupyter]: <https://jupyter.org>

[^JupyterApp]: <https://github.com/jupyterlab/jupyterlab-desktop>


## Package Installation and Setup on National Computing infrastructure

Finally, some users might want to install __pe*ga*susQC__ on the **NCIS**. The following notes might prove useful.

	a. Setting up the python environment in NCI

Using the instructions in this [link](https://opus.nci.org.au/display/Help/3.3+Using+a+custom+Python+virtual+environment+in+JupyterLab), set up a `python3` environment in your specific location. Due to memory/space limitations for each user, it is suggested that you make this custom python environment in your project directory in NCI. The environment will have a folder with the same name in your working directory after it is built:

For the QAQC python scripts, we need to install these extra `python3` packages within our new python environment:
- &nbsp; Xarray (version 2022.12.0 or higher)
- &nbsp; Rioxarray
- &nbsp; netCDF4
- &nbsp; Ipympl
- &nbsp; Colorcet
- &nbsp; Filebrowser
- &nbsp; Matplotlib(version 3.6.2) 

In addition to the above packages, two customised scripts called `graphics.py` and `colors.py` can be downloaded from [here](https://github.com/jobar8/graphics). These packages shoud be placed in `site_packages` available in your python environment.

Also make sure that you have all the following packages as well - if anything missing:

- &nbsp; H5py
- &nbsp; Verde
- &nbsp; Pooch

Finally, in order to run the current QAQC python scripts on NCI, the following steps need to be taken:

- &nbsp; Analyst need to have an NCI account
- &nbsp; Login to https://ood.nci.org.au/pun/sys/dashboard
- &nbsp; Choose Virtual Desktop 
- &nbsp; Except project that needs to be changed to your active project the rest of settings should be accepted as default 
- &nbsp; Launch the VDI app
- &nbsp; From the menu choose the Terminal 
- &nbsp; Navigate to your working directory 
- &nbsp; Load `python3` package available in NCI
``` 
[user@ood-vn30 user]$ module load python3/3.9.2
	Loading python3/3.9.2
	 Loading requirement: intel-mkl/2020.4.304
```
In your working directory you will have two active subfolders:
- &nbsp; Directory that contains your code and notebooks
- &nbsp; Directory that contains your specific python environment
Load your specific `python3` environment for instance loading environment called QA 	
```
[user@ood-vn30 user]$ . QA/bin/activate
     (QA) [user@ood-vn30 user]$
```
- &nbsp;  Run jupyter notebook in your VDI using the command:
```
Jupyter notebook
```

At this stage, you will be redirected to a web browser that you can load the specific notebook and quality check your airborne gravity dataset. 
