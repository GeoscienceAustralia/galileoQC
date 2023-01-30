# AirGravQC
AirGravQC is a python package to check the quality of the airborne gravity dataset delivered by suppliers.

**Package Installation and Setup on National Computing infrastructure**

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
