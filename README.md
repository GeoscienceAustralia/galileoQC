# AirGravQC
A python package to quality check the airborne gravity dataset

**Package Installation and Setup on National Computing infrastructure**

	a. Setting up the python environment in NCI

Using the instructions in this [link](https://opus.nci.org.au/display/Help/3.3+Using+a+custom+Python+virtual+environment+in+JupyterLab), and set up a `python3` environment in your specific location. For memory/space limitations for each user, it is suggested that you make this custom python environment in your project directory. The environment will have a folder with the same name in your working directory after it is built:

For the QAQC python scripts, we need to install these extra python3 packages within our new python environment:

	-Xarray (version 2022.12.0 or higher)
	-Rioxarray
	-netCDF4
	-Ipympl
	-Colorcet
	-Filebrowser
	-Matplotlib(version 3.6.2) 

In addition to the above packages, two customised scripts called `graphics.py` and `colors.py` can be downloaded from [here](https://github.com/jobar8/graphics). These packages shoud be placed both in the subfolder whizz as well as site_packages available in your python environment:

/path to the working directory in NCI/QA/lib/python3.9/site-packages

After these steps make sure that you have all the following packages as well - if any missing:

	-Numpy
 
	-H5py

	-Verde

	-Pooch



	b. Running notebook in NCI
In order to run the current QAQC python scripts developed by Dr. Mark Dransfield on NCI, the following steps need to be taken:

	- Analyst need to have an NCI account
	- Login to https://ood.nci.org.au/pun/sys/dashboard
	- Choose Virtual Desktop 
	- Except project that needs to be changed to your active project the rest of settings should be accepted as default 
	- Launch the VDI app
	- From the menu choose the Terminal 
	- Navigate to your working directory 
	- Load python3 package available in NCI 

		` [user@ood-vn30 user]$ module load python3/3.9.2
		  Loading python3/3.9.2
 		    Loading requirement: intel-mkl/2020.4.304`

	- In your working directory you will have two active subfolders:

		- Directory that contains your code and notebooks
 
		- Directory that contains your specific python environment
 
	- Load your specific python3 environment for instance loading environment call QA
		
		` [user@ood-vn30 user]$ . QA/bin/activate
                  (QA) [user@ood-vn30 user]$`

	- Run jupyter notebook in your VDI using the command:
	
		`Jupyter notebook` 

At this stage, you will be redirected to browser that you can load the specific notebook to test on your airborne gravity dataset. 
 
   
