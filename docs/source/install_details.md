(installs-target)=
# Install Details

This document describes the steps for installing pegasusQC in a variety of operating systems. The steps include those for installing python and setting up the environment. 

There are many ways to accomplish these steps. Here one way is documented, together with a little explanation. Feel free to use a different way if you have the expertise (or at least can find the knowledge on the web).

These instructions use the python venv environment manager and pip installer.

It is important to have a single virtual environment set up, with pagasusQC and jupyterlab, installed in that environment. Throughout this installation guide, the environment used for pegasusQC is named `whizz`. If you use a different name, just replace `whizz` with that name.

For each operating system, the steps are:

1. create a python virtual environment;
1. download and install the pegasusQC package;
1. install jupyterlab.

For Windows installations, you may need to install python first. This can be done in several ways; the `winget` method is given here.

*At the time of writing, pegasusQC is in a private repo and this installation will only work if you have access to that repo.*

## Windows 10 and Windows 11 pip install

Go into Settings -> Apps -> Advanced app settings -> App execution aliases, and turn off the App Installers for python.exe and python3.exe.

The following commands must all be entered from the `cmd` window (not the Powershell).

Firstly install Python, then create a virtual environment `whizz` (you can choose another name if you like but it is `whizz` in the documentation), and then activate the virtual environment.

```bash
$ winget install Python.Python.3.13
$ python -m venv whizz
$ whizz\Scripts\activate
```

The prompt shows that you are in the whizz environment, ready to install `pegasusQC`. Currently, pegasusQC is in a private repository. While this remains the case, we will proceed by downloading the repo using a web browser (one needs permission here), and then installing locally.

Do this from your web browser by proceeding to the `pegasusQC` repo; click on Code (green button); then Download ZIP. Once downloaded, extract the contents to your HOME directory, or other directory of choice.

Now use pip to install it.

```bash
(whizz) $ pip install .\name_of_directory_containing_setup_py
```

And then install jupyterlab.

```bash
(whizz) $ pip install jupyterlab
```

Now everything is installed and you should be able to run pegasusQC. To see how this is done, run some of the tutorials from jupyterlab. This should open up in your default browser (a 64-bit browser is strongly recommended).

```bash
(whizz) $ jupyter lab
```

## Linux pip install (macos, fedora, debian and ubuntu)

First, you will need a working version of `python3`, as well as the `pip` installer and the `venv` virtual environment manager. These should be readily available on any linux system so we assume they are ready to go. If not, a search on the web will provide instructions for installing them (some brief notes are included further below). The installations are different for different operating systems. On macos, this needs some care because there is usually an older version of python already installed but the python install instructions will help you here.

Here, we create an environment called `whizz` within which we will install and run pegasusQC. (You don't have to call the environment `whizz` but that is the name used in all the documentation.)

Open a terminal window and enter the following commands. Note that this documentation does not yet specify version numbers of software.

Create the whizz environment and activate it.

```bash
$ python3 -m venv whizz
$ source whizz/bin/activate
```

The prompt shows that you are in the `whizz` environment, ready to install `pegasusQC`. Currently, `pegasusQC` is in a private repository. While this remains the case, we will proceed by downloading the repo using a web browser (one needs permission here), and then installing locally.

Do this from your web browser by proceeding to the `pegasusQC` repo; click on Code (green button); then Download ZIP. Once downloaded, extract the contents to your HOME directory, or other directory of choice.

Now use pip to install it.

```bash
(whizz) $ pip install ./name_of_directory_containing_setup_py
```

Now you need a way of running `pegasusQC`. Recommended for this purpose is `jupyterlab`. This is easily installed:

```bash
(whizz) $ pip install jupyterlab
```

Now everything is installed and you should be able to run `pegasusQC`. To see how this is done, run some of the tutorials from jupyterlab:

```bash
(whizz) $ jupyter lab
```

The Jupyterlab web app should open in your browser.

### Preparing linux

The following commands might be useful in preparing your system for the `pegasusQC` and jupyterlab installations.

__DEBIAN__

```bash
$ sudo apt update && sudo apt upgrade
$ sudo apt install python3.11-venv

```

__UBUNTU__

```bash
$ sudo apt update && sudo apt upgrade
$ sudo apt install python3.10-venv
```

__FEDORA__

`venv` is installed with python3.

```bash
$ sudo dnf upgrade
```

## Running pegasusQC

The Jupyterlab web app should look something like this screenshot:

```{figure} jupyterlab_start.png
:alt: 
:align: center
:scale: 25

Start page for pegasusQC in jupyterlab.
```

Navigate to pegasusQC-main/docs/source/tutorials and run one or more of the tutorial `.ipynb` files, or start a new notebook from the `File` menu. Here is the Odd - Even Analysis notebook from the tutorials.

```{figure} jupyterlab_tutorial.png
:alt: 
:align: center
:scale: 25

Running the pegasusQC Odd - Even Analysis tutorial.
```

## Upgrading pegasusQC

Because we are installing from a downloaded copy of the `pegasusQC` repository, upgrading to the latest version involves a few extra steps. First, uninstall `pegasusQC` from the whizz environment:

```powershell
(whizz) > python -m pip uninstall pegasusQC
```

Then follow the steps above to download the current version and install it.

## Uninstalls

### Remove pegasusQC and the whole environment

If you have a problem with the installation and want to get rid of everything and start again (or if you just want to remove `pegasusQC` completely), you can remove the `whizz` environment. (**This is not advisable if you have installed other packages in the environment that you want to use**).

To remove only `pegasusQC` is eaily done in any of our tested operating systems. Navigate to the working directory with `whizz` activated and use `pip`.

```bash
(whizz) > pip uninstall pegasusQC
```

If you want to uninstall the `whizz` environment and all of its contents:

__Windows__

```bash
(whizz) > deactivate
> rmdir /s /q whizz
```

__linux__ and __macos__

```bash
(whizz) $ deactivate
(base) $ sudo rm -rf whizz
```

You may be asked to confirm the remove. Accept.

This removes `pegasusQC`, `jupyterlab`, and all other packages in the `whizz` environment. To recover, you will need to create a new environment, and install `pegasusQC` and `jupyterlab` in the environment, following the installation steps.



