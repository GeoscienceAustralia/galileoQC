# Getting Started

The __galileoQC__ package is run from a Jupyterlab notebook. You should create a virtual working environment (called `whizz` in this documentation) with python venv.

## Installation

To use __galileoQC__, first install it into your `whizz` environment using pip (conda installations of __galileoQC__ can work but conda has not been tested as fully as pip):

```bash
(whizz) $ pip install galileoQC
```

Then start up a jupyterlab notebook. For more information about jupyterlab, see [^Jupyter] . The Jupyterlab-Desktop app [^JupyterApp] also works.

For more information on installation see [Install Details](#installs-target).

```{toctree}
:maxdepth: 2
:hidden:
install_details.md
```

## Session Setup

In a jupyterlab notebook, set up your session with the following commands.

Firstly import the necessary Python packages ...

```python
from pathlib import Path # useful for file names
%matplotlib widget # needed to, for example, zoom into plots
```

... then import __galileoQC__ ... 

```python
import galileoQC as qc
```

... and finally set the path to your field and plan data. This might look something like this:

```python
# Setup the path to the acquired data, ...
data_root = r'.source/tutorials/'
dx = Path(data_root + r'CanobieData/Canobie.xyz')
dh = dx.with_suffix('.hdf5')

# ..., and plan data.
plan_root = data_root
px = Path(plan_root + r'CanobieData/CanobiePlan.xyz')
ph = px.with_suffix('.hdf5')
```

The variables `dx`, `dh`, `px` and `ph` are often used to represent the four key data filenames. Of course, you can use any variable names for your projects but it is worth noting these because you will see them occasionally in the documentation.

If you want to see how to use __galileoQC__, go to [Tutorials](#tutorials-target) and follow along with the examples. You can also download the python tutorial notebooks and example data from [github](https://github.com/GeoscienceAustralia/galileoQC/tree/main/docs/source/tutorials) and experiment with the functions yourself.

Now you are ready to QC your data!

[^Jupyter]: <https://jupyter.org>

[^JupyterApp]: <https://github.com/jupyterlab/jupyterlab-desktop>
