# Getting Started

## Installation

To use AirGravQC, first install it into your `whizz` environment using pip or conda:

```bash
(whizz) $ pip install AirGravQC
```

```bash
(whizz) $ conda install AirGravQC
```

Then start up a Jupyter-lab notebook. For more information about this, see [^Jupyter] . The Jupyter-lab app [^JupyterApp] is recommended.

:::{tip}
Installation not yet working!

While AirGravQC remains on the GA private repository, it can be installed with the following steps:

Clone the repository using Git:

```bash
git clone https://github.com/GeoscienceAustralia/AirGravQC.git
```
Navigate to the cloned directory:

```bash
cd AirGravQC
```
Install the package via setup.py:

```bash
sudo python setup.py install
```
:::

## Session Setup

In a Jupyter-lab notebook, set up your session with the following commands.

Firstly import the necessary Python packages ...

```python
from pathlib import Path # useful for file names
%matplotlib widget # needed to, for example, zoom into plots
```

... then import AirGravQC ... 

```python
import AirGravQC as qc
```

... and finally set the path to your field and plan data:

```python
# Setup the path to the acquired data, ...
data_root = "/Users/username/Documents/GitHub/AirGravQC/examples/"
dx = Path(data_root + r'AGG/Canobie/20211130.xyz')
dh = dx.with_suffix('.hdf5')

# ..., and plan data.
plan_root = data_root
px = Path(plan_root + r'AGG/Canobie/902212_1.xyz')
ph = px.with_suffix('.hdf5')
```

The variables `dx`, `dh`, `px` and `ph` are used throughout the documentation to represent the four key data filenames. Of course, you can use any variable names for your projects but it is worth noting these because you will see them often in the documentation.

Now you are ready to QC the data!

Go to [Methods](#methods-target) and choose the QC methods you need, or go to [Tutorials](#tutorials-target) and follow along one of the examples.

[^Jupyter]: <https://jupyter.org>

[^JupyterApp]: <https://github.com/jupyterlab/jupyterlab-desktop>
