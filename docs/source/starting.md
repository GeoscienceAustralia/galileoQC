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
:::

## Session Setup

In a Jupyter-lab notebook, set up your session with the following commands. These should be interspersed with useful commentary.

Firstly import the necessary Python packages (actually not all of these are necessary). TODO - CLEAN UP!!

```python
# to set the path to your copy of `AirGravQC` - might not be needed.
import sys
# In case you change the s/w. If you do, let us know! Maybe not needed
import importlib as im 

import numpy as np # perhaps not needed
import xarray as xr # perhaps not needed
import matplotlib.pyplot as plt # needed for plots
%matplotlib widget # needed to interact, eg zoom, with plots

from pathlib import Path # useful for file names
```

```python
# replace the text in "" with the path to your copy of `AirGravQC`
local_docs = "/Users/markdransfield/"
src_path = local_docs + "Documents/GitHub/AirGravQC/src"
sys.path.insert(0, src_path)

# Now import the `AirGravQC` modules
import pointfiles as mhd
import qualityAnalysis as qc
import whizzPlot as wp
import gridfiles as erm
```

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

Now you are ready to QC the data!

Go to `Methods` and choose the QC methods you need, or go to `Cookbooks` and follow along one of the examples.

[^Jupyter]: <https://jupyter.org>

[^JupyterApp]: <https://github.com/jupyterlab/jupyterlab-desktop>
