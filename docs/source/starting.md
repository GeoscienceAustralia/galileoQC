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

In a Jupyter-lab notebook, set up your session with the following commands.

Firstly import the necessary Python packages ...

```python
# to set the path to your copy of `AirGravQC` - might not be needed.
import sys
from pathlib import Path # useful for file names
```

... then set up the path to AirGravQC (this should not be needed when the package is released), and import the needed modules ... 

```python
# replace the text in "" with the path to your copy of `AirGravQC`
local_docs = "/Users/markdransfield/"
src_path = local_docs + "Documents/GitHub/AirGravQC/AirGravQC"
sys.path.insert(0, src_path)

# Now import the `AirGravQC` modules
import pointfiles as mhd
import qualityAnalysis as qc
import whizzPlot as wp
import gridfiles as erm
%matplotlib widget # needed to interact, eg zoom, with plots
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

Go to `Methods` and choose the QC methods you need, or go to `Cookbooks` and follow along one of the examples.

[^Jupyter]: <https://jupyter.org>

[^JupyterApp]: <https://github.com/jupyterlab/jupyterlab-desktop>
