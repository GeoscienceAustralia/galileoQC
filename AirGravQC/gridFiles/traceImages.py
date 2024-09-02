from AirGravQC.gridFiles.xdImage import xdImage
from AirGravQC.gridFiles.gridfiles import gridfile_to_xa
import colorcet as cc


def traceImages(file1, file2, file3, plotTitle):
    """
    Plots an image of the trace of the three diagonal components of the tensor.

    Parameters
    ----------
    file1 : String
        A grid file containing the first diagonal tensor component. May be either an `ERS` or `NC` grid file.
    file2 : String
        A grid file containing the second diagonal tensor component. May be either an `ERS` or `NC` grid file.
    file3 : String
        A grid file containing the third diagonal tensor component. May be either an `ERS` or `NC` grid file.
    plotTitle : String
        A title for the plot of the trace image.

    Returns
    -------
    None.

    """
    (x1, _) = gridfile_to_xa(file1, bandout=0)
    (x2, _) = gridfile_to_xa(file2, bandout=0)
    (x3, _) = gridfile_to_xa(file3, bandout=0)
    trace = x1 + x2 + x3
    xdImage(trace, plotTitle, colormap=cc.m_CET_L9, cmap_norm='no', hs=False, azdeg=90)

    # n1, e1, z1 = ers.read_ers_image(file1)
    # n2, e2, z2 = ers.read_ers_image(file2)
    # n3, e3, z3 = ers.read_ers_image(file3)
    # trace = z1 + z2 + z3
    # graphicsShaded(n1, e1, trace, plotTitle,  hs=False, colormap=cc.m_CET_L9, cmap_norm='no', azdeg=90)
