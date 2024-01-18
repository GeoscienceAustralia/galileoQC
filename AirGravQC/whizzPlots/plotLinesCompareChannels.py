import numpy as np
import matplotlib.pyplot as plt
import h5py
import xarray as xr
import verde as vd
import pooch

import AirGravQC.config as config
# import AirGravQC.qc.qualityAnalysis as qc
import AirGravQC.whizzFiles.pointfiles as gw
import AirGravQC.whizzPlots.whizzPlot as wpl
import AirGravQC.utility.utility as util
import matplotlib.ticker as tkr
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica Neue']})

groupName = config.groupName


def plotLinesCompareChannels(whizzFile, flightLines, x, channels, xOffset=True, mean_remove=False):
    """
    For the given flightLine in the whizzFile, plot both channel1 and channel2
    against x.

    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    flightLine : String
        A flightline, e.g. '1000110.0'.
    x : String
        The name of the independent variable for the plot.
    channel : String
        The name of the first channel or field to plot.
    channel : String
        The name of the second channel or field to plot.
    xOffset : Bool, optional
        If True, map x to x - x[0] before plotting. The default is True.
    mean_remove : Bool, optional
        If True, the y data will have their means subtracted before plotting. The default is False.

    Returns
    -------
    None.

    """
    
    filename = str(whizzFile)
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if x == '':
            x = f[groupName]['CoordinateFrame'].attrs['TimeChannel']
        plotPreTitle = wpl.make_plot_title(f[groupName])
        xDel = 0.0
        
        for line in flightLines:
            xData = gw.getLineData(g[line], x)
            xUnits = gw.getLineDataUnits(g[line], x)
            if xUnits == '':
                xlabelstr = f'{x}'
            else:
                xlabelstr = f'{x} [{xUnits}]'
            if xOffset:
                xData = xData - xData[0]
            fig = plt.figure(figsize=(6,9))
            thou_format = tkr.FuncFormatter(util._space_thou)
            ax = fig.add_subplot(1,1,1)
            for channel in channels:
                yData = gw.getLineData(g[line], channel)
                yUnits = gw.getLineDataUnits(g[line], channel)
                if yUnits == '':
                    ylabelstr = f'{channel}'
                else:
                    ylabelstr = f'{channel} [{yUnits}]'
                ax.plot(xData, yData, label=ylabelstr, lw=0.3)

            plotTitle = f'{plotPreTitle}, L{line}'
            ax.legend(fontsize=6)
            ax.xaxis.set_major_formatter(thou_format)
            ax.yaxis.set_major_formatter(thou_format)
            plt.xlabel(xlabelstr, fontsize = 6)
            plt.title(plotTitle, fontsize = 8)
            plt.grid(True)
            for label in ax.get_xticklabels(): label.set_fontsize(6)
            for label in ax.get_yticklabels(): label.set_fontsize(6)
            plt.show()























