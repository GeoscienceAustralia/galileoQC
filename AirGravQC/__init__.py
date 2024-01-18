"""
A library of geophysical QA/QC functions for airborne gravimetry, gravity gradiometry surveys
(with a few QC functions for airborne magnetic surveys).

"""

<<<<<<< Updated upstream
=======
__version__ = "0.0.0"

# Make it easy for user to call any function `myfun` as `qc.myfun`.
# These lines access all functions in AirGravQC that are intended for direct access.
# Intended useage:
#     > import AirGravQC as qc
#     > qc.myfun()
from AirGravQC.whizzFiles.pointfiles import (xyzToHDF, updateCoordFrame, updateProject, 
    reportSampling, reportFlights, reportWhizz, updateChannelAttributes, updateLineAttributes,
    asegToHDF)
from AirGravQC.whizzPlots.whizzPlot import *
# from AirGravQC.qc.qualityAnalysis import *

# 1D Plotting
from AirGravQC.whizzPlots.plotChannelsLines import plotChannelsLines
from AirGravQC.whizzPlots.plotLinesCompareChannels import plotLinesCompareChannels
from AirGravQC.whizzPlots.linesMap import linesMap
from AirGravQC.whizzPlots.plotLinesOnGroundStns import plotLinesOnGroundStns

# 2D Plotting
from AirGravQC.gridFiles.read_ers import *
from AirGravQC.gridFiles.gridfiles import (image_pygmt, diff_n_image, display_grid, ers_to_netcdf4)
from AirGravQC.gridFiles.graphicsShaded import graphicsShaded
from AirGravQC.gridFiles.grid_n_image import grid_n_image
from AirGravQC.gridFiles.xdImage import xdImage
from AirGravQC.gridFiles.xarray_to_grid import xarray_to_grid
from AirGravQC.gridFiles.whizz_to_xarray import whizz_to_xarray

# Mag
from AirGravQC.qc.checkDiurnal import checkDiurnal
from AirGravQC.qc.checkBasemag import checkBasemag
from AirGravQC.qc.checkTCDiff4 import checkTCDiff4

# AGG (Falcon) checks
from AirGravQC.qc.diffNoiseVturb import diffNoiseVturb
from AirGravQC.qc.diffGravVturb import diffGravVturb

# FTG and AGG checks
from AirGravQC.qc.checkHighFreq import checkHighFreq
from AirGravQC.qc.checkPhase import checkPhase

# FTG checks
from AirGravQC.qc.checkFrobenius import checkFrobenius
from AirGravQC.qc.ilsNoiseVturb import ilsNoiseVturb
from AirGravQC.qc.checkInlineSum import checkInlineSum


# Navigation and position checks
from AirGravQC.qc.checkClearance import (checkSafeClearance, checkClearance, checkDrape)
from AirGravQC.qc.checkIntersection import checkIntersection
from AirGravQC.qc.checkSpeeds import checkSpeeds
from AirGravQC.qc.checkOverlaps import checkOverlaps
from AirGravQC.qc.checkXYPlan import checkXYPlan
from AirGravQC.qc.checkVertPlan import checkVertPlan
from AirGravQC.qc.checkLineLengths import checkLineLengths
from AirGravQC.qc.checkGNSS import checkGNSS
from AirGravQC.qc.checkHeading import checkHeading

# General checks
from AirGravQC.qc.checkErsHeaders import checkErsHeaders
from AirGravQC.qc.psdChannelDiff import psdChannelDiff
from AirGravQC.qc.allChanStats import allChanStats
from AirGravQC.qc.checkDiffStats import checkDiffStats
from AirGravQC.qc.vertAccStats import vertAccStats
from AirGravQC.qc.checkGaps import checkGaps
from AirGravQC.qc.checkConstantSlope import checkConstantSlope
from AirGravQC.qc.checkSpikes import checkSpikes

# Airborne gravimetry checks
from AirGravQC.qc.diffGroundGrid import diffGroundGrid
from AirGravQC.qc.checkRepeatLines import checkRepeatLines
from AirGravQC.qc.checkFreeAirCorr import checkFreeAirCorr
from AirGravQC.qc.checkEotvosCorr import checkEotvosCorr
from AirGravQC.qc.checkLatCorr import checkLatCorr
from AirGravQC.qc.checkAtmosEffect import checkAtmosEffect

from AirGravQC.graphics import *

>>>>>>> Stashed changes
