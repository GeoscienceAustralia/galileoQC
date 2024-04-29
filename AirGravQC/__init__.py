"""
A library of geophysical QA/QC functions for line data.

References:
    C. Aiken, J. M. Brozena, B. Coakley, D. Dater, G. Flanagan, R. Forsberg, 
    J. Kellogg, R. Kucks, X. Li, A. Mainville, R. Morin, M. Pilkington, 
    D. Plouff, D. Roman, J. Urrutia-Fucugauchi, M. V ́eronneau, and D. Winester. 
    New standards for reducing gravity data: The North American gravity 
    database. Geophysics, 70(4):J25, 2005.
    
    C. Jekeli, Theoretical fundamentals of airborne gradiometry. In Airborne 
    Gravity for Geodesy Summer School, 23-27 May, 2016.
"""

__version__ = "0.0.0"

# Make it easy for user to call any function `myfun` as `aqc.myfun`.
# These lines access all functions and classes in all .py files
# in the named directories. The user can then do:
#     import AirGravQC as qc
#     qc.myfun()
from AirGravQC.whizzFiles.pointfiles import (updateCoordFrame, updateProject, 
    updateChannelAttributes)
from AirGravQC.whizzPlots.whizzPlot import *
from AirGravQC.gridFiles.gridfiles import (image_pygmt, diff_n_image, display_grid)
from AirGravQC.gridFiles.read_ers import *
# from AirGravQC.qc.qualityAnalysis import * Commented out to see if we actually ever use any of these!!

# Data management
from AirGravQC.whizzFiles.asegReader import asegToHDF
from AirGravQC.whizzFiles.xyzReader import xyzToHDF
from AirGravQC.whizzFiles.updateLineAttributes import updateLineAttributes
from AirGravQC.whizzFiles.retrieveData import (getWhizzData, getLineData, getChannelAttrs, getLineXChannel)
from AirGravQC.whizzFiles.reportData import (reportWhizz, reportFlights, reportSampling)

# 1D Plotting
from AirGravQC.whizzPlots.plotChannelsLines import plotChannelsLines
from AirGravQC.whizzPlots.plotLinesCompareChannels import plotLinesCompareChannels
from AirGravQC.whizzPlots.linesMap import linesMap
from AirGravQC.whizzPlots.psdLineChannels import psdLineChannels
from AirGravQC.whizzPlots.plotBoxWhisker import plotBoxWhisker
from AirGravQC.whizzPlots.plotLinesOnGroundStns import plotLinesOnGroundStns

# 2D Plotting
from AirGravQC.gridFiles.graphicsShaded import graphicsShaded
from AirGravQC.gridFiles.grid_n_image import grid_n_image
from AirGravQC.gridFiles.xdImage import xdImage
from AirGravQC.gridFiles.xarray_to_grid import xarray_to_grid
from AirGravQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from AirGravQC.gridFiles.graphics import imshow_hs # ??? TODO Is this needed here?

# Mag
from AirGravQC.qc.checkDiurnal import checkDiurnal
from AirGravQC.qc.checkBasemag import checkBasemag
from AirGravQC.qc.checkTCDiff4 import checkTCDiff4

# AGG (Falcon) checks
from AirGravQC.qc.diffNoiseVturb import diffNoiseVturb
from AirGravQC.qc.diffGravVturb import diffGravVturb
from AirGravQC.qc.checkPhase import checkPhase

# FTG and AGG checks
from AirGravQC.qc.checkHighFreq import checkHighFreq

# FTG checks
from AirGravQC.qc.checkFrobenius import checkFrobenius
from AirGravQC.qc.ilsNoiseVturb import ilsNoiseVturb
from AirGravQC.qc.checkInlineSum import checkInlineSum


# Navigation and position checks
from AirGravQC.qc.checkClearance import (checkSafeClearance, checkClearance, checkDrape)
from AirGravQC.qc.checkIntersection import checkIntersection
from AirGravQC.qc.checkRMSIntersection import checkRMSIntersection
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
from AirGravQC.qc.psdChannelDiff import psdChannelGain
from AirGravQC.qc.allChanStats import allChanStats
from AirGravQC.qc.diffChanStats import diffChanStats
from AirGravQC.qc.checkGaps import checkGaps
from AirGravQC.qc.checkSpikes import checkSpikes
from AirGravQC.qc.checkConstantSlope import checkConstantSlope

# Airborne gravimetry checks
from AirGravQC.qc.diffGroundGrid import diffGroundGrid
from AirGravQC.qc.checkRepeatLines import checkRepeatLines
from AirGravQC.qc.checkFreeAirCorr import checkFreeAirCorr
from AirGravQC.qc.checkEotvosCorr import checkEotvosCorr
from AirGravQC.qc.checkLatCorr import checkLatCorr
from AirGravQC.qc.checkAtmosEffect import checkAtmosEffect


