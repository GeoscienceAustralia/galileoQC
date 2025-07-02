"""
A library of geophysical QA/QC functions for line data.

References:
    W. J. Hinze, C. Aiken, J. M. Brozena, B. Coakley, D. Dater, G. Flanagan, 
    R. Forsberg, T. G. Hildenbrand, G. R. Keller, J. Kellogg, R. Kucks, X. Li, 
    A. Mainville, R. Morin, M. Pilkington, D. Plouff, D. Ravat, D. Roman, 
    J. Urrutia-Fucugauchi, M. V ́eronneau, M. Webring, and D. Winester. 
    New standards for reducing gravity data: The North American gravity database. 
    Geophysics, 70(4):J25, 2005.
    
    C. Jekeli, Theoretical fundamentals of airborne gradiometry. In Airborne 
    Gravity for Geodesy Summer School, 23-27 May, 2016.
"""

# __version__ = "0.0.0"

# Make it easy for user to call any function `myfun` as `aqc.myfun`.
# These lines access all functions and classes in all .py files
# in the named directories. The user can then do:
#     import AirGravQC as qc
#     qc.myfun()
from AirGravQC.whizzFiles.pointfiles import (updateCoordFrame, updateProject, 
    updateChannelAttributes)
from AirGravQC.qualitycontrol.qualityAnalysis import *
from AirGravQC.whizzPlots.whizzPlot import *
from AirGravQC.gridFiles.gridfiles import (gridfile_to_xa)
from AirGravQC.gridFiles.read_ers import *

# Data management
from AirGravQC.whizzFiles.gdfReader import asegToHDF
from AirGravQC.whizzFiles.xyzReader import xyzToHDF
from AirGravQC.whizzFiles.updateLineAttributes import updateLineAttributes
from AirGravQC.whizzFiles.retrieveData import (getWhizzData, getLineData, getChannelAttrs, getLineXChannel)
from AirGravQC.whizzFiles.reportData import (reportWhizz, reportFlights, reportSampling, reportChannels, reportLines)
from AirGravQC.whizzFiles.updateLineSpacing import updateLineSpacing
from AirGravQC.whizzFiles.updateOddOrEven import updateOddOrEven
from AirGravQC.gridFiles.gridToNc4 import gridToNc4

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
from AirGravQC.gridFiles.display_grid import display_grid
from AirGravQC.gridFiles.imageStats import imageStats
from AirGravQC.gridFiles.traceImages import traceImages
from AirGravQC.gridFiles.checkTCratio import checkTCratio
from AirGravQC.gridFiles.xdImage import xdImage
from AirGravQC.gridFiles.diff_n_image import diff_n_image
from AirGravQC.gridFiles.xarray_to_grid import xarray_to_grid
from AirGravQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from AirGravQC.gridFiles.oddeven import updateLineTracks
from AirGravQC.gridFiles.oddeven import oddevenlines, altsample_grid
from AirGravQC.gridFiles.imageAllInDir import imageAllInDir

# AGG (Falcon) checks
from AirGravQC.qualitycontrol.diffNoiseVturb import diffNoiseVturb
from AirGravQC.qualitycontrol.diffGravVturb import diffGravVturb
from AirGravQC.qualitycontrol.checkPhase import checkPhase

# FTG and AGG checks
from AirGravQC.qualitycontrol.checkHighFreq import checkHighFreq

# FTG checks
from AirGravQC.qualitycontrol.checkFrobenius import checkFrobenius
from AirGravQC.qualitycontrol.ilsNoiseAnalysis import ilsNoiseAnalysis
from AirGravQC.qualitycontrol.checkInlineSum import checkInlineSum


# Navigation and position checks
from AirGravQC.qualitycontrol.checkClearance import (checkSafeClearance, checkClearance, checkDrape)
from AirGravQC.qualitycontrol.checkIntersection import checkIntersection
from AirGravQC.qualitycontrol.checkRMSIntersection import checkRMSIntersection
from AirGravQC.qualitycontrol.checkSpeeds import checkSpeeds
from AirGravQC.qualitycontrol.checkOverlaps import checkOverlaps
from AirGravQC.qualitycontrol.checkXYPlan import checkXYPlan
from AirGravQC.qualitycontrol.checkVertPlan import checkVertPlan
from AirGravQC.qualitycontrol.checkLineLengths import checkLineLengths
from AirGravQC.qualitycontrol.checkGNSS import checkGNSS
from AirGravQC.qualitycontrol.checkHeading import checkHeading

# General checks
from AirGravQC.qualitycontrol.checkErsHeaders import checkErsHeaders
from AirGravQC.qualitycontrol.psdChannelDiff import psdChannelDiff
from AirGravQC.qualitycontrol.psdChannelDiff import psdChannelGain
from AirGravQC.qualitycontrol.allChanStats import allChanStats
from AirGravQC.qualitycontrol.diffChanStats import diffChanStats
from AirGravQC.qualitycontrol.checkGaps import checkGaps
from AirGravQC.qualitycontrol.checkSpikes import checkSpikes
from AirGravQC.qualitycontrol.checkConstantSlope import checkConstantSlope

# Airborne gravimetry checks
from AirGravQC.qualitycontrol.diffGroundGrid import diffGroundGrid
from AirGravQC.qualitycontrol.checkRepeatLines import checkRepeatLines
from AirGravQC.qualitycontrol.checkFreeAirCorr import checkFreeAirCorr
from AirGravQC.qualitycontrol.checkEotvosCorr import checkEotvosCorr
from AirGravQC.qualitycontrol.checkLatCorr import checkLatCorr
from AirGravQC.qualitycontrol.checkAtmosEffect import checkAtmosEffect


