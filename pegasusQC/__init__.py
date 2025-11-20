#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A library of geophysical QA/QC functions for airborne gravity (and magnetic) survey data.

Author: Mark Helm Dransfield

Created: 2023

License: CC BY-SA
"""

# Make it easy for user to call any function `myfun` as `aqc.myfun`.
# These lines access all functions and classes in all .py files
# in the named directories. The user can then do:
#     import pegasusQC as qc
#     qc.myfun()
from pegasusQC.whizzFiles.pointfiles import (updateCoordFrame, updateProject)
from pegasusQC.qualitycontrol.qualityAnalysis import *
from pegasusQC.whizzPlots.whizzPlot import *
from pegasusQC.gridFiles.grid_to_xarray import gridfile_to_xa
from pegasusQC.gridFiles.read_ers import *

# Data management
from pegasusQC.whizzFiles.gdfReader import asegToHDF
from pegasusQC.whizzFiles.xyzReader import xyzToHDF
from pegasusQC.whizzFiles.colReader import colsToHDF
from pegasusQC.whizzFiles.updateLineAttributes import updateLineAttributes
from pegasusQC.whizzFiles.updateChannelAttributes import updateChannelAttributes
from pegasusQC.whizzFiles.retrieveData import (getWhizzData, getLineData, getChannelAttrs, getLineXChannel)
from pegasusQC.whizzFiles.reportData import (reportWhizz, reportFlights, reportSampling, reportChannels, reportLines)
from pegasusQC.whizzFiles.updateLineSpacing import updateLineSpacing
from pegasusQC.whizzFiles.updateLineVariety import updateLineVariety
from pegasusQC.whizzFiles.updateOddOrEven import updateOddOrEven
from pegasusQC.gridFiles.gridToNc4 import gridToNc4

# 1D Plotting
from pegasusQC.whizzPlots.plotChannelsLines import plotChannelsLines
from pegasusQC.whizzPlots.plotLinesCompareChannels import plotLinesCompareChannels
from pegasusQC.whizzPlots.linesMap import linesMap
from pegasusQC.whizzPlots.psdLineChannels import psdLineChannels
from pegasusQC.whizzPlots.plotBoxWhisker import plotBoxWhisker
from pegasusQC.whizzPlots.plotLinesOnGroundStns import plotLinesOnGroundStns

# 2D Plotting
from pegasusQC.gridFiles.graphicsShaded import graphicsShaded
from pegasusQC.gridFiles.grid_n_image import grid_n_image
from pegasusQC.gridFiles.minc import minc
from pegasusQC.gridFiles.display_grid import display_grid
from pegasusQC.gridFiles.imageStats import imageStats
from pegasusQC.gridFiles.traceImages import traceImages
from pegasusQC.gridFiles.checkTCratio import checkTCratio
from pegasusQC.gridFiles.xdImage import xdImage
from pegasusQC.gridFiles.diff_n_image import diff_n_image
from pegasusQC.gridFiles.xarray_to_grid import xarray_to_grid
from pegasusQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from pegasusQC.gridFiles.oddeven import updateLineTracks
from pegasusQC.gridFiles.oddeven import oddevenlines, altsample_grid
from pegasusQC.gridFiles.imageAllInDir import imageAllInDir

# AGG (Falcon) checks
from pegasusQC.qualitycontrol.diffNoiseVturb import diffNoiseVturb
from pegasusQC.qualitycontrol.diffGravVturb import diffGravVturb
from pegasusQC.qualitycontrol.checkPhase import checkPhase

# FTG and AGG checks
from pegasusQC.qualitycontrol.checkHighFreq import checkHighFreq

# FTG checks
from pegasusQC.qualitycontrol.checkFrobenius import checkFrobenius
from pegasusQC.qualitycontrol.ilsNoiseAnalysis import ilsNoiseAnalysis

# Aeromag checks
from pegasusQC.qualitycontrol.checkDiurnal import checkDiurnal
from pegasusQC.qualitycontrol.checkDifference import checkDiff


# Navigation and position checks
from pegasusQC.qualitycontrol.checkClearance import (checkSafeClearance, checkClearance, checkDrape)
from pegasusQC.qualitycontrol.checkIntersection import checkIntersection
from pegasusQC.qualitycontrol.checkSpeeds import checkSpeeds
from pegasusQC.qualitycontrol.checkOverlaps import checkOverlaps
from pegasusQC.qualitycontrol.checkXYPlan import checkXYPlan
from pegasusQC.qualitycontrol.checkVertPlan import checkVertPlan
from pegasusQC.qualitycontrol.checkLineLengths import checkLineLengths
from pegasusQC.qualitycontrol.checkGNSS import checkGNSS
from pegasusQC.qualitycontrol.checkHeading import checkHeading

# General checks
from pegasusQC.qualitycontrol.checkErsHeaders import checkErsHeaders
from pegasusQC.qualitycontrol.psdChannelDiff import psdChannelDiff
from pegasusQC.qualitycontrol.psdChannelDiff import psdChannelGain
from pegasusQC.qualitycontrol.psdChannelDiff import psdChannel
from pegasusQC.qualitycontrol.allChanStats import allChanStats
from pegasusQC.qualitycontrol.diffChanStats import diffChanStats
from pegasusQC.qualitycontrol.checkGaps import checkGaps
from pegasusQC.qualitycontrol.checkSpikes import checkSpikes
from pegasusQC.qualitycontrol.checkConstantSlope import checkConstantSlope

# Airborne gravimetry checks
from pegasusQC.qualitycontrol.diffGroundGrid import diffGroundGrid
from pegasusQC.qualitycontrol.checkRepeatLines import checkRepeatLines
from pegasusQC.qualitycontrol.checkFreeAirCorr import checkFreeAirCorr
from pegasusQC.qualitycontrol.checkEotvosCorr import checkEotvosCorr
from pegasusQC.qualitycontrol.checkLatCorr import checkLatCorr
from pegasusQC.qualitycontrol.checkAtmosEffect import checkAtmosEffect

# Data Transforms
from pegasusQC.transforms.craig_transform import craig_transform
