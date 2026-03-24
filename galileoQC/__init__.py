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
#     import galileoQC as qc
#     qc.myfun()
from galileoQC.whizzFiles.pointfiles import (updateCoordFrame, updateProject)
from galileoQC.qualitycontrol.qualityAnalysis import *
from galileoQC.whizzPlots.whizzPlot import *
from galileoQC.gridFiles.grid_to_xarray import gridfile_to_xa
from galileoQC.gridFiles.read_ers import *

# Data management
from galileoQC.whizzFiles.gdfReader import asegToHDF
from galileoQC.whizzFiles.xyzReader import xyzToHDF
from galileoQC.whizzFiles.colReader import colsToHDF
from galileoQC.whizzFiles.updateLineAttributes import updateLineAttributes
from galileoQC.whizzFiles.updateChannelAttributes import updateChannelAttributes
from galileoQC.whizzFiles.retrieveData import (getWhizzData, getLineData, getChannelAttrs, getLineXChannel)
from galileoQC.whizzFiles.reportData import (reportWhizz, reportFlights, reportSampling, reportChannels, reportLines)
from galileoQC.whizzFiles.updateLineSpacing import updateLineSpacing
from galileoQC.whizzFiles.updateLineVariety import updateLineVariety
from galileoQC.whizzFiles.updateOddOrEven import updateOddOrEven
from galileoQC.gridFiles.gridToNc4 import gridToNc4

# 1D Plotting
from galileoQC.whizzPlots.plotChannelsLines import plotChannelsLines
from galileoQC.whizzPlots.plotLinesCompareChannels import plotLinesCompareChannels
from galileoQC.whizzPlots.linesMap import linesMap
from galileoQC.whizzPlots.psdLineChannels import psdLineChannels
from galileoQC.whizzPlots.plotBoxWhisker import plotBoxWhisker
from galileoQC.whizzPlots.plotLinesOnGroundStns import plotLinesOnGroundStns

# 2D Plotting
from galileoQC.gridFiles.graphicsShaded import graphicsShaded
from galileoQC.gridFiles.grid_n_image import grid_n_image
from galileoQC.gridFiles.minc import minc
from galileoQC.gridFiles.display_grid import display_grid
from galileoQC.gridFiles.imageStats import imageStats
from galileoQC.gridFiles.traceImages import traceImages
from galileoQC.gridFiles.checkTCratio import checkTCratio
from galileoQC.gridFiles.xdImage import xdImage
from galileoQC.gridFiles.diff_n_image import diff_n_image
from galileoQC.gridFiles.xarray_to_grid import xarray_to_grid
from galileoQC.gridFiles.whizz_to_xarray import whizz_to_xarray
from galileoQC.gridFiles.oddeven import updateLineTracks
from galileoQC.gridFiles.oddeven import oddevenlines, altsample_grid
from galileoQC.gridFiles.imageAllInDir import imageAllInDir

# AGG (Falcon) checks
from galileoQC.qualitycontrol.diffNoiseVturb import diffNoiseVturb
from galileoQC.qualitycontrol.diffGravVturb import diffGravVturb
from galileoQC.qualitycontrol.checkPhase import checkPhase

# FTG and AGG checks
from galileoQC.qualitycontrol.checkHighFreq import checkHighFreq

# FTG checks
from galileoQC.qualitycontrol.checkFrobenius import checkFrobenius
from galileoQC.qualitycontrol.ilsNoiseAnalysis import ilsNoiseAnalysis

# Aeromag checks
from galileoQC.qualitycontrol.checkDiurnal import checkDiurnal
from galileoQC.qualitycontrol.checkDifference import checkDiff


# Navigation and position checks
from galileoQC.qualitycontrol.checkClearance import (checkSafeClearance, checkClearance, checkDrape)
from galileoQC.qualitycontrol.checkIntersection import checkIntersection
from galileoQC.qualitycontrol.checkSpeeds import checkSpeeds
from galileoQC.qualitycontrol.checkOverlaps import checkOverlaps
from galileoQC.qualitycontrol.checkXYPlan import checkXYPlan
from galileoQC.qualitycontrol.checkVertPlan import checkVertPlan
from galileoQC.qualitycontrol.checkLineLengths import checkLineLengths
from galileoQC.qualitycontrol.checkGNSS import checkGNSS
from galileoQC.qualitycontrol.checkHeading import checkHeading

# General checks
from galileoQC.qualitycontrol.checkErsHeaders import checkErsHeaders
from galileoQC.qualitycontrol.psdChannelDiff import psdChannelDiff
from galileoQC.qualitycontrol.psdChannelDiff import psdChannelGain
from galileoQC.qualitycontrol.psdChannelDiff import psdChannel
from galileoQC.qualitycontrol.allChanStats import allChanStats
from galileoQC.qualitycontrol.diffChanStats import diffChanStats
from galileoQC.qualitycontrol.checkGaps import checkGaps
from galileoQC.qualitycontrol.checkSpikes import checkSpikes
from galileoQC.qualitycontrol.checkConstantSlope import checkConstantSlope

# Airborne gravimetry checks
from galileoQC.qualitycontrol.diffGroundGrid import diffGroundGrid
from galileoQC.qualitycontrol.checkRepeatLines import checkRepeatLines
from galileoQC.qualitycontrol.checkFreeAirCorr import checkFreeAirCorr
from galileoQC.qualitycontrol.checkEotvosCorr import checkEotvosCorr
from galileoQC.qualitycontrol.checkLatCorr import checkLatCorr
from galileoQC.qualitycontrol.checkAtmosEffect import checkAtmosEffect

# Data Transforms
from galileoQC.transforms.craig_transform import craig_transform
from galileoQC.transforms.conform import conform
from galileoQC.transforms.conform import conform_from_file
