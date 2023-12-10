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
from AirGravQC.whizzFiles.pointfiles import (xyzToHDF, updateCoordFrame, updateProject, 
    reportSampling, reportFlights, reportWhizz, updateChannelAttributes, updateLineAttributes,
    asegToHDF)
from AirGravQC.whizzPlots.whizzPlot import *
from AirGravQC.gridFiles.gridfiles import (image_pygmt, graphicsShaded, grid_n_image, diff_n_image,
    display_grid)
from AirGravQC.gridFiles.read_ers import *
from AirGravQC.qc.qualityAnalysis import *
from AirGravQC.qc.qcagg import *
from AirGravQC.qc.qcftg import *
from AirGravQC.qc.qcgrv import *
from AirGravQC.qc.qcmag import *
from AirGravQC.qc.checkClearance import (checkSafeClearance, checkClearance, checkDrape)
from AirGravQC.qc.checkIntersection import checkIntersection
from AirGravQC.qc.checkSpeeds import checkSpeeds
from AirGravQC.qc.checkOverlaps import checkOverlaps
from AirGravQC.qc.checkXYPlan import checkXYPlan
from AirGravQC.qc.checkVertPlan import checkVertPlan
from AirGravQC.qc.checkLineLengths import checkLineLengths
from AirGravQC.qc.checkGNSS import checkGNSS
from AirGravQC.qc.checkLineLengths import checkLineLengths
# from AirGravQC.qc.qcnav import *
from AirGravQC.graphics import *

