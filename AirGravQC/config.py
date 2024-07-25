#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 14 19:08:29 2021

@author: markdransfield

Set the global constants for AirGravQC:

missingReal = -99999.99

missingInt = -999

missingDate = ""

Set the Whizz version number and default project name:

groupName = 'Whizz Version 1.0'

projectName = 'Unknown Project'
"""

missingReal = -99999.99
missingInt = -999
missingDate = ""

# Get Project Name, reference [default: ""], and rest of header info.
groupName = 'Whizz Version 1.0'
projectName = 'Unknown Project'

# Just ensuring that the plot titles and labels have a uniform look that is clean and clear.
from matplotlib import rc
rc('figure', max_open_warning=50)
rc('font',**{'family':'sans-serif'})

SMALL_SIZE = 6
MEDIUM_SIZE = 8
BIGGER_SIZE = 10

rc('font', size=SMALL_SIZE)          # controls default text sizes
rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
