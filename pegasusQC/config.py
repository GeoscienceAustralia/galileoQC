#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set the global constants (defaults for various functions) for pegasusQC.
This file can be edited to change plot fonts and default colormap and
default project name, for example.
Author: Mark Helm Dransfield
Created: Sat Aug 14 19:08:29 2021
License: CC BY-SA
"""

missingReal = -99999.99
missingInt = -999
missingDate = ""

# Set the Whizz version number and default project name.
groupName = 'Whizz Version 1.0'
projectName = 'Unknown Project'

# Just ensuring that the plot titles and labels have a uniform look that is clean and clear.
import colorcet as cc
from matplotlib import rc
rc('figure', max_open_warning=0)
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

qc_colormap = "cet_linear_bgyw_20_98_c66" # cc.m_CET_L9


