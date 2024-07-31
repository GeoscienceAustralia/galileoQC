import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
from matplotlib.ticker import StrMethodFormatter

import AirGravQC.config as config
import AirGravQC.utility.utility as util

groupName = config.groupName


def checkLineSeparation(whizzfile, nominalsep, lines=[], measX='', measY='', allowance=200.0, maxDistance=0, maxCounter=0, known='', plot_flag=False, verbose=False):
    """
    Reports exceedances of survey line separation for an airborne survey Whizz database.
    
    For each survey line, its two neighbours are found and the perpendicular separation vector
    calculated. If the absolute difference between this separation and the nominal separation
    exceeds the allowance for maxDistance consecutive metres, or maxCounter consecutive
    fiducials, then an out-of-specification exceedance is reported for that line.

    Parameters
    ----------
    whizzfile : String or pathlib Path
        Name of a HDF5 Whizz file, including path and extension, with the survey
        measured data.
    nominalsep : Float
        The nominal (planned) line separation in metres.
    lines : String list, optional.
        The line numbers to be checked. Default is all lines in the whizzFile.
    measX : String, optional
        The name of the geoWhizz field or channel containing the measured x position. The
        default is to read the xChannel field name from the Coordinate Frame.
    measY : String, optional
        The name of the geoWhizz field or channel containing the measured y position. The
        default is to read the yChannel field name from the Coordinate Frame.
    allowance : Float, optional
        The allowed horizontal separation of a survey line from its neighbours in metres.
        The default is 200.0.
    maxDistance : Float, optional
        The maximum number of consecutive metres for which an exceedance
        greater than allowance is permitted. If 0, then the constraint is
        ignored. The default is 1000.0.
    maxFiducials : Int, optional
        The maximum number of consecutive fids for which an exceedance
        greater than allowance is permitted. If 0, then the constraint is
        ignored. The default is 0.
    known : String, optional
        If present, the name of the channel containing the "known error" flag.
        This is reported against any error so that known errors can be distinguished
        from unknown errors.
    plot_flag : Bool, optional
        If True, plot exceedances for each failed line. Default False.
    verbose : Bool, optional
        If True, verbose reporting is given which is annoying if there are many errors.
        Default False.

    Returns
    -------
    None.

    """
    # start_x = 0.0
    # end_x = 0.0
    # start_y = 0.0
    # end_y = 0.0

    if True:
        print('checkLineSeparation() is not yet available. It awaits an algorithm that is acceptably fast.')
        return
    else:
        # exceedances_known = False
        # this_exc_known = False
        # number_known = 0

        # measFile = str(whizzfile)
        

        # with h5py.File(measFile, 'r+') as fm:
        #     linesgroup = fm[groupName]['Lines']
        #     if measX == '':
        #         measX = fm[groupName]['CoordinateFrame'].attrs['XChannel']
        #     if measY == '':
        #         measY = fm[groupName]['CoordinateFrame'].attrs['YChannel']
        #     numLines = len(linesgroup.items())

        #     message = ''
        #     num_lines_exceeded = 0
        #     total_num_excs = 0
        #     num_lines_unplanned = 0

        #     if lines == []:
        #         lines = linesgroup.keys()

        #     for line in lines:
        #         adjacentlines = findadjacentlines(line, linesgroup)
        #         for aline in adjacentlines:
        #             separation = calcseparation(line, aline)
        #             excess_sep = allowance - np.abs(separation - nominalsep)

        #             # no exceedances
        #             if np.max(excess_sep) < 0.0:
        #                 continue
        #             elif 
        #             if util._exceedance_fail(num_fids_in_exceedance, len_exceedance, maxCounter, maxDistance):
        return




