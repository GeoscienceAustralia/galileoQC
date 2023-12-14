import numpy as np
import h5py
import matplotlib.pyplot as plt
# from scipy.signal import butter, lfilter

import AirGravQC.config as config
import AirGravQC.whizzFiles.pointfiles as gw
import AirGravQC.utility.utility as util

groupName = config.groupName
                            

def checkTCDiff4(whizzFile, TCDiff4='', rawMag='', limit = 0.02, nSamples = 3000, plotAll = False):
    """
    Checks the total magnetic field fourth difference channel in a whizzFile
    against the specification that the peak to peak variation over a set
    number of samples must not exceed some peak value. If `TCDiff4` is not
    available, then the rawMag channel may be used.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    TCDiff4 : String, optional
        The name of the channel in whizzFile containing the 4th difference mag data.
        Default is '', in which case rawMag is used.
    rawMag : String, optional
        The name of the channel in whizzFile containing the raw magnetic data data.
        Default is ''; if both rawMag and TCDiff4 are '', then an error is reported.
    limit : Float, optional
        The maximum allowed peak to peak variation. Default = 0.02
    nSamples : Integer, optional
        The number of samples (moving window) over which the test is applied.
        Default = 3000
    plotAll : Bool, optional
        If True, all plots are generated.

    Returns
    -------
    None

    """
    
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        _groupMagDiff4(g, TCDiff4=TCDiff4, rawMag=rawMag, limit = limit, nSamples = nSamples, plotAll = plotAll)
        
        
def _groupMagDiff4(group, TCDiff4='', rawMag='', limit = 0.02, nSamples = 3000, plotAll = False):
    """
    Checks the total magnetic field fourth difference channel in an HDF Group
    against the specification that the peak to peak variation over a set
    number of samples must not exceed some peak value. If `TCDiff4` is not
    available, then the rawMag channel may be used.

    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    TCDiff4 : String, optional
        The name of the channel in whizzFile containing the 4th difference mag data.
        Default is '', in which case rawMag is used.
    rawMag : String, optional
        The name of the channel in whizzFile containing the raw magnetic data data.
        Default is ''; if both rawMag and TCDiff4 are '', then an error is reported.
    limit : Float, optional
        The maximum allowed peak to peak variation. Default = 0.02
    nSamples : Integer, optional
        The number of samples (moving window) over which the test is applied.
        Default = 3000
    plotAll : Bool, optional
        If True, all plots are generated.

    Returns
    -------
    None

    """
    # check 4th difference noise.
    for line in group.keys():
        if TCDiff4 == '':
            if rawMag == '':
                print('ERROR - no rawmag or 4th difference channel name supplied.')
            else:
                mag = gw.getLineData(group[line], rawMag) # np.array(group[line][rawMag])
                md4 = np.diff(mag, n=4)
                data = np.append(np.append(md4[0:2],md4),md4[-3:-1])
        else:
            data = gw.getLineData(group[line], TCDiff4)
        rangeTooHigh = False
        plotTitle = line + ' ' + TCDiff4 + ' Range'
        
        rangeTooHigh, numExceedances = util._failsDeviation(data, limit, nSamples)
        
        if rangeTooHigh:
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            ax.plot(data)
            plt.title(plotTitle)
            plt.grid(True)
            plt.show()
            print(line, ' ', numExceedances, ' > ', nSamples, ' - FAIL')
          
        elif numExceedances > 0:
            print(line, ' ', numExceedances, ' < ', nSamples, ' - PASS')
            if plotAll:
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                ax.plot(data)
                plt.title(plotTitle)
                plt.grid(True)
                plt.show()
            
        else:
            print(line, ' - PASS')


