#=========================
#
# MAG
#
#=========================

import numpy as np
import h5py
import matplotlib.pyplot as plt
# from scipy.signal import butter, lfilter

import AirGravQC.config as config
import AirGravQC.whizzFiles.pointfiles as gw
import AirGravQC.utility.utility as util

groupName = config.groupName

def checkDiurnal(whizzFile, basemag, rangeLimit = 5.0, nSamples = 3000, diff4Limit = 0.5):
    # TODO: add check for singleValueExceedance()
    
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        _reportAllDiurnal(g, basemag, rangeLimit=rangeLimit, nSamples=nSamples, diff4Limit=diff4Limit)
                    

def checkLineDiurnal(whizzFile, line, basemag, rangeLimit = 5.0, nSamples = 3000, diff4Limit = 0.5):
    # TODO: add check for singleValueExceedance()
    # NOT USED??
    
    filename = str(whizzFile)
    with h5py.File(filename, 'r') as f:
        g = f[groupName]

        _reportLineDiurnal(g, line, basemag, rangeLimit, nSamples, diff4Limit)


def _reportAllDiurnal(group, basemag, rangeLimit = 5.0, nSamples = 3000, diff4Limit = 0.5):
    # check 4th difference noise.
    for line in group.keys():
        _reportLineDiurnal(group, line, basemag, rangeLimit, nSamples, diff4Limit)
        
        
def _reportLineDiurnal(group, line, basemag, rangeLimit = 5.0, nSamples = 3000, diff4Limit = 0.5):
    # get data, check we have sufficient.
    diurnalExceeded = False
    failedSample = 0
    bigExtremum = 0.0
    data = gw.getLineData(group[line], basemag)
    data = data[np.logical_not(np.isnan(data))]
        
    if nSamples > len(data):
        print(f'\n  Short line: {len(data)} < {nSamples}.')
        nSamples = len(data)
    if (nSamples % 2) == 0:
        nSamples = nSamples - 1
    nSam = (nSamples - 1) // 2

    for ii in range(nSam, len(data)-nSam):
        localData = data[ii-nSam:ii+nSam]
        localSlope = (localData[-1] - localData[0]) / nSamples
        # deviation = np.zeros(localData.shape)
        deviation = localData - localSlope * range(0, len(localData)) - localData[0]
        # for jj in range(0, len(localData)):
        #     deviation[jj] = localData[jj] - localData[0] - localSlope * jj
        extremum = np.max(deviation) if np.max(deviation) > -np.min(deviation) else -np.min(deviation)
        if extremum > rangeLimit:
            diurnalExceeded = True
            if extremum > bigExtremum:
                bigExtremum = extremum
                failedSample = ii
            
    if diurnalExceeded:
        print(f'\n  Diurnal for {basemag} at sample number {failedSample} diverges from chord by {bigExtremum:.2f}, exceeding {rangeLimit:.1f} - FAIL')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(data)
        plotTitle = f'Line {line} Channel {basemag}: reaches {bigExtremum:.2f} at {failedSample}, exceeding {rangeLimit} - FAIL'
        plt.title(plotTitle, fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(6)
        for label in ax.get_yticklabels(): label.set_fontsize(6)
        fig.tight_layout()
        plt.show()


def diurnal(filename, lines = [], name = 'Basemag', rangeLimit = 5.0, nSamples = 3000, showPlot=False):
    """
    QC Standard: May not exceed rangeLimit nT over nSamples.

    From each individual observation:
    
    Examine nSamples / 2 into the past and nSamples / 2 into the future.
    Interpolate linearly.
    Measure the deviation from the interpolated line.
    The deviation may not be larger than rangeLimit.

    Parameters
    ----------
    name : String
        Name of the channel containg 'basemag' data.
    basemag : Numpy array
        the basemag data to be checked.
    rangeLimit : Float, optional
        The largest allowable deviation from a straight chord (usu nT). The default is 5.0.
    nSamples : Int, optional
        The number of samples of 'basemag' for the windowing for the QC check 
        (at 10 Hz, 3000 samples = 5 min) . The default is 3000.
    showPlot : Bool, optional
        if showPlot and ~ok: plot the results.

    Returns
    -------
    ok : Bool
        True if 'basemag' passed the test.
    report : String
        "" if ok == True, else an error description.

    """
    finalStatus = True
    report = ''
    with h5py.File(filename, 'r') as f:
        g = f[groupName]
        if lines == []:
            lines = g.keys()
        
        for line in lines:
            basemag = gw.getLineData(g[line], name)
            lineStatus = True
            if nSamples > len(basemag):
                nSamples = len(basemag)
            if (nSamples % 2) == 0:
                nSamples = nSamples - 1
            nSam = (nSamples - 1) // 2
            for ii in range(nSam, len(basemag)-nSam):
                localData = basemag[ii-nSam:ii+nSam]
                localSlope = (localData[-1] - localData[0]) / nSamples
                # deviation = np.zeros(localData.shape)
                deviation = localData - localSlope * range(0, len(localData)) - localData[0]
                # for jj in range(0, len(localData)):
                #     deviation[jj] = localData[jj] - localData[0] - localSlope * jj
                extremum = np.max(deviation) if np.max(deviation) > -np.min(deviation) else -np.min(deviation)
                if extremum > rangeLimit:
                    lineStatus = False
                    finalStatus = False
                    print(f'\n  Line {line} Diurnal for {name} reaches {extremum:.1f}, exceeding {rangeLimit:.1f} - FAIL')
                    report += f'\n  Line {line} Diurnal for {name} reaches {extremum:.1f}, exceeding {rangeLimit:.1f} - FAIL'
                    break
                
            print(lineStatus)   
            if lineStatus == False and showPlot == True:
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                ax.plot(basemag, lw=0.5)
                plotTitle = f'Line {line} Channel {name}: reaches {extremum:.1f} at {ii}, exceeding {rangeLimit:.1f} - FAIL'
                plt.title(plotTitle, fontsize = 8)
                plt.grid(True)
                for label in ax.get_xticklabels(): label.set_fontsize(6)
                for label in ax.get_yticklabels(): label.set_fontsize(6)
                fig.tight_layout()
                plt.show()
                
    return finalStatus, report
            