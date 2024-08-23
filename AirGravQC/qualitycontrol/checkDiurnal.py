#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check magnetic diurnal data against specification.
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

import AirGravQC.config as config
import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.utility.utility as util

groupName = config.groupName


def checkDiurnal(whizzFile, basemag, lines=[], rangeLimit = 5.0, nSamples = 3000, diff4Limit = 0.5, plot_flag=False):
    # TODO: add check for singleValueExceedance()
    
    filename = str(whizzFile)
    report = ''
    num_failed_lines = 0

    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if lines == []:
            lines = g.keys()
        numLines = len(lines)
        for line in lines:
            diurnalExceeded = False
            failedSample = 0
            bigExtremum = 0.0
            data = rd.getLineData(g[line], basemag)
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
                report += f'\n  Diurnal for {basemag} at sample number {failedSample} diverges from chord by {bigExtremum:.2f}, exceeding {rangeLimit:.1f} - FAIL'
                num_failed_lines += 1
                if plot_flag:
                    fig = plt.figure()
                    ax = fig.add_subplot(1,1,1)
                    ax.plot(data)
                    plotTitle = f'Line {line} Channel {basemag}: reaches {bigExtremum:.2f} at {failedSample}, exceeding {rangeLimit} - FAIL'
                    plt.title(plotTitle, fontsize = 8)
                    plt.grid(True)
                    for label in ax.get_xticklabels(): label.set_fontsize(6)
                    for label in ax.get_yticklabels(): label.set_fontsize(6)
                    fig.tight_layout()

    print(f'  Checked {numLines} lines, {num_failed_lines} failed.\n')
    print(report)
    if plot_flag and num_failed_lines > 0:
        plt.show()


