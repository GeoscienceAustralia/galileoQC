#=========================
#
# AGG
#
#=========================

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

import AirGravQC.config as config
import AirGravQC.whizzFiles.pointfiles as gw

groupName = config.groupName

def diffGravVturb(whizzFile, turbulence, aD, bD, error_spec=5.0, low_cut=0.001, measX='', measY=''):
    """
    For a Falcon AGG. For each line, reports the gD difference noise,
    stdev(aD-bD)/2, and plots this as a scatter plot against the mean turbulence.
    All lines with difference noise greater than `error_spec` are reported.
    
    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    turbulence : String
        The name of the channel containing the turbulence field.
    aD : String
        The name of the channel containing the A complement gravity data.
    bD : String
        The name of the channel containing the B complement gravity data.
    error_spec : Float, optional
        The value above which the difference noise is excessive and is
        reported. Default 5.0.
    low_cut : Float, optional
        The low frequency cut-off frequency (in Hz) for the band-pass filtering
        applied before differencing. Default 0.001 (ie 1 mHz).

    Returns
    -------
    None.

    """
    filename = str(whizzFile)

    num_lines_failed = 0
    report = 'diffGravVturb() estimates the noise in (A+B)/2 as the stdev(A-B)/2.\n'
    period = 1.0 / low_cut
    wavelength = 60.0 * period / 1000.0 # km
    report += f'The input data are band-pass filtered at [{period}, 1.0] sec or [{wavelength}, 0.06] km at 60m/s.\n'
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        numLines = len(g.items())
        turbMean = np.zeros((numLines,))
        errmean = np.zeros((numLines,))
        lineNo = np.zeros((numLines,))
        count = 0

        line_strs = list(g.keys())
        flightLine = line_strs[0]
        turb_units = g[flightLine][turbulence].attrs['Units']
        err_units = g[flightLine][aD].attrs['Units']

        if measX == '':
            measX = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if measY == '':
            measY = f[groupName]['CoordinateFrame'].attrs['YChannel']


        failed_lines = 0
        for line in g.keys():
            xM = gw.getLineData(g[line], measX)#np.array(g[line][measX])
            yM = gw.getLineData(g[line], measY)#np.array(g[line][measY])
            line_length = _displacement2(xM[0], xM[-1], yM[0], yM[-1]) / 1000.0

            turb = gw.getLineData(g[line], turbulence)#np.array(g[line][turbulence])
            A_d = gw.getLineData(g[line], aD)#np.array(g[line][aD])
            B_d = gw.getLineData(g[line], bD)#np.array(g[line][bD])
            idx = np.where(~np.isnan(A_d + B_d))
            Ad = _butter_bandpass_filter(A_d[idx], low_cut, 1.0, 8.0, order=3)
            Bd = _butter_bandpass_filter(B_d[idx], low_cut, 1.0, 8.0, order=3)
            turb_clean = turb[idx]
            err_data = (Ad - Bd)/2.0
            err_data = err_data - np.mean(err_data)

            turbMean[count] = np.mean(turb_clean)
            errmean[count] = np.std(err_data)
            
            if errmean[count] > error_spec:
                report += f'{line} (length {line_length:6.1f} km) fails with noise {errmean[count]:.1f} > {error_spec}.\n'
                num_lines_failed += 1

            count += 1

        fig = plt.figure()
        fig.suptitle(f'Noise vs Turbulence - {projName}', fontsize=12)
        fig.subplots_adjust(top=0.85)
        ax = fig.add_subplot(1,1,1)
        neplot, = ax.plot(turbMean, errmean, 'bo', label='$E_{D}$')
        for ii in range(0, len(turbMean)):
            ax.text(turbMean[ii], errmean[ii], f'{line_strs[ii]}', fontsize=8)
        plt.ylabel(f'Difference Noise [{err_units}]', fontsize = 8)
        plt.xlabel(f'Turbulence [{turb_units}]', fontsize = 8)
        plt.grid(True)
        for label in ax.get_xticklabels(): label.set_fontsize(8)
        for label in ax.get_yticklabels(): label.set_fontsize(8)
        fig.tight_layout()
        plt.show()
        report = f'{num_lines_failed} failed of {count} lines total.' + report
    print(report)
        

def diffNoiseVturb(whizzFile, turbulence, lines=[], aNE='', aUV='', bNE='', bUV='', eNE='', eUV='', error_spec=5.0, labelLines=False):
    """
    For a Falcon AGG. For each line, reports the mean NE and UV difference noise,
    and plots these as a scatter plot against the mean turbulence.
    
    Parameters
    ----------
    whizzFile : String or pathlib.PosixPath
        Name of a HDF5 Whizz file, including path and extension.
    turbulence : String
        The name of the channel containing the turbulence field.
    lines : String list, optional
        The line numbers to be checked. Default is all lines in the whizzFile.
    aNE : String
        The name of the channel containing the A_NE field. If '', then the
        eNE and eUV channels are used, otherwise, the aNE, aUV, bNE
        and bUV channels are used.
    aUV : String
        The name of the channel containing the A_UV field. If '', then the
        eNE and eUV channels are used, otherwise, the aNE, aUV, bNE
        and bUV channels are used.
    bNE : String
        The name of the channel containing the B_NE field. If '', then the
        eNE and eUV channels are used, otherwise, the aNE, aUV, bNE
        and bUV channels are used.
    bUV : String
        The name of the channel containing the B_UV field. If '', then the
        eNE and eUV channels are used, otherwise, the aNE, aUV, bNE
        and bUV channels are used.
    eNE : String, optional
        The name of the channel containing the NE difference error. If '',
        then the aNE, aUV, bNE and bUV channels are used.
    eUV : String
        The name of the channel containing the UV difference error. If '',
        then the aNE, aUV, bNE and bUV channels are used.
    error_spec : Float, optional
        The value above which the difference noise is excessive and is reported.
        Default 5.0.
    labelLines : Bool, optional
        if True, label (with the line number) all points on the plot where the
        line failed the specification. Defaults to False

    Returns
    -------
    None.

    """
    filename = str(whizzFile)
    if aNE=='' and aUV=='' and bNE=='' and bUV=='':
        need_calc = False
    elif eNE=='' and eUV=='':
        need_calc = True
    else:
        print('ERROR - must specify either all four raw channels or both error channels')
        return

    report = ''
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        numLines = len(g.items())
        turbMean = np.zeros((numLines,))
        errNEmean = np.zeros((numLines,))
        errUVmean = np.zeros((numLines,))
        lineNo = np.chararray((numLines,))
        count = 0
        flightLine = list(g.keys())[0]

        turb_units = g[flightLine][turbulence].attrs['Units']
        if need_calc:
            err_units = g[flightLine][aNE].attrs['Units']
        else:
            err_units = g[flightLine][eNE].attrs['Units']

        labelx = []
        labely = []
        labelt = []
        failed_lines = 0

        if lines == []:
            lines = g.keys()
        for line in lines:
            turb = gw.getLineData(g[line], turbulence)#np.array(g[line][turbulence])
            if need_calc:
                A_ne = gw.getLineData(g[line], aNE)#np.array(g[line][aNE])
                A_uv = gw.getLineData(g[line], aUV)#np.array(g[line][aUV])
                B_ne = gw.getLineData(g[line], bNE)#np.array(g[line][bNE])
                B_uv = gw.getLineData(g[line], bUV)#np.array(g[line][bUV])
                idx = np.where(~np.isnan(A_ne + A_uv + B_ne + B_uv))
                Ane = A_ne[idx]
                Auv = A_uv[idx]
                Bne = B_ne[idx]
                Buv = B_uv[idx]
                turb_clean = turb[idx]
    
                errNE_data = (Ane - Bne)/np.sqrt(8)
                errUV_data = (Auv - Buv)/np.sqrt(8)
            else:
                E_ne = gw.getLineData(g[line], eNE)#np.array(g[line][eNE])
                E_uv = gw.getLineData(g[line], eUV)#np.array(g[line][eUV])
                idx = np.where(~np.isnan(E_ne + E_uv))
                errNE_data = E_ne[idx]
                errUV_data = E_uv[idx]
                turb_clean = turb[idx]

            turbMean[count] = np.mean(turb_clean)
            errNEmean[count] = np.std(errNE_data)
            errUVmean[count] = np.std(errUV_data)
            avge_noise = 0.5 * (errNEmean[count] + errUVmean[count])
            
            if 0.5 * (errNEmean[count] + errUVmean[count]) > error_spec:
                report += f'{line} fails with noise {avge_noise:.2f} > {error_spec}, mean turbulence = {turbMean[count]:.2f}.\n'
                if labelLines:
                    labelx.append(turbMean[count])
                    labely.append(max(errNEmean[count], errUVmean[count]))
                    labelt.append(line)
                failed_lines += 1
            count += 1

        fig = plt.figure()
        fig.suptitle(f'Noise vs Turbulence - {projName}', fontsize=12)
        fig.subplots_adjust(top=0.85)
        ax = fig.add_subplot(1,1,1)
        ax.vlines(turbMean, errNEmean, errUVmean, 'k', lw=0.3)
        neplot, = ax.plot(turbMean, errNEmean, 'bo', label='$E_{ne}$')
        uvplot, = ax.plot(turbMean, errUVmean, 'go', label='$E_{uv}$')
        if labelLines:
            for ii in range(failed_lines):
                plt.text(labelx[ii], labely[ii], labelt[ii], va='top', ha='right', size=6.0)
        plt.ylabel(f'Difference Noise [{err_units}]', fontsize = 8)
        plt.xlabel(f'Turbulence [{turb_units}]', fontsize = 8)
        plt.grid(True)
        ax.legend(handles=[neplot, uvplot])
        for label in ax.get_xticklabels(): label.set_fontsize(8)
        for label in ax.get_yticklabels(): label.set_fontsize(8)
        fig.tight_layout()
        plt.show()
    print(report)
        
