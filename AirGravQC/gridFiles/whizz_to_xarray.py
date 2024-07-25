import numpy as np
from pathlib import Path
import xarray as xr
import h5py

import AirGravQC.whizzFiles.retrieveData as rd
import AirGravQC.config as config

groupName = config.groupName
projectName = config.projectName


def whizz_to_xarray(whizz_file, z_chan, n_chan='', e_chan='', lines=[], remove_mean=False, diff_one=False, skipcontrols=False, controls=[]):
    """
    Return a point-located xArray Dataset of (northing, easting, z), over the `fiducials` dimension,
    from a whizz_file.

    Parameters
    ----------
    whizzFile : Path or String
        The Path to, or String name of, the whizz file in HDF5 format.
    z_chan : String
        The name of the channel in `whizz_file` to be imaged.
    n_chan : String, optional
        The name of the channel in `whizz_file` containing the northings (y).
        Default "" causes the name of the "YChannel" in `whizz_file` to be used.
    e_chan : String, optional
        The name of the channel in `whizz_file` containing the eastings (x).
        Default "" causes the name of the "XChannel" in `whizz_file` to be used.
    lines : String list, optional
        The line numbers to be checked. Default is all lines in the whizzFile.
    remove_mean : Bool, optional
        If True, the mean is subtracted from each survey line of data before
        writing to `my_dataset`. Default False.
    diff_one : Bool, optional
        If True, the first difference along each survey line of data is
        written to `my_dataset`. Default False.
    controls : String list, optional
        The line numbers to be skipped if skipcontrols=True. Default is empty.
    skipcontrols : Bool, optional
        If True, control lines are excluded. Default False.

    Returns
    -------
    my_dataset : xArray Dataset
        Contains the data. If `whizz_to_xarray()` is unable to return data,
        it returns an empty xArray Dataset (test with `len(aa.attrs) == 0`).

    """
    filename = str(whizz_file)
    
    with h5py.File(filename, 'r') as f:
        lines_group = f[groupName]['Lines']
        projName = f[groupName].attrs['ProjectName']
        if e_chan == '':
            e_chan = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if n_chan == '':
            n_chan = f[groupName]['CoordinateFrame'].attrs['YChannel']

        if not remove_mean and not diff_one and (z_chan == e_chan or z_chan == n_chan):
            print(f'Cannot process {z_chan}, same as {e_chan} or {n_chan}.')
            return xr.Dataset()
        totalNumFids = 0

        if lines == []:
            lines = list(lines_group.keys())

        if skipcontrols and len(controls) > 0:
            for ctrl in controls:
                lines.remove(ctrl)

        for line in lines:
            xData = rd.getLineData(lines_group[line], e_chan)
            totalNumFids += len(xData)
        print(f'{len(lines)} lines; total number of fids in whizz file = {totalNumFids}.')

        # initialise an xarray to take the data, with name and units
        fiducials = np.arange(0, totalNumFids)
        data = np.zeros((totalNumFids,))
        eastings = np.zeros((totalNumFids,))
        northings = np.zeros((totalNumFids,))

        if remove_mean and diff_one:
            xr_zchan = f'D1_MR_{z_chan}'
        elif remove_mean:
            xr_zchan = f'MR_{z_chan}'
        elif diff_one:
            xr_zchan = f'D1_{z_chan}'
        else:
            xr_zchan = f'{z_chan}'

        my_dataset = xr.Dataset({
            e_chan: xr.DataArray(
                data = eastings,
                coords={'fiducials': fiducials}, 
                dims = ['fiducials'],
                attrs = {
                    'units': 'm'
                }
            ),
            n_chan: xr.DataArray(
                data = northings,
                coords={'fiducials': fiducials}, 
                dims = ['fiducials'],
                attrs = {
                    'units': 'm'
                }
            ),
            xr_zchan: xr.DataArray(
                data = np.zeros((totalNumFids,)),
                coords={'fiducials': fiducials}, 
                dims = ['fiducials'],
                attrs = {
                    'units': '-'
                }
            )},
            attrs = {
            'author': 'Mark Dransfield',
            'x_channel': e_chan,
            'y_channel': n_chan,
            'z_channel': xr_zchan
            }
        )

        sfid = 0
        efid = 0
        for line in lines:
            sfid = efid
            xData = rd.getLineData(lines_group[line], e_chan)
            yData = rd.getLineData(lines_group[line], n_chan)
            zData = rd.getLineData(lines_group[line], z_chan)
            efid += len(yData)
            if remove_mean:
                zData = zData - np.mean(zData)
            if diff_one:
                zData = np.append(np.diff(zData), zData[-1]-zData[-2])

            my_dataset[e_chan][sfid:efid] = xData
            my_dataset[n_chan][sfid:efid] = yData
            my_dataset[xr_zchan][sfid:efid] = zData

            my_dataset[e_chan].attrs['units'] = rd.getChannelAttrs(lines_group[line], e_chan)
            my_dataset[n_chan].attrs['units'] = rd.getChannelAttrs(lines_group[line], n_chan)
            my_dataset[xr_zchan].attrs['units'] = rd.getChannelAttrs(lines_group[line], z_chan)
            if remove_mean and diff_one:
                my_dataset.attrs['title'] = f'{z_chan} (mr) (d1)'
            elif remove_mean:
                my_dataset.attrs['title'] = f'{z_chan} (mr)'
            elif diff_one:
                my_dataset.attrs['title'] = f'{z_chan} (d1)'
            else:
                my_dataset.attrs['title'] = z_chan
            
        print(f'    {my_dataset.attrs["title"]}: min = {my_dataset[xr_zchan].data.min()}, max = {my_dataset[xr_zchan].data.max()}.')
    return my_dataset

