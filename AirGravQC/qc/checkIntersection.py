# import necessary modules
import numpy as np
import h5py
import matplotlib.pyplot as plt

import AirGravQC.config as config

groupName = config.groupName
        

def checkIntersection(whizzFile, controls=[], travs=[], xChannel='', yChannel='', zChannel='', max_allowed_deltaZ=10.0, plot_flag=False):
    """
    Checks that the values of the data in `zChannel` at the intersection of traverse and control
    lines are different by no more than the maximum allowed delta.
    
    Parameters
    ----------
    whizzFile : HDF5 Whizz file pathlib Path
        The pathlib Path to the Whizz HDF5 file containing the survey line data.
    controls : [String]
        A list of control flightlines, e.g. ['1000110.0', '1000210.0', '1000310.0'].
    travs : [String], optional
        A list of traverse flightlines, e.g. ['1000110.0', '1000210.0', '1000310.0']. Defaults
        to all flight lines in the `whizzFile`.
    xChannel : String, optional
        The name of the geoWhizz field or channel containing the measured x positions. The
        default is to read the xChannel field name from the Coordinate Frame.
    yChannel : String, optional
        The name of the geoWhizz field or channel containing the measured y positions. The
        default is to read the yChannel field name from the Coordinate Frame.
    zChannel : String, optional
        The name of the geoWhizz field or channel containing the data to be tested. The
        default is to read the AltitudeChannel field name from the Coordinate Frame.
    max_allowed_deltaZ : Float, optional
        The maximum allowed difference in `zChannel` between the traverse and control lines
        at each intersection point. Defaults to 10.0.
    plot_flag : Bool, optional
        If True, plot a map of each pair of intersecting lines where the `zChannel`
        values differ by more than `max_allowed_deltaZ`. Default False.

    Returns
    -------
    None

    """
    data_is_good = True
    report = ''
    filename = str(whizzFile)
    num_intersections_checked = 0
    num_failed_intersections = 0
    
    with h5py.File(filename, 'r') as f:
        g = f[groupName]['Lines']
        if xChannel == '':
            xChannel = f[groupName]['CoordinateFrame'].attrs['XChannel']
        if yChannel == '':
            yChannel = f[groupName]['CoordinateFrame'].attrs['YChannel']
        if zChannel == '':
            zChannel = f[groupName]['CoordinateFrame'].attrs['AltitudeChannel']
        all_lines = list(g.keys())
        alltravs, allcontrols = controls_lessthan_1000(all_lines)
        if controls == []:
            controls = allcontrols
        if travs == []:
            lines = alltravs
        else:
            lines = travs
        #if the channel has an attribute 'Units'
        dd = g[lines[0]][zChannel]
        z_units = ''
        if 'Units' in dd.attrs.keys():
            z_units = dd.attrs['Units']

        for linec in controls:
            x_ctrl = np.array(g[linec][xChannel])
            y_ctrl = np.array(g[linec][yChannel])
            z_ctrl = np.array(g[linec][zChannel])

            bear_ctrl = _calc_bearing(x_ctrl, y_ctrl)
            (y_ctrl1, x_ctrl1) = _rotateCoords(x_ctrl-x_ctrl[0], y_ctrl-y_ctrl[0], -bear_ctrl)
            for linet in lines:
                # if linet == linec, then it is a control line, not a traverse.
                # TODO: compare the PlannedLine for linet and linec rather than the lines themselves,
                #       since we don't want to compare different segments of the same line!
                if linet == linec:
                    continue
                x_trav = np.array(g[linet][xChannel])
                y_trav = np.array(g[linet][yChannel])
                z_trav = np.array(g[linet][zChannel])
                (y_trav1, x_trav1) = _rotateCoords(x_trav-x_ctrl[0], y_trav-y_ctrl[0], -bear_ctrl)
                # if _lines_cross(x_ctrl, y_ctrl, x_trav, y_trav):
                if _intersect(x_ctrl[0], y_ctrl[0], x_ctrl[-1], y_ctrl[-1], x_trav[0], y_trav[0], x_trav[-1], y_trav[-1]):
                    # print(f'bearings: {bearingt}, {bearingt} -- {np.abs(np.cos(bearingc - bearingt))}')
                    dh = _intersection_height(x_trav1, y_trav1, z_trav, x_ctrl1, y_ctrl1, z_ctrl, bear_ctrl)
                    num_intersections_checked += 1
                    if dh > max_allowed_deltaZ:
                        num_failed_intersections += 1
                        # print(dh)
                        bc_deg = bear_ctrl * 180.0 / np.pi
                        report += f'\n  {linet} : {linec} [bearing={bc_deg:.1f}] intersection {zChannel} difference = {dh:.1f} > {max_allowed_deltaZ:.1f}'
                        if z_units != '':
                            report += ' ' + z_units
                        report += '.'
                        data_is_good = False
                        if plot_flag:
                            fig = plt.figure()
                            fig.suptitle(f'Title {linet}', fontsize=10)
                            fig.subplots_adjust(top=0.85)
                            
                            ax = fig.add_subplot(1,2,1)
                            ax.plot(x_trav, y_trav, x_ctrl, y_ctrl)
                            plt.ylabel('y_trav', fontsize = 6)
                            plt.grid(True)
                            for label in ax.get_xticklabels(): label.set_fontsize(6)
                            for label in ax.get_yticklabels(): label.set_fontsize(6)

                            ax = fig.add_subplot(1,2,2)
                            ax.plot(x_trav1, y_trav1, x_ctrl1, y_ctrl1)
                            plt.ylabel('y_ctrl', fontsize = 6)
                            plt.grid(True)
                            for label in ax.get_xticklabels(): label.set_fontsize(6)
                            for label in ax.get_yticklabels(): label.set_fontsize(6)
                            plt.show()
                # else:
                #     report += f'\n  {linet} : {linec} un-tested since not perpendicular.'
    if data_is_good:
        report += f'All {num_intersections_checked} intersection {zChannel} differences were less than {max_allowed_deltaZ:.1f}'
        if z_units != '':
            report += ' ' + z_units
        report += '.'
    else:
        tmpstr = f'Of {num_intersections_checked} intersections checked'
        tmpstr += f', {num_failed_intersections} exceeded the '
        tmpstr += f'{max_allowed_deltaZ} allowed {zChannel} difference.\n'
        report = tmpstr + report
    print(report)
    return


def _ccw(x1, y1, x2, y2, x3, y3):
    return (y3-y1)*(x2-x1) > (y2-y1)*(x3-x1)


def _intersect(cx1, cy1, cx2, cy2, tx1, ty1, tx2, ty2):
    return _ccw(cx1, cy1, tx1, ty1, tx2, ty2) != _ccw(cx2, cy2, tx1, ty1, tx2, ty2) and _ccw(cx1, cy1, cx2, cy2, tx1, ty1) != _ccw(cx1, cy1,cx2, cy2, tx2, ty2)


def _intersection_height(x_trav, y_trav, z_trav, x_ctrl, y_ctrl, z_ctrl, bearingc):
    """
    """

    y = np.abs(y_trav - _mean_1std(y_ctrl))
    it_arr = np.where(y == y.min())
    it = it_arr[0]
    if len(it) > 1:
        print('\nBug in _intersection_height - this intersection height cant be calculated')
        print(it, type(it), len(it))
        return 0.0

    x = np.abs(x_ctrl - x_trav[it])
    ic_arr = np.where(x == x.min())
    ic = ic_arr[0]

    # print(it, z_trav[it], ic, z_ctrl[ic])
    return np.abs(z_trav[it] - z_ctrl[ic])[0]


def _lines_cross(x_ctrl, y_ctrl, x_trav, y_trav):
    """
    Assumes x,y coordinates rotated so that y_ctrl is approximately constant (y-axis
    approximately parallel to the control line).
    """

    # The traverse is 'north' or 'south' of the control line
    if y_trav.min() > y_ctrl.max() or y_trav.max() < y_ctrl.min():
        return False

    # We don't allow shallow angle crossovers (and won't count lines that were supposed to be parallel!).
    min_cosine = 0.1
    bear_ctrl = _calc_bearing(x_ctrl, y_ctrl)
    bear_trav = _calc_bearing(x_trav, y_trav)
    if np.abs(np.cos(bear_ctrl - bear_trav)) < min_cosine:
        return False

    # Find the index to the traverse sample whose y coordinate is closest to the mean control y coordinate
    y = np.abs(y_trav - _mean_1std(y_ctrl))
    it_arr = np.where(y == y.min())
    if np.size(it_arr):
        itc = it_arr[0]
    else:
        return False

    # Now can check if the traverse is 'east' or 'west' of the control line.
    if x_ctrl.min() > x_trav[itc] or x_ctrl.max() < x_trav[itc]:
        return False

    return True


def _calc_bearing(x, y):
    """
    arctan(mean(diff(x) / mean(diff(y))))
    """
    return np.arctan2(_mean_1std(np.diff(x)), _mean_1std(np.diff(y)))


def _mean_1std(x):
    """
    Calculate the mean of the values in x that fall in the range of +/- stdev(x).
    This is a simplistic "mean excluding outliers".
    """
    mean1 = np.mean(x)
    std1 = np.std(x)
    idx = np.argwhere(np.logical_and(x > (mean1 - std1), x < (mean1 + std1)))
    return np.mean(x[idx])


def controls_lessthan_1000(all_lines):
    ctrl_strs = []
    trav_strs = []
    for line in all_lines:
        if float(line) < 999.99:
            ctrl_strs.append(line)
        else:
            trav_strs.append(line)
    return trav_strs, ctrl_strs


