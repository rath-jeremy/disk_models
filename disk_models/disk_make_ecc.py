from __future__ import division
import numpy as np; import h5py
import pandas as pd; import os
import scipy.optimize as so
import shutil

## Helper Functions ##
def E_from_M(M, e):
    func = lambda E: E - e * np.sin(E) - M
    return so.brentq(func, 0, 2 * np.pi)

def intersection(lst1, lst2): 
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return np.array(lst3)

## Convert Axisymmetric Disk to an Eccentric Disk ##
# Filename is the location of a disk h5py object.
# ecc and peri are functions of r.
# length is the halflength of a side of the box
def make_disk_eccentric(filename, ecc, per, length = 110, taper = True):
    # Make New File #
    ext = filename.rfind('.')
    output = filename[:ext] + '_eccentric' + filename[ext:]
    if os.path.isfile(output):
        os.remove(output)
    shutil.copyfile(filename, output)
    
    # Load Data #
    root = h5py.File(output, 'r+')
    data = root[root.keys()[1]]
    keys = data.keys()
    
    # Find Disk Elements #
    pos_d = np.array(data[keys[0]])
    mag_p = (pos_d[:, 0] - length)**2 + (pos_d[:, 1] - length)**2
    vel_d = np.array(data[keys[4]])
    mag_v = vel_d[:, 0]**2 + vel_d[:, 1]**2
    
    # Cut in Velocity #
    rdx_1 = np.where(mag_v < 10**(-5))[0]
    rdx_2 = np.where(mag_p > (length / 2)**2)[0]
    #rdx_2 = np.where(mag_p > 0)[0]
    redux = intersection(rdx_1, rdx_2)
    
    # Shift Positions #
    center = [length, length, 0]
    x_axis = pos_d[:, 0] - center[0]
    y_axis = pos_d[:, 1] - center[0]
    
    # Convert to Polar Coordinates #
    phi = np.arctan2(y_axis, x_axis)
    a = np.sqrt(x_axis**2 + y_axis**2)
    E = np.zeros(a.shape)
    
    # Taper e-profile #
    if taper:
        out_r = np.amin(a[redux])
        out_e = ecc(out_r)
        tap_r = out_r * (1 - out_e)
        eccen = ecc
        ecc = lambda r: (1 / np.pi * np.arctan(1 * (tap_r - r)) + 0.5) * eccen(r)
    
    # Make peri at Disk Edge 0 #
    #per = lambda r: peri(r) - peri(np.amin(a[redux]))
    
    # Convert to Orbital Elements #
    M = (phi - per(a)) % (2 * np.pi)
    e = ecc(a)
    for i in range(len(E)):
        E[i] = E_from_M(M[i], e[i])
    cf = (np.cos(E) - e) / (1 - e * np.cos(E))
    sf = (np.sqrt(1 - e**2) * np.sin(E)) / (1 - e * np.cos(E))
    cw = np.cos(per(a))
    sw = np.sin(per(a))
    ct = sf * cw + cf * sw
    st = cf * cw - sf * sw
    r = a * (1 - e**2) / (1 + e * cf)
    v = np.sqrt(2/r - 1/a)
    
    # Computer New Coordinates - Relative to Peri #
    tmp_px = r * cf
    tmp_py = r * sf
    tmp_vx = -sf / np.sqrt(a * (1 - e**2))
    tmp_vy = +(e + cf) / np.sqrt(a * (1 - e**2))
    
    # Compute New Coordinates - Absolute #
    new_px = tmp_px * cw - tmp_py * sw + center[0]
    new_py = tmp_px * sw + tmp_py * cw + center[1]
    new_vx = tmp_vx * cw - tmp_vy * sw
    new_vy = tmp_vx * sw + tmp_vy * cw
    
    # Change Exterior Grid #
    #out_r = np.amin(a[redux])
    #out_e = ecc(out_r)
    #out_w = per(out_r)
    
    new_px[redux] = pos_d[redux, 0]
    new_py[redux] = pos_d[redux, 1]
    new_vx[redux] = vel_d[redux, 0]
    new_vy[redux] = vel_d[redux, 1]
    
    # Alter Disk #
    data[keys[0]][:, 0] = new_px
    data[keys[0]][:, 1] = new_py
    data[keys[4]][:, 0] = new_vx
    data[keys[4]][:, 1] = new_vy
    root.close()
    
    return output