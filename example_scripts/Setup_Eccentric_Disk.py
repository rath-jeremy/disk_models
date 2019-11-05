from __future__ import division
import matplotlib.pyplot as plt
import numpy as np; import h5py
import disk_data_analysis as da
import os

## Filenames ##
#name = "disk_qb1.00_eb0.10_alpha0.01_h0.10.hdf5"
#name = "disk_qb0.00_eb0.10_alpha0.01_h0.10.hdf5"
name = "disk_qb0.00_eb0.00_alpha0.01_h0.01.hdf5"

## Use Absolute Path ##
name = os.path.join(os.path.abspath(os.path.dirname(__file__)), name)

## Disk Profiles ##
def e_profile(r, out = 50):
    #return r**2 * np.e**(-r)
    return 0.2 * r/r

def w_profile(r):
    return np.pi / 2 * r/r

def a_from_r_v(pos, vel):
    r1 = np.sqrt(pos[:, 0]**2 + pos[:, 1]**2 + pos[:, 2]**2)
    v2 = vel[:, 0]**2 + vel[:, 1]**2 + vel[:, 2]**2
    return (2 / r1 - v2)**(-1)

## Eccentricify ##
import disk_models.disk_make_ecc as dme
outp = dme.make_disk_eccentric(name, e_profile, w_profile)

## Look at New Disk ##
oldd = h5py.File(name, 'r')
newd = h5py.File(outp, 'r')
oldp = oldd[oldd.keys()[1]]['Coordinates'] - np.array([110, 110, 0])
newp = newd[newd.keys()[1]]['Coordinates'] - np.array([110, 110, 0])
oldv = oldd[oldd.keys()[1]]['Velocities']
newv = newd[newd.keys()[1]]['Velocities']
oldm = oldd[oldd.keys()[1]]['Masses']
newm = newd[newd.keys()[1]]['Masses']

## Computer Eccentricity Vectors ##
olde = da.circumbinary.disk_orbital_elements.compute_disk_eccentriciy_vector(oldp, oldv)#, (oldm[:])**(-1))
newe = da.circumbinary.disk_orbital_elements.compute_disk_eccentriciy_vector(newp, newv)#, (newm[:])**(-1))

## Plot Data ##
# Plot Disk Positions #
plt.figure(1)
plt.subplot(121); plt.gca().set_aspect(1); plt.title("Old Disk")
plt.scatter(oldp[:, 0], oldp[:, 1], s = 0.1, edgecolor = None, c = 'b')
plt.subplot(122); plt.gca().set_aspect(1); plt.title("New Disk")
plt.scatter(newp[:, 0], newp[:, 1], s = 0.1, edgecolor = None, c = 'b')

# Plot Disk Eccentricity Vector #
plt.figure(2)
plt.subplot(121); plt.gca().set_aspect(1); plt.title("Old Disk")
plt.scatter(olde[:, 0], olde[:, 1], s = 0.1, edgecolor = None, c = 'b')
plt.xlim(-1, 1); plt.ylim(-1, 1)
plt.subplot(122); plt.gca().set_aspect(1); plt.title("New Disk")
plt.scatter(newe[:, 0], newe[:, 1], s = 0.1, edgecolor = None, c = 'b')
plt.xlim(-1, 1); plt.ylim(-1, 1)

# Plot Disk Eccentricity Profile #
plt.figure(3)
plt.subplot(121); plt.title("Old Disk")
plt.scatter(a_from_r_v(oldp, oldv), np.sqrt(olde[:, 0]**2 + olde[:, 1]**2), s = 0.1, edgecolor = None, c = 'b')
plt.subplot(122); plt.title("New Disk")
plt.scatter(a_from_r_v(newp, newv), np.sqrt(newe[:, 0]**2 + newe[:, 1]**2), s = 0.1, edgecolor = None, c = 'b')

plt.show()