from __future__ import division
import matplotlib.pyplot as plt
import numpy as np; import h5py
import disk_data_analysis as da
import os

## Filenames ##
#name = "disk_qb0.00_eb0.00_alpha0.10_h0.01.hdf5"
name = "disk_qb0.00_eb0.00_alpha0.10_h0.01_eccentric.hdf5"

## Load Data ##
name = os.path.join(os.path.abspath(os.path.dirname(__file__)), name)
file = h5py.File(name, 'r+')
Pids = file['PartType0']['ParticleIDs'][...]
pind = np.where(Pids >= 0)[0]
pids = Pids[pind]
nind = np.where(Pids < 0)[0]
nids = Pids[nind]

## Fix Negative Indices ##
if len(nids) != len(set(nids)):
    print "Fixing negative indices:"
    s_nind = np.sort(nind)
    print s_nind
    file['PartType0']['ParticleIDs'][s_nind
    ] = np.arange(len(nind)) + 1
    print "Done."

## Fix Positive Indices #
if len(pids) != len(set(pids)):
    print "Fixing positive indices:"
    s_pind = np.sort(pind)
    file['PartType0']['ParticleIDs'][s_pind] = np.arange(len(pind)) + 1
    print "Done."

## Save ##
#print len(file['PartType0']['ParticleIDs'][...]), len(set(file['PartType0']['ParticleIDs'][...]))
file.close()