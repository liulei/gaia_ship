#!/usr/bin/env python

# simple script to plot the iso fitting result. 

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

_, diso =   np.load('Z_python3.npy', allow_pickle=True)

def plot(g, b_r, name_fit_iso):

    plt.clf()
    fig =   plt.figure()
# isochrone
    ax  =   fig.add_subplot(111)

#    fit     =   np.load(name_fit_iso, allow_pickle=True, encoding='latin1').flat[0]
    fit     =   np.load(name_fit_iso, allow_pickle=True).flat[0]
    idx_age =   fit['idx_age']
    idx_Z   =   fit['idx_Z'] 
    age     =   fit['age']
    Z       =   fit['Z']
    dg  =   fit['shift'][0]
    db_r=   fit['shift'][1]

#    gmax    =   np.max(g)
#    gmin    =   np.min(g)
#    cmax    =   np.max(b_r)
#    cmin    =   np.min(b_r)

    ms  =   2
    plt.plot(b_r, g, '.', c = 'steelblue', ms = ms)

    iso =   diso[idx_Z]['iso'][idx_age]
    b_r0    =   iso['b'] - iso['r'] - db_r
    g0      =   iso['g'] - dg
    plt.plot(b_r0, g0, 'k.', ms = ms)
#    plt.xlim(cmin - 0.3, cmax + 0.3)
#    plt.ylim(18.9, gmin - 2.0)

    plt.xlabel('$G_{\\rm BP}-G_{\\rm RP}$ [mag]')
    plt.ylabel('$G$ [mag]')

    plt.show()
