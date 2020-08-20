#!/usr/bin/env python

import os, sys
#from mpi4py import MPI
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky, Galactic, ICRS, match_coordinates_sky
import scutil

# $GAIA/select/50/mag17


def read_cat_fof(name):
    return np.load(name)

def main():
    
    cat_K   =   scutil.read_cat('../../../catalog.dat')
    cat_fof =   read_cat_fof('../ginfos_merge.npy')

#    print 'Total clusters in K13: %d' % (len(cat_K))

    id0 =   np.arange(len(cat_K), dtype = int)
    ids =   np.where(np.abs(cat_K['b']) < 25.)[0]
    cat_K =   cat_K[ids]
    id0 =   id0[ids]
#    print 'Total cluster in K13 with |b| < 25 degree: %d' % (len(cat_K))

    c_K     =   Galactic(l = cat_K['l'] * u.degree, b = cat_K['b'] * u.degree)
    c_fof   =   Galactic(l = cat_fof[:, 1] * u.degree, b = cat_fof[:, 2] * u.degree)

    ids_K, sep2d, _   =   match_coordinates_sky(c_fof, c_K)
#    print ids_K
#    print sep2d.degree

    r_match =   sep2d.degree

    r_K     =   cat_K['r2'][ids_K]
    id0     =   id0[ids_K]
    r_fof   =   cat_fof[:, 3]

    b0  =   r_match < r_fof
    b1  =   r_match < r_K

    ids =   np.where(np.logical_and(b0, b1))[0]
#    ids =   np.where(np.logical_not(np.logical_and(b0, b1)))[0]
#    ids =   np.arange(len(b0))

    cls_sc  =   np.loadtxt('../sc_info.txt')[:, -1].astype(int)

#    cls_sc  =   np.zeros(len(cat_fof), dtype = int)
#    for icls in range(1, 5):
#        arr =   np.loadtxt('c%1d/c%1d.txt' % (icls, icls)) 
#        n   =   len(arr)
#        idc =   arr[:, 0].astype(int)
#        cls_sc[idc] =   icls
    
    print '# idx_in_fof    separation    r_fof    r_K13'
    for id in ids:
        print '%04d\t%.3f\t%.3f\t%.3f\t%04d\t%1d' % \
            (id, r_match[id], r_fof[id], r_K[id], id0[id], cls_sc[id])

    cls =   cls_sc[ids]
    l1  =   len(np.where(cls == 1)[0])
    l2  =   len(np.where(cls == 2)[0])
    l3  =   len(np.where(cls == 3)[0])
    print '# %d %d %d' % (l1, l2, l3)

if __name__ == '__main__':
    main()
