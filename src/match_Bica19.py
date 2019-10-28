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
    
    cat   =   scutil.read_B19('../../../Bica.txt')
    cat_fof =   read_cat_fof('../ginfos_merge.npy')

#    print 'Total clusters in K13: %d' % (len(cat_K))

    id0 =   np.arange(len(cat), dtype = int)
    ids =   np.where(np.abs(cat['b']) < 25.)[0]
    cat =   cat[ids]
    id0 =   id0[ids]
    print ('# Total cluster in B19 with |b| < 25 degree: %d' % (len(cat)))

    c     =   Galactic(l = cat['l'] * u.degree, b = cat['b'] * u.degree)
    c_fof   =   Galactic(l = cat_fof[:, 1] * u.degree, b = cat_fof[:, 2] * u.degree)

    ids, sep2d, _   =   match_coordinates_sky(c_fof, c)
#    print ids_K
#    print sep2d.degree

    r_match =   sep2d.degree

    r       =   cat['r'][ids]
    id0     =   id0[ids]
    r_fof   =   cat_fof[:, 3]

    b0  =   r_match < r_fof
    b1  =   r_match < r

    ids =   np.where(np.logical_and(b0, b1))[0]
#    ids =   np.where(np.logical_not(np.logical_and(b0, b1)))[0]

#    ids =   np.where(np.logical_or(b0, b1))[0]
#    ids =   np.where(np.logical_not(np.logical_or(b0, b1)))[0]

    cls_sc  =   np.loadtxt('../sc_info.txt')[:, -1].astype(int)

    print ('# idx_in_fof    separation    r_fof    r_B19')
    for id in ids:
        print ('%04d\t%.3f\t%.3f\t%.3f\t%04d\t%1d' % \
            (id, r_match[id], r_fof[id], r[id], id0[id], cls_sc[id]))

    cls =   cls_sc[ids]
    l1  =   len(np.where(cls == 1)[0])
    l2  =   len(np.where(cls == 2)[0])
    l3  =   len(np.where(cls == 3)[0])
    print ('# %d %d %d' % (l1, l2, l3))

if __name__ == '__main__':
    main()
