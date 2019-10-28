#!/usr/bin/env python

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import isochrone

# $GAIA/select/50
# input: ginfos_merge.npy, fof_scXXXX.npy, mag17/fit_iso_scXXXX.npy

def classify_r_dist(n17, r_dist, r_hrd, d2):
    
    if n17 >= 50 and r_dist < 3 and r_hrd < 0.1 and d2 < 0.05:
        return 1
    elif n17 >= 50 and r_dist < 3 and r_hrd < 0.1:
        return 2
    else:
        return 3

age_cut =   5E6
#age_cut =   1E7
# age in yr
# using age as cut
def classify(n17, age, r_hrd, d2):
    
    if n17 >= 50 and r_hrd < 0.1 and age >= age_cut and d2 < 0.05:
        return 1
    elif n17 >= 50 and r_hrd < 0.1 and age >= age_cut:
        return 2
    else:
        return 3

def classify_no50(n17, age, r_hrd, d2):
    
    if r_hrd < 0.1 and age >= age_cut and d2 < 0.05:
        return 1
    elif r_hrd < 0.1 and age >= age_cut:
        return 2
    else:
        return 3

def remove_nan(arr):
    
    idb =   np.isnan(arr['mag_bp'])
    idr =   np.isnan(arr['mag_rp'])

    ids =   np.where(np.logical_not(np.logical_or(idb, idr)))[0]

    return arr[ids]

def calc_ratio(x, y):

    x   =   x - np.average(x)
    y   =   y - np.average(y)
    xx  =   np.average(x * x)
    xy  =   np.average(x * y)
    yy  =   np.average(y * y)
    m   =   np.array([[xx, xy], [xy, yy]])
    w, v    =   np.linalg.eig(m)

    w   =   np.sort(np.abs(w))
    return w[0] / w[1]

def extract_fit_iso_ratio(id_sc):

    name    =   'mag17/fit_iso_sc%04d.npy' % (id_sc)
    d =   np.load(name, allow_pickle=True).item()

    name    =   'fof_sc%04d.npy' % (id_sc)
    arr =   np.load(name)
    n0  =   len(arr)
    arr =   remove_nan(arr)
    ids =   np.where(arr['mag_g'] < 17.)[0]
    arr =   arr[ids]
    n1  =   len(arr)
    g   =   arr['mag_g']
    b   =   arr['mag_bp']
    r   =   arr['mag_rp']
    b_r =   b - r
    ratio   =   calc_ratio(b_r, g)

    return [id_sc, n0, n1, d['d2'], ratio, d['Z'], d['age'], d['shift'][0], d['shift'][1]]

def extract():

    plx =   np.load('ginfos_merge.npy')[:, -2]

    nsc =   2443
   
    f   =   open('sc_info.txt', 'w')
    f.write('# id_sc ntot n17 d2 r_hrd Z age dg db_r d_plx d_g r_dist class\n')
    for i in range(nsc):
        print ('sc %d ...' % (i))
        l   =   extract_fit_iso_ratio(i)

        d_plx   =   1. / plx[i] * 1E3
        dg      =   l[7]
        d_g     =   np.power(10., -dg / 5. + 1.)
        r_dist  =   d_g / d_plx
        
        n17     =   int(l[2])
        d2      =   l[3]
        r_hrd   =   l[4]

        age     =   l[6]

        cls =   classify(n17, age, r_hrd, d2)

        f.write('%5d\t%5d\t%5d\t%.5f\t\t%.5f\t\t%.2f\t\t%.3f\t\t%.3f\t\t%.3f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%1d\n' % \
                (l[0], l[1], l[2], l[3], l[4], \
                 np.log10(l[5]/0.0152), l[6]/1E9, l[7], l[8], \
                 d_plx, d_g, r_dist, cls))
    f.close()

if __name__ == '__main__':
    extract()
