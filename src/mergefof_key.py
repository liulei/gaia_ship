#!/usr/bin/env python

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import scutil
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic, ICRS

partitions  =   np.load('gaia_partition.npy')

nseg        =   200
npartition  =   len(partitions)
nmin        =   50
fmin        =   0.5

def is_merge(g0, k0, g1, k1):

    n0  =   int(g0[0])
    n1  =   int(g1[0])
    
    assert n0 == len(k0)
    assert n1 == len(k1)

    kist    =   np.intersect1d(k0, k1)
    ni  =   len(kist)
    if ni > fmin * min(n0, n1):
        return True
    return False

def is_overlap(g0, g1):

    fac =   1.0
    
    p0  =   Galactic(l = g0[1] * u.degree, b = g0[2] * u.degree)
    p1  =   Galactic(l = g1[1] * u.degree, b = g1[2] * u.degree)
    dr  =   p0.separation(p1).degree
    if dr > (g0[3] + g1[3]) * fac:
        return False
    
    dpm =   np.sqrt((g0[4] - g1[4]) ** 2 + (g0[5] - g1[5]) ** 2)
    if dpm > (g0[6] + g1[6]) * fac:
        return False

    dplx    =   np.abs(g0[7] - g1[7])
    if dplx > (g0[8] + g1[8]) * fac:
        return False

    return True 

def calc_avg_r(v0, v1, r0, r1, n0, n1):
    
    n   =   n0 + n1
    v   =   (v0 * n0 + v1 * n1) / n
    vmin    =   np.min([v0 - r0, v1 - r1])
    vmax    =   np.max([v0 + r0, v1 + r1])
    r   =   np.max([np.abs(v - vmin), np.abs(v - vmax)])
    return v, r

def merge_group(g0, k0, g1, k1):

    n0, n1  =   g0[0], g1[0]

    dec, rdec   =   calc_avg_r(g0[2], g1[2], g0[3], g1[3], n0, n1)
    c   =   1. / np.cos(dec / 180. * np.pi)
    if np.abs(g0[1] - g1[1]) > 180:
        if g1[1] > g0[1]:
            g1[1]   -=  360.
        else:
            g1[1]   +=  360.
    ra, rra     =   calc_avg_r(g0[1], g1[1], g0[3] * c, g1[3] * c, n0, n1)
    if ra > 360:
        ra  -=  360.
    elif ra < 0:
        ra  +=  360.
    rra /=  c
    r   =   np.max([rra, rdec])
    
    pmra, rpmra     =   calc_avg_r(g0[4], g1[4], g0[6], g1[6], n0, n1)
    pmdec, rpmdec   =   calc_avg_r(g0[5], g1[5], g0[6], g1[6], n0, n1)
    rpm =   np.max([rpmra, rpmdec]) 
    
    plx, rplx   =   calc_avg_r(g0[7], g1[7], g0[8], g1[8], n0, n1)

    ku  =   np.union1d(k0, k1)
    n   =   len(ku)

    return [n, ra, dec, r, pmra, pmdec, rpm, plx, rplx], ku

def print_group(g):
    print 'Len %d, ra %.2f, dec %.2f, r %.2f, pmra %.2f, pmdec %.2f, rpm %.2f, plx %.2f, dplx %.2f' % (g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8])
 
def insert_to_group(gl, kl, g, k):

    n   =   len(gl)
    if n == 0:
        gl.append(g)
        kl.append(k)
        return
    for i in range(n):
#        if is_overlap(gl[i], g):
        if is_merge(gl[i], kl[i], g, k):
            print 'Group %d:' % (i)
            print_group(gl[i])
            print_group(g)
            gl[i], kl[i]   =   merge_group(gl[i], kl[i], g, k)
            print_group(gl[i])
            print ''
            return
    gl.append(g)
    kl.append(k)

def gen_ginfo(arr):

    l   =   arr['param'][:, 0]
    b   =   arr['param'][:, 1]
    plx =   arr['param'][:, 2]
    pmra    =   arr['param'][:, 3]
    pmdec   =   arr['param'][:, 4]

# correct for 0 ~ 360 jump
    if np.max(l) - np.min(l) > 180.:
        ids =   np.where(l > 180.)[0]
        l[ids]  -=  360.

    l0      =   np.average(l)
    b0      =   np.average(b)
    pmra0   =   np.average(pmra)
    pmdec0  =   np.average(pmdec)
    plx0    =   np.average(plx)
            
    g0  =   Galactic(l = l0 * u.degree, b = b0 * u.degree)
    gs  =   Galactic(l = l * u.degree, b = b * u.degree)

    rs  =   g0.separation(gs).degree
    r_max   =   np.max(rs)
            
    dpm     =   np.sqrt((pmra - pmra0) ** 2 + (pmdec - pmdec0) ** 2)
    dpm_max =   np.max(dpm) 
            
    dplx_max=   np.max(np.abs(plx - plx0))
            
    tpl =   [len(arr), l0, b0, r_max, pmra0, pmdec0, dpm_max, plx0, dplx_max]
    return tpl

def update_gk(k):
    
    segs    =   k >> 32
    idssel   =   k - (segs << 32)
    
    l_arr       =   []
    l_keyseg    =   []
    for i in range(nseg):
        ids =   np.where(segs == i)[0]
        if len(ids) == 0:
            continue
        arr =   np.load('../sel/sel%04d.npy' % (i))[idssel[ids]]
        l_arr.append(arr)
        keyseg  =   ((i << 32) + arr['idx']).astype('i8')
        l_keyseg.append(keyseg)
    arr =   np.concatenate(l_arr, axis = 0)
    keyseg  =   np.concatenate(l_keyseg, axis = 0)
    ginfo   =   gen_ginfo(arr) 
    return ginfo, keyseg
    
def main():

    gl  =   []
    kl  =   []

    for i in range(npartition):
        gs   =   np.load('ginfos_p%04d.npy' % (i))
        if len(gs) == 0:
            continue
        ks  =   np.load('keys_sel_p%04d.npy' % (i))
        for k in range(len(gs)):
            if len(ks[k]) < nmin:
                continue
            insert_to_group(gl, kl, gs[k], ks[k])
            if (i + 1) % 100 == 0:
                print '%d groups inserted.' % (i + 1)
                print '' 

    nc0  =   len(gl)
    assert nc0 == len(kl)

    while True:

        gl0 =   []
        kl0 =   []
        for k in range(len(kl)):
            insert_to_group(gl0, kl0, gl[k], kl[k])
        nc  =   len(gl0)
        assert nc == len(kl0)
        print '##### nc: %d, nc0: %d' % (nc, nc0)
        gl  =   gl0
        kl  =   kl0
        if nc0 == nc:
            break
        else:
            nc0 =   nc

    ksegs   =   []
    for i in range(nc):
        print 'update sc %d...' % (i)
        gl[i], kseg   =   update_gk(kl[i])
        ksegs.append(kseg)

    print 'total groups after merge: %d.' % (len(gl))
    np.save('ginfos_merge.npy', gl)
    np.save('keys_seg_cluster.npy', ksegs)
        
if __name__ == '__main__':
    main()
