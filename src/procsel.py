#!/usr/bin/env python

import os, sys
from mpi4py import MPI
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky, Galactic

nseg        =   200
prefix_seg  =   '..'

comm    =   MPI.COMM_WORLD
rank    =   comm.Get_rank()
size    =   comm.Get_size()

tag_req     =   1
tag_task    =   2

dtype_select    =   np.dtype([('idx', 'i4'), ('param', 'f8', 5)])

def select(id_seg):

    seg_name    =   '%s/gaia_seg%04d.npy' % (prefix_seg, id_seg)
    arr =   np.load(seg_name)

    mas2deg =   1E-3 / 3600.0

    ns  =   len(arr)
    id0 =   np.arange(ns, dtype = int)

# remove pmnan
    b_pmra  =   np.isfinite(arr['pmra'])
    b_pmdec =   np.isfinite(arr['pmdec'])
    ids =   np.where(np.logical_and(b_pmra, b_pmdec))[0]
    arr =   arr[ids]
    id0 =   id0[ids]

# select by mag_g
    ids =   np.where(arr['mag_g'] < 18)     
    arr =   arr[ids]
    id0 =   id0[ids]

# select by parallax (0.2 mas - 7 mas)
    b0  =   arr['parallax'] > 0.2
    b1  =   arr['parallax'] < 7.0
    ids =   np.where(np.logical_and(b0, b1))[0]
    arr =   arr[ids]
    id0 =   id0[ids]

# select by pm < 30 mas/yr
    b0  =   np.abs(arr['pmra']) < 30
    b1  =   np.abs(arr['pmdec']) < 30
    ids =   np.where(np.logical_and(b0, b1))[0]
    arr =   arr[ids]
    id0 =   id0[ids]

# select by b:
    
    ra  =   arr['ra']
    dec =   arr['dec']

    pmra    =   arr['pmra']
    pmdec   =   arr['pmdec']

    dt  =   -15.5 # from J2015.5 to J2000
# from mas to degree, and project to ra
    ra  +=  dt * pmra / (3600.0E3) / np.cos(dec * np.pi / 180.0)
    dec +=  dt * pmdec / (3600.0E3)

    c   =   SkyCoord(ra = ra * u.degree, dec = dec * u.degree, frame = 'icrs')
    
    c_G =   c.transform_to(Galactic)
    l   =   c_G.l.degree
    b   =   c_G.b.degree

    ids =   np.where(np.abs(b) < 25.0)[0]

    sel =   np.zeros(len(ids), dtype = dtype_select)

    sel['idx'][:]       =   id0[ids]
    sel['param'][:, 0]  =   l[ids]
    sel['param'][:, 1]  =   b[ids]
    sel['param'][:, 2]  =   arr['parallax'][ids]
    sel['param'][:, 3]  =   pmra[ids]
    sel['param'][:, 4]  =   pmdec[ids]

    np.save('sel%04d.npy' % (id_seg), sel)

def assign():

    seg0    =   int(sys.argv[1])
    seg1    =   int(sys.argv[2])
    count_seg   =   seg0

    while count_seg < seg1:
        rank =  comm.recv(source = MPI.ANY_SOURCE, tag = tag_req)
        comm.send((count_seg, ), dest = rank, tag = tag_task)
        print 'Seg %d, send to proc %d.' % (count_seg, rank)
        count_seg   +=  1

    count_calc  =   0
    ncalc   =   size - 1
    while count_calc < ncalc:
        rank =  comm.recv(source = MPI.ANY_SOURCE, tag = tag_req)
        comm.send((-1,), dest = rank, tag = tag_task) 
        count_calc  +=  1
    
def calc():

    while True:
        comm.send(rank, dest = 0, tag = tag_req)
        (id_seg, ) =   comm.recv(source = 0, tag = tag_task)
        if id_seg < 0:
            break
        select(id_seg)

def main():

    if len(sys.argv) < 3:
        if rank == 0:
            print './procselect.py seg0 seg1'
        sys.exit(0)

    if size < 2:
        if rank == 0:
            print './procselect.py: at least 2 procs are required!'
        sys.exit(0)

    if rank == 0:
        assign()
    else:
        calc()

def main_serial():
    
    select(100)

if __name__ == '__main__':
    main()
#    main_serial()
