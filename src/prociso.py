#!/usr/bin/env python

import os, sys
import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib import rc
from mpi4py import MPI
from sklearn.neighbors import KDTree
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic
from isochrone import *

nseg    =   200
comm    =   MPI.COMM_WORLD
rank    =   comm.Get_rank()
size    =   comm.Get_size()

tag_req     =   1
tag_task    =   2

iso = ISO()
iso.load_npy('../../Z.npy')

def fit_iso(id_sc):

    arr =   np.load('fof_sc%04d.npy' % (id_sc))
    arr =   remove_nan_b_r(arr)

    g   =   arr['mag_g']
    b   =   arr['mag_bp']
    r   =   arr['mag_rp']
    b_r =   b - r
    
    d_fit   =   iso.fit_age_Z(g, b_r)
    np.save('fit_iso_sc%04d.npy' % (id_sc), d_fit)

def assign():

    task0    =   int(sys.argv[1])
    task1    =   int(sys.argv[2])
    count_task   =   task0

    while count_task < task1:
        rank =  comm.recv(source = MPI.ANY_SOURCE, tag = tag_req)
        comm.send((count_task, ), dest = rank, tag = tag_task)
        print 'Seg %d, send to proc %d.' % (count_task, rank)
        count_task   +=  1

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
        fit_iso(id_seg)

def main():

    if len(sys.argv) < 3:
        if rank == 0:
            print './prociso.py sc0 sc1'
        sys.exit(0)

    if size < 2:
        if rank == 0:
            print './prociso.py: at least 2 procs are required!'
        sys.exit(0)

    if rank == 0:
        assign()
    else:
        calc()

def main_serial():
    fit_iso(0)

if __name__ == '__main__':
    main()
#    main_serial()
