#!/usr/bin/env python

import os, sys
import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib import rc
from mpi4py import MPI
from sklearn.neighbors import KDTree
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic
import scutil

nseg    =   200
comm    =   MPI.COMM_WORLD
rank    =   comm.Get_rank()
size    =   comm.Get_size()

tag_req     =   1
tag_task    =   2

keys_sc =   'keys_seg_cluster.npy'
prefix  =   '/data/ll/gaia'

def keyseg2sc(id_sc):
    
    keys    =   np.load(keys_sc)[id_sc]

    segs    =   keys >> 32
    idsseg  =   keys - (segs << 32)
    
    l_arr   =   []
    for i in range(nseg):
        ids =   np.where(segs == i)[0]
        if len(ids) == 0:
            continue
        print 'sc %d, open seg %d' % (id_sc, i)
        arr =   np.load('%s/gaia_seg%04d.npy' % (prefix, i))[idsseg[ids]]
        l_arr.append(arr)

    arr =   np.concatenate(l_arr, axis = 0)
    np.save('fof_sc%04d.npy' % (id_sc), arr)
    
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
        keyseg2sc(id_seg)

def main():

    if len(sys.argv) < 3:
        if rank == 0:
            print './prockeysel2sc.py sc1 sc2'
        sys.exit(0)

    if size < 2:
        if rank == 0:
            print './prockeysel2sc.py: at least 2 procs are required!'
        sys.exit(0)

    if rank == 0:
        assign()
    else:
        calc()

def main_serial():
#    fof(0)
    keyseg2sc(0)

if __name__ == '__main__':
    main()
#    main_serial()
