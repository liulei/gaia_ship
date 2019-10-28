#!/usr/bin/env python

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from mpi4py import MPI
import scutil

nseg    =   200
comm    =   MPI.COMM_WORLD
rank    =   comm.Get_rank()
size    =   comm.Get_size()

tag_req     =   1
tag_task    =   2

l0, l1  =   -180.,  180.    # in degree
b0, b1  =   -25.,   25.     # in degree
plx0, plx1  =   0.2, 7.0    # in mas
d0, d1  =   1./plx1, 1./plx0     # in kpc
l_sc    =   0.02    # typical scale of star cluster, in kpc
sigma_plx   =   0.2 # parallax dispersion, in mas

def select_range(inds, vs, v0, v1):
    
    b0 =   vs > v0
    b1 =   vs < v1
    ids =   np.where(np.logical_and(b0, b1))[0]
    return inds[ids]

def select_stars_per_part_in_seg(arr, ptn):
    
    l   =   arr['param'][:, 0]
    b   =   arr['param'][:, 1]
    plx =   arr['param'][:, 2]
    
    ntot    =   len(arr)
    inds    =   np.arange(ntot, dtype = int)
    
    plxl    =   ptn[4] - sigma_plx
    if plxl < plx0:
        plxl    =   plx0
    plxh    =   ptn[5] + sigma_plx
    if plxh > plx1:
        plxh    =   plx1
    inds0   =   select_range(inds, plx, plxl, plxh)

    d_min   =   1. / ptn[5]
    r_max   =   l_sc / d_min * 0.5 / np.pi * 180.
    bl      =   ptn[2] - r_max
    if bl < b0:
        bl  =   b0
    bh      =   ptn[3] + r_max
    if bh > b1:
        bh  =   b1
    inds1       =   select_range(inds0, b[inds0], bl, bh)

    b_abs_max   =   np.maximum(np.abs(ptn[2]), np.abs(ptn[3]))
    r_max_l     =   r_max / np.cos(b_abs_max / 180. * np.pi)
    ll  =   ptn[0] - r_max_l
    lh  =   ptn[1] + r_max_l
    l   =   l[inds1]
    if ll < -180.:
        inds20 =   select_range(inds1, l, ll + 360., 180.)
        inds21 =   select_range(inds1, l, -180, lh)
        inds2  =   np.concatenate((inds20, inds21), axis = 0)
    elif lh > 180.:
        inds20 =   select_range(inds1, l, lh, 180.)
        inds21 =   select_range(inds1, l, -180, lh - 360.)
        inds2  =   np.concatenate((inds20, inds21), axis = 0)
    else:
        inds2  =   select_range(inds1, l, ll, lh) 

    return inds2

def select_stars_in_seg(arr, partition_table):

    npart   =   len(partition_table)
    stars_in_seg  =   []
    for id_part in range(npart):
        partition   =   partition_table[id_part]
        stars_in_seg.append(select_stars_per_part_in_seg(arr, partition)) 

    return stars_in_seg

def select_stars(id_seg):

    partition_table    =   np.load('gaia_partition.npy')
    npart   =   len(partition_table)
    prefix  =   '.'
    
    stars_in_segs  =   []
    for id_part in range(npart):
        stars_in_segs.append([])

    print 'seg %04d ...' % (id_seg)
    arr =   np.load('%s/sel%04d.npy' % (prefix, id_seg)) 
    ids =   np.where(arr['param'][:, 0] > 180.)[0]
    arr['param'][ids, 0]    -=  360.
    stars_in_seg = select_stars_in_seg(arr, partition_table)
    
    np.save('stars_in_seg%04d.npy' % (id_seg), stars_in_seg)

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
        (id_task, ) =   comm.recv(source = 0, tag = tag_task)
        if id_task < 0:
            break
        select_stars(id_task)
#        fof(id_task)

def main():

    if len(sys.argv) < 3:
        if rank == 0:
            print './procpartition.py seg0 seg1'
        sys.exit(0)

    if size < 2:
        if rank == 0:
            print './procseg.py: at least 2 procs are required!'
        sys.exit(0)

    if rank == 0:
        assign()
    else:
        calc()

if __name__ == '__main__':
    main()
