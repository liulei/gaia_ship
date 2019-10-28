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

D2R         =   np.pi / 180.
n_star_min  =   20

#l0, l1  =   -180.,  180.    # in degree
#b0, b1  =   -25.,   25.     # in degree
#plx0, plx1  =   0.2, 7.0    # in mas
#d0, d1  =   1./plx1, 1./plx0     # in kpc
#l_sc    =   0.02    # typical scale of star cluster, in kpc
#sigma_plx   =   0.2 # parallax dispersion, in mas

prefix_sel  =   '../select'

def load_stars(id_task):

    l_arr    =   []
    l_key    =   []
    for i in range(nseg):
        name_ids    =   '%s/stars_in_seg%04d.npy' % (prefix_sel, i)
        ids     =   np.load(name_ids)[id_task]
        if len(ids) == 0:
            continue
        name_sel    =   '%s/sel%04d.npy' % (prefix_sel, i)
        arr     =   np.load(name_sel)[ids]
        ns      =   len(arr)
        print 'Partition %d, load seg %d, %d stars.' % \
                (id_task, i, ns)
        l_arr.append(arr)
        key =   (np.ones(len(ids), dtype = 'i8') << 32) * i
        key +=  ids
        l_key.append(key)
        
    return np.concatenate(l_arr, axis = 0), np.concatenate(l_key, axis = 0)

def norm(arr):

    v0  =   np.min(arr)
    v1  =   np.max(arr)

    return (arr - v0) / (v1 - v0)

class Group(object):

    def __init__(self, arr, keys):

        self.nstar  =   len(arr)
        print 'stars: %d' % (self.nstar)
        self.l      =   arr['param'][:, 0]
        self.b      =   arr['param'][:, 1]
        self.plx    =   arr['param'][:, 2]
        self.pmra   =   arr['param'][:, 3]
        self.pmdec  =   arr['param'][:, 4]
        self.arr    =   arr
        self.keys   =   keys

        if np.max(self.l) - np.min(self.l) > 180.:
            ids =   np.where(self.l > 180.)[0]
            self.l[ids] -=  360.
    
        wra     =   np.cos(np.median(self.b) * D2R)
        self.w  =   np.array([wra, 1., 0.5, 1., 1.])
        self.w  /=  np.average(self.w)
        w   =   self.w
#        V   =   w[0] * w[1] * w[2] * w[3] * w[4]

        self.b_fof   =   np.power(1.0 / self.nstar, 1. / 5.) * 0.2
        print 'linking length: %f' % (self.b_fof)

        self.groups  =   []
        self.ngroup  =   0
        self.s2g     =   np.zeros(self.nstar, dtype = int) - 1
        self.nassigned  =   0

    def build_tree(self):

        x0  =   norm(self.l)[:, np.newaxis]
        x1  =   norm(self.b)[:, np.newaxis]
        x2  =   norm(self.plx)[:, np.newaxis]
        x3  =   norm(self.pmra)[:, np.newaxis]
        x4  =   norm(self.pmdec)[:, np.newaxis]

        self.X   =   np.concatenate((x0, x1, x2, x3, x4), axis = 1)
        self.X   *=  self.w

        tree    =   KDTree(self.X, leaf_size = 2)

        return tree

    def merge_group(self, gid0, gid1):
        for k in self.groups[gid1]:
            self.s2g[k]    =   gid0
        self.groups[gid0]   +=  self.groups[gid1]
        self.groups[gid1]   =   []

    def update_neighbor(self, gid0, ids):
        gid =   gid0
        for k in ids:
            if self.s2g[k]  < 0:
                self.s2g[k] = gid
                self.groups[gid].append(k)
                self.nassigned  +=  1
            elif self.s2g[k] != gid:
                gid1    =   self.s2g[k]
                if gid > gid1:
                    gid, gid1  =   gid1, gid
                self.merge_group(gid, gid1)

    def clear_group(self):

        while len(self.groups[-1]) == 0:
            del self.groups[-1]
            self.ngroup -=  1

        gid =   0
        while gid < self.ngroup - 1:
            if len(self.groups[gid]) == 0:
                self.merge_group(gid, self.ngroup - 1)
                while len(self.groups[-1]) == 0:
                    del self.groups[-1]
                    self.ngroup -=  1
            gid +=  1

    def fof(self):

        tree    =   self.build_tree()
        nloop   =   0
        while self.nassigned < self.nstar:
        
            for i in range(self.nstar):
                if self.s2g[i] > 0:
                    continue
                ids =   tree.query_radius([self.X[i]], r = self.b_fof)[0]
                gid =   -1
                for k in ids:
                    if self.s2g[k] >= 0:
                        gid =   self.s2g[i]
                        break
                if gid < 0:
                    self.groups.append([])
                    gid =   self.ngroup
                    self.ngroup +=  1
            
                self.update_neighbor(gid, ids)

            nloop   +=  1
            print 'loop %d: %d stars assigned, %d groups' % \
                (nloop, self.nassigned, self.ngroup)
        
        self.clear_group()
        print '%d groups constructed.' % (self.ngroup)

        ginfo0  =   []
        lens    =   []
        l_key   =   []
        for grp in self.groups:
            if len(grp) < n_star_min:
                continue

            lens.append(len(grp))

            l       =   self.l[grp]
            b       =   self.b[grp]
            plx     =   self.plx[grp]
            pmra    =   self.pmra[grp]
            pmdec   =   self.pmdec[grp]
            l_key.append(self.keys[grp])

# 0 to 360 jump of l has been corrected.  
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
            
#            tpl =   (l0, b0, r_max, pmra0, pmdec0, dpm_max, plx0, dplx_max, grp, seg) 
            tpl =   [len(grp), l0, b0, r_max, pmra0, pmdec0, dpm_max, plx0, dplx_max]
            ginfo0.append(tpl)

        ids =   np.argsort(lens)[::-1]
        l_key1  =   []
        ginfos  =   []
        for id in ids:
            ginfos.append(ginfo0[id])
            l_key1.append(l_key[id])

        return ginfos, l_key1
 
def fof(id_task):
    
    if os.path.exists('ginfos_p%04d.npy' % (id_task)):
        return
    arr, keys =   load_stars(id_task)
    g   =   Group(arr, keys) 
    ginfos, l_keys  =   g.fof()
    np.save('ginfos_p%04d.npy' % (id_task), ginfos)
    np.save('keys_sel_p%04d.npy' % (id_task), l_keys)
     
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
#        select_stars(id_seg)
        fof(id_seg)

def main():

    if len(sys.argv) < 3:
        if rank == 0:
            print './procfof.py partition0 partition1'
        sys.exit(0)

    if size < 2:
        if rank == 0:
            print './procfof.py: at least 2 procs are required!'
        sys.exit(0)

    if rank == 0:
        assign()
    else:
        calc()

def main_serial():
    fof(0)

if __name__ == '__main__':
#    main()
    main_serial()
