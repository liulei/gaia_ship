#!/usr/bin/env python

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import scutil

l0, l1  =   -180.,  180.    # in degree
b0, b1  =   -25.,   25.     # in degree
plx0, plx1  =   0.2, 7.0    # in mas
d0, d1  =   1./plx1, 1./plx0     # in kpc
l_sc    =   0.02    # typical scale of star cluster, in kpc
sigma_plx   =   0.2 # parallax dispersion, in mas

nseg_plx   =   8
nseg_b     =   8
nseg_l     =   64

def load_sel(frac = 0.01, prefix = 'select/sel'):

    n0_tot  =   0
    n1_tot  =   0
    l   =   []
    for id_seg in range(200):
#    for id_seg in range(10):
        
        name_sel    =   '%s/sel%04d.npy' % (prefix, id_seg)
#        print '%s ...' % (name_sel)
        a   =   np.load(name_sel)
        n0  =   len(a)
        n1  =   int(n0 * frac)
        n0_tot  +=  n0
        n1_tot  +=  n1
        a   =   np.random.choice(a, size = n1) 
        l.append(a) 

    print 'total: %d, selected: %d' % (n0_tot, n1_tot)
    return np.concatenate(l)

class Node(object):

    def __init__(self):
        self.inds   =   []
        self.vs     =   []
        self.v0     =   -np.inf
        self.v1     =   np.inf
        self.node0  =   None
        self.node1  =   None

def split(n, f):

    if not f(n):
        return

#    print 'node %.3f - %.3f splitted, v0: %.3f, v1: %.3f.' % \
#            (n.v0, n.v1, n.vs[0], n.vs[-1])

    v       =   (n.v0 + n.v1) * 0.5
    ids     =   np.where(n.vs < v)[0]

    if len(ids) == 0:
        iv  =   -1
    else:
        iv      =   ids[-1]

    n0      =   Node()
    n0.inds =   n.inds[:iv+1].copy()
    n0.vs   =   n.vs[:iv+1].copy()
    n0.v0   =   n.v0
    n0.v1   =   v

    n1      =   Node()
    n1.inds =   n.inds[iv+1:].copy()
    n1.vs   =   n.vs[iv+1:].copy()
    n1.v0   =   v
    n1.v1   =   n.v1

    n.inds   =   [] 
    n.vs     =   []
    if len(n0.inds) > 0:
        n.node0  =   n0
        split(n.node0, f)
    if len(n1.inds) > 0:
        n.node1  =   n1
        split(n.node1, f)

def print_plx_node(n):
    plxl, plxh      =   n.v0, n.v1
    dl, dh          =   1./n.v1, 1./n.v0
    print '%.3f - %.3f kpc, %.3f - %.3f mas, %d stars' % \
            (dl, dh, plxl, plxh, len(n.inds))

def print_b_node(n):
    bl, bh      =   n.v0, n.v1
    print '%.3f - %.3f degree, %d stars' % \
            (bl, bh, len(n.inds))

def print_l_node(n):
    ll, lh      =   n.v0, n.v1
    print '%.3f - %.3f degree, %d stars' % \
            (ll, lh, len(n.inds))

def extract_leaf(n, l):

    if n is None:
        return
    
    if len(n.inds) > 0:
#        print_node(n)
        l.append(n) 
    else:
        extract_leaf(n.node0, l)
        extract_leaf(n.node1, l)

def l_partition(inds, ls, plx_max, b_abs_max):

    ids     =   np.argsort(ls)
    inds    =   inds[ids]
    ls      =   ls[ids]

    node    =   Node()
    node.inds   =   inds
    node.vs     =   ls
    node.v0     =   l0
    node.v1     =   l1

    n_tot   =   len(inds) 
    n_min   =   n_tot / nseg_l
    n_min   =   np.maximum(n_min, n_min_g)

    d_min   =   1. / plx_max

# estimate maximum angular radius of sc, in degree
    r_max   =   l_sc / d_min * 0.5 / np.pi * 180.

    fac =   np.cos(b_abs_max / 180. * np.pi) 
    def f(n):

        dl    =   n.v1 - n.v0
        if      len(n.inds) > n_min \
            and dl * fac > r_max * 2:
            return True
        return False
            
    split(node, f)
    
    l_node  =   []
    extract_leaf(node, l_node)

    n_l =   len(l_node)
#    print 'total node: %d, total particles: %d' % (n_l, n_tot)
#    print 'maximum sc angluar size: %.3f degree' % (r_max * 2)
    for i in range(n_l):
        n   =   l_node[i]
#        print_l_node(n)

    return l_node

def b_partition(inds, bs, plx_max):

    ids =   np.argsort(bs)
    inds    =   inds[ids]
    bs      =   bs[ids]

    node    =   Node()
    node.inds   =   inds
    node.vs     =   bs
    node.v0     =   b0
    node.v1     =   b1

    n_tot   =   len(inds) 
    n_min   =   n_tot / nseg_b
    n_min   =   np.maximum(n_min, n_min_g)

    d_min   =   1. / plx_max

# estimate maximum angular radium of sc, in degree
    r_max   =   l_sc / d_min * 0.5 / np.pi * 180.

    def f(n):

        db    =   n.v1 - n.v0
        if      len(n.inds) > n_min \
            and db > r_max * 2:
            return True
        return False
            
    split(node, f)
    
    l_node  =   []
    extract_leaf(node, l_node)

    n_l =   len(l_node)
#    print 'total node: %d, total particles: %d' % (n_l, n_tot)
#    print 'maximum sc angluar size: %.3f degree' % (r_max * 2)
    for i in range(n_l):
        n   =   l_node[i]
#        print_b_node(n)

    return l_node
 
def plx_partition(inds, plxs):

    ids =   np.argsort(plxs)
    inds    =   inds[ids]
    plxs    =   plxs[ids]

    node    =   Node()
    node.inds   =   inds
    node.vs     =   plxs
    node.v0     =   plx0
    node.v1     =   plx1

    n_tot   =   len(inds) 
    n_min   =   n_tot / nseg_plx
    n_min   =   np.maximum(n_min, n_min_g)

    def f(n):

        dplx    =   n.v1 - n.v0
        dd      =   1. / n.v0 - 1. / n.v1
        if      len(n.inds) > n_min \
            and dplx > sigma_plx * 2 \
            and dd > l_sc * 0.5:
            return True
        return False
            
    split(node, f)
    
    l_node  =   []
    extract_leaf(node, l_node)

    n_l =   len(l_node)
    print 'total node: %d, total particles: %d' % (n_l, n_tot)
    for i in range(n_l):
        n   =   l_node[i]
        print_plx_node(n)

    return l_node
        
def main():

    p       =   load_sel(frac = 1)['param']

    l       =   p[:, 0]
    b       =   p[:, 1]
    plx     =   p[:, 2]

# original input to the tree
    ids     =   np.where(l >= 180.)[0]
    l[ids]  -=  360.
    ntot    =   len(l)

    global n_min_g  
    n_min_g =   ntot / (nseg_l * nseg_b * nseg_plx)

    inds    =   np.arange(ntot, dtype = int)

    print ''
    print 'total particles: %d' % (ntot)
    nodes   =   []
    nodes_plx   =   plx_partition(inds, plx)
    d2r =   np.pi / 180.
    for n0 in nodes_plx:
        nodes_b = b_partition(n0.inds, b[n0.inds], n0.v1) 
        for n1 in nodes_b:
            b_abs_max   =   np.maximum(np.abs(n1.v0), np.abs(n1.v1))
            nodes_l = l_partition(n1.inds, l[n1.inds], n0.v1, b_abs_max)
            for n2 in nodes_l:
                nodes.append([n2.v0, n2.v1, n1.v0, n1.v1, n0.v0, n0.v1])
                print 'l: %.3f - %.3f deg, b: %.3f - %.3f deg, plx: %.3f - %.3f mas, %d particles' % (n2.v0, n2.v1, n1.v0, n1.v1, n0.v0, n0.v1, len(n2.inds))
    np.save('gaia_partition.npy', nodes)

if __name__ == '__main__':
    main()
