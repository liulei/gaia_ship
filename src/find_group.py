#!/usr/bin/env python

import os, sys
#from mpi4py import MPI
import numpy as np
import pandas as pd
import scutil
from sklearn.neighbors import KDTree
from matplotlib import pyplot as plt
from matplotlib import rc

def wrap_at(l):
    
    ids =   np.where(l > 180.)[0]
    l[ids]  -=  360.
    return l

def new_flag_func(row):
    if row.K13 < 0 and row.CG18 < 0 and row.B19 < 0:
        return True
    return False

class Group(object):

    def __init__(self, df):

        self.df =   df
        self.n  =   df.shape[0]

        self.b_fof   =   100.0 # 100 pc

        self.groups  =   []
        self.ngroup  =   0
        self.s2g     =   np.zeros(self.n, dtype = int) - 1
        self.nassigned  =   0

    def build_tree(self):
        x1  =   self.df.x[:, np.newaxis]
        x2  =   self.df.y[:, np.newaxis]
        x3  =   self.df.z[:, np.newaxis]

        self.X   =   np.concatenate((x1, x2, x3), axis = 1)

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
        while self.nassigned < self.n:
        
            for i in range(self.n):
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
            print('loop %d: %d stars assigned, %d groups' % \
                (nloop, self.nassigned, self.ngroup))
        
        self.clear_group()
        print('%d groups constructed.' % (self.ngroup))

        df  =   self.df
        df.index    =   np.arange(df.shape[0])

        fp  =   open('group/sc_groups.txt', 'w')

        npair  =   0
        for grp in self.groups:
            if len(grp) < 2:
                continue
            fp.write('%4d  %d\n' % (npair, len(grp)))
#            self.df.iloc[grp].to_csv('pair/sc_group%04d.txt' % (np), sep = ' ', index = False)
            f   =   open('group/sc_group%04d.txt' % (npair), 'w')
            f.write('id l l_err b b_err r plx plx_err pmra pmra_err pmdec pmdec_err n age Z class K13 CG18 B19 x y z\n')
            for i in grp:
                f.write("%4d  %7.3f  %3.3f  %7.3f  %3.3f  %7.3f  %7.3f  %3.3f  %7.3f  %3.3f  %7.3f  %3.3f  %4d  %7.3f  %7.2f  %2d  %4d  %4d  %4d  %7.3f  %7.3f  %7.3f\n" % \
                (df['id'][i], df.l[i], df.l_err[i], df.b[i], df.b_err[i], 
                 df.r[i], df.plx[i], df.plx_err[i], \
                 df.pmra[i], df.pmra_err[i], df.pmdec[i], df.pmdec_err[i], \
                 df.n[i], df.age[i], df.Z[i], df['class'][i], \
                 df.K13[i], df.CG18[i], df.B19[i], df.x[i], df.y[i], df.z[i]))
            f.close()
            npair  +=  1
        fp.close()
        print('%d groups (len >= 2) constructed.' % (npair))

   
def main():

    df  =   pd.read_csv('match/cat_all.txt', delim_whitespace=True, header = 0)

#    df_sc  =   pd.read_csv('sc_info.txt', delim_whitespace=True, header = -1, \
#            comment = '#')
    
    df['new_flag']  =   df.apply(new_flag_func, axis = 1)

    df  =   df[df['class'] == 1]
    df1 =   df[df['new_flag']]
    df0 =   df[np.logical_not(df['new_flag'])]

    l   =   df.l
    b   =   df.b

    d   =   1E3/df.plx  # in pc

    z   =   d * np.sin(b / 180. * np.pi)
    x   =   -d * np.cos(b / 180. * np.pi) * np.cos(l / 180. * np.pi)
    y   =   d * np.cos(b / 180. * np.pi) * np.sin(l / 180. * np.pi)

    df['x'] =   x
    df['y'] =   y
    df['z'] =   z

    g   =   Group(df)
    g.fof()

if __name__ == '__main__':
    main()
