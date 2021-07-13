#!/usr/bin/env python

import numpy as np
from sklearn.neighbors import KDTree
from scipy.optimize import minimize

dtype_iso   =   np.dtype([  ('Mass',    'f4'), \
                            ('logTe',   'f4'), \
                            ('g',       'f4'), \
                            ('b',       'f4'), \
                            ('r',       'f4')])
EPS     =   1E-3

# Bounds for g, b_r fit, default
bounds  =   None    

# Bounds for b_r only fit
#bounds  =   ((-EPS, EPS), (-3.0, 3.0))

def remove_nan_b_r(arr):
    
    idb =   np.isnan(arr['mag_bp'])
    idr =   np.isnan(arr['mag_rp'])

    ids =   np.where(np.logical_not(np.logical_or(idb, idr)))[0]

    return arr[ids]

def calc_weight(x0, x1, dx0, dx1):
    
    x00 =   np.min(x0)
    x01 =   np.max(x0)
    x10 =   np.min(x1)
    x11 =   np.max(x1)

    n0  =   int(np.ceil((x01- x00) / dx0))
    n1  =   int(np.ceil((x11- x10) / dx1))

    grid    =   np.zeros((n0, n1), dtype = int)

    nx  =   len(x0)
    for i in range(nx):
        
        ix0 =   int((x0[i] - x00) / dx0)
        ix1 =   int((x1[i] - x10) / dx1)
        grid[ix0, ix1]  +=  1

    w   =   np.zeros(nx, dtype = float)
    for i in range(nx):
        
        ix0 =   int((x0[i] - x00) / dx0)
        ix1 =   int((x1[i] - x10) / dx1)
        w[i]=   1. / grid[ix0, ix1]

    return w
    
def iso_ref(g0, b_r0):

    ids =   np.argsort(g0)
    n   =   len(ids)
    i0  =   int(n * 0.1)
    i1  =   int(n * 0.9)
    ids =   ids[i0:i1]
    g   =   g0[ids]
    b_r =   b_r0[ids]

    ids =   np.where(g > -10)[0]
    g   =   g[ids]
    b_r =   b_r[ids]

    i   =   np.argmin(b_r)
    return np.array([g[i], b_r[i]])

def mag_ref(g, b_r):
    
    cs  =   [g, b_r]
    n   =   len(g)
    cm  =   []
    for c in cs:
        cm.append(np.sort(c)[n / 2])
    return np.array(cm)

def mag_ref3(g, b, r):

    cs  =   [g, b, r]
    n   =   len(g)
    cm  =   []
    for c in cs:
        cm.append(np.sort(c)[n / 2])
    return np.array(cm)

def sc_ref(g, b_r):

    ids =   np.argsort(b_r)
    n   =   len(ids)
    i   =   ids[int(n * 0.1)]
    return np.array([g[i], b_r[i]])

class ISO(object):
    
    def __init__(self):

        self.prefix =   ''
        self.l_Z    =   []
        self.n_Z    =   -1
        self.d      =   {}

    def load_dat(self, prefix):

        self.prefix =   prefix
        table_Z     =   np.loadtxt(prefix + '/table_Z.dat')
        self.l_Z    =   table_Z[:, 2]
        self.n_Z    =   len(self.l_Z)

        for idx_Z in range(self.n_Z):
            print ('%d, Z: %.7f' % (idx_Z, self.l_Z[idx_Z]))
            self.d[idx_Z]    =   self.dat2age_iso_tree(idx_Z)

    def save_npy(self, name):
        np.save(name, (self.l_Z, self.d)) 

    def load_npy(self, name):
        self.l_Z, self.d    =   np.load(name, allow_pickle=True)
        self.n_Z    =   len(self.l_Z)

    def gen_iso_tree(self, iso):

        iso =   iso[:, np.newaxis]
    
        X   =   np.concatenate((iso['g'], iso['b'] - iso['r']), axis = 1)
        tree    =   KDTree(X, leaf_size = 2)

        return tree

    def gen_iso_tree3(self, iso):

        iso =   iso[:, np.newaxis]
    
        X   =   np.concatenate((iso['g'], iso['b'], iso['r']), axis = 1)
        tree    =   KDTree(X, leaf_size = 2)

        return tree


    def dat2age_iso_tree(self, idx_Z):

        l_age   =   []
        l_iso   =   []
        l_tree  =   []
        l_tree3 =   []
   
        name_dat    =   '%s/%d.dat' % (self.prefix, idx_Z)
        arr =   np.loadtxt(name_dat, comments = '#', dtype = float)

        i0      =   0
        ages    =   arr[:, 1]
        age0    =   ages[0]
        for i in range(len(arr)):
            if ages[i] - age0 > 1.0:

                n   =   i - i0        

                iso =   np.zeros(n, dtype = dtype_iso)
                iso['Mass']     =   arr[i0:i, 3]
                iso['logTe']    =   arr[i0:i, 5]
                iso['g']        =   arr[i0:i, 23]
                iso['b']        =   arr[i0:i, 24]
                iso['r']        =   arr[i0:i, 25]

#                iso1    =   interpolate_iso(iso)
                iso1    =   iso
                l_iso.append(iso1)
                l_age.append(age0)
                age0    =   ages[i] 
                i0      =   i
                tree    =   self.gen_iso_tree(iso1)
#                tree3   =   self.gen_iso_tree3(iso1)
                l_tree.append(tree)
#                l_tree3.append(tree3)

        d   =   {}
        d['age']    =   l_age
        d['iso']    =   l_iso
        d['tree']   =   l_tree
#        d['tree3']   =   l_tree3
        return d

    def fit_HRD3(self, idx_Z, idx_iso, X0):

        iso =   self.d[idx_Z]['iso'][idx_iso]
        tree3=   self.d[idx_Z]['tree3'][idx_iso]

        x0  =   mag_ref3(iso['g'], iso['b'], iso['r']) \
                - mag_ref3(X0[:, 0], X0[:, 1], X0[:, 2])

        def f(x):
            X       =   X0.copy()
            X[:, 0] +=  x[0]
            X[:, 1] +=  x[1]
            X[:, 2] +=  x[2]
            d, _    =   tree3.query(X)
            d       =   d.reshape((-1))

            d      =   np.sort(d)
            n       =   len(d)
#            return np.average(d[:int(n * 0.9)] ** 2)
            return np.average(d[:int(n * 0.7)] ** 2)

            return np.sum(d * d) / len(d)

        res =   minimize(f, x0, method = 'Nelder-Mead')
        x  =   res.x

        return f(x), x

    def fit_HRD(self, idx_Z, idx_iso, X0):

        iso =   self.d[idx_Z]['iso'][idx_iso]
        tree=   self.d[idx_Z]['tree'][idx_iso]

        x0  =   iso_ref(iso['g'], iso['b'] - iso['r']) \
                - sc_ref(X0[:, 0], X0[:, 1])
        
#        x0  =   mag_ref(iso['g'], iso['b'] - iso['r']) \
#                - mag_ref(X0[:, 0], X0[:, 1])

        def f(x):
            X       =   X0.copy()
            X[:, 0] +=  x[0]
            X[:, 1] +=  x[1]
            d, _    =   tree.query(X)
            d       =   d.reshape((-1))

            d      =   np.sort(d)
            n       =   len(d)
#            return np.average(d[:int(n * 0.9)] ** 2)
            return np.average(d ** 2, weights = self.w)

            return np.sum(d * d) / len(d)

        if bounds is not None:
            x0[0]   =   0.0

        res =   minimize(f, x0, method = 'Nelder-Mead', bounds = bounds)
        x  =   res.x

        return f(x), x

    def fit_age3(self, idx_Z, X):

        l_tree3    =   self.d[idx_Z]['tree3']

        n   =   len(l_tree3) 
        dn  =   5
        i   =   dn // 2

# Loop 1:
#        print 'Loop 1:'
        l_i     =   []
        l_d2    =   []
        while i < n:
            l_i.append(i)
            d2, shift   =   self.fit_HRD3(idx_Z, i, X)
#            print 'Z %f, age %f, d2 %f' % (self.l_Z[idx_Z], self.d[idx_Z]['age'][i]/1E6, d2)
            l_d2.append(d2)
            i   +=  dn
        i_min   =   l_i[np.argmin(l_d2)]
        
# Loop 2:
#        print 'Loop 2:'
        l_d2    =   []
        l_shift =   []
        i0  =   i_min - dn * 2
        if i0 < 0:
            i0  =   0
        i1  =   i_min + dn * 2
        if i1 > n:
            i1  =   n
        l_i =   range(i0, i1)
        print ('i_min %d, i_max %d' % (np.min(l_i), np.max(l_i)))
        for i in l_i:
            d2, shift  =   self.fit_HRD3(idx_Z, i, X)
#            print 'Z %f, age %f, d2 %f' % (self.l_Z[idx_Z], self.d[idx_Z]['age'][i]/1E6, d2)
            l_d2.append(d2)
            l_shift.append(shift)
        idx     =   np.argmin(l_d2)
        i_min   =   l_i[idx]
        
        d_fit   =   {}
        d_fit['age']    =   self.d[idx_Z]['age'][i_min]
        d_fit['idx_age']=   i_min
        d_fit['d2']     =   l_d2[idx]
        d_fit['shift']  =   l_shift[idx]

        print ('Z = %f, age = %f Myr, d2 = %f, d[g, b, r] = [%f, %f, %f]' % \
            (self.l_Z[idx_Z], d_fit['age'] / 1E6, d_fit['d2'], \
             d_fit['shift'][0], d_fit['shift'][1], d_fit['shift'][2]))

        return d_fit


    def fit_age(self, idx_Z, X):

        l_tree    =   self.d[idx_Z]['tree']

        n   =   len(l_tree) 
        dn  =   5
        i   =   dn // 2

# Loop 1:
#        print 'Loop 1:'
        l_i     =   []
        l_d2    =   []
        while i < n:
            l_i.append(i)
            d2, shift   =   self.fit_HRD(idx_Z, i, X)
#            print 'Z %f, age %f, d2 %f' % (self.l_Z[idx_Z], self.d[idx_Z]['age'][i]/1E6, d2)
            l_d2.append(d2)
            i   +=  dn
        i_min   =   l_i[np.argmin(l_d2)]
        
# Loop 2:
#        print 'Loop 2:'
        l_d2    =   []
        l_shift =   []
        i0  =   i_min - dn * 2
        if i0 < 0:
            i0  =   0
        i1  =   i_min + dn * 2
        if i1 > n:
            i1  =   n
        l_i =   range(i0, i1)
        print ('i_min %d, i_max %d' % (np.min(l_i), np.max(l_i)))
        for i in l_i:
            d2, shift  =   self.fit_HRD(idx_Z, i, X)
#            print 'Z %f, age %f, d2 %f' % (self.l_Z[idx_Z], self.d[idx_Z]['age'][i]/1E6, d2)
            l_d2.append(d2)
            l_shift.append(shift)
        idx     =   np.argmin(l_d2)
        i_min   =   l_i[idx]
        
        d_fit   =   {}
        d_fit['age']    =   self.d[idx_Z]['age'][i_min]
        d_fit['idx_age']=   i_min
        d_fit['d2']     =   l_d2[idx]
        d_fit['shift']  =   l_shift[idx]

        print ('Z = %f, age = %f Myr, d2 = %f, d[g, d-r] = [%f, %f]' % \
            (self.l_Z[idx_Z], d_fit['age'] / 1E6, d_fit['d2'], \
             d_fit['shift'][0], d_fit['shift'][1]))

        return d_fit

    def fit_age_Z3(self, g, b, r):

        Z_table =   [8, 7, 9, 6, 10, 5, 4, 3, 2, 1, 0]

        X   =   np.zeros((len(g), 3), dtype = float)
        X[:, 0] =   g[:]
        X[:, 1] =   b[:]
        X[:, 2] =   r[:]

        d2_min  =   1E9
        n_postfit   =   0
        d_fit_age =   {}
#        for idx_Z in range(self.n_Z):
        for idx_Z in Z_table:
            d_fit_age[idx_Z]  =   self.fit_age3(idx_Z, X)

            d2  =   d_fit_age[idx_Z]['d2']
            if d2_min > d2:
                d2_min  =   d2
                n_postfit   =   0
            else:
                n_postfit   +=  1
            if n_postfit >= 2:
                break 
        
        d2_min  =   1E9
        idx_min =   -1
        shift   =   []
        age     =   -1.0

        d_this  =   {}
        for idx_Z, d in d_fit_age.items():
            if d['d2'] < d2_min:
                d2_min  =   d['d2']
                idx_min =   idx_Z
                d_this['shift']   =   d['shift']
                d_this['age']     =   d['age']
                d_this['idx_age']   =   d['idx_age']
                d_this['idx_Z']     =   idx_Z

        d_this['d2']    =   d2_min
        d_this['Z']     =   self.l_Z[idx_min]

        return d_this

# This is the default entry for isochrone fitting
    def fit_age_Z(self, g, b_r):

        ns  =   len(g)
# For accuracy, select those with g < 17 for fitting
        ids =   np.where(g < 17)[0]
        g   =   g[ids]
        b_r =   b_r[ids]
        print('total: %d, mag_g < 17: %d' % (ns, len(g)))

# For efficiency, set maximum star number to 1000, you may change it
        nmax    =   1000
        ns  =   len(g)
        if ns > nmax:
            ids =   np.random.choice(ns, nmax)
            g   =   g[ids]
            b_r =   b_r[ids]

# Weight are set to 1 for every star in current version
        self.w  =   np.ones(len(g))

# Sequence of metallicty for fitting. Starting with solar metallicity.
# Discarded feature
#        Z_table =   [8, 7, 9, 6, 10, 5, 4, 3, 2, 1, 0]

        X   =   np.zeros((len(g), 2), dtype = float)
        X[:, 0] =   g[:]
        X[:, 1] =   b_r[:]

        d2_min  =   1E9

# n_posfit is not used any more. 
        n_postfit   =   0
        d_fit_age =   {}
        for idx_Z in range(self.n_Z):
#        for idx_Z in Z_table:
            d_fit_age[idx_Z]  =   self.fit_age(idx_Z, X)

            d2  =   d_fit_age[idx_Z]['d2']
            if d2_min > d2:
                d2_min  =   d2
                n_postfit   =   0
            else:
                n_postfit   +=  1
# n_postfit is discarded.
#            if n_postfit >= 2:
#                break 
        
        d2_min  =   1E9
        idx_min =   -1
        shift   =   []
        age     =   -1.0

        d_this  =   {}
        for idx_Z, d in d_fit_age.items():
            if d['d2'] < d2_min:
                d2_min  =   d['d2']
                idx_min =   idx_Z
                d_this['shift']   =   d['shift']
                d_this['age']     =   d['age']
                d_this['idx_age']   =   d['idx_age']
                d_this['idx_Z']     =   idx_Z

        d_this['d2']    =   d2_min
        d_this['Z']     =   self.l_Z[idx_min]

# Label for failed fitting, overlook it.
        if len(g) < 20:
            d_this['d2']    =   len(g)

# Fitting params are stored in this dict
        return d_this

def fit3(id_sc):
    
    arr =   np.load('fof_sc%04d.npy' % (id_sc))
    arr =   remove_nan_b_r(arr)

    g   =   arr['mag_g']
    b   =   arr['mag_bp']
    r   =   arr['mag_rp']
    
    iso =   ISO()
    iso.load_npy('Z.npy')
    d_fit   =   iso.fit_age_Z3(g, b, r)
    print (d_fit)

def fit(id_sc):
    
    arr =   np.load('fof_sc%04d.npy' % (id_sc))
    arr =   remove_nan_b_r(arr)

    g   =   arr['mag_g']
    b   =   arr['mag_bp']
    r   =   arr['mag_rp']
    b_r =   b - r
    
    iso =   ISO()
    iso.load_npy('Z_python3.npy')
    d_fit   =   iso.fit_age_Z(g, b_r)
    print (d_fit)
        
def main_save_npy():
    iso =   ISO()
    iso.load_dat('../data/isochrone')
    iso.save_npy('Z_python3.npy')

if __name__ == '__main__':
    main_save_npy()
#    fit(0)

