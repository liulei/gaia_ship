#!/usr/bin/env python

import os, sys
import numpy as np
from matplotlib import pyplot as plt

def loghist(vals, x0, x1, nb):
    
    lx0 =   np.log10(x0)
    lx1 =   np.log10(x1)
    dlx =   (lx1 - lx0) / nb
    count   =   np.zeros(nb, dtype = float) 
    
    lvs =   np.log10(vals)
    for lv in lvs:
        i   =   np.int((lv - lx0) / dlx)
        if i < 0:
#            i   =   0
            continue
        if i > nb - 1:
            continue
        count[i]    +=  1.0
    lx  =   lx0 + (np.arange(nb) + 0.5) * dlx 
    return count, np.power(10, lx)

# K13
dtype_cat = np.dtype([ \
                ('name',    'a17'), \
                ('ra',      'f8'), \
                ('dec',     'f8'), \
                ('l',       'f8'), \
                ('b',       'f8'), \
                ('r2',      'f8'), \
                ('pmra',    'f8'), \
                ('pmdec',   'f8'), \
                ('d',       'f8'), \
                ('age',     'f8'), \
                ('Z',       'f8')])

def partition(xmin, xmax, nx, f, nb):

    dx  =   (xmax - xmin) / nx
    x   =   np.arange(nx) * dx
    y   =   f(x, x + dx)
    print (y)

    ycum    =   np.cumsum(y)
    ytot    =   ycum[-1]
    
    ypart   =   ytot / nb 

    print ('ytot: %f, ypart: %f' % (ytot, ypart))
    
    y0  =   0.0
    k   =   0
    yb  =   []
    kb  =   []
    for i in range(nb):
        y1  =   y0 + ypart
        y0  =   y1
        while k < nx - 1:
            if ycum[k + 1] > y1:
                yb.append(ycum[k])
                kb.append(k)
                print ('i: %d, yb: %f, kb: %d' % (i, yb[-1], kb[-1]))
                break
            k   +=  1
    return kb

def read_cat(name_cat):
    
    f   =   open(name_cat)
    lines   =   f.readlines()
    f.close()

    nsc =   len(lines)
    cat =   np.zeros(nsc, dtype = dtype_cat)
    
    for i, l in enumerate(lines):
        cat[i]['name']  =   l[5:22]
        cat[i]['ra']    =   float(l[25:34])
        cat[i]['dec']   =   float(l[34:42])
        cat[i]['l']     =   float(l[42:50])
        cat[i]['b']     =   float(l[50:58])
        cat[i]['r2']    =   float(l[72:79])
        cat[i]['pmra']  =   float(l[79:86])
        cat[i]['pmdec'] =   float(l[86:93])
        cat[i]['d']     =   float(l[142:150])
        cat[i]['age']   =   float(l[185:192])
        cat[i]['Z']     =   float(l[261:269])
        
#    print 'total sc in %s: %d' % (name_cat, len(cat))
    return cat

dtype_CG18 = np.dtype([ \
                ('l',       'f8'), \
                ('b',       'f8'), \
                ('ra',      'f8'), \
                ('dec',     'f8'), \
                ('r50',     'f8'), \
                ('pmra',    'f8'), \
                ('pmdec',   'f8'), \
                ('plx',     'f8'), \
                ('age',     'f8'), \
                ('name',    'a17')])

def read_CG18b(name_cat):

    f   =   open(name_cat)
    lines   =   f.readlines()
    f.close()

    nsc =   len(lines)
    cat =   np.zeros(nsc, dtype = dtype_CG18)
    
    for i, l in enumerate(lines):

        w   =   l.split()

        cat[i]['name']  =   w[2]
        cat[i]['ra']    =   float(w[3])
        cat[i]['dec']   =   float(w[4])
       
        cat[i]['l']     =   float(w[5])
        cat[i]['b']     =   float(w[6])

        cat[i]['r50']   =   float(w[7])

        cat[i]['pmra']  =   float(w[9])
        cat[i]['pmdec'] =   float(w[11])

        cat[i]['plx']   =   float(w[13])

    print ('total sc in %s: %d' % (name_cat, len(cat)))
    return cat

def read_CG18(name_cat):
    
    f   =   open(name_cat)
    lines   =   f.readlines()
    f.close()

    nsc =   len(lines)
    cat =   np.zeros(nsc, dtype = dtype_CG18)
    
    for i, l in enumerate(lines):

        w   =   l.split()

        cat[i]['name']  =   l[24:24+17]
        cat[i]['ra']    =   float(l[42:49])
        cat[i]['dec']   =   float(l[50:57])
       
        cat[i]['l']     =   float(l[58:65])
        cat[i]['b']     =   float(l[66:73])

        cat[i]['r50']   =   float(l[75:80])

        cat[i]['pmra']  =   float(w[13])
        cat[i]['pmdec'] =   float(w[14])

        cat[i]['plx']   =   float(w[15])

    print ('total sc in %s: %d' % (name_cat, len(cat)))

    return cat

# Bossini et al. 2019
dtype_B19 = np.dtype([ \
                ('ra',      'f8'), \
                ('dec',     'f8'), \
                ('age',     'f8')])

def read_Bossini19(name_cat):
    
    f   =   open(name_cat)
    lines   =   f.readlines()
    f.close()

    nsc =   len(lines)
    cat =   np.zeros(nsc, dtype = dtype_B19)
    
    for i, l in enumerate(lines):
        w   =   l.split()
        cat[i]['ra']    =   float(w[3])
        cat[i]['dec']   =   float(w[4])
        cat[i]['age']   =   float(w[5])

#    print 'total sc in %s: %d' % (name_cat, len(cat))
    return cat

dtype_Bica19 = np.dtype([ \
                ('l',       'f8'), \
                ('b',       'f8'), \
                ('ra',      'f8'), \
                ('dec',     'f8'), \
                ('rmaj',    'f8'), \
                ('rmin',    'f8'), \
                ('r',       'f8')])

def read_B19(name_cat):
    
    f   =   open(name_cat)
    lines   =   f.readlines()
    f.close()

    nsc =   len(lines)
    cat =   np.zeros(nsc, dtype = dtype_Bica19)
    
    for i, l in enumerate(lines):

        h   =   int(l[14:16])
        m   =   int(l[17:19])
        s   =   int(l[20:22])
        cat[i]['ra']    =   (h + m / 60. + s / 3600.) * 15.

        d   =   int(l[24:26])
        m   =   int(l[27:29])
        s   =   int(l[30:32])
        cat[i]['dec']   =   d + m / 60. + s / 3600.
       
        cat[i]['l']     =   float(l[0:6])
        cat[i]['b']     =   float(l[7:13])

        cat[i]['rmaj']  =   float(l[33:40]) / 60. * 0.5
        cat[i]['rmin']  =   float(l[41:48]) / 60. * 0.5

        cat[i]['r']     =   max(cat[i]['rmaj'], cat[i]['rmin'])

#    print 'total sc in %s: %d' % (name_cat, len(cat))
    return cat

def main():

    cat =   read_cat('catalog.dat')
    age =   cat['age']
    Z   =   cat['Z']
    
    plt.hist(age, bins = 100)    
#    Z0  =   Z[np.where(Z < 10)[0]]
#    plt.hist(Z0, bins = 50)    
    plt.show()
    
def main_CG18():

    cat =   read_CG18('CG18.txt')
#    print (cat[0:10])

def main_B19():

    cat =   read_B19('Bica.txt')
#    print (cat[0:10]

if __name__ == '__main__':
    main_CG18()

