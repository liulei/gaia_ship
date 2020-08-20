#!/usr/bin/env python

import os, sys
#from mpi4py import MPI
import numpy as np
import pandas as pd
from astropy import units as u
import scutil
from matplotlib import pyplot as plt
from matplotlib import rc

# $GAIA/select/50/mag17

prefix  =   '.'

nsc =   2443

def gen_match_arr(cat):
    
    m   =   np.ones(nsc, dtype = int) * (-1)

    idfof   =   cat[:, 0].astype(int)
    idcat   =   cat[:, -2].astype(int)

    m[idfof]    =   idcat

    return m

def main():

    cat1  =   np.loadtxt('%s/and/match_K13.txt' % (prefix))
    cat2  =   np.loadtxt('%s/and/match_CG18.txt' % (prefix))
    cat3  =   np.loadtxt('%s/and/match_B19.txt' % (prefix))

#    cat1  =   np.loadtxt('%s/plus/match_K13.txt' % (prefix))
#    cat2  =   np.loadtxt('%s/plus/match_CG18.txt' % (prefix))
#    cat3  =   np.loadtxt('%s/plus/match_B19.txt' % (prefix))

    m1  =   gen_match_arr(cat1)
    m2  =   gen_match_arr(cat2)
    m3  =   gen_match_arr(cat3)

    g0   =   np.load('%s/ginfos_merge.npy' % (prefix))
    gsig  =   np.load('%s/ginfos_merge_sigma.npy' % (prefix))

    sc_info0 = np.loadtxt('%s/sc_info.txt' % (prefix)) 
    cls0 =   sc_info0[:, -1].astype(int)
    
    l      =    g0[:, 1]
    b      =    g0[:, 2]
    plx    =    g0[:, -2]
    r       =   g0[:, 3]
    pmra   =    g0[:, 4]
    pmdec  =    g0[:, 5]

    age    =   sc_info0[:, 6]
    n      =   sc_info0[:, 1].astype(int)
    Z      =   sc_info0[:, 5]

    f   =   open('cat_all.tex', 'w') 
    fmt_04f =   "%d & %.3f $\pm$ %.3f & %.3f $\pm$ %.3f & %.3f & %.3f $\pm$ %.3f & %.3f $\pm$ %.3f & %.3f $\pm$ %.3f & %d & %.4f $\pm$ %0.4f & %.3f & %1d & %4d & %4d & %4d\\\\\n"
    fmt_02f =   "%d & %.3f $\pm$ %.3f & %.3f $\pm$ %.3f & %.3f & %.3f $\pm$ %.3f & %.3f $\pm$ %.3f & %.3f $\pm$ %.3f & %d & %.2f $\pm$ %0.2f & %.3f & %1d & %4d & %4d & %4d\\\\\n"

    for i in range(len(g0)):

        sig =   gsig[i, :]
        if  age[i]  < 0.016:
            fmt =   fmt_04f
        else:
            fmt =   fmt_02f
        f.write( fmt % \
                (i, l[i], sig[0], b[i], sig[1], r[i], plx[i], sig[4], \
                 pmra[i], sig[2], pmdec[i], sig[3], \
                 n[i], age[i], age[i] * 0.06, Z[i], cls0[i], m1[i], m2[i], m3[i]))
        if i == 9:
            break
    f.close()

    f   =   open('cat_new.tex', 'w') 
    fmt_04f =   "%d & %.3f $\pm$ %.3f & %.3f $\pm$ %.3f & %.3f & %.3f $\pm$ %.3f & %.3f $\pm$ %.3f & %.3f $\pm$ %.3f & %d & %.4f $\pm$ %0.4f & %.3f\\\\\n"
    fmt_02f =   "%d & %.3f $\pm$ %.3f & %.3f $\pm$ %.3f & %.3f & %.3f $\pm$ %.3f & %.3f $\pm$ %.3f & %.3f $\pm$ %.3f & %d & %.2f $\pm$ %0.2f & %.3f\\\\\n"

    iout    =   0
    for i in range(len(g0)):

        if cls0[i] != 1:
            continue
        if not (m1[i] == -1 and m2[i] == -1 and m3[i] == -1):
            continue

        sig =   gsig[i, :]
        if  age[i]  < 0.016:
            fmt =   fmt_04f
        else:
            fmt =   fmt_02f
        f.write( fmt % \
                (i, l[i], sig[0], b[i], sig[1], r[i], plx[i], sig[4], \
                 pmra[i], sig[2], pmdec[i], sig[3], \
                 n[i], age[i], age[i] * 0.06, Z[i]))
        iout    +=  1
        if iout == 10:
            break
    f.close()

    nn  =   np.zeros(3, dtype = int)
    f   =   open('cat_all.txt', 'w')
    f1  =   open('cat_new.txt', 'w')
    f.write('id l l_err b b_err r plx plx_err pmra pmra_err pmdec pmdec_err n age age_err Z class K13 CG18 B19\n')
    f1.write('id l l_err b b_err r plx plx_err pmra pmra_err pmdec pmdec_err n age age_err Z\n')
    for i in range(len(g0)):
        sig =   gsig[i, :]
        f.write("%4d  %7.3f  %3.3f  %7.3f  %3.3f  %7.3f  %7.3f  %3.3f  %7.3f  %3.3f  %7.3f  %3.3f  %4d  %7.3f  %7.4f  %7.2f  %2d  %4d  %4d  %4d\n" % \
                (i, l[i], sig[0], b[i], sig[1], r[i], plx[i], sig[4], \
                 pmra[i], sig[2], pmdec[i], sig[3], \
                 n[i], age[i], age[i] * 0.06, Z[i], cls0[i], m1[i], m2[i], m3[i]))
        if m1[i] == -1 and m2[i] == -1 and m3[i] == -1:
#        if m1[i] == -1 and m2[i] == -1:
            nn[cls0[i]-1] +=  1
            if cls0[i] == 1:
                f1.write("%4d  %7.3f  %3.3f  %7.3f  %3.3f  %7.3f  %7.3f  %3.3f  %7.3f  %3.3f  %7.3f  %3.3f  %4d  %7.3f  %7.4f  %7.2f\n" % \
                (i, l[i], sig[0], b[i], sig[1], r[i], plx[i], sig[4], \
                 pmra[i], sig[2], pmdec[i], sig[3], n[i], age[i], age[i] * 0.06, Z[i]))
    f.close()
    f1.close()

    print(nn)

if __name__ == '__main__':
    main()
