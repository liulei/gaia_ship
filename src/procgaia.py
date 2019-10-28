#!/usr/bin/env python

import os, sys
from mpi4py import MPI
import numpy as np
import pandas as pd

nseg    =   200

csv_list_filename = '1.txt'

prefix  =   '/home/liulei/program/gaia/ali-gaia_dr2/gaia'

comm    =   MPI.COMM_WORLD
rank    =   comm.Get_rank()
size    =   comm.Get_size()

tag_req     =   1
tag_task    =   2

dtype_gaia  =   np.dtype([  \
                            ('source_id',   'i8'), \
                            ('ra',          'f8'), \
                            ('dec',         'f8'), \
                            ('pmra',        'f8'), \
                            ('pmdec',       'f8'), \
                            ('parallax',    'f8'), \
                            ('mag_g',       'f4'), \
                            ('mag_bp',      'f4'), \
                            ('mag_rp',      'f4'), \
                            ('rv',          'f8'), \
                            ('ra_err',      'f4'), \
                            ('dec_err',     'f4'), \
                            ('pmra_err',    'f4'), \
                            ('pmdec_err',   'f4'), \
                            ('parallax_err','f4'), \
                            ('rv_err',      'f4'), \
                            ('id_file',     'i4'), \
                            ('id_line',     'i4')])

dtype_gaia_err  =   np.dtype([  \
                            ('source_id',   'i8'), \
                            ('ra_err',      'f4'), \
                            ('dec_err',     'f4'), \
                            ('pmra_err',    'f4'), \
                            ('pmdec_err',   'f4'), \
                            ('parallax_err','f4'), \
                            ('rv_err',      'f4'), \
                            ])

def read_csv_list_file(fname):
    f   =   open(fname, 'r')
    lines   =   f.readlines()
    f.close()

    l   =   []
    for line in lines:
        l.append(line.split('/')[-1].rstrip())

    return l

def read_csv(id_file, gz_name):
    
    csv_name    =   gz_name[0:-3]
    full_gz_name    =   '%s/%s' % (prefix, gz_name)
    cmd =   'gunzip < %s > %s' % (full_gz_name, csv_name)
    if os.system(cmd) != 0:
        print 'read_csv(): failed to unzip %s!' % (full_gz_name)
        return np.zeros(0, dtype = dtype_gaia)
    
    df  =   pd.read_csv(full_gz_name, header = 0)
    nrow    =   df.shape[0]

    arr =   np.zeros(nrow, dtype = dtype_gaia)
    arr['source_id']    =   df['source_id']
    arr['ra']           =   df['ra']
    arr['dec']          =   df['dec']
    arr['pmra']         =   df['pmra']
    arr['pmdec']        =   df['pmdec']
    arr['parallax']     =   df['parallax']
    arr['mag_g']        =   df['phot_g_mean_mag']
    arr['mag_bp']       =   df['phot_bp_mean_mag']
    arr['mag_rp']       =   df['phot_rp_mean_mag']
    arr['rv']           =   df['radial_velocity']
    arr['id_file']      =   id_file
    arr['id_line']      =   np.arange(nrow)

#    arr =   np.zeros(nrow, dtype = dtype_gaia_err)
#    arr['source_id']    =   df['source_id']
    arr['ra_err']       =   df['ra_error']
    arr['dec_err']      =   df['dec_error']
    arr['pmra_err']     =   df['pmra_error']
    arr['pmdec_err']    =   df['pmdec_error']
    arr['parallax_err'] =   df['parallax_error']
    arr['rv_err']       =   df['radial_velocity_error']
#
    cmd =   'rm %s' % (csv_name)
    os.system(cmd)

    return arr

def assign():

#    l   =   read_csv_list_file('ali-gaia_dr2_source.txt')
    l   =   read_csv_list_file(csv_list_filename)
    l   =   comm.bcast(l, root = 0) 

    seg0    =   int(sys.argv[1])
    seg1    =   int(sys.argv[2])
    count_seg   =   seg0

    nfile   =   len(l)
    print 'Total csv files in %s:   %d' % (csv_list_filename, nfile)
    nfile_per_seg   =   np.ceil(np.float(nfile) / nseg).astype(int)
    print 'Total segs:              %d' % (nseg)
    print 'Total files per seg:     %d' % (nfile_per_seg)
#    while count_seg < nseg:
    while count_seg < seg1:
        rank =  comm.recv(source = MPI.ANY_SOURCE, tag = tag_req)
        id_file0    =   count_seg * nfile_per_seg
        id_file1    =   id_file0 + nfile_per_seg
        if id_file1 > nfile:
            id_file1 = nfile
        comm.send((count_seg, id_file0, id_file1), dest = rank, tag = tag_task)
        print 'Seg %d (line %d - %d), send to proc %d.' % \
                (count_seg, id_file0, id_file1, rank)
        count_seg   +=  1

    count_calc  =   0
    ncalc   =   size - 1
    while count_calc < ncalc:
        rank =  comm.recv(source = MPI.ANY_SOURCE, tag = tag_req)
        comm.send((-1, -1, -1), dest = rank, tag = tag_task) 
        count_calc  +=  1
    
def calc():
    l   =   None
    l   =   comm.bcast(l, root = 0)
#    print 'rank %d, total csv files: %d' % (rank, len(l))

#    arr =   np.zeros(0, dtype = dtype_gaia)
    while True:
        comm.send(rank, dest = 0, tag = tag_req)

# [id_file0, id_file1):
        (id_seg, id_file0, id_file1) =   comm.recv(source = 0, tag = tag_task)
        if id_seg < 0:
            break
            
        seg_name    =   'gaia_seg%04d.npy' % (id_seg)
        if os.path.exists(seg_name):
            continue

        arr =   []
        for id_file in range(id_file0, id_file1):
            tmp =   read_csv(id_file, l[id_file])
#            tmp =   np.zeros(1, dtype = dtype_gaia)
            arr.append(tmp)
#            arr =   np.concatenate((arr, tmp), axis = 0)

        arr =   np.concatenate(arr, axis = 0)
        np.save(seg_name, arr)

def main():

    if len(sys.argv) < 3:
        if rank == 0:
            print './procgaia.py seg0 seg1'
        sys.exit(0)

    if size < 2:
        if rank == 0:
            print './procgaia.py: at least 2 procs are required!'
        sys.exit(0)

    if rank == 0:
        assign()
    else:
        calc()

def main_serial():
    
    l   =   read_csv_list_file(csv_list_filename)
    for i in range(len(l)):
        print '%s...' % (l[i])
        arr =   read_csv(i, l[i])
        np.save('line%d.npy' % i, arr)

if __name__ == '__main__':
    main()
#    main_serial()
