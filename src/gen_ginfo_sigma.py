#!/usr/bin/env python

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic

# $GAIA/rev1

prefix  =   '../select/50'
nsc =   2443

sigma_l =   []

for i in range(nsc):
    print('sc %d ...' % (i))
    arr =   np.load('%s/fof_sc%04d.npy' % (prefix, i))

    ra  =   arr['ra']
    dec =   arr['dec']

    pmra    =   arr['pmra']
    pmdec   =   arr['pmdec']

    plx     =   arr['parallax']

    dt  =   -15.5 # from J2015.5 to J2000
# from mas to degree, and project to ra
    ra  +=  dt * pmra / (3600.0E3) / np.cos(dec * np.pi / 180.0)
    dec +=  dt * pmdec / (3600.0E3)

    c   =   SkyCoord(ra = ra * u.degree, dec = dec * u.degree, frame = 'icrs')
    
    c_G =   c.transform_to(Galactic)
    l   =   c_G.l.degree
    b   =   c_G.b.degree

    # correct for 0 ~ 360 jump
    if np.max(l) - np.min(l) > 180.:
        ids =   np.where(l > 180.)[0]
        l[ids]  -=  360.

    l0      =   np.average(l)
    b0      =   np.average(b)
    pmra0   =   np.average(pmra)
    pmdec0  =   np.average(pmdec)
    plx0    =   np.average(plx)

    sigma_l.append([np.std(l), np.std(b), np.std(pmra), np.std(pmdec), \
            np.std(plx), \
            np.average(arr['ra_err']), np.average(arr['dec_err']), \
            np.average(arr['pmra_err']), np.average(arr['pmdec_err']), \
            np.average(arr['parallax_err'])]) 

np.save('ginfos_merge_sigma.npy', sigma_l)
