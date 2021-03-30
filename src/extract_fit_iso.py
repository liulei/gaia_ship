#!/usr/bin/env python

import numpy as np
import pandas as pd
import os

src_dir =  '../select/50' 

Zsol    =   0.0152

nsc =   2443

f   =   open('fit_iso.txt', 'w')
f.write('id d_g d_b-r age Z\n')
for i in range(nsc):
    print('sc %d...' % (i))
    d   =   np.load('%s/fit_iso_sc%04d.npy' % (src_dir, i), \
            allow_pickle=True, encoding='latin1').flat[0]
    s   =   d['shift']
    f.write('%4d  %7.3f  %7.3f  %7.3f  %7.2f\n' % \
            (i, s[0], s[1], d['age']/1E9, np.log10(d['Z']/Zsol)))
f.close()
