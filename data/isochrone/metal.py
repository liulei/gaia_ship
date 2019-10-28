#!/usr/bin/env python

import numpy as np

lzl  =   -2.0
lzh  =   0.5

zsol    =   0.0152

lz  =   lzl
dlz =   0.25
i   =   0
print '# idx    log10(Z/Zsol)   Z'
while lz <= lzh:
    
    print '%d\t\t%.2f\t\t%.7f' % (i, lz, zsol * np.power(10, lz))
    i   +=  1
    lz  +=  dlz
    
