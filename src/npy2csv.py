#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys

if len(sys.argv) < 2:
    print('Usage: ./npy2csv.py fof_scXXXX.npy')
    sys.exit(0)

src =   sys.argv[1]
dst =   src[:-4] + '.csv'
arr =   np.load(src)
df  =   pd.DataFrame(arr)
df.to_csv(dst, index=False)

