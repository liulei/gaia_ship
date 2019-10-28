#!/bin/bash

#mpirun -f host_gaia -np 121 ../procpartition.py 0 200
#mpirun -f host_gaia -np 121 ../procfof.py 0 4170
#./mergefof_key.py
#mpirun -f host_gaia -np 121 ../prockeyseg2sc.py 0 2443
mpirun -f host_gaia -np 121 ../prociso.py 0 2443
