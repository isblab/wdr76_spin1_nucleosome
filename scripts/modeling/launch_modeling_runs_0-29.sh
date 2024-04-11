#! /bin/bash

for runid in `seq 0 29` ; do mpirun -np 8 $imp python modeling.py prod $runid & done
