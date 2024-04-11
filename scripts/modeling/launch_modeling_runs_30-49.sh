#! /bin/bash

for runid in `seq 30 49` ; do mpirun -np 8 $imp python modeling.py prod $runid & done
