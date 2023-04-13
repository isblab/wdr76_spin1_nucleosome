# /bin/bash

for runid in `seq 31 49` ; do mpirun -np 4 $imp python /home/shreyas/Projects/washburn/wdr76_spin1_nucleosome/scripts/modeling.py prod $runid & done
