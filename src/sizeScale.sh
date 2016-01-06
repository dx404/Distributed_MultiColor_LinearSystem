#!/bin/bash

mpiexec -np 10 ./GS_OOP_v8_noPipe.out 2 2 2 1 ../data/256.txt 32

####   run_mycode    ####
#BSUB -n 10
#BSUB –R “span[hosts=1]”
./mycode
##### end of run_mycode ####

for d in $(seq 64 32 256);do mpiexec -np 10 ./GS_OOP_v8_noPipe.out 2 2 2 1 ../data/$d.txt 32; done
