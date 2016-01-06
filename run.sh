#!/bin/bash
# module unload mvapich_intel
# module unload mvapich_pgi
# module unload mvapich_gcc
# module load mvapich2_gcc/4.8.1

version=v8_noPipe
steps=64
mpicxx -O3 -std=c++1y ./src/GS_OOP_$version.cpp -o ./debug/$version.out
mv output/* recycle/
bsub -n 10 \
	-q hour -x \
	 -e ./output/err.%J \
	 -o ./output/out.%J \
	mpirun ./debug/$version.out 2 2 2 1 ./data/256_256_256.txt $steps
