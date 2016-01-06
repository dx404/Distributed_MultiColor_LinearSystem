#!/bin/bash
queue=hour
d0=$1
test=$2 #003
d1=$(($d0+1))
np0=$(($d0*$d0*$d0))
np1=$(($d0*$d0*$d1))
np2=$(($d0*$d1*$d1))

#echo $d0 $d0 $d0 $np0 512_0${d0}0${d0}0${d0}-$queue-${test}.txt
#echo $d0 $d0 $d1 $np1 512_0${d0}0${d0}0${d1}-$queue-${test}.txt
#echo $d0 $d1 $d1 $np2 512_0${d0}0${d1}0${d1}-$queue-${test}.txt 

#bsub -M 20 -Ip -n $(($np0+2)) -q $queue mpiexec -np $(($np0+1)) ./GS_OOP_v8_noPipe.out $d0 $d0 $d0 1 ../data/512.txt 32 > 512_0${d0}0${d0}0${d0}-$queue-T${test}.txt &
bsub -M 20 -Ip -n $(($np1+2)) -q $queue mpiexec -np $(($np1+1)) ./GS_OOP_v8_noPipe.out $d0 $d0 $d1 1 ../data/512.txt 32 > 512_0${d0}0${d0}0${d1}-$queue-T${test}.txt &
bsub -M 20 -Ip -n $(($np2+2)) -q $queue mpiexec -np $(($np2+1)) ./GS_OOP_v8_noPipe.out $d0 $d1 $d1 1 ../data/512.txt 32 > 512_0${d0}0${d1}0${d1}-$queue-T${test}.txt &
