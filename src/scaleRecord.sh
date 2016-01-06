#!/bin/bash
$id=$1
for i in $(seq 62 62 310)
do 
        q=$(($i/62))
        q000=$(($q*$q*$q))
        q001=$(($q*$q*($q+1)))
        q011=$(($q*($q+1)*($q+1)))
        bsub -M 20 -Ip -n $(($q000+2)) -q hour mpiexec -np $(($q000+1)) ./GS_OOP_v8_noPipe.out $q $q $q 1 ../data/$i-$i-$i.txt 32 > scale$id--$i-$i-$i--$q-$q-$q.result &
        bsub -M 20 -Ip -n $(($q001+2)) -q hour mpiexec -np $(($q001+1)) ./GS_OOP_v8_noPipe.out $q $q $(($q+1)) 1 ../data/$i-$i-$(($i+62)).txt 32 > scale$id--$i-$i-$(($i+62))--$q-$q-$(($q+1)).result &
        bsub -M 20 -Ip -n $(($q011+2)) -q hour mpiexec -np $(($q011+1)) ./GS_OOP_v8_noPipe.out $q $(($q+1)) $(($q+1)) 1 ../data/$i-$(($i+62))-$(($i+62)).txt 32 > scale$id--$i-$(($i+62))-$(($i+62))--$q-$(($q+1))-$(($q+1)).result & 
done
