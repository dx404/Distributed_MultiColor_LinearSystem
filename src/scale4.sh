#!/bin/bash
steps=32
for i in $(seq 62 62 248)
do 
        q=$(($i/62))
        q000=$(($q*$q*$q))
        q001=$(($q*$q*($q+1)))
        q011=$(($q*($q+1)*($q+1)))
        bsub -M 20 -Ip -n $(($q000+2)) -q hour mpiexec -np $(($q000+1)) ./v8_Isend.out $q $q $q 1 ../data/$i-$i-$i.txt $steps > scale4--$i-$i-$i--$q-$q-$q.result &
        bsub -M 20 -Ip -n $(($q001+2)) -q hour mpiexec -np $(($q001+1)) ./v8_Isend.out $q $q $(($q+1)) 1 ../data/$i-$i-$(($i+62)).txt $steps > scale4--$i-$i-$(($i+62))--$q-$q-$(($q+1)).result &
        bsub -M 20 -Ip -n $(($q011+2)) -q hour mpiexec -np $(($q011+1)) ./v8_Isend.out $q $(($q+1)) $(($q+1)) 1 ../data/$i-$(($i+62))-$(($i+62)).txt $steps > scale4--$i-$(($i+62))-$(($i+62))--$q-$(($q+1))-$(($q+1)).result & 
done
