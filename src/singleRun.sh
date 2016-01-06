#!/bin/bash
for w in $(seq 1 192)
do
	./singleCore.out ../data/256_256_256.txt $w
	./singleCore.out ../data/256_256_256.txt $w
	./singleCore.out ../data/256_256_256.txt $w
done
