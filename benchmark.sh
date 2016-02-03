#!/bin/bash
for config in "64 0 nbody.ppm 100000" "128 0 nbody.ppm 100000" "256 0 nbody.ppm 100000" "100 0 nbody.ppm 100" "1000 0 nbody.ppm 100" "10000 0 nbody.ppm 100"; do
	echo $config >> results
	echo seq >> results
	for i in {1..3}; do
		prun -v -1 -np 1 nbody/nbody-seq $config 2>&1 | grep seconds | cut -f4 -d" " >> results
	done
	for i in 1 2 4 8 16; do
		echo par $i >> results
		for j in {1..3}; do
			prun -v -1 -np $i $PRUN_ETC/prun-openmpi nbody/nbody-par $config 2>&1 | grep seconds | cut -f4 -d" " >> results
			wait
		done
	done
done