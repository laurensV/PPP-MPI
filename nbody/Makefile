nbody-seq: nbody-seq.c
	gcc -Wall -O3 -o nbody-seq nbody-seq.c -lm

nbody-par: nbody-par.c
	mpicc -Wall -O2 -o nbody-par nbody-par.c -lm

clean:
	rm -f *.o nbody-seq nbody-par *~ *core
