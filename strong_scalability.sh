#!/bin/sh

cd time
rm *.out
cd ..

cd bin

for i in 2 4
do
    mpirun -np ${i} tree_mpi.x >> ../time/mpi.out
done

for i in 2 4
do
    export OMP_NUM_THREADS=${i} 
    ./tree_omp.x >> ../time/omp.out
done

for i in 2 4
do
    export OMP_NUM_THREADS=${i}
    ./tree_serial.x >> ../time/serial.out
done