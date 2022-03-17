#!/bin/sh

module load openmpi-4.1.1+gnu-9.3.0

make all src=mpi
make all src=omp
make all src=serial

cd bin

declare -a nprocs=(2 4 8 16 24 32 48)

for i in "${nprocs[@]}"
do
    echo "executing mpi with: ${i}"
    mpirun -np ${i} tree_mpi.x ../datasets/float/dataset.csv >> ../time/strong/mpi.out
done

for i in "${nprocs[@]}"
do
    echo "executing omp with: ${i}"
    export OMP_NUM_THREADS=${i} 
    ./tree_omp.x ../datasets/float/dataset.csv >> ../time/strong/omp.out 
done

for i in "${nprocs[@]}"
do
    echo "executing serial with: ${i}"
    export OMP_NUM_THREADS=${i}
    ./tree_serial.x ../datasets/float/dataset.csv >> ../time/strong/serial.out
done
