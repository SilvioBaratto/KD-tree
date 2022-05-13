#!/bin/sh

module load openmpi-4.1.1+gnu-9.3.0

make src=mpi

cd time/strong/

rm mpi.csv

cd ..

cd ..

cd bin

OMP_NUM_THREADS=1

declare -a nprocs=(2 4 8 16 24 32 48)
#declare -a nprocs=(2 4 8 16 24)

for i in "${nprocs[@]}"
do
    echo "executing mpi with: ${i}"
    mpirun -np ${i} tree_mpi.x ../datasets/float/dataset.csv >> ../time/strong/mpi.csv
done
