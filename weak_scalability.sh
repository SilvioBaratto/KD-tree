#!/bin/sh

module load openmpi-4.1.1+gnu-9.3.0

make src=mpi

cd time/weak/
rm mpi.csv

cd ..
cd ..

cd bin

declare -a nprocs=(1 2 4 8 16)

for i in "${nprocs[@]}"
do
    echo "executing mpi with: ${i}"
    mpirun -np ${i} tree_mpi.x ../datasets/float/dataset.csv >> ../time/weak/mpi.csv
done
