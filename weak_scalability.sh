#!/bin/sh

module load openmpi-4.1.1+gnu-9.3.0

make all src=mpi
make all src=omp
make all src=serial

cd datasets/float/

c++ write_data.cpp -o write_data.x
echo "writing data"
./write_data.x
echo "finish writing datasets"

cd ..
cd ..

cd bin

declare -a nprocs=(1 2 4 8 16)

for i in "${nprocs[@]}"
do
    echo "executing mpi with: ${i}"
    mpirun -np ${i} tree_mpi.x ../datasets/float/dataset.csv >> ../time/weak/mpi.out
done

for i in "${nprocs[@]}"
do
    echo "executing omp with: ${i}"
    export OMP_NUM_THREADS=${i} 
    ./tree_omp.x ../datasets/float/dataset.csv >> ../time/weak/omp.out 
done

for i in "${nprocs[@]}"
do
    echo "executing serial with: ${i}"
    export OMP_NUM_THREADS=${i}
    ./tree_serial.x ../datasets/float/dataset.csv >> ../time/weak/serial.out
done

cd ..
cd datasets/float/
rm *.csv