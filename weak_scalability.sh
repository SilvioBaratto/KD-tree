#!/bin/sh

#PBS -q dssc_gpu
#PBS -l nodes=1:ppn=48
#PBS -l walltime=01:00:00

cd $PBS_O_WORKDIR

cd Assignment2/KD-tree/

module load openmpi-4.1.1+gnu-9.3.0

cd time
rm -d weak 
mkdir weak 
cd ..

module load openmpi-4.1.1+gnu-9.3.0

export OMP_PLACES=cores
export OMP_PROC_BIND=true
export OMP_WAIT_POLICY=active
export OMP_DYNAMIC=false

cd bin

for i in 2 4 8 16 32
do
    mpirun -np ${i} --map-by socket tree_mpi.x >> ../time/weak/mpi.out
done

for i in 2 4 8 16 32
do
    export OMP_NUM_THREADS=${i} 
    ./tree_omp.x >> ../time/weak/omp.out
done

for i in 2 4 8 16 32
do
    export OMP_NUM_THREADS=${i}
    ./tree_serial.x >> ../time/weak/serial.out
done