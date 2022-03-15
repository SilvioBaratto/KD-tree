#!/bin/sh

#PBS -q dssc_gpu
#PBS -l nodes=1:ppn=48
#PBS -l walltime=01:00:00

cd $PBS_O_WORKDIR

cd Assignment2/KD-tree/

module load openmpi-4.1.1+gnu-9.3.0

cd time
rm -rf strong 
mkdir strong
cd ..

export OMP_PLACES=cores
export OMP_PROC_BIND=true
export OMP_WAIT_POLICY=active
export OMP_DYNAMIC=false

cd bin

for i in {2..48}
do
    mpirun -np ${i} --map-by socket tree_mpi.x >> ../time/strong/mpi.out
done

for i in {2..48}
do
    export OMP_NUM_THREADS=${i} 
    ./tree_omp.x >> ../time/strong/omp.out
done

for i in {2..48}
do
    export OMP_NUM_THREADS=${i}
    ./tree_serial.x >> ../time/strong/serial.out
done
