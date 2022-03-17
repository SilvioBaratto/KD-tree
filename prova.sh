#!/bin/sh

#PBS -q dssc_gpu
#PBS -l nodes=1:ppn=48
#PBS -l walltime=01:00:00

#cd $PBS_O_WORKDIR

module load openmpi-4.1.1+gnu-9.3.0

cd time
mkdir prova

