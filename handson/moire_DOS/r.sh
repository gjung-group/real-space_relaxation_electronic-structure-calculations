#!/bin/sh
##PBS -l nodes=2:mem128:ppn=20
##PBS -l nodes=2:mem512:ppn=20
#PBS -l nodes=1:ppn=20
#PBS -q mem128
#PBS -N SpectralFunction

##ulimit -s unlimited
#NPROCS=`wc -l < $PBS_NODEFILE`

hostname
date

#cd $PBS_O_WORKDIR

date

export OMP_NUM_THREADS=20
#export KMP_STACKSIZE=10g

mpirun -genv I_MPI_DEBUG 5 -np 1 /home/nleconte//github/JeilJungGroupCodes/lanczosKuboCode/bin/grabnes Gendata.in > job.out

