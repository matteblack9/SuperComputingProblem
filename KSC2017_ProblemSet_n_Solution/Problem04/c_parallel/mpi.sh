#!/bin/bash
#$ -V
#$ -cwd
#$ -N Prob4_Cmpi
#$ -pe mpi_8cpu 32
#$ -q ksc2017@tachyon0013,ksc2017@tachyon0014,ksc2017@tachyon0015,ksc2017@tachyon0016
#$ -R yes
#$ -o $JOB_NAME.$JOB_ID.out -j y
#$ -l h_rt=48:00:00
mpirun -machinefile $TMPDIR/machines -np $NSLOTS ls
