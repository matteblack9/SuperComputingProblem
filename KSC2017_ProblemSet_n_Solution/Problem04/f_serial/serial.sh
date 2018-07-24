#!/bin/bash
#$ -V
#$ -cwd
#$ -N Prob4_Fserial
#$ -pe mpi_1cpu 1
#$ -q ksc2017@tachyon0013,ksc2017@tachyon0014,ksc2017@tachyon0015,ksc2017@tachyon0016
#$ -R yes
#$ -o $JOB_NAME.$JOB_ID.out -j y
#$ -l h_rt=48:00:00
time ./f_serial.ex
