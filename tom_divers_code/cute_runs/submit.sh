#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=0:59:00
#PBS -N CUTE
#PBS -o stderr/CUTE.o
#PBS -e stderr/CUTE.e
#PBS -m abe
##PBS -M td5g14@soton.ac.uk
#PBS -V
#PBS -q fast.q
cd /users/tdivers/ELGs_eBOSS/CUTE-master/CUTE_box

export OMP_NUM_THREADS=8
./CUTE_box input/Galaxy_Santi_M1_10.00_DLM0.00.txt


