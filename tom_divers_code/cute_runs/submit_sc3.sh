#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=100:59:00
#PBS -N CUTE
#PBS -o stderr/CUTE.o
#PBS -e stderr/CUTE.e
#PBS -m abe
##PBS -M td5g14@soton.ac.uk
#PBS -V
#PBS -q sciama3.q
cd /users/tdivers/ELGs_eBOSS/CUTE-master/CUTE_box

export OMP_NUM_THREADS=16
./CUTE_box input/parameter.txt










