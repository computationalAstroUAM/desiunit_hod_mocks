#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=0:59:00
#PBS -N HOD
#PBS -o stderr/HOD.o
#PBS -e stderr/HOD.e
#PBS -m abe
##PBS -M santiago.avila@port.ac.uk
#PBS -V
#PBS -q fast.q
cd /users/savila/ELGs_eBOSS/HOD_Santi

./HOD 13.0 0.5 11.0
