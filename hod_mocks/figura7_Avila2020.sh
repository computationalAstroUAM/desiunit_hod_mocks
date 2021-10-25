#!/bin/bash

#Mocks Table 3 from Avila, 2020

#La chi squared son del orden de 10-20

#       COMPILER        mu      Ac      As  alphav beta  K   vt vtdisp

./HOD_NFW_V14_more_BVG 10.944 0.00843 0.017237 1 -2 1 0 0      #Mock 1
./HOD_NFW_V14_more_BVG 11.069 0.01213 0.021702 1 0 0.25 0 0      #Mock 2
./HOD_NFW_V14_more_BVG 11.759 0.08799 0.057523 1.5 0 1 0 0      #Mock 3
./HOD_NFW_V14_more_BVG 11.759 0.08799 0.057523 1 0 1 500 200     #Mock 4
./HOD_NFW_V14_more_BVG 11.240 0.02003 0.029333 1 0 0.4 0 0      #Mock 6
./HOD_NFW_V14_more_BVG 11.264 0.02153 0.030546 1 0 0.25 0 0      #Mock 13
./HOD_NFW_V14_more_BVG 10.500 0.00636 0.007098 1 0 0.25 0 0      #Mock 16
#./HOD_NFW_V14_more_BVG 11.069 0.01213 0.021702 1     0  0.25  0   0      #Mock 11 reference

#Mock 2 is equal to 11 in order to compute wp

#Best fit to data (xi2 = 10.9):

./HOD_NFW_V14_more_BVG 11.143 0.01510 0.024783 1 0.1 0.15 0 0      #Mock9 

