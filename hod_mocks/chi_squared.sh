#!/bin/bash

#CHI SQUARED MOCKS




#1. VARYING FSAT (1.5 - 0.6) AND BETA (-0.2,0.2)

#1.0.  Calibration with Avila, 2020. We choose mock 0 in order to compute wp chi squared and compare with the paper

#./HOD_NFW_V14_more_BVG_product 11.738 0.083107 0.056817 1 0 1 0 0
#chi_wp^^2 paper: 21.7 ; chi_wp^2 mia: 20.34

#1.1.  with number density n_eBOS -> clustering not changes significatively, only more noise than considering 10*n_eBOSS.

#./HOD_NFW_V14_more_BVG_product 11.877 0.012306 0.0057432 1 -0.08 1 0 0

#1.2.  with number density 10*n_eBOSS

#fsat = 0
#./HOD_NFW_V14_more_BVG_product 12.138 0.26041 0 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.138 0.26041 0 1 -0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.138 0.26041 0 1 -0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.138 0.26041 0 1 -0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.138 0.26041 0 1 -0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.138 0.26041 0 1 0.0 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.138 0.26041 0 1 0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.138 0.26041 0 1 0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.138 0.26041 0 1 0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.138 0.26041 0 1 0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.138 0.26041 0 1 0.2 1 0 0


#fsat = 0.05
#./HOD_NFW_V14_more_BVG_product 12.057 0.20588 0.03232 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.057 0.20588 0.03232 1 -0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.057 0.20588 0.03232 1 -0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.057 0.20588 0.03232 1 -0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.057 0.20588 0.03232 1 -0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.057 0.20588 0.03232 1 0.0 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.057 0.20588 0.03232 1 0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.057 0.20588 0.03232 1 0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.057 0.20588 0.03232 1 0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.057 0.20588 0.03232 1 0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 12.057 0.20588 0.03232 1 0.2 1 0 0


#fsat = 0.10
#./HOD_NFW_V14_more_BVG_product 11.970 0.16028 0.050109 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.970 0.16028 0.050109 1 -0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.970 0.16028 0.050109 1 -0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.970 0.16028 0.050109 1 -0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.970 0.16028 0.050109 1 -0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.970 0.16028 0.050109 1 0.0 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.970 0.16028 0.050109 1 0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.970 0.16028 0.050109 1 0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.970 0.16028 0.050109 1 0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.970 0.16028 0.050109 1 0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.970 0.16028 0.050109 1 0.2 1 0 0


#fsat = 0.15
#./HOD_NFW_V14_more_BVG_product 11.877 0.12306 0.057432 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.877 0.12306 0.057432 1 -0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.877 0.12306 0.057432 1 -0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.877 0.12306 0.057432 1 -0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.877 0.12306 0.057432 1 -0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.877 0.12306 0.057432 1 0.0 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.877 0.12306 0.057432 1 0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.877 0.12306 0.057432 1 0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.877 0.12306 0.057432 1 0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.877 0.12306 0.057432 1 0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.877 0.12306 0.057432 1 0.2 1 0 0

#fsat = 0.2
#./HOD_NFW_V14_more_BVG_product 11.779 0.093298 0.057967 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.779 0.093298 0.057967 1 -0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.779 0.093298 0.057967 1 -0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.779 0.093298 0.057967 1 -0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.779 0.093298 0.057967 1 -0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.779 0.093298 0.057967 1 0.0 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.779 0.093298 0.057967 1 0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.779 0.093298 0.057967 1 0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.779 0.093298 0.057967 1 0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.779 0.093298 0.057967 1 0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.779 0.093298 0.057967 1 0.2 1 0 0

#fsat = 0.25
#./HOD_NFW_V14_more_BVG_product 11.675 0.06954 0.054137 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.675 0.06954 0.054137 1 -0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.675 0.06954 0.054137 1 -0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.675 0.06954 0.054137 1 -0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.675 0.06954 0.054137 1 -0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.675 0.06954 0.054137 1 0.0 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.675 0.06954 0.054137 1 0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.675 0.06954 0.054137 1 0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.675 0.06954 0.054137 1 0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.675 0.06954 0.054137 1 0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.675 0.06954 0.054137 1 0.2 1 0 0

#fsat = 0.3
#./HOD_NFW_V14_more_BVG_product 11.566 0.051077 0.048096 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.566 0.051077 0.048096 1 -0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.566 0.051077 0.048096 1 -0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.566 0.051077 0.048096 1 -0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.566 0.051077 0.048096 1 -0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.566 0.051077 0.048096 1 0.0 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.566 0.051077 0.048096 1 0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.566 0.051077 0.048096 1 0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.566 0.051077 0.048096 1 0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.566 0.051077 0.048096 1 0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.566 0.051077 0.048096 1 0.2 1 0 0

#fsat = 0.35
#./HOD_NFW_V14_more_BVG_product 11.452 0.036988 0.041141 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.452 0.036988 0.041141 1 -0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.452 0.036988 0.041141 1 -0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.452 0.036988 0.041141 1 -0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.452 0.036988 0.041141 1 -0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.452 0.036988 0.041141 1 0.0 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.452 0.036988 0.041141 1 0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.452 0.036988 0.041141 1 0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.452 0.036988 0.041141 1 0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.452 0.036988 0.041141 1 0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.452 0.036988 0.041141 1 0.2 1 0 0

#fsat = 0.4
#./HOD_NFW_V14_more_BVG_product 11.335 0.026464 0.034333 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.335 0.026464 0.034333 1 -0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.335 0.026464 0.034333 1 -0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.335 0.026464 0.034333 1 -0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.335 0.026464 0.034333 1 -0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.335 0.026464 0.034333 1 0.0 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.335 0.026464 0.034333 1 0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.335 0.026464 0.034333 1 0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.335 0.026464 0.034333 1 0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.335 0.026464 0.034333 1 0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.335 0.026464 0.034333 1 0.2 1 0 0

#fsat = 0.45
#./HOD_NFW_V14_more_BVG_product 11.215 0.018646 0.028087 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.215 0.018646 0.028087 1 -0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.215 0.018646 0.028087 1 -0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.215 0.018646 0.028087 1 -0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.215 0.018646 0.028087 1 -0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.215 0.018646 0.028087 1 0.0 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.215 0.018646 0.028087 1 0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.215 0.018646 0.028087 1 0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.215 0.018646 0.028087 1 0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.215 0.018646 0.028087 1 0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.215 0.018646 0.028087 1 0.2 1 0 0

#fsat = 0.5
#./HOD_NFW_V14_more_BVG_product 11.094 0.013067 0.022717 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.094 0.013067 0.022717 1 -0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.094 0.013067 0.022717 1 -0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.094 0.013067 0.022717 1 -0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.094 0.013067 0.022717 1 -0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.094 0.013067 0.022717 1 0.0 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.094 0.013067 0.022717 1 0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.094 0.013067 0.022717 1 0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.094 0.013067 0.022717 1 0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.094 0.013067 0.022717 1 0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 11.094 0.013067 0.022717 1 0.2 1 0 0

#fsat = 0.55
#./HOD_NFW_V14_more_BVG_product 10.969 0.009056 0.018062 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.969 0.009056 0.018062 1 -0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.969 0.009056 0.018062 1 -0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.969 0.009056 0.018062 1 -0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.969 0.009056 0.018062 1 -0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.969 0.009056 0.018062 1 0.0 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.969 0.009056 0.018062 1 0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.969 0.009056 0.018062 1 0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.969 0.009056 0.018062 1 0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.969 0.009056 0.018062 1 0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.969 0.009056 0.018062 1 0.2 1 0 0

#fsat = 0.6
#./HOD_NFW_V14_more_BVG_product 10.841 0.006236 0.014182 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.841 0.006236 0.014182 1 -0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.841 0.006236 0.014182 1 -0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.841 0.006236 0.014182 1 -0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.841 0.006236 0.014182 1 -0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.841 0.006236 0.014182 1 0.0 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.841 0.006236 0.014182 1 0.04 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.841 0.006236 0.014182 1 0.08 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.841 0.006236 0.014182 1 0.12 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.841 0.006236 0.014182 1 0.16 1 0 0
#./HOD_NFW_V14_more_BVG_product 10.841 0.006236 0.014182 1 0.2 1 0 0


#POSSIBLE ERRORS IN THE GALAXY BIAS OBTENTION

#1. Considering 1% of particles in the calculation of halo bias:

#considering n = 10n_eboss
#fsat = 0.6
#./HOD_NFW_V14_more_BVG_product_particles 11.090 0.010219 0.026822 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product_particles 11.090 0.010219 0.026822 1 0.0 1 0 0
#fsat = 0.15
#./HOD_NFW_V14_more_BVG_product_particles 12.017 0.168016 0.085443 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product_particles 12.017 0.168016 0.085443 1 0.0 1 0 0

#Considering n = n_eBOSS

#fsat = 0.6
#./HOD_NFW_V14_more_BVG_product_particles 11.090 0.001022 0.002682 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product_particles 11.090 0.001022 0.002682 1 0.0 1 0 0
#fsat = 0.15
#./HOD_NFW_V14_more_BVG_product_particles 12.017 0.016802 0.008544 1 -0.2 1 0 0
#./HOD_NFW_V14_more_BVG_product_particles 12.017 0.016802 0.008544 1 0.0 1 0 0

#extending mass bins

#fsat = 0.6, n = 10n_eboss b = 1.327
./HOD_NFW_V14_more_BVG_product_particles 10.803 0.005786 0.012801 1 -0.2 1 0 0

#fsat = 0.6, n = 10n_eboss b = 1.359
./HOD_NFW_V14_more_BVG_product_particles 10.913 0.007197 0.016950 1 0 1 0 0 #poisson
./HOD_NFW_V14_more_BVG_product_particles 10.913 0.007197 0.016950 1 -0.2 1 0 0

#has de product incorporated in order to not have problems with the gamma divisions, the halo catalog that we upload has not subhalos and has a minimum mass of 10**10.32, 10**10.5 and 10**11, the halo bias used to compute fsat,Ac,As is computed comparing with particles 0.5%

#unit_001
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 10.900 0.007016 0.016400 1 -0.2 1 0 0  #Mh  no below limit
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 10.900 0.007016 0.016400 1 -0.2 1 0 0  #Mh > 10.32
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 10.900 0.007016 0.016400 1 -0.2 1 0 0  #Mh > 10.5
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 10.500 0.083751 0.006811 1 -0.2 1 0 0  #Mh > 11
10.652000 5.620524e-02 9.201159e-03 #Mh > 11 (fsat = 0.58!)

#inv_phase_001
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 10.889000 0.006856 0.015941 1 -0.2 1 0 0  #Mh no below limit and Mh > 10.32
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 10.892000 0.006899 0.016063 1 -0.2 1 0 0  #Mh > 10.5
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 10.528000 0.083798 0.006994 1 -0.2 1 0 0  #Mh > 11  (for fsat 0.58!)



#21 particles, mass range 10.4 - 15, bias = 1.359, number density = 10n_eboss, order 5 polynomial taking into account weights
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 12.192 0.294535 0 1 0 1 0 0 #poisson, fsat = 0
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 11.821 0.102201 0.064595 1 0 1 0 0 #poisson, fsat = 0.20
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 11.386 0.029494 0.039046 1 0 1 0 0 #poisson, fsat = 0.40
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 10.888 0.006851 0.015899 1 0 1 0 0 #poisson, fsat = 0.60

#21 particles, mass range 10.4 - 15, bias = 1.363, number density = 7n_eboss, order 5 polynomial taking into account weights
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 11.829 0.072890 0.046270 1 0 1 0 0 #poisson, fsat = 0.20
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 11.397 0.021170 0.028150 1 0 1 0 0 #poisson, fsat = 0.40
./HOD_NFW_V14_more_BVG_product_nosubhalos_particles_minimumhalomass 10.901 0.004907 0.011507 1 0 1 0 0 #poisson, fsat = 0.60



