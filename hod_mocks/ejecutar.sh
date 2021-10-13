#!/bin/bash





# Para n = 0.0025 y kmax = 0.3  #ESTOS NO SON CORRECTOS
./HOD_NFW_V14 10.782 0.00637 0.01253 1 -2 1 0 0 
./HOD_NFW_V14 10.907 0.009017 0.01575 1 0 0.25 0 0  
./HOD_NFW_V14 11.621 0.06861 0.04434 1.5 0 1 0 0
./HOD_NFW_V14 11.621 0.06861 0.04434 1 0 1 500 200 
./HOD_NFW_V14 10.983 0.01117 0.01807 1 0.10 0.15 0 0







#LOS BUENOS
# Para n = 0.0025 y kmax =0.1 and polynomial grade 6

./HOD_NFW_V14_more 10.954 0.008899 0.01954 1 -2 1 0 0
./HOD_NFW_V14_more 11.079 0.01286 0.02471 1 0 0.25 0 0
./HOD_NFW_V14_more 11.723 0.08596 0.05916 1.5 0 1 0 0
./HOD_NFW_V14_more 11.723 0.08596 0.05916 1 0 1 500 200 
./HOD_NFW_V14_more 11.153 0.01603 0.02829 1 0.10 0.15 0 0

              # COMPROBACION NEAREST INTEGER, BINOMIAL, POISSON, NEGATIVE BINOMIAL, we only change beta and beta2 parameters
nohup ./HOD_NFW_V14_more_binomial 10.954 0.008899 0.01954 1 -2 1 0 0 0   #we have one number more, the last one, which represents beta2, a parameter from the binomial distribution
nohup ./HOD_NFW_V14_more_binomial 10.954 0.008899 0.01954 1 0 1 0 0 0.1
nohup ./HOD_NFW_V14_more_binomial 10.954 0.008899 0.01954 1 0 1 0 0 0
nohup ./HOD_NFW_V14_more_binomial 10.954 0.008899 0.01954 1 0.1 1 0 0 0


              # ESTO ES PARA COMPROBAR ONA COSA CON ALPHA_V, VAN APARTE
./HOD_NFW_V14_more 11.723 0.08596 0.05916 1 0 1 0 0
./HOD_NFW_V14_more 11.723 0.08596 0.05916 0.2 0 1 0 0

# Para n = 0.0020 y kmax = 0.1 and polynomial grade 6

./HOD_NFW_V14_more 10.954 0.007119 0.01563 1 -2 1 0 0
./HOD_NFW_V14_more 11.079 0.01029 0.01977 1 0 0.25 0 0
./HOD_NFW_V14_more 11.723 0.06877 0.04733 1.5 0 1 0 0
./HOD_NFW_V14_more 11.723 0.06877 0.04733 1 0 1 500 200
./HOD_NFW_V14_more 11.153 0.01283 0.02263 1 0.10 0.15 0 0

# Para n = 0.0005 y kmax = 0.1 and polynomial grade 6

./HOD_NFW_V14_more 10.954 0.001780 0.003909 1 -2 1 0 0
./HOD_NFW_V14_more 11.079 0.002571 0.004943 1 0 0.25 0 0
./HOD_NFW_V14_more 11.723 0.01719 0.01183 1.5 0 1 0 0
./HOD_NFW_V14_more 11.723 0.01719 0.01183 1 0 1 500 200
./HOD_NFW_V14_more 11.153 0.003207 0.005658 1 0.10 0.15 0 0

# A modo de comprobacion utilizao un model con fsat=0

./HOD_NFW_V14_more 12.079 0.2437 0.0 1 -2 1 0 0



