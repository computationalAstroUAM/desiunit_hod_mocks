*write_input* This code writes the input for cute (the 2PCF code used here to estimate the bias). This can be run using qsub.sh + run.sh.

*check_cute.py* to check CUTE's outputs and plot them.

*meanxi_bias.py* to calculate the mean two point correlation function from all the boxes, for the statistical errors.

*test_integration.py* to compare integration methods.

-------------Santi Sept 2019
PK 

/users/savila/ELGs_eBOSS/ELGs_LSS_V5/data_raw/OuterRim_z0.865_new_xi_from_T.dat 

# cosmological parameters
omega_M= 0.2648
omega_L= 0.7352
w= -1

# binning
log_bin= 1
#log_bin= 0
#n_logint= 8
n_logint= 10

#dim1_max= 200.
#dim1_nbin= 40

#before, esto estaba usando hasta ahora
#dim1_max= 100.
#dim1_nbin= 18
#Esto  es lo que estoy usando últimamente, porque es mas fácil par comparar con lo que me pasa faizan:
dim1_max= 177.827941
dim1_nbin= 33


Para qué escalas usar, creo que lo mejor es ver hasta donde se ajusta bien el fit en los halos, mirando un bin que tenga bastantes objetos. 

Pero a ojo diría que de 10 o 15 hasta 60 o 70 MPc/h

-----
Aquí están los 'mocks oficiales’, que usó Arnaud (y alguno más), por si quieres mirar lo de las velocidades etc 


/mnt/lustre/eboss/ELGs_mocks_OuterRim/ELG_V4_mocks/mocks
-------------------------------------------------------------


