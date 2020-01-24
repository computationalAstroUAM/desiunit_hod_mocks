*dprop.py* Generate new catalogue files with: x, y ,z (Mpc/h), vx, vy, vz (comoving peculiar in km/s), lmass (np.log10(fof_halo_mass)), cen (-1=sat, 0=cen), dvx, dvy, dvz, dx, dy, dz, tag (fof_halo_tag)


# OuterRim boxes ##############
/mnt/lustre/eboss/OuterRim/OuterRim_sim/ascii/OuterRim_STEP266_z0.865/subvols27/OuterRim_STEP266_fofproperties_***.txt

About 68 million haloes, about 6GB

------------------------------
De Santi:
Código de HOD: /users/savila/ELGs_eBOSS/ELGs_LSS_V7/HOD/*c

Para fittear los paámetro de HOD: /users/savila/ELGs_eBOSS/ELGs_LSS_V7/scripts/fit_HOD1.py

# Mocks ###########################

NFW:
x, y ,z (Mpc/h), vx, vy, vz (comoving peculiar in km/s), lmass (np.log10(fof_halo_mass)), cen (-1=sat, 0=cen), tag (fof_halo_tag)

Particles:
x, y ,z (Mpc/h), vx, vy, vz (comoving peculiar in km/s), lmass (np.log10(fof_halo_mass)), cen (-1=sat, 0=cen), dvx, dvy, dvz, dx, dy, dz, tag (fof_halo_tag)

The mocks contain about 2million galaxies.

Path mocks: /mnt/lustre/savila/HOD_NFW/output_V1/galaxies_1000Mpc_V1.0more_NFW_mu*

Lista de mocks que se han utilizado para los plots del paper so far (y que habría que incluir en las gráficas que faltando de PDF, density profile, vel prof. (Pongo sólo el primero: 000)

HODs/fsat/  galaxies_1000Mpc_V1.0more_NFW_mu11.154_Ac0.0193_As0.02990_vfact1.00_beta0.000_K1.00_vt0pm0_mock000.dat
HODs/fsat/  galaxies_1000Mpc_V1.0more_NFW_mu11.820_Ac0.1291_As0.06275_vfact1.00_beta0.000_K1.00_vt0pm0_mock000.dat
HODs/fsat/  galaxies_1000Mpc_V1.0more_NFW_mu12.088_Ac0.1970_As0.00000_vfact1.00_beta0.000_K1.00_vt0pm0_mock000.dat
position_prof/NFW_K/  galaxies_1000Mpc_V1.0more_NFW_mu11.515_Ac0.0537_As0.05301_vfact1.00_beta0.000_K0.40_vt0pm0_mock000.dat
position_prof/NFW_K/  galaxies_1000Mpc_V1.0more_NFW_mu11.515_Ac0.0537_As0.05301_vfact1.00_beta0.000_K0.60_vt0pm0_mock000.dat
position_prof/NFW_K/  galaxies_1000Mpc_V1.0more_NFW_mu11.515_Ac0.0537_As0.05301_vfact1.00_beta0.000_K0.80_vt0pm0_mock000.dat
prob_distribution/neg_bin_beta0100/  galaxies_1000Mpc_V1.0more_NFW_mu11.515_Ac0.0537_As0.05301_vfact1.00_beta0.100_K1.00_vt0pm0_mock000.dat
prob_distribution/neg_bin_beta0200/  galaxies_1000Mpc_V1.0more_NFW_mu11.515_Ac0.0537_As0.05301_vfact1.00_beta0.200_K1.00_vt0pm0_mock000.dat
prob_distribution/next_int/  galaxies_1000Mpc_V1.0more_NFW_mu11.515_Ac0.0537_As0.05301_vfact1.00_beta-2.000_K1.00_vt0pm0_mock000.dat
prob_distribution/poisson/  galaxies_1000Mpc_V1.0more_NFW_mu11.515_Ac0.0537_As0.05301_vfact1.00_beta0.000_K1.00_vt0pm0_mock000.dat. #este es el “default"
velocity_prof/viralth_velbias/galaxies_1000Mpc_V1.0more_NFW_mu11.515_Ac0.0537_As0.05301_vfact0.20_beta0.000_K1.00_vt0pm0_mock000.dat
velocity_prof/viralth_velbias/galaxies_1000Mpc_V1.0more_NFW_mu11.515_Ac0.0537_As0.05301_vfact0.60_beta0.000_K1.00_vt0pm0_mock000.dat
velocity_prof/viralth_velbias/galaxies_1000Mpc_V1.0more_NFW_mu11.515_Ac0.0537_As0.05301_vfact1.40_beta0.000_K1.00_vt0pm0_mock000.dat
velocity_prof/viralth_velbias_vt/galaxies_1000Mpc_V1.0more_NFW_mu11.515_Ac0.0537_As0.05301_vfact0.20_beta0.000_K1.00_vt500pm200_mock000.dat
velocity_prof/viralth_velbias_vt/galaxies_1000Mpc_V1.0more_NFW_mu11.515_Ac0.0537_As0.05301_vfact0.60_beta0.000_K1.00_vt500pm200_mock000.dat
velocity_prof/viralth_velbias_vt/galaxies_1000Mpc_V1.0more_NFW_mu11.515_Ac0.0537_As0.05301_vfact1.40_beta0.000_K1.00_vt500pm200_mock000.dat

Y luego el best fit sería este: 
bestfits/ galaxies_1000Mpc_V1.0more_NFW_mu11.356_Ac0.0342_As0.04229_vfact1.00_beta0.110_K1.00_vt0pm0_mock000.dat

Ahora he creado los que contienen partículas: 

(K)            /mnt/lustre/savila/HOD_part/output_V1/galaxies_1000Mpc_V1.0more_part_mu11.515_Ac0.0054_As0.00530_vfact1.00_beta0.000_K?.?0_vt0pm0_mock???.dat
(Alpha)        /mnt/lustre/savila/HOD_part/output_V1/galaxies_1000Mpc_V1.0more_part_mu11.515_Ac0.0054_As0.00530_vfact?.?0_beta0.000_K1.00_vt500pm200_mock???.dat
(Alpha + vt)   /mnt/lustre/savila/HOD_part/output_V1/galaxies_1000Mpc_V1.0more_part_mu11.515_Ac0.0054_As0.00530_vfact?.?0_beta0.000_K1.00_vt0pm0_mock???.dat 


