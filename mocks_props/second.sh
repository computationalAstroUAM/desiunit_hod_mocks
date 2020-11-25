#!/bin/bash


echo Calculates the  bias, number density of galaxies and the fraction of satellites of each one of the galaxy mocks generated.


#kmax = 0.3 polynomial grade 5


#nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/galaxies_1000Mpc_V1.4_NFW_mu10.782_Ac0.0064_As0.01253_vfact1.00_beta-2.000_K1.00_vt0pm0_mock.dat > n_b_fsat_mocks/b_fsat_n_mock1.txt &

#nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/galaxies_1000Mpc_V1.4_NFW_mu10.907_Ac0.0090_As0.01575_vfact1.00_beta0.000_K0.25_vt0pm0_mock.dat > n_b_fsat_mocks/b_fsat_n_mock2.txt &

#nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/galaxies_1000Mpc_V1.4_NFW_mu11.621_Ac0.0686_As0.04434_vfact1.50_beta0.000_K1.00_vt0pm0_mock.dat > n_b_fsat_mocks/b_fsat_n_mock3.txt &

#nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/galaxies_1000Mpc_V1.4_NFW_mu11.621_Ac0.0686_As0.04434_vfact1.00_beta0.000_K1.00_vt-500pm200_mock.dat > n_b_fsat_mocks/b_fsat_n_mock4.txt &

#nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/galaxies_1000Mpc_V1.4_NFW_mu10.983_Ac0.0112_As0.01807_vfact1.00_beta0.100_K0.15_vt0pm0_mock.dat > n_b_fsat_mocks/b_fsat_n_mock9.txt &


#kmax = 0.1  polynomial grade 6
# n = 25e-4

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0089_As0.01954_vfact1.00_beta-2.000_K1.00_vt0pm0_mock.dat > n_b_fsat_mocks/k_01_mock1.txt &

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.079_Ac0.0129_As0.02471_vfact1.00_beta0.000_K0.25_vt0pm0_mock.dat >  n_b_fsat_mocks/k_01_mock2.txt &

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0860_As0.05916_vfact1.50_beta0.000_K1.00_vt0pm0_mock.dat > n_b_fsat_mocks/k_01_mock3.txt &

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0860_As0.05916_vfact1.00_beta0.000_K1.00_vt500pm200_mock.dat > n_b_fsat_mocks/k_01_mock4.txt &

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_25e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.153_Ac0.0160_As0.02829_vfact1.00_beta0.100_K0.15_vt0pm0_mock.dat > n_b_fsat_mocks/k_01_mock9.txt &


# n = 20e-4

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0071_As0.01563_vfact1.00_beta-2.000_K1.00_vt0pm0_mock.dat > n_b_fsat_mocks/k_01_mock1_n_20e-4.txt &

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.079_Ac0.0103_As0.01977_vfact1.00_beta0.000_K0.25_vt0pm0_mock.dat >  n_b_fsat_mocks/k_01_mock2_n_20e-4.txt &

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0688_As0.04733_vfact1.50_beta0.000_K1.00_vt0pm0_mock.dat > n_b_fsat_mocks/k_01_mock3_n_20e-4.txt &

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0688_As0.04733_vfact1.00_beta0.000_K1.00_vt500pm200_mock.dat > n_b_fsat_mocks/k_01_mock4_n_20e-4.txt &

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_20e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.153_Ac0.0128_As0.02263_vfact1.00_beta0.100_K0.15_vt0pm0_mock.dat > n_b_fsat_mocks/k_01_mock9_n_20e-4.txt &


# n = 5e-4

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu10.954_Ac0.0018_As0.00391_vfact1.00_beta-2.000_K1.00_vt0pm0_mock.dat > n_b_fsat_mocks/k_01_mock1_n_5e-4.txt &

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.079_Ac0.0026_As0.00494_vfact1.00_beta0.000_K0.25_vt0pm0_mock.dat > n_b_fsat_mocks/k_01_mock2_n_5e-4.txt &

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0172_As0.01183_vfact1.50_beta0.000_K1.00_vt0pm0_mock.dat > n_b_fsat_mocks/k_01_mock3_n_5e-4.txt &

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.723_Ac0.0172_As0.01183_vfact1.00_beta0.000_K1.00_vt500pm200_mock.dat > n_b_fsat_mocks/k_01_mock4_n_5e-4.txt &

nohup python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/n_5e-4/galaxies_1000Mpc_V1.4more_NFW_mu11.153_Ac0.0032_As0.00566_vfact1.00_beta0.100_K0.15_vt0pm0_mock.dat > n_b_fsat_mocks/k_01_mock9_n_5e-4.txt &


#prueba con fsat=0

python3 bias_fsat_numdensity_galaxies.py hod_mocks/output_V1/galaxies_1000Mpc_V1.4more_NFW_mu12.079_Ac0.2437_As0.00000_vfact1.00_beta-2.000_K1.00_vt0pm0_mock.dat

