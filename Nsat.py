import numpy as np


mock = np.loadtxt('../DESI_outputs/output_V1/mocks_fig7_Av2020/galaxies_1000Mpc_V1.4more_NFW_mu10.500_Ac0.0064_As0.00710_vfact1.00_beta0.000_K0.25_vt0pm0_BVG.dat')
Nsat_mock = mock[:,7]


Max = 6
Nsat = np.zeros((Max-1))
for i in Nsat_mock:
    for j in range(1,Max):
        if i == j:
            Nsat[j-1] += 1/j

sel = np.where(Nsat_mock != 0)
Total_sats = len(Nsat_mock[sel])

Nsat = Nsat/Total_sats

print(Nsat)

