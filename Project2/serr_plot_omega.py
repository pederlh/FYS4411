import numpy as np
import matplotlib.pyplot as plt

fz = 14

omegas = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0])
analytical = omegas/2
energies_anasig = np.array([0.0507685005656783,0.100697660495258,0.149912547266176,0.201181188107038,0.249760897657481,0.300380701749821,0.349857579776897,0.400248609121509,0.449938577026867,0.500093005316653,1.00001243917618,1.49976896099163])
energies_noanasig = np.array([0.252463070184053,0.257988334921609,0.162629118370204,0.207559415411015,0.253001021159003,0.339655117821471,0.372382559323085,0.400198808194329,0.45266476900845,0.500047478621961,1.25196355735102,2.50438706844833
])
rel_diff_ana = np.abs((analytical-energies_anasig)/analytical)
rel_diff_noana = np.abs((analytical-energies_noanasig)/analytical)


plt.plot(omegas, analytical, label="Analytical solution")
plt.plot(omegas,energies_anasig, label="Numerical $\sigma = 1/\sqrt{\omega}$")
plt.plot(omegas, energies_noanasig, label="Numerical $\sigma = 1.0$")
plt.title("$E_L$ vs. HO-frequency (1 particle, 1D)",fontsize = fz)
plt.xlabel("$\omega$",fontsize = fz)
plt.ylabel("$E_L$",fontsize = fz)
plt.legend()
plt.show()

plt.title("Relative error in numerical solution",fontsize = fz)
plt.plot(omegas,np.log(rel_diff_ana), label = "$\sigma = 1/\sqrt{\omega}$")
plt.plot(omegas, np.log(rel_diff_noana), label = "$\sigma = 1.0$")
plt.xlabel("$\omega$", fontsize = fz)
plt.ylabel("Relative error", fontsize = fz)
plt.legend()
plt.show()

0.252463070184053,0.257988334921609,0.162629118370204,0.207559415411015,0.253001021159003,0.339655117821471,0.372382559323085,0.400198808194329,0.45266476900845,0.500047478621961,1.25196355735102,2.50438706844833
