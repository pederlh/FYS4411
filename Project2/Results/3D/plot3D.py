import numpy as np
import os
import matplotlib.pyplot as plt


mean = np.load("mean_3D.npy")
omega = np.load("omega_3D.npy")
std = np.load("std_3D.npy")
print(std)

taut_1_over_omega = np.array([4, 20, 54.7386, 115.299, 208.803])
taut_E = np.array([0.6250, 0.1750, 0.0822, 0.0477, 0.0311 ])*2


omega = 1/omega

plt.errorbar(omega,mean,yerr=std, fmt = "or",capsize=5, elinewidth=1,markeredgewidth=1, label = "Hast")
plt.plot(taut_1_over_omega,taut_E, label = "Taut")
plt.xlabel("Omega")
plt.ylabel("EL")
plt.legend()
plt.show()

#omega, mean, std = zip(*sorted(zip(omega,mean,std)))
