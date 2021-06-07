import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(font_scale=1.3, rc={'legend.facecolor': 'White', 'legend.framealpha': 0.5, 'lines.markersize':5})

mean = np.load("mean_3D.npy")
omega = np.load("omega_3D.npy")
std = np.load("std_3D.npy")
print("omega: ", omega)
print("omega inv.: ", np.reciprocal(omega))


taut_1_over_omega = np.array([4, 20, 54.7386, 115.299, 208.803])
taut_E = np.array([0.6250, 0.1750, 0.0822, 0.0477, 0.0311 ])*2


omega = 1/omega

plt.errorbar(omega,mean,yerr=std, fmt = "or",capsize=5, elinewidth=1,markeredgewidth=1, label = "Numerical")
plt.plot(taut_1_over_omega,taut_E, "-o", label = "Analytical (Taut)")
plt.xlabel(r"$1\,/\,\omega$")
plt.ylabel(r"$\langle E_L \rangle$")
plt.legend()
plt.tight_layout()
plt.savefig("3D_interaction_vs_omega.pdf")
plt.show()

#omega, mean, std = zip(*sorted(zip(omega,mean,std)))
