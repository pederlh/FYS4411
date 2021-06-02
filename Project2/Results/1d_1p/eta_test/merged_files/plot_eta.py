import numpy as np
import matplotlib.pyplot as plt


etas = np.load("etas.npy")
its = np.load("its.npy")
means = np.load("means.npy")
stds = np.load("stds.npy")
print(means)

analytical = np.ones(len(means))*0.5

etas,its,means,stds = zip(*sorted(zip(etas,its,means,stds)))

plt.plot(etas,analytical, label = "Analytical")
plt.errorbar(etas,means,yerr=stds, fmt = "or",capsize=5, elinewidth=1,markeredgewidth=1, label = "Numerical")
plt.ylim(0.4998,0.5002)
plt.legend()
plt.show()

plt.plot(etas,its)
plt.show()
