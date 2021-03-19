import numpy as np
import matplotlib.pyplot as plt


filename = "Interaction_One_body_density_N_10_stringID_0_alpha_0.497933.txt"
OBDs = np.loadtxt(filename)

filename_I = "Interaction_One_body_density_N_10_stringID_0_alpha_0.497818.txt"
OBDs_I = np.loadtxt(filename_I)

r = np.linspace(0,5,len(OBDs))

plt.plot(r,OBDs, label = "HEAVY interaction")
plt.plot(r,OBDs_I, label = "Interaction")
plt.legend()
plt.show()
