import numpy as np
import matplotlib.pyplot as plt


filename = "One_body_density_N_10_stringID_0_alpha_0.500000.txt"
OBDs = np.loadtxt(filename)

filename_I = "Interaction_One_body_density_N_10_stringID_0_alpha_0.498882.txt"
OBDs_I = np.loadtxt(filename_I)

r = np.linspace(0,5,len(OBDs))

plt.plot(r,OBDs, label = "No interaction")
plt.plot(r,OBDs_I, label = "Interaction")
plt.legend()
plt.show()
