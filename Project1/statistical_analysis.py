import numpy as np
import matplotlib.pyplot as plt

"""Requires that you have run Solver w/GD to produce data """


def Bootstrap_Local_Energy(B,filename):

    energies = []

    with open(filename,"r") as infile:
        lines = infile.readlines()
        for line in lines:
            vals = line.split()
            energies.append(float(vals[0]))


    energies = np.array(energies)


    mean = np.mean(energies)
    energy_squared = energies*energies;
    mean_2 = np.mean(energy_squared)

    reg_M = mean
    reg_V = mean_2 - mean*mean

    #reg_M =0
    #reg_V = 0
    k = 10000;

    boot_M = np.zeros(k)
    boot_V = np.zeros(k)
    for b in range(k):
        idx = np.random.randint(0, len(energies), size = B);
        energies_copy = energies[idx]
        mean = np.mean(energies_copy)
        energy_squared = energies_copy*energies_copy;
        mean_2 = np.mean(energy_squared)
        boot_M[b] = mean
        boot_V[b]= mean_2 - mean*mean

    boot_M = np.mean(boot_M)
    boot_V = np.mean(boot_V)




    return boot_M, boot_V, reg_M, reg_V


N = [10]
alphas = []
with open("alpha_values_GD.txt", "r") as afile:
    lines = afile.readlines()
    for line in lines:
        vals = line.split()
        alphas.append(float(vals[0]))

alphas = np.array(alphas)

local_energy_b = np.zeros(len(alphas))
variance_b = np.zeros(len(alphas))

local_energy = np.zeros(len(alphas))
variance = np.zeros(len(alphas))

N = 10
for i in range(len(alphas)):
    readfile = str(N)+"_part_alpha_"+ str(alphas[i]) + "_E_L_samples.txt"
    local_energy_b[i], variance_b[i], local_energy[i], variance[i] = Bootstrap_Local_Energy(9800, readfile);




expectation_E = 15   #When N = 10
rel_error = np.zeros(len(alphas))
for i in range(len(alphas)):
    rel_error[i] = np.abs((expectation_E-local_energy_b[i])/expectation_E)



fig, (ax1, ax2) = plt.subplots(1, 2)

#ax1.plot(alphas, rel_error, "*", label = "Relative error from expectation")
ax1.plot(alphas,local_energy_b,"*", label="Bootstrapped values")
ax1.plot(alphas,local_energy,"*", label = "non B")
ax1.set_title("Local Energy")
ax1.set_xlabel("Alpha")
ax1.set_ylabel("Local energy")
ax1.legend()
ax1.invert_xaxis()

ax2.plot(alphas, variance_b,"*", label="Bootstrapped values")
ax2.plot(alphas, variance,"*", label = "non B")
ax2.set_title("Variance")
ax2.invert_xaxis()
ax2.set_xlabel("Alpha")
ax2.legend()
ax2.set_ylabel("Variance")

fig.tight_layout()
plt.show()
