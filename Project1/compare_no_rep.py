import numpy as np
import matplotlib.pyplot as plt

N = 24
path = "./Results/1b_low_MC/spherical_HO_"
dims = ["1D_", "2D_", "3D_"]
methods = ["ana_","num_"]
particles = [1,10,100,500]
"""
num_e_1d = []
num_e_2d = []
num_e_3d = []

ana_e_1d = []
ana_e_2d = []
ana_e_3d = []

num_var_1d = []
num_var_2d = []
num_var_3d = []

ana_var_1d = []
ana_var_2d = []
ana_var_3d = []

t_num_1d=[]
t_num_2d=[]
t_num_3d=[]

t_ana_1d=[]
t_ana_2d=[]
t_ana_3d=[]

#1D-ana
for p in particles:
    infilename = path + dims[0] + methods[0] + "N_" + str(p) + "_MC_100000.txt"
    with open(infilename, "r") as infile:
        for i, line in enumerate(infile):
            if i == 9:
                val  = line.split()
                ana_e_1d.append(float(val[1]))
                ana_var_1d.append(float(val[2]))
            if i == 12:
                val  = line.split()
                t_ana_1d.append(float(val[1]))

#1D-num
for p in particles:
    infilename = path + dims[0] + methods[1] + "N_" + str(p) + "_MC_100000.txt"
    with open(infilename, "r") as infile:
        for i, line in enumerate(infile):
            if i == 9:
                val  = line.split()
                num_e_1d.append(float(val[1]))
                num_var_1d.append(float(val[2]))
            if i == 12:
                val  = line.split()
                t_num_1d.append(float(val[1]))

#2D-ana
for p in particles:
    infilename = path + dims[1] + methods[0] + "N_" + str(p) + "_MC_100000.txt"
    with open(infilename, "r") as infile:
        for i, line in enumerate(infile):
            if i == 9:
                val  = line.split()
                ana_e_2d.append(float(val[1]))
                ana_var_2d.append(float(val[2]))
            if i == 12:
                val  = line.split()
                t_ana_2d.append(float(val[1]))

#2D-num
for p in particles:
    infilename = path + dims[1] + methods[1] + "N_" + str(p) + "_MC_100000.txt"
    with open(infilename, "r") as infile:
        for i, line in enumerate(infile):
            if i == 9:
                val  = line.split()
                num_e_2d.append(float(val[1]))
                num_var_2d.append(float(val[2]))
            if i == 12:
                val  = line.split()
                t_num_2d.append(float(val[1]))

#3D-ana
for p in particles:
    infilename = path + dims[2] + methods[0] + "N_" + str(p) + "_MC_100000.txt"
    with open(infilename, "r") as infile:
        for i, line in enumerate(infile):
            if i == 9:
                val  = line.split()
                ana_e_3d.append(float(val[1]))
                ana_var_3d.append(float(val[2]))
            if i == 12:
                val  = line.split()
                t_ana_3d.append(float(val[1]))

#3D-num
for p in particles:
    infilename = path + dims[2] + methods[1] + "N_" + str(p) + "_MC_100000.txt"
    with open(infilename, "r") as infile:
        for i, line in enumerate(infile):
            if i == 9:
                val  = line.split()
                num_e_3d.append(float(val[1]))
                num_var_3d.append(float(val[2]))
            if i == 12:
                val  = line.split()
                t_num_3d.append(float(val[1]))


particles =np.array(particles)
e_1d = (1/2)*particles
e_2d = 1*particles
e_3d = (3/2)*particles

num_e_1d = abs((np.array(num_e_1d)-e_1d)/e_1d)
num_e_2d = abs((np.array(num_e_2d)-e_2d)/e_2d)
num_e_3d = abs((np.array(num_e_3d)-e_3d)/e_3d)

ana_e_1d = abs((np.array(ana_e_1d)-e_1d)/e_1d)
ana_e_2d = abs((np.array(ana_e_2d)-e_2d)/e_2d)
ana_e_3d = abs((np.array(ana_e_3d)-e_3d)/e_3d)

num_var_1d = np.array(num_var_1d)
num_var_2d = np.array(num_var_2d)
num_var_3d = np.array(num_var_3d)

ana_var_1d = np.array(ana_var_1d)
ana_var_2d = np.array(ana_var_2d)
ana_var_3d = np.array(ana_var_3d)

t_num_1d = np.array(t_num_1d)
t_num_2d = np.array(t_num_2d)
t_num_3d = np.array(t_num_3d)

t_ana_1d = np.array(t_ana_1d)
t_ana_2d = np.array(t_ana_2d)
t_ana_3d = np.array(t_ana_3d)

t_rel_1d = t_num_1d/t_ana_1d
t_rel_2d = t_num_2d/t_ana_2d
t_rel_3d = t_num_3d/t_ana_3d


#Plot 1D energy
plt.plot(particles,num_e_1d,"*", label= "Brute force")
plt.plot(particles,ana_e_1d,"*", label= "Analytical")
plt.legend()
plt.xlabel("Number of particles")
plt.ylabel("Relative error")
plt.show()


#Plot 2D energy
plt.plot(particles,num_e_2d,"*", label= "Brute force")
plt.plot(particles,ana_e_2d,"*", label= "Analytical")
plt.legend()
plt.xlabel("Number of particles")
plt.ylabel("Relative error")
plt.show()


#Plot 3D energy
plt.plot(particles,num_e_3d,"*", label= "Brute force")
plt.plot(particles,ana_e_3d,"*", label= "Analytical")
plt.legend()
plt.xlabel("Number of particles")
plt.ylabel("Relative error")
plt.show()

#Plot 1D time
plt.plot(particles,t_rel_1d ,label= "1D")
plt.plot(particles,t_rel_2d ,label= "2D")
plt.plot(particles,t_rel_3d ,label= "3D")
plt.legend()
plt.xlabel("Number of particles")
plt.ylabel(r"$\frac{t_{num}}{t_{ana}}$", fontsize = 20)
plt.show()


#Plot 2D time
plt.plot(particles,t_rel_2d,"*", label= "Time difference between brute force and analytical")
plt.legend()
plt.xlabel("Number of particles")
plt.ylabel(r"$\frac{t_{num}}{t_{ana}}$", fontsize = 20)
plt.show()

#Plot 3D time
plt.plot(particles,t_rel_3d,"*", label= "Time difference between brute force and analytical")
plt.legend()
plt.xlabel("Number of particles")
plt.ylabel(r"$\frac{t_{num}}{t_{ana}}$", fontsize = 20)
plt.show()
"""



#Relativ feil i energi som funksjon av alpha med N=10 particler i 3D
alphas = []
energies = []
variance = []

filename = "./Results/1b_used/spherical_HO_3D_num_N_10_MC_10000000.txt"
with open(filename,"r") as file:
    lines = file.readlines()[1:-2]
    for line in lines:
        val = line.split()
        alphas.append(float(val[0]))
        energies.append(float(val[1]))
        variance.append(float(val[2]))

alphas = np.array(alphas)
energies = np.array(energies)
variance = np.array(variance)

e_ana = (3/2)*10

rel_e = abs((energies-e_ana)/e_ana)
min_idx = np.where(rel_e == np.min(rel_e))
plt.plot(alphas, rel_e, label = "Relative error in energy")
plt.plot(alphas[min_idx], rel_e[min_idx], "*", label = "Minimum relative error")
plt.legend()
plt.xlabel("Alphas")
plt.ylabel("Relative error")
plt.show()

min_idx = np.where(variance == np.min(variance))
plt.plot(alphas,variance, label = "Variance")
plt.plot(alphas[min_idx],variance[min_idx], "*",label = "Minimum variance")
plt.legend()
plt.xlabel("Alphas")
plt.ylabel("Variance")
plt.show()
