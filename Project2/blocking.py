import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import sys

plt.style.use('seaborn')
sns.set(font_scale=1.4)
""" Script that performs blocking on local energy data sets """

#type = sys.argv[1]

def block(x):
    # preliminaries
    n = len(x)
    d = int(np.log2(n))
    s, gamma = np.zeros(d), np.zeros(d)
    mu = np.mean(x)

    # estimate the auto-covariance and variances
    # for each blocking transformation
    for i in np.arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*np.sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = np.var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])

    # generate the test observator M_k from the theorem
    M = (np.cumsum( ((gamma/s)**2*2**np.arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =np.array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in np.arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    return mu, s[k]/2**(d-k)

"""
data = np.load("EnergySamples_GD_N_2_D_2_H_2_eta_0.010000_MC_65536_sigma_1.000000_ID_0_interaction_1_omega_1.000000_its_70_.npy")
print("Files has been read, start blocking...")
(mean, var) = block(data)
std = np.sqrt(var)
print("Statistical values:\n")
print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))

"""


"""
samples = []

if type == "par":
    num_threads = 2;
    for i in range(num_threads):
        filename = "EnergySamples_ID_" + str(i) + ".txt"
        with open(filename) as infile:
            lines = infile.readlines()
            for line in lines:
                word = line.split()
                samples.append(float(word[0]))

if type == "ser":
    filename = "EnergySamples.txt"
    with open(filename) as infile:
        lines = infile.readlines()
        for line in lines:
            word = line.split()
            samples.append(float(word[0]))



mc_test =[]
means_test = []
stds_test = []

mc_NC_test =[]
means_NC_test = []
stds_NC_test = []

mc_chi =[]
means_chi = []
stds_chi = []

mc_NC_chi =[]
means_NC_chi = []
stds_NC_chi = []

files = os.listdir()
for samples in files:
    if "TEST_MC" in samples:
        if "NOTCONVERGED" in samples:
            words = samples.split("_")
            mc_NC_test.append(int(words[2]))
            data = np.load(samples)
            print("Files has been read, start blocking...")
            (mean, var) = block(data)
            std = np.sqrt(var)
            print("Statistical values:\n")
            print(samples)
            print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))
            means_NC_test.append(mean)
            stds_NC_test.append(std)
        else:
            words = samples.split("_")
            mc_test.append(int(words[2]))
            data = np.load(samples)
            print("Files has been read, start blocking...")
            (mean, var) = block(data)
            std = np.sqrt(var)
            print("Statistical values:\n")
            print(samples)
            print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))
            means_test.append(mean)
            stds_test.append(std)
    elif "CHI_MC" in samples:
        if "NOTCONVERGED" in samples:
            words = samples.split("_")
            mc_NC_chi.append(int(words[2]))
            data = np.load(samples)
            print("Files has been read, start blocking...")
            (mean, var) = block(data)
            std = np.sqrt(var)
            print("Statistical values:\n")
            print(samples)
            print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))
            means_NC_chi.append(mean)
            stds_NC_chi.append(std)
        else:
            words = samples.split("_")
            mc_chi.append(int(words[2]))
            data = np.load(samples)
            print("Files has been read, start blocking...")
            (mean, var) = block(data)
            std = np.sqrt(var)
            print("Statistical values:\n")
            print(samples)
            print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))
            means_chi.append(mean)
            stds_chi.append(std)

mc_test = np.array(mc_test)
means_test = np.array(means_test)
stds_test = np.array(stds_test)

np.savetxt("mc_test.txt",mc_test)
np.savetxt("means_test.txt",means_test)
np.savetxt("stds_test.txt",stds_test)

mc_chi = np.array(mc_chi)
means_chi = np.array(means_chi)
stds_chi = np.array(stds_chi)

np.savetxt("mc_chi.txt",mc_chi)
np.savetxt("means_chi.txt",means_chi)
np.savetxt("stds_chi.txt",stds_chi)



mc_NC_test = np.array(mc_NC_test)
means_NC_test = np.array(means_NC_test)
stds_NC_test = np.array(stds_NC_test)

np.savetxt("mc_NC_test.txt",mc_NC_test)
np.savetxt("means_NC_test.txt",means_NC_test)
np.savetxt("stds_NC_test.txt",stds_NC_test)

mc_chi = np.array(mc_chi)
means_chi = np.array(means_chi)
stds_chi = np.array(stds_chi)

np.savetxt("mc_NC_chi.txt",mc_NC_chi)
np.savetxt("means_NC_chi.txt",means_NC_chi)
np.savetxt("stds_NC_chi.txt",stds_NC_chi)


H = []
Hw = []

M = []
Mw = []

all = os.listdir()
all_new = []
for file in all:
    if "HAST_" in file:
        w = file.split("_")
        H.append(file)
        Hw.append(float(w[4]))

for file in all:
    if "HAST_" not in file:
        all_new.append(file)

all = all_new


for file in all:
    if "MC_" in file:
        w = file.split("_")
        M.append(file)
        Mw.append(float(w[3]))


H2 =[]
M2 =[]


for h in H:
    data = np.load(h)
    (mean, var) = block(data)
    H2.append(mean)

for m in M:
    data = np.load(m)
    (mean, var) = block(data)
    M2.append(mean)

H2 = np.array(H2)
M2 = np.array(M2)
Hw = np.array(Hw)
Mw = np.array(Mw)

np.save("H2",H2)
np.save("M2",M2)
np.save("Hw",Hw)
np.save("Mw",Mw)

"""
