import numpy as np
import matplotlib.pyplot as plt
import math

mc_chi = np.loadtxt("mc_chi.txt")
means_chi = np.loadtxt("means_chi.txt")
stds_chi = np.loadtxt("stds_chi.txt")

mc_test = np.loadtxt("mc_test.txt")
means_test = np.loadtxt("means_test.txt")
stds_test = np.loadtxt("stds_test.txt")

mc_NC_chi = np.loadtxt("mc_NC_chi.txt")
means_NC_chi = np.loadtxt("means_NC_chi.txt")
stds_NC_chi = np.loadtxt("stds_NC_chi.txt")

mc_NC_test = np.loadtxt("mc_NC_test.txt")
means_NC_test = np.loadtxt("means_NC_test.txt")
stds_NC_test = np.loadtxt("stds_NC_test.txt")

nc = np.size(mc_chi)
nt = np.size(mc_test)
nnc = np.size(mc_NC_chi)
nnt = len(mc_NC_test)

mc_c = np.zeros(nc)
means_c =  np.zeros(nc)
stds_c =  np.zeros(nc)

mc_t =  np.zeros(nt)
means_t =  np.zeros(nt)
stds_t =  np.zeros(nt)

mc_NC_c = np.zeros(nnc)
means_NC_c =  np.zeros(nnc)
stds_NC_c =  np.zeros(nnc)

mc_NC_t =  np.zeros(nnt)
means_NC_t =  np.zeros(nnt)
stds_NC_t =  np.zeros(nnt)

if nc == 0:
    print("nc = 0!!!")

if nnc == 0:
    print("ncc = 0!!!")

if nt == 0:
    print("nt = 0!!!")

if nnt == 0:
    print("nt = 0!!!")

if nc == 1:
    mc_c[0] = math.log2(int(mc_chi)*4)
    means_c[0] = float(means_chi)
    stds_c[0] = float(stds_chi)
else:
    for i in range(nc):
        mc_c[i] = math.log2(int(mc_chi[i])*4)
        means_c[i] = float(means_chi[i])
        stds_c[i] = float(stds_chi[i])

if nt == 1:
    mc_t[0] = math.log2(int(mc_test)*4)
    means_t[0] = float(means_test)
    stds_t[0] = float(stds_test)
else:
    for i in range(nt):
        mc_t[i] = math.log2(int(mc_test[i])*4)
        means_t[i] = float(means_test[i])
        stds_t[i] = float(stds_test[i])

if nnc == 1:
    mc_NC_c[0] = math.log2(int(mc_NC_chi)*4)
    means_NC_c[0] = float(means_NC_chi)
    stds_NC_c[0] = float(stds_NC_chi)
else:
    for i in range(nnc):
        mc_NC_c[i] = math.log2(int(mc_NC_chi[i])*4)
        means_NC_c[i] = float(means_NC_chi[i])
        stds_NC_c[i] = float(stds_NC_chi[i])


if nnt == 1:
    mc_NC_t[0] = math.log2(int(mc_NC_test)*4)
    means_NC_t[0] = float(means_NC_test)
    stds_NC_t[0] = float(stds_NC_test)
else:
    for i in range(nnt):
        mc_NC_t[i] = math.log2(int(mc_NC_test[i])*4)
        means_NC_t[i] = float(means_NC_test[i])
        stds_NC_t[i] = float(stds_NC_test[i])


plt.errorbar(mc_c,means_c,yerr=stds_c, fmt = "or",capsize=5, elinewidth=1,markeredgewidth=1, label = "Chi")
plt.errorbar(mc_t,means_t,yerr=stds_t, fmt = "ob",capsize=5, elinewidth=1,markeredgewidth=1, label = "Test")
plt.plot(mc_NC_c,means_NC_c,"rx", label = "Chi not converged")
plt.plot(mc_NC_t,means_NC_t,"bx", label = "Test not converged")
plt.axhline(y=3.0, color='k', linestyle='-', label = "Analytical")
plt.legend()
plt.xlabel("log2(MC-cycles)")
plt.ylabel("<EL>")
plt.show()
