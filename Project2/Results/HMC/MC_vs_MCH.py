import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(font_scale=1.3, rc={'legend.facecolor': 'White', 'legend.framealpha': 0.5, 'lines.markersize':5})

"""
all_files = os.listdir()
H_files =[]
M_files =[]
for file in all_files:
    if "HAST" in file:
        H_files.append(file)

    if "BF" in file:
        M_files.append(file)
H_sort = []
M_sort = []

h = len(H_files)
m = len(M_files)


for i in range(h):
    w1 = H_files[i].split("_")
    for j in range(i+1,h):
        w2 = H_files[j].split("_")
        if w1[4] == w2[4]:
            H_sort.append(H_files[i])
            H_sort.append(H_files[j])


for i in range(m):
    w1 = M_files[i].split("_")
    for j in range(i+1,m):
        w2 = M_files[j].split("_")
        if w1[4] == w2[4]:
            M_sort.append(M_files[i])
            M_sort.append(M_files[j])

H_w = []
H_m = []

M_w = []
M_m = []

j=0
for i in range(0,h,2):
    data =[]

    w = H_sort[i].split("_")
    # Reading data from file1
    v1 = np.load(H_sort[i])
    v2 = np.load(H_sort[i+1])

    vec = np.concatenate((v1,v2))
    mean = np.mean(vec)

    H_w.append(float(w[4]))
    H_m.append(mean)

    j+=1

j=0
for i in range(0,h,2):
    data =[]

    w = M_sort[i].split("_")
    # Reading data from file1
    v1 = np.load(M_sort[i])
    v2 = np.load(M_sort[i+1])

    vec = np.concatenate((v1,v2))
    mean = np.mean(vec)

    M_w.append(float(w[4]))
    M_m.append(mean)

    j+=1


H_w, H_m = zip(*sorted(zip(H_w, H_m)))
M_w, M_m = zip(*sorted(zip(M_w, M_m)))
"""
H_w = np.load("H_w.npy")
H_m = np.load("H_m.npy")
M_w = np.load("M_w.npy")
M_m = np.load("M_m.npy")


H_w = H_w[2:]
H_m = H_m[2:]
M_w = M_w[2:]
M_m = M_m[2:]

Hw = np.load("Hw.npy")
H2 = np.load("H2.npy")
Mw = np.load("Mw.npy")
M2 = np.load("M2.npy")

Hw,H2 = zip(*sorted(zip(Hw,H2)))
Mw,M2 = zip(*sorted(zip(Mw,M2)))

Hw = np.array(Hw)
H2 =np.array(H2)
Mw =np.array(Mw)
M2 =np.array(M2)

Hw = Hw[2:]
H2 = H2[2:]
Mw = Mw[2:]
M2 = M2[2:]

taut_1_over_omega = np.array([4, 20, 54.7386, 115.299, 208.803])
taut_E_half = np.array([0.6250, 0.1750, 0.0822, 0.0477, 0.0311 ])


# plt.figure(1)
# plt.plot(1/H_w,H_m/2,":o" ,label ="Metropolis-Hastings")
# plt.plot(1/M_w,M_m/2,":o" ,label ="Metropolis")
# plt.plot(taut_1_over_omega, taut_E_half,"-x" ,label ="Theoretical value (Taut)")

# plt.xlabel(r"$1\, / \, \omega$")
# plt.ylabel(r"$\langle E_L \rangle\, / \, 2$")
# plt.legend()
# plt.tight_layout()


plt.figure(2)
plt.plot(1/Hw,H2,":o", label ="Metropolis-Hastings")
plt.plot(1/Mw,M2,":o", label ="Metropolis")
plt.plot(taut_1_over_omega, taut_E_half*2,"-*", color="k", alpha=0.7, ms=8,label ="Theoretical value (Taut)")

plt.xlabel(r"$1\, / \, \omega$")
plt.ylabel(r"$\langle E_L \rangle$")
plt.legend()
plt.tight_layout()
plt.savefig("Brute_force_Hastings_Taut_vs_omega.pdf")
plt.show()
