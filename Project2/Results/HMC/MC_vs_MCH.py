import numpy as np
import os
import matplotlib.pyplot as plt
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


plt.plot(1/H_w,H_m/2,"-o" ,label ="Hastings")
plt.plot(1/M_w,M_m/2,"-o" ,label ="MC")
plt.xlabel("1/w")
plt.ylabel("El/2")
plt.legend()
plt.show()
