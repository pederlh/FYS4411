import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import sys

plt.style.use('seaborn')
sns.set(font_scale=1.4)
""" Script that performs blocking on local energy data sets """


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

files_rand = os.listdir()
files = []
for f in files_rand:
    if not f == "block_etas.py" and not f == "plot_eta.py":
        files.append(f)


eta = []
its = np.zeros(len(files))
means = np.zeros(len(files))
stds = np.zeros(len(files))

i = 0
for fil in files:
    words = fil.split("_")
    eta.append(words[1])
    its[i] = float(words[-2])
    samples = np.loadtxt(fil)
    data = np.array(samples)
    print("Files has been read, start blocking...")
    (mean, var) = block(data)
    std = np.sqrt(var)
    print("Statistical values:\n")
    print("N=2")
    print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))

    means[i] = mean
    stds[i] = std
    i+=1

e_save = []
for e in eta:
    elem = ""
    words = list(e)
    n = len(words)
    elem += words[0] + "."
    for i in range(1,n):
        elem += words[i]
    e_save.append(float(elem))

e_save = np.array(e_save)

np.save("etas.npy",e_save)
np.save("its.npy",its)
np.save("means.npy",means)
np.save("stds.npy",stds)
