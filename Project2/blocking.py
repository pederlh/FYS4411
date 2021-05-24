import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import sys

plt.style.use('seaborn')
sns.set(font_scale=1.4)
""" Script that performs blocking on local energy data sets """

type = sys.argv[1]

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


data = np.array(samples)
print("Files has been read, start blocking...")
(mean, var) = block(data)
std = np.sqrt(var)
print("Statistical values:\n")
print("N=2")
print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))
