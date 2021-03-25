import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
import sys


filename = sys.argv[1]
"""
def Bootstrap(B,filename):

    energies = np.loadtxt(filename)
    energies = np.delete(energies,0)   #delete first line as it is the time used, not data
    n = len(energies)

    mean = np.mean(energies)
    energy_squared = energies*energies;
    mean_2 = np.mean(energy_squared)

    reg_M = mean
    reg_V = mean_2 - mean*mean

    k = 10000;

    boot_M = np.zeros(B)
    boot_V = np.zeros(B)
    for b in range(B):
        idx = np.random.randint(0, n, size = n)
        energies_copy = energies[idx]
        mean_b = np.mean(energies_copy)
        energy_squared_copy = energies_copy*energies_copy
        mean_b_2 = np.mean(energy_squared_copy)
        boot_M[b] = mean_b
        boot_V[b]= mean_b_2 - mean_b*mean_b

    boot_Mean = np.mean(boot_M)
    boot_Var = np.mean(boot_V)

    hist, bins = np.histogram(boot_M, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.show()

    return boot_Mean, boot_Var, reg_M, reg_V

"""
""" The blocking code, based on the article of Marius Jonsson  """
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


x = np.loadtxt(filename)
x = np.delete(x,0)   #delete first line as it is the time used, not data
(mean, var) = block(x)
std = np.sqrt(var)
data ={'Mean':[mean], 'STDev':[std]}
frame = pd.DataFrame(data,index=['Values'])
print(frame)


#Bootstrap(1000,filename)
