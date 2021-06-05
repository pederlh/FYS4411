#Python script for running a test instance of the RBM class
import numpy as np
import os, sys


print("Note: This program will produce data and perform bootstrapping on samples before \ndeleting the dataset.")

""" Compilation """

mac_linux_prompt = input("Run on macOS ('m') or Linux ('l')? \n")
os.chdir("./Src/")
if mac_linux_prompt == 'l':
    os.system("g++  -o main.out BoltzmannMachine.cpp main.cpp -std=c++11 -fopenmp -larmadillo")
elif mac_linux_prompt == 'm':
    os.system("c++  -o  main.out -O3 BoltzmannMachine.cpp main.cpp -std=c++11 -larmadillo -lomp -Xpreprocessor -fopenmp")
else:
    print("Input not recognized. Aborting.")
    sys.exit(1)


""" Choose task """
task_prompt = input("Which task to run? Choices are '1' (1p, 1D), '2' (2p, 2D)\n")
if task_prompt not in ['1' , '2']:
    print("Input not recognized. Aborting.")
    sys.exit(1)

interaction = 0
N_particles = 0
D_dimentions = 0
if task_prompt == "2":
    N_particles = 2
    D_dimentions = 2
    inter = input("Interaction? yes [y] or no [n]\n")
    if inter == "y":
        interaction = 1



if task_prompt == "1":
    N_particles = 1
    D_dimentions = 1

omega = 1
MC_cycles = 2**18
H_layers = 2
T_threads = 1
learning_rate = 0.035
type_sampling = 0  # 1 for hastings, 0 for brute force
optimizer = 1      # 1 for ADAM, 0 for vanilla GD

os.system("./main.out " + str(N_particles) + " " + str(D_dimentions) + " " + str(MC_cycles) + " " + str(H_layers) + " " + str(T_threads) + " " + str(learning_rate) + " " + str(type_sampling) + " " + str(interaction) + " " + str(optimizer) + " " + str(omega))

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


all = os.listdir()
for f in all:
    if "joined.txt" in f:
        os.remove(f)

for f in all:
    if "EnergySamples" in f:
        samples = np.loadtxt(f)
        print(samples)
        print("File has been read, start blocking...")
        (mean,var) = block(samples)
        std = np.sqrt(var)
        print("Statistical values:\n")
        print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))


for f in all:
    if "EnergySamples" in f:
        os.remove(f)
