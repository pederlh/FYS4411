import numpy as np
import matplotlib.pyplot as plt
import os
"""
ADAM = []
A_eta = []
A_its = []

GD = []
G_eta = []
G_its =[]

all = os.listdir()
all.remove("adamGD.py")

for f in all:
    if "ADAM" in f:
        w = f.split("_")
        ADAM.append(f)
        A_eta.append(float(w[9]))
        A_its.append(int(w[-2]))
    if "GD" in f:
        w = f.split("_")
        GD.append(f)
        G_eta.append(float(w[9]))
        G_its.append(int(w[-2]))


A_eta,ADAM,A_its = zip(*sorted(zip(A_eta,ADAM,A_its)))
G_eta,GD,G_its = zip(*sorted(zip(G_eta,GD,G_its)))

G_eta = np.array(G_eta)
A_eta = np.array(A_eta)

G_its = np.array(G_its)
A_its = np.array(A_its)

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

am =[]
ast = []

for a in ADAM:
    data = np.load(a)
    (mean, var) = block(data)
    std = np.sqrt(var)
    am.append(mean)
    ast.append(std)


gm =[]
gst = []

for g in GD:
    data = np.load(g)
    (mean, var) = block(data)
    std = np.sqrt(var)
    gm.append(mean)
    gst.append(std)

am = np.array(am)
ast = np.array(ast)
gm = np.array(gm)
gst = np.array(gst)


np.save("am",am)
np.save("ast",ast)
np.save("gm",gm)
np.save("gst",gst)

np.save("git", G_its)
np.save("ait", A_its)
np.save("get", G_eta)
np.save("aet", A_eta)


"""
"""
ADAM = []
A_eta = []
A_its = []

GD = []
G_eta = []
G_its =[]

all = os.listdir()
all.remove("adamGD.py")

for f in all:
    if "ADAM" in f:
        w = f.split("_")
        ADAM.append(f)
        A_eta.append(float(w[9]))
        A_its.append(int(w[-2]))
    if "GD" in f:
        w = f.split("_")
        GD.append(f)
        G_eta.append(float(w[9]))
        G_its.append(int(w[-2]))


A_eta,ADAM,A_its = zip(*sorted(zip(A_eta,ADAM,A_its)))
G_eta,GD,G_its = zip(*sorted(zip(G_eta,GD,G_its)))

G_eta = np.array(G_eta)
A_eta = np.array(A_eta)

G_its = np.array(G_its)
A_its = np.array(A_its)

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

am =[]
ast = []

for a in ADAM:
    data = np.load(a)
    (mean, var) = block(data)
    std = np.sqrt(var)
    am.append(mean)
    ast.append(std)


gm =[]
gst = []

for g in GD:
    data = np.load(g)
    (mean, var) = block(data)
    std = np.sqrt(var)
    gm.append(mean)
    gst.append(std)

am = np.array(am)
ast = np.array(ast)
gm = np.array(gm)
gst = np.array(gst)


np.save("am14",am)
np.save("ast14",ast)
np.save("gm14",gm)
np.save("gst14",gst)

np.save("git14", G_its)
np.save("ait14", A_its)
np.save("get14", G_eta)
np.save("aet14", A_eta)

"""
#Tolerance = 9e-4
ADAM = np.load("am.npy")
GD = np.load("gm.npy")

ADAM_std = np.load("ast.npy")
GD_std = np.load("gst.npy")

ADAM_its = np.load("ait.npy")
GD_its = np.load("git.npy")

ADAM_eta = np.load("aet.npy")
GD_eta = np.load("get.npy")


plt.errorbar(ADAM_eta,ADAM,yerr=ADAM_std, fmt = "ob",capsize=5, elinewidth=1,markeredgewidth=1, label = "ADAM")
plt.errorbar(GD_eta,GD,yerr=GD_std, fmt = "or",capsize=5, elinewidth=1,markeredgewidth=1, label = "GD")
plt.axhline(2.0,label = "Analytical")
plt.xlabel(r"$\eta$")
plt.ylabel(r"$\langle E_L \rangle$")
plt.legend()
plt.show()


plt.plot(ADAM_eta,ADAM_its, label = "ADAM")
plt.plot(GD_eta,GD_its, label = "GD")
plt.xlabel(r"$\eta$")
plt.ylabel("Iterations")
plt.legend()
plt.show()


# Tolerance = 1e-4
ADAM = np.load("am14.npy")
GD = np.load("gm14.npy")

ADAM_std = np.load("ast14.npy")
GD_std = np.load("gst14.npy")

ADAM_its = np.load("ait14.npy")
GD_its = np.load("git14.npy")

ADAM_eta = np.load("aet14.npy")
GD_eta = np.load("get14.npy")


plt.errorbar(ADAM_eta,ADAM,yerr=ADAM_std, fmt = ":ob",capsize=5, elinewidth=1,markeredgewidth=1, label = "ADAM")
plt.errorbar(GD_eta,GD,yerr=GD_std, fmt = ":or",capsize=5, elinewidth=1,markeredgewidth=1, label = "GD")
plt.axhline(2.0,label = "Analytical")
plt.xlabel(r"$\eta$")
plt.ylabel(r"$\langle E_L \rangle$")
plt.legend()
plt.show()


plt.plot(ADAM_eta,ADAM_its,"-o", label = "ADAM")
plt.plot(GD_eta,GD_its, "-o",label = "GD")
plt.xlabel(r"$\eta$")
plt.ylabel("Iterations")
plt.legend()
plt.show()
