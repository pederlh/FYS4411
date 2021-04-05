import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns

plt.style.use('seaborn')
sns.set(font_scale=1.4)
""" Script that both creates plots and performs blocking on local energy data sets """

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


def make_plots(task):
    #Plots local energy and std for gridsearch brute force MC
    if task == "b":
        N = [1,10,100,500]
        MC = 1e6
        path = "./Results/1b_simple_noninteracting/"

        t_avg_a = np.zeros(len(N))
        t_avg_n = np.zeros(len(N))

        t_std_a = np.zeros(len(N))
        t_std_n = np.zeros(len(N))

        count = 0

        alphas = []
        energies_a = []
        variances_a = []
        times_a = []
        energies_n = []
        variances_n = []
        times_n = []
        for n in N:
            mc = int(MC/n)
            filename_num = path + "spherical_HO_3" + "D_num" +"_N_"+ str(n) + "_MC_"+ str(mc) +"_stringID_0" + ".txt"
            filename_ana = path + "spherical_HO_3" + "D_ana" +"_N_"+ str(n) + "_MC_"+ str(mc) +"_stringID_0" + ".txt"

            with open(filename_ana, "r") as infile:
                lines = infile.readlines()
                for line in lines[1:]:
                    vals = line.split()
                    alphas.append(float(vals[0]))
                    energies_a.append(float(vals[1]))
                    variances_a.append(float(vals[2]))
                    times_a.append(float(vals[3]))

            with open(filename_num, "r") as infile:
                lines = infile.readlines()
                for line in lines[1:]:
                    vals = line.split()
                    energies_n.append(float(vals[1]))
                    variances_n.append(float(vals[2])) 
                    times_n.append(float(vals[3]))

            analytical_E = (3/2)*n

            n_str = "%i"%n
            plotname = "./Results/Plots/b/N_" + n_str +".pdf"
            plt.figure(n)

            # # Dumb standard error
            std_a = np.sqrt(np.asarray(variances_a)) / np.sqrt(MC)
            std_n = np.sqrt(np.asarray(variances_n)) / np.sqrt(MC)
            plt.errorbar(alphas, energies_a,yerr=std_a, fmt= "or",capsize=5, elinewidth=1, label = "Analytical derivative", markeredgewidth=1)
            plt.errorbar(alphas, energies_n,yerr=std_n, fmt= "ok",capsize=5, elinewidth=1, label = "Numerical derivative", markeredgewidth=1)

            # Blocking standard error
            # std_a_block = np.sqrt(block(np.asarray(energies_a))[1])
            # std_n_block = np.sqrt(block(np.asarray(energies_n))[1])
            # plt.errorbar(alphas, energies_a,yerr=std_a_block, fmt= "or",capsize=5, elinewidth=1, label = "Analytical derivative", markeredgewidth=1)
            # plt.errorbar(alphas, energies_n,yerr=std_n_block, fmt= "ok",capsize=5, elinewidth=1, label = "Numerical derivative", markeredgewidth=1)
            

            plt.axhline(y=analytical_E, color='mediumblue', linestyle=':', label = "Exact solution")
            plt.xlabel(r"$\alpha$")
            plt.ylabel(r"$\langle E_L \rangle$")
            plt.xticks() # plt.xticks(fontsize=fsz)
            plt.yticks()
            plt.legend()
            plt.tight_layout()
            plt.savefig(plotname)
            plt.grid('--')
            plt.show()

            t_avg_a[count] = np.mean(times_a)
            t_avg_n[count] = np.mean(times_n)

            t_std_a[count] = np.std(times_a)
            t_std_n[count] = np.std(times_n)



            # Reset for next value of N
            alphas = []
            energies_a = []
            variances_a = []
            times_a = []
            energies_n = []
            variances_n = []
            times_n = []

            count += 1
        """
        plt.plot(N,t_avg_a,"*" ,label = "Avg. CPU time w/analytical EL")
        plt.plot(N,t_avg_n,"*" ,label = "Avg. CPU time w/numerical EL")
        plt.xlabel("Number of particles")
        plt.ylabel("time [s]")
        plt.legend()
        plt.show()
        """

        print("avg time num:")
        print(t_avg_n)
        print("std time num:")
        print(t_std_n, "\n")
        print("avg time ana:")
        print(t_avg_a)
        print("std time ana:")
        print(t_std_a, "\n")

    #Plots local energy and std for gridsearch metropolis-hastings MC
    if task == "c":
        N = [1,10,100,500]
        MC = 1e6
        path = "./Results/1c_implementing_importance_sampling/"

        t_avg_a = np.zeros(len(N))
        t_avg_n = np.zeros(len(N))

        t_std_a = np.zeros(len(N))
        t_std_n = np.zeros(len(N))

        count = 0
        alphas = []
        energies_a = []
        variances_a = []
        times_a = []
        energies_n = []
        variances_n = []
        times_n = []
        for n in N:
            mc = int(MC/n)
            filename_num = path + "importance_spherical_HO_3" + "D_num" +"_N_"+ str(n) + "_MC_"+ str(mc) +"_stringID_0" + ".txt"
            filename_ana = path + "importance_spherical_HO_3" + "D_ana" +"_N_"+ str(n) + "_MC_"+ str(mc) +"_stringID_0" + ".txt"

            with open(filename_ana, "r") as infile:
                lines = infile.readlines()
                for line in lines[1:]:
                    vals = line.split()
                    alphas.append(float(vals[0]))
                    energies_a.append(float(vals[1]))
                    variances_a.append(float(vals[2]))
                    times_a.append(float(vals[3]))

            with open(filename_num, "r") as infile:
                lines = infile.readlines()
                for line in lines[1:]:
                    vals = line.split()
                    energies_n.append(float(vals[1]))
                    variances_n.append(float(vals[2]))
                    times_n.append(float(vals[3]))

            analytical_E = (3/2)*n

            n_str = "%i" % n
            plotname = "./Results/Plots/c/N_" + n_str + ".pdf"
            plt.figure(n)

            # Dumb standard error
            std_a = np.sqrt(np.asarray(variances_a)) / np.sqrt(MC)
            std_n = np.sqrt(np.asarray(variances_n)) / np.sqrt(MC)
            plt.errorbar(alphas, energies_a,yerr=std_a, fmt= "or",capsize=5, elinewidth=1, label = "Analytical derivative", markeredgewidth=1)
            plt.errorbar(alphas, energies_n,yerr=std_n, fmt= "ok",capsize=5, elinewidth=1, label = "Numerical derivative", markeredgewidth=1)

            # Blocking standard error
            # std_a_block = np.sqrt(block(np.asarray(energies_a))[1])
            # std_n_block = np.sqrt(block(np.asarray(energies_n))[1])
            # plt.errorbar(alphas, energies_a,yerr=std_a_block, fmt= "or",capsize=5, elinewidth=1, label = "Analytical derivative", markeredgewidth=1)
            # plt.errorbar(alphas, energies_n,yerr=std_n_block, fmt= "ok",capsize=5, elinewidth=1, label = "Numerical derivative", markeredgewidth=1)


            plt.axhline(y=analytical_E, color='mediumblue', linestyle=':', label = "Exact solution")
            plt.xlabel(r"$\alpha$")
            plt.ylabel(r"$\langle E_L \rangle$")

            plt.xticks()
            plt.yticks()
            plt.legend()
            plt.tight_layout()
            plt.savefig(plotname)
            plt.grid('--')
            plt.show()

            t_avg_a[count] = np.mean(times_a)
            t_avg_n[count] = np.mean(times_n)

            t_std_a[count] = np.std(times_a)
            t_std_n[count] = np.std(times_n)

            alphas = []
            energies_a = []
            variances_a = []
            times_a = []
            energies_n = []
            variances_n = []
            times_n = []

            count += 1

        """
        plt.plot(N,t_avg_a,"*" ,label = "Avg. CPU time w/analytical EL")
        plt.plot(N,t_avg_n,"*" ,label = "Avg. CPU time w/numerical EL")
        plt.xlabel("Number of particles")
        plt.ylabel("time [s]")
        plt.legend()
        plt.show()
        """

        print("avg time num:")
        print("")
        print(t_avg_n)
        print("")
        print("")
        print("avg time ana:")
        print("")
        print(t_avg_a)

    #Plots gradient descent evolution for non-interacting bosons
    if task == "d":
        path = "./Results/GD_alphas/"
        os.chdir(path)
        GD = os.listdir()
        GD.sort()

        numbers = []

        for files in GD:
            for words in files.split("_"):
                if words.isdigit():
                    numbers.append(int(words))



        numbers, GD = zip(*sorted(zip(numbers, GD)))

        GD_2 = np.loadtxt(GD[0])
        GD_16 = np.loadtxt(GD[1])
        GD_32 =  np.loadtxt(GD[2])
        GD_64 = np.loadtxt(GD[3])
        GD_128 = np.loadtxt(GD[4])
        alphas= [GD_2,GD_16,GD_32,GD_64,GD_128]

        it_2 = np.linspace(0,len(GD_2),len(GD_2))
        it_16 = np.linspace(0,len(GD_16),len(GD_16))
        it_32 = np.linspace(0,len(GD_32),len(GD_32))
        it_64 = np.linspace(0,len(GD_64),len(GD_64))
        it_128 = np.linspace(0,len(GD_128),len(GD_128))
        its = [it_2,it_16,it_32,it_64,it_128]

        colors = ["r","b","g","c","m"]
        for i in range(len(its)):
            plt.plot(its[i], alphas[i], "-o", markersize=4, label = "$N$ = %i" %numbers[i], color = colors[i])
            plt.plot(its[i][-1], alphas[i][-1],"*", markersize=10, markeredgecolor='k', markeredgewidth=0.5, color = colors[i], zorder=10)


        os.chdir("../")
        plotname = "../Results/Plots/d/GD_evolution.pdf"
        plt.ylabel(r"$\alpha$")
        plt.xlabel("Number of iterations")
        # plt.title("Gradient descent for non-interacting bosons")
        plt.legend()
        plt.xticks()
        plt.yticks()
        plt.tight_layout()
        plt.savefig(plotname)
        plt.show()

    if task == "e":
        print("Blocking results from GD without repulsion")

        path = "./Results/1e_implementing_gradient_descent_and_blocking/"
        os.chdir(path)
        runs = os.listdir()

        t2 = []
        t16 = []
        t32 = []
        t64 = []
        t128 = []

        N32 = [s for s in runs if "32_N_" in s]
        runs = [s for s in runs if "32_N_" not in s]
        N64 = [s for s in runs if "64_N_" in s]
        N2 = [s for s in runs if "2_N_" in s]
        N16 = [s for s in runs if "16_N_" in s]
        N128 = [s for s in runs if "128_N_" in s]

        thread0 = np.loadtxt(N2[0])
        t2.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N2[1])
        t2.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N2[2])
        t2.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N2[3])
        t2.append(thread3[0])
        thread3 = thread3[1:]

        joined_N2 = np.concatenate((thread0,thread1,thread2,thread3))
        joined_N2_mean = np.mean(joined_N2)
        joined_N2_std = np.std(joined_N2)

        thread0 = np.loadtxt(N16[0])
        t16.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N16[1])
        t16.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N16[2])
        t16.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N16[3])
        t16.append(thread3[0])
        thread3 = thread3[1:]

        joined_N16 = np.concatenate((thread0,thread1,thread2,thread3))
        joined_N16_mean = np.mean(joined_N16)
        joined_N16_std = np.std(joined_N16)

        thread0 = np.loadtxt(N32[0])
        t32.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N32[1])
        t32.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N32[2])
        t32.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N32[3])
        t32.append(thread3[0])
        thread3 = thread3[1:]

        joined_N32 = np.concatenate((thread0,thread1,thread2,thread3))
        joined_N32_mean = np.mean(joined_N32)
        joined_N32_std = np.std(joined_N32)

        thread0 = np.loadtxt(N64[0])
        t64.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N64[1])
        t64.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N64[2])
        t64.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N64[3])
        t64.append(thread3[0])
        thread3 = thread3[1:]

        joined_N64 = np.concatenate((thread0,thread1,thread2,thread3))
        joined_N64_mean = np.mean(joined_N64)
        joined_N64_std = np.std(joined_N64)

        thread0 = np.loadtxt(N128[0])
        t128.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N128[1])
        t128.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N128[2])
        t128.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N128[3])
        t128.append(thread3[0])
        thread3 = thread3[1:]

        joined_N128 = np.concatenate((thread0,thread1,thread2,thread3))
        joined_N128_mean = np.mean(joined_N128)
        joined_N128_std = np.std(joined_N128)

        print("Files have been read, start blocking...")
        (mean, var) = block(joined_N2)
        std = np.sqrt(var)
        """
        data ={'Mean':[mean], 'STDev':[std]}
        frame_2 = pd.DataFrame(data,index=['Values'])
        index = frame_2.index
        index.name = "N = 2"
        print(frame_2)
        """

        print("Statistical values:\n")
        print("N=2")
        print("Naive:     Mean = %.7f        STDev = %.7f" % (joined_N2_mean, joined_N2_std))
        print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))


        (mean, var) = block(joined_N16)
        std = np.sqrt(var)
        """
        data ={'Mean':[mean], 'STDev':[std]}
        frame_16 = pd.DataFrame(data,index=['Values'])
        index = frame_16.index
        index.name = "N = 16"
        print(frame_16)
        """

        print("\n N=16")
        print("Naive:     Mean = %.7f        STDev = %.7f" % (joined_N16_mean, joined_N16_std))
        print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))

        (mean, var) = block(joined_N32)
        std = np.sqrt(var)
        """
        data ={'Mean':[mean], 'STDev':[std]}
        frame_32 = pd.DataFrame(data,index=['Values'])
        index = frame_32.index
        index.name = "N = 32"
        print(frame_32)
        """
        print("\n N=32")
        print("Naive:     Mean = %.7f        STDev = %.7f" % (joined_N32_mean, joined_N32_std))
        print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))

        (mean, var) = block(joined_N64)
        std = np.sqrt(var)
        """
        data ={'Mean':[mean], 'STDev':[std]}
        frame_64 = pd.DataFrame(data,index=['Values'])
        index = frame_64.index
        index.name = "N = 64"
        print(frame_64)
        """

        print("\n N=64")
        print("Naive:     Mean = %.7f        STDev = %.7f" % (joined_N64_mean, joined_N64_std))
        print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))

        (mean, var) = block(joined_N128)
        std = np.sqrt(var)


        print("\n N=128")
        print("Naive:     Mean = %.7f        STDev = %.7f" % (joined_N128_mean, joined_N128_std))
        print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))
        """
        data ={'Mean':[mean], 'STDev':[std]}
        frame_128 = pd.DataFrame(data,index=['Values'])
        index = frame_128.index
        index.name = "N = 128"
        print(frame_128)
        """

    if task == "g":

        print("Blocking results from GD with repulsion")

        path = "./Results/1g_implementing_repulsion/"
        os.chdir(path)
        runs = os.listdir()

        t2 = []
        t16 = []
        t32 = []
        t64 = []

        N32 = [s for s in runs if "32_N_" in s]
        runs = [s for s in runs if "32_N_" not in s]
        N64 = [s for s in runs if "64_N_" in s]
        N2 = [s for s in runs if "2_N_" in s]
        N16 = [s for s in runs if "16_N_" in s]

        thread0 = np.loadtxt(N2[0])
        t2.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N2[1])
        t2.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N2[2])
        t2.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N2[3])
        t2.append(thread3[0])
        thread3 = thread3[1:]

        joined_N2 = np.concatenate((thread0,thread1,thread2,thread3))
        joined_N2_mean = np.mean(joined_N2)
        joined_N2_std = np.std(joined_N2)

        thread0 = np.loadtxt(N16[0])
        t16.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N16[1])
        t16.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N16[2])
        t16.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N16[3])
        t16.append(thread3[0])
        thread3 = thread3[1:]

        joined_N16 = np.concatenate((thread0,thread1,thread2,thread3))
        joined_N16_mean = np.mean(joined_N16)
        joined_N16_std = np.std(joined_N16)

        thread0 = np.loadtxt(N32[0])
        t32.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N32[1])
        t32.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N32[2])
        t32.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N32[3])
        t32.append(thread3[0])
        thread3 = thread3[1:]

        joined_N32 = np.concatenate((thread0,thread1,thread2,thread3))
        joined_N32_mean = np.mean(joined_N32)
        joined_N32_std = np.std(joined_N32)

        thread0 = np.loadtxt(N64[0])
        t64.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N64[1])
        t64.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N64[2])
        t64.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N64[3])
        t64.append(thread3[0])
        thread3 = thread3[1:]

        joined_N64 = np.concatenate((thread0,thread1,thread2,thread3))
        joined_N64_mean = np.mean(joined_N64)
        joined_N64_std = np.std(joined_N64)


        print("Files has been read, start blocking...")
        (mean, var) = block(joined_N2)
        std = np.sqrt(var)
        """
        data ={'Mean':[mean], 'STDev':[std]}
        frame_2 = pd.DataFrame(data,index=['Values'])
        index = frame_2.index
        index.name = "N = 2"
        print(frame_2)
        """

        print("Statistical values:\n")
        print("N=2")
        print("Naive:     Mean = %.7f        STDev = %.7f" % (joined_N2_mean, joined_N2_std))
        print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))


        (mean, var) = block(joined_N16)
        std = np.sqrt(var)
        """
        data ={'Mean':[mean], 'STDev':[std]}
        frame_16 = pd.DataFrame(data,index=['Values'])
        index = frame_16.index
        index.name = "N = 16"
        print(frame_16)
        """


        print("Statistical values:\n")
        print("N=16")
        print("Naive:     Mean = %.7f        STDev = %.7f" % (joined_N16_mean, joined_N16_std))
        print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))

        (mean, var) = block(joined_N32)
        std = np.sqrt(var)
        """
        data ={'Mean':[mean], 'STDev':[std]}
        frame_32 = pd.DataFrame(data,index=['Values'])
        index = frame_32.index
        index.name = "N = 32"
        print(frame_32)
        """

        print("Statistical values:\n")
        print("N=32")
        print("Naive:     Mean = %.7f        STDev = %.7f" % (joined_N32_mean, joined_N32_std))
        print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))

        (mean, var) = block(joined_N64)
        std = np.sqrt(var)
        """
        data ={'Mean':[mean], 'STDev':[std]}
        frame_64 = pd.DataFrame(data,index=['Values'])
        index = frame_64.index
        index.name = "N = 64"
        print(frame_64)
        """

        print("Statistical values:\n")
        print("N=64")
        print("Naive:     Mean = %.7f        STDev = %.7f" % (joined_N64_mean, joined_N64_std))
        print("Blocking:  Mean = %.7f        STDev = %.7f" % (mean,std))

    if task == "h":

        path = "./Results/1h_one_body_densities/"
        os.chdir(path)
        OBD_list = os.listdir()

        OBD_list_I = [s for s in OBD_list if "Interaction" in s]
        OBD_list = [s for s in OBD_list if "Interaction" not in s]

        OBD_N2 = [s for s in OBD_list if "_N_2_" in s]
        OBD_N16 = [s for s in OBD_list if "_N_16_" in s]
        OBD_N32 = [s for s in OBD_list if "_N_32_" in s]
        OBD_N64 = [s for s in OBD_list if "_N_64_" in s]
        OBD_N128 = [s for s in OBD_list if "_N_128_" in s]

        OBD_N2_I = [s for s in OBD_list_I if "_N_2_" in s]
        OBD_N16_I = [s for s in OBD_list_I if "_N_16_" in s]
        OBD_N32_I = [s for s in OBD_list_I if "_N_32_" in s]
        OBD_N64_I = [s for s in OBD_list_I if "_N_64_" in s]

        thread0 = np.loadtxt(OBD_N2[0])
        thread1 = np.loadtxt(OBD_N2[1])
        thread2 = np.loadtxt(OBD_N2[2])
        thread3 = np.loadtxt(OBD_N2[3])

        joined_OBD_N2 = (thread0 + thread1 + thread2 + thread3)/4

        thread0 = np.loadtxt(OBD_N16[0])
        thread1 = np.loadtxt(OBD_N16[1])
        thread2 = np.loadtxt(OBD_N16[2])
        thread3 = np.loadtxt(OBD_N16[3])

        joined_OBD_N16 = (thread0 + thread1 + thread2 + thread3)/4

        thread0 = np.loadtxt(OBD_N32[0])
        thread1 = np.loadtxt(OBD_N32[1])
        thread2 = np.loadtxt(OBD_N32[2])
        thread3 = np.loadtxt(OBD_N32[3])

        joined_OBD_N32 = (thread0 + thread1 + thread2 + thread3)/4

        thread0 = np.loadtxt(OBD_N64[0])
        thread1 = np.loadtxt(OBD_N64[1])
        thread2 = np.loadtxt(OBD_N64[2])
        thread3 = np.loadtxt(OBD_N64[3])

        joined_OBD_N64 = (thread0 + thread1 + thread2 + thread3)/4

        thread0 = np.loadtxt(OBD_N128[0])
        thread1 = np.loadtxt(OBD_N128[1])
        thread2 = np.loadtxt(OBD_N128[2])
        thread3 = np.loadtxt(OBD_N128[3])

        joined_OBD_N128 = (thread0 + thread1 + thread2 + thread3)/4

        thread0 = np.loadtxt(OBD_N2_I[0])
        thread1 = np.loadtxt(OBD_N2_I[1])
        thread2 = np.loadtxt(OBD_N2_I[2])
        thread3 = np.loadtxt(OBD_N2_I[3])

        joined_OBD_N2_I = (thread0 + thread1 + thread2 + thread3)/4

        thread0 = np.loadtxt(OBD_N16_I[0])
        thread1 = np.loadtxt(OBD_N16_I[1])
        thread2 = np.loadtxt(OBD_N16_I[2])
        thread3 = np.loadtxt(OBD_N16_I[3])

        joined_OBD_N16_I = (thread0 + thread1 + thread2 + thread3)/4

        thread0 = np.loadtxt(OBD_N32_I[0])
        thread1 = np.loadtxt(OBD_N32_I[1])
        thread2 = np.loadtxt(OBD_N32_I[2])
        thread3 = np.loadtxt(OBD_N32_I[3])

        joined_OBD_N32_I = (thread0 + thread1 + thread2 + thread3)/4

        thread0 = np.loadtxt(OBD_N64_I[0])
        thread1 = np.loadtxt(OBD_N64_I[1])
        thread2 = np.loadtxt(OBD_N64_I[2])
        thread3 = np.loadtxt(OBD_N64_I[3])

        joined_OBD_N64_I = (thread0 + thread1 + thread2 + thread3)/4


        r = np.linspace(0,8,50)
        r_ana = np.linspace(0,8,500)
        os.chdir("../")
        path = "../Results/Plots/h/"

        # N = 2
        plotname = "N_2.pdf"
        A = 2/np.trapz(joined_OBD_N2*r**2,r)
        B = 2/np.trapz(joined_OBD_N2_I*r**2,r)
        plt.plot(r[:20],A*joined_OBD_N2[:20], "-o", ms=5, label="No interaction")
        plt.plot(r[:20],B*joined_OBD_N2_I[:20], "-o", ms=5, label="Interaction")

        rho_solution = 2 * (np.sqrt(np.pi)**(-3) * 4 * np.pi * np.exp(-r_ana**2))
        plt.plot(r_ana[:200],rho_solution[:200], label = "Analytical solution")

        plt.xlabel(r"$r$")
        plt.ylabel(r"$\rho(r)$")
        plt.legend()
        plt.xticks()
        plt.yticks()
        plt.savefig(path+plotname)
        plt.show()

        # N = 16
        plotname = "N_16.pdf"
        A = 16/np.trapz(joined_OBD_N16*r**2,r)
        B = 16/np.trapz(joined_OBD_N16_I*r**2,r)
        plt.plot(r[:20],A*joined_OBD_N16[:20], "-o", ms=5, label="No interaction")
        plt.plot(r[:20],B*joined_OBD_N16_I[:20], "-o", ms=5, label="Interaction")

        rho_solution = 16 * (np.sqrt(np.pi)**(-3) * 4 * np.pi * np.exp(-r_ana**2))
        plt.plot(r_ana[:200],rho_solution[:200], label = "Analytical solution")

        plt.xlabel(r"$r$")
        plt.ylabel(r"$\rho(r)$")
        plt.legend()
        plt.xticks()
        plt.yticks()
        plt.savefig(path+plotname)
        plt.show()


        # TEST NORMALIZATION: N = 16
        print(f"{np.trapz(r_ana**2 * rho_solution,r_ana) :3.15f}") # Should be 16
        print(f"{np.trapz(r**2*A*joined_OBD_N16,r) :3.15f}")    # Should be 16
        print(f"{np.trapz(r**2*B*joined_OBD_N16_I,r) :3.15f}")  # Should be 16 
        
        plt.plot(r[:20],A*joined_OBD_N16[:20]*r[:20]**2, "-o", ms=5, label="No interaction")
        plt.plot(r[:20],B*joined_OBD_N16_I[:20]*r[:20]**2, "-o", ms=5, label="Interaction")

        rho_solution = 16 * (np.sqrt(np.pi)**(-3) * 4 * np.pi * np.exp(-r_ana**2))
        plt.plot(r_ana[:200],rho_solution[:200]*r_ana[:200]**2, label = "Analytical solution")

        plt.xlabel(r"$r$")
        plt.ylabel(r"$r^2 \cdot \rho(r)$")
        plt.legend()
        plt.xticks()
        plt.yticks()
        plt.savefig(path + "test_normalization_N_16.pdf")
        plt.show()

        plotname = "N_32.pdf"
        A = 32/np.trapz(joined_OBD_N32*r**2,r)
        B = 32/np.trapz(joined_OBD_N32_I*r**2,r)
        plt.plot(r[:20],A*joined_OBD_N32[:20], "-o", ms=5, label="No interaction")
        plt.plot(r[:20],B*joined_OBD_N32_I[:20], "-o", ms=5, label="Interaction")

        rho_solution = 32 * (np.sqrt(np.pi)**(-3) * 4 * np.pi * np.exp(-r_ana**2))
        plt.plot(r_ana[:200],rho_solution[:200], label = "Analytical solution")


        plt.xlabel(r"$r$")
        plt.ylabel(r"$\rho(r)$")
        plt.legend()
        plt.xticks()
        plt.yticks()
        plt.savefig(path+plotname)
        plt.show()


        plotname = "N_64.pdf"
        A = 64/np.trapz(joined_OBD_N64*r**2,r)
        B = 64/np.trapz(joined_OBD_N64_I*r**2,r)
        plt.plot(r[:20],A*joined_OBD_N64[:20], "-o", ms=5, label="No interaction")
        plt.plot(r[:20],B*joined_OBD_N64_I[:20], "-o", ms=5, label="Interaction")

        rho_solution = 64 * (np.sqrt(np.pi)**(-3) * 4 * np.pi * np.exp(-r_ana**2))
        plt.plot(r_ana[:200],rho_solution[:200], label = "Analytical solution")


        plt.xlabel(r"$r$")
        plt.ylabel(r"$\rho(r)$")
        plt.legend()
        plt.xticks()
        plt.yticks()
        plt.savefig(path+plotname)
        plt.show()

        plotname = "N_128.pdf"
        A = 128/np.trapz(joined_OBD_N128*r**2,r)
        plt.plot(r[:20],A*joined_OBD_N128[:20],"-o", color="b", ms=5, label="No interaction")

        rho_solution = 128 * (np.sqrt(np.pi)**(-3) * 4 * np.pi * np.exp(-r_ana**2))
        plt.plot(r_ana[:200],rho_solution[:200],"g", label = "Analytical solution")

        plt.xlabel(r"$r$")
        plt.ylabel(r"$\rho(r)$")
        plt.legend()
        plt.xticks()
        plt.yticks()
        plt.savefig(path+plotname)
        plt.show()

    if task =="delta_t":

        path = "./Results/delta_t_comparisons/"
        #dts = [1.0,5.0,0.1,0.5,0.01,0.05,0.001,0.005,0.0001,0.0005,0.00001,0.00005]
        dts= ["5.0","1.0","0.5","0.1","0.05","0.01","0.005","0.001","0.0005","0.0001","0.00005","0.00001"]

        energy_045 = np.zeros(len(dts))     # alpha = 0.45
        energy_050 = np.zeros(len(dts))     # alpha = 0.50
        energy_055 = np.zeros(len(dts))     # alpha = 0.55

        variance_045 = np.zeros(len(dts))   
        variance_050 = np.zeros(len(dts))   
        variance_055 = np.zeros(len(dts))   

        outfile_ana = path + "delta_t_comparrison_ana.txt"

        for i in range(len(dts)):
            filename_ana = path + "dt_run_" + dts[i] + "_ana.txt"

            with open(filename_ana,"r") as infile:
                lines = infile.readlines()

                energy_045[i] = float(lines[4].split()[1])
                energy_050[i] = float(lines[5].split()[1])
                energy_055[i] = float(lines[6].split()[1])

                variance_045[i] = float(lines[4].split()[2])
                variance_050[i] = float(lines[5].split()[2])
                variance_055[i] = float(lines[6].split()[2])

        # Ghetto variance
        with open(outfile_ana,"w") as outfile:
            outfile.write("dt - E_045 - var_045 - E_050 - var_055 - E_055 - var_055 \n")
            for i in range(len(dts)):
                outfile.write(dts[i]+ " & " + str(energy_045[i]) + " & " + str(variance_045[i]) + " & " + str(energy_050[i]) \
                                    + " & " + str(variance_050[i]) + " & " + str(energy_055[i]) + " & " + str(variance_055[i]) + " \\\\ \n")

        """ Below code not working - FIX"""
        # # Block std. error
        # std_045_block = std_050_block = std_055_block = np.zeros(len(dts))
        # for i in range(len(dts)):
        #     std_045_block[i] = np.sqrt(block(energy_045)[1])
        #     std_050_block[i] = np.sqrt(block(energy_050)[1])
        #     std_055_block[i] = np.sqrt(block(energy_055)[1])

        # with open(outfile_ana,"w") as outfile:
        #     outfile.write("dt - E_045 - std_block_045 - E_050 - std_block_055 - E_055 - std_block_055 \n")
        #     for i in range(len(dts)):
        #         outfile.write(dts[i] + " & " + str(energy_045[i]) + " & " + str(std_045_block[i]) + " & " + str(energy_050[i]) \
        #                             + " & " + str(std_050_block[i]) + " & " + str(energy_055[i]) + " & " + str(std_055_block[i]) + " \\\\ \n")



if __name__ == "__main__":
    print("Possible commands: 'b' (simplest), 'c' (importance sampling), 'd' (gradient descent evolution), \n" \
        + "'e' (no-repulsion), 'g' (repulsion), 'h' (one body densities) and 'delta_t' (find optimal delta t)")
    task = input("Which task to make plots for? ")
    assert task in ["b", "c", "d", "e", "g", "h", "delta_t"], "Input not recognized."
    make_plots(task)
