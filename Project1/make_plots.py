import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

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

            std_a = np.sqrt(variances_a)
            std_n = np.sqrt(variances_n)

            fsz = 14
            n_str = "%i"%n
            plotname = "./Results/Plots/b/N_" + n_str +".pdf"
            plt.title("Number of particles = %i" % n, fontsize = fsz)
            plt.errorbar(alphas, energies_a,yerr=std_a, fmt= "or",capsize=5, elinewidth=1, label = "Local energy (analytical)")
            plt.errorbar(alphas, energies_n,yerr=std_n, fmt= "ok",capsize=5, elinewidth=1, label = "Local energy (numerical)")
            plt.axhline(y=analytical_E, color='b', linestyle='-', label = "Local energy (exact)")
            plt.xlabel("Alpha", fontsize = fsz )
            plt.ylabel("Local energy", fontsize = fsz)
            plt.xticks(fontsize=fsz )
            plt.yticks(fontsize= fsz)
            plt.legend(fontsize = fsz-2)
            plt.savefig(plotname)
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
        print(t_avg_n)
        print("std time num:")
        print(t_std_n, "\n")
        print("avg time ana:")
        print(t_avg_a)
        print("std time ana:")
        print(t_std_a, "\n")

    #Plots local energy and std for gridsearch mtropolis-hastings MC
    if task == "c":
        N = [1,10,100,500]
        MC = 1e6
        string_ID = 0
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

            std_a = np.sqrt(variances_a)
            std_n = np.sqrt(variances_n)

            fsz = 14
            n_str = "%i" % n
            plotname = "./Results/Plots/c/N_" + n_str + ".pdf"
            plt.title("Number of particles = %i" % n, fontsize = fsz)
            plt.errorbar(alphas, energies_a,yerr=std_a, fmt= "or",capsize=5, elinewidth=1, label = "Local energy (analytical)")
            plt.errorbar(alphas, energies_n,yerr=std_n, fmt= "ok",capsize=5, elinewidth=1, label = "Local energy (numerical)")
            plt.axhline(y=analytical_E, color='b', linestyle='-', label = "Local energy (exact)")
            plt.xlabel("alpha",fontsize = fsz)
            plt.ylabel("Local energy",fontsize = fsz)
            plt.xticks(fontsize=fsz )
            plt.yticks(fontsize= fsz)
            plt.legend(fontsize = fsz-2)
            plt.savefig(plotname)
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
            plt.plot(its[i], alphas[i],label = "N = %i" %numbers[i], color = colors[i])
            plt.plot(its[i][-1], alphas[i][-1],"*", color = colors[i])


        os.chdir("../")
        plotname = "../Results/Plots/d/GD_evolution.pdf"
        fsz = 14
        plt.ylabel("Alpha value", fontsize = fsz)
        plt.xlabel("Number of iterations", fontsize = fsz)
        plt.title("Gradient descent for non-interacting bosons", fontsize = fsz)
        plt.legend(fontsize = fsz)
        plt.xticks(fontsize=fsz )
        plt.yticks(fontsize= fsz)
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

        print("Files has been read, start blocking...")
        (mean, var) = block(joined_N2)
        std = np.sqrt(var)
        data ={'Mean':[mean], 'STDev':[std]}
        frame_2 = pd.DataFrame(data,index=['Values'])
        index = frame_2.index
        index.name = "N = 2"
        print(frame_2)

        (mean, var) = block(joined_N16)
        std = np.sqrt(var)
        data ={'Mean':[mean], 'STDev':[std]}
        frame_16 = pd.DataFrame(data,index=['Values'])
        index = frame_16.index
        index.name = "N = 16"
        print(frame_16)

        (mean, var) = block(joined_N32)
        std = np.sqrt(var)
        data ={'Mean':[mean], 'STDev':[std]}
        frame_32 = pd.DataFrame(data,index=['Values'])
        index = frame_32.index
        index.name = "N = 32"
        print(frame_32)

        (mean, var) = block(joined_N64)
        std = np.sqrt(var)
        data ={'Mean':[mean], 'STDev':[std]}
        frame_64 = pd.DataFrame(data,index=['Values'])
        index = frame_64.index
        index.name = "N = 64"
        print(frame_64)

        (mean, var) = block(joined_N128)
        std = np.sqrt(var)
        data ={'Mean':[mean], 'STDev':[std]}
        frame_128 = pd.DataFrame(data,index=['Values'])
        index = frame_128.index
        index.name = "N = 128"
        print(frame_128)

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


        print("Files has been read, start blocking...")
        (mean, var) = block(joined_N2)
        std = np.sqrt(var)
        data ={'Mean':[mean], 'STDev':[std]}
        frame_2 = pd.DataFrame(data,index=['Values'])
        index = frame_2.index
        index.name = "N = 2"
        print(frame_2)

        (mean, var) = block(joined_N16)
        std = np.sqrt(var)
        data ={'Mean':[mean], 'STDev':[std]}
        frame_16 = pd.DataFrame(data,index=['Values'])
        index = frame_16.index
        index.name = "N = 16"
        print(frame_16)

        (mean, var) = block(joined_N32)
        std = np.sqrt(var)
        data ={'Mean':[mean], 'STDev':[std]}
        frame_32 = pd.DataFrame(data,index=['Values'])
        index = frame_32.index
        index.name = "N = 32"
        print(frame_32)

        (mean, var) = block(joined_N64)
        std = np.sqrt(var)
        data ={'Mean':[mean], 'STDev':[std]}
        frame_64 = pd.DataFrame(data,index=['Values'])
        index = frame_64.index
        index.name = "N = 64"
        print(frame_64)

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
        fsz = 14
        os.chdir("../")
        path = "../Results/Plots/h/"

        plt.title("N = 2",fontsize = fsz)
        plotname = "N_2.pdf"
        A = 2/np.trapz(joined_OBD_N2*r**2,r)
        B = 2/np.trapz(joined_OBD_N2_I*r**2,r)
        plt.plot(r[:20],A*joined_OBD_N2[:20],label="No interaction")
        plt.plot(r[:20],B*joined_OBD_N2_I[:20],label="Interaction")

        better_ana_rho = 2 * (np.sqrt(np.pi) ** (-3) * 4 * np.pi * np.exp(-r ** 2))
        plt.plot(r[:20],better_ana_rho[:20], label = "Analytical rho")

        plt.xlabel("r",fontsize = fsz)
        plt.ylabel("$\\rho(r)$",fontsize = fsz)
        plt.legend(fontsize = fsz)
        plt.xticks(fontsize=fsz )
        plt.yticks(fontsize= fsz)
        plt.savefig(path+plotname)
        plt.show()


        plt.title("N = 16",fontsize = fsz)
        plotname = "N_16.pdf"
        A = 16/np.trapz(joined_OBD_N16*r**2,r)
        B = 16/np.trapz(joined_OBD_N16_I*r**2,r)
        plt.plot(r[:20],A*joined_OBD_N16[:20],label="No interaction")
        plt.plot(r[:20],B*joined_OBD_N16_I[:20],label="Interaction")

        better_ana_rho = 16 * (np.sqrt(np.pi) ** (-3) * 4 * np.pi * np.exp(-r ** 2))
        plt.plot(r[:20],better_ana_rho[:20], label = "Analytical rho")

        plt.xlabel("r",fontsize = fsz)
        plt.ylabel("$\\rho(r)$",fontsize = fsz)
        plt.legend(fontsize = fsz)
        plt.xticks(fontsize=fsz )
        plt.yticks(fontsize= fsz)
        plt.savefig(path+plotname)
        plt.show()

        plt.title("N = 32",fontsize = fsz)
        plotname = "N_32.pdf"
        A = 32/np.trapz(joined_OBD_N32*r**2,r)
        B = 32/np.trapz(joined_OBD_N32_I*r**2,r)
        plt.plot(r[:20],A*joined_OBD_N32[:20],label="No interaction")
        plt.plot(r[:20],B*joined_OBD_N32_I[:20],label="Interaction")

        better_ana_rho = 32 * (np.sqrt(np.pi) ** (-3) * 4 * np.pi * np.exp(-r ** 2))
        plt.plot(r[:20],better_ana_rho[:20], label = "Analytical rho")

        plt.xlabel("r",fontsize = fsz)
        plt.ylabel("$\\rho(r)$",fontsize = fsz)
        plt.legend(fontsize = fsz)
        plt.xticks(fontsize=fsz )
        plt.yticks(fontsize= fsz)
        plt.savefig(path+plotname)
        plt.show()


        plt.title("N = 64",fontsize = fsz)
        plotname = "N_64.pdf"
        A = 64/np.trapz(joined_OBD_N64*r**2,r)
        B = 64/np.trapz(joined_OBD_N64_I*r**2,r)
        plt.plot(r[:20],A*joined_OBD_N64[:20],label="No interaction")
        plt.plot(r[:20],B*joined_OBD_N64_I[:20],label="Interaction")

        better_ana_rho = 64 * (np.sqrt(np.pi) ** (-3) * 4 * np.pi * np.exp(-r ** 2))
        plt.plot(r[:20],better_ana_rho[:20], label = "Analytical rho")

        plt.xlabel("r",fontsize = fsz)
        plt.ylabel("$\\rho(r)$",fontsize = fsz)
        plt.legend(fontsize = fsz)
        plt.xticks(fontsize=fsz )
        plt.yticks(fontsize= fsz)
        plt.savefig(path+plotname)
        plt.show()

        plt.title("N = 128",fontsize = fsz)
        plotname = "N_128.pdf"
        A = 128/np.trapz(joined_OBD_N128*r**2,r)
        plt.plot(r[:20],A*joined_OBD_N128[:20],"b",label="No interaction")

        better_ana_rho = 128 * (np.sqrt(np.pi) ** (-3) * 4 * np.pi * np.exp(-r ** 2))
        plt.plot(r[:20],better_ana_rho[:20],"g", label = "Analytical rho")

        plt.xlabel("r",fontsize = fsz)
        plt.ylabel("$\\rho(r)$",fontsize = fsz)
        plt.legend(fontsize = fsz)
        plt.xticks(fontsize=fsz )
        plt.yticks(fontsize= fsz)
        plt.savefig(path+plotname)
        plt.show()

    if task =="delta_t":

        path = "./Results/data_for_report/"
        #dts = [1.0,5.0,0.1,0.5,0.01,0.05,0.001,0.005,0.0001,0.0005,0.00001,0.00005]
        dts= ["5.0","1.0","0.5","0.1","0.05","0.01","0.005","0.001","0.0005","0.0001","0.00005","0.00001"]
        variance = np.zeros(len(dts))
        time = np.zeros(len(dts))
        energy = np.zeros(len(dts))

        outfile_num = path + "delta_t_comparrison_num.txt"
        outfile_ana = path + "delta_t_comparrison_ana.txt"

        variance_n = np.zeros(len(dts))
        time_n = np.zeros(len(dts))
        energy_n = np.zeros(len(dts))

        for i in range(len(dts)):
            filename_num = path + "dt_run_" + dts[i] + "_num.txt"
            filename_ana = path + "dt_run_" + dts[i] + "_ana.txt"

            with open(filename_num,"r") as infile:
                lines = infile.readlines()
                vals = lines[5].split()
                energy_n[i] = float(vals[1])
                variance_n[i] = float(vals[2])
                time_n[i] = float(vals[3])

            with open(filename_ana,"r") as infile:
                lines = infile.readlines()
                vals = lines[5].split()
                energy[i] = float(vals[1])
                variance[i] = float(vals[2])
                time[i] = float(vals[3])


        with open(outfile_num,"w") as outfile:
            for i in range(len(dts)):
                outfile.write(dts[i]+ " " + str(energy_n[i]) + " " + str(np.sqrt(variance_n[i])) + " " + str(time_n[i])+"\n")

        with open(outfile_ana,"w") as outfile:
            for i in range(len(dts)):
                outfile.write(dts[i]+ " " + str(energy[i]) + " " + str(np.sqrt(variance[i])) + " " + str(time[i])+"\n")


if __name__ == "__main__":
    print("Possible commands: 'b' (simplest), 'c' (importance sampling), 'd' (gradient descent evolution), \n" \
        + "'e' (no-repulsion), 'g' (repulsion), 'h' (one body densities) and 'delta_t' (find optimal delta t)")
    task = input("Which task to make plots for? ")
    assert task in ["b", "c", "d", "e", "g", "h", "deltat"], "Input not recognized."
    make_plots(task)
