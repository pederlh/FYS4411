import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd


def make_plots(task):
    if task == "b":
        N = [1,10,50,100]
        MC = 1000
        string_ID = 0
        path = "./Results/1b_simple_noninteracting/"

        alphas = []
        energies_a = []
        variances_a = []
        times_a = []
        energies_n = []
        variances_n = []
        times_n = []
        for n in N:
            filename_num = path + "spherical_HO_3" + "D_num" +"_N_"+ str(n) + "_MC_"+ str(MC) +"_stringID_0" + ".txt";
            filename_ana = path + "spherical_HO_3" + "D_ana" +"_N_"+ str(n) + "_MC_"+ str(MC) +"_stringID_0" + ".txt";

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

            plt.title("Number of particles = %i" % n)
            plt.errorbar(alphas, energies_a,yerr=std_a, fmt= "or",capsize=5, elinewidth=1, label = "Local energy (analytical)")
            plt.errorbar(alphas, energies_n,yerr=std_n, fmt= "ok",capsize=5, elinewidth=1, label = "Local energy (numerical)")
            plt.axhline(y=analytical_E, color='b', linestyle='-', label = "Local energy (exact)")
            plt.xlabel("alpha")
            plt.ylabel("Local energy")
            plt.legend()
            plt.show()


            alphas = []
            energies_a = []
            variances_a = []
            times_a = []
            energies_n = []
            variances_n = []
            times_n = []

    if task == "c":
        N = [1,10,100,500]
        MC = 10000
        string_ID = 0
        path = "./Results/1c_implementing_importance_sampling/"

        alphas = []
        energies_a = []
        variances_a = []
        times_a = []
        energies_n = []
        variances_n = []
        times_n = []
        for n in N:
            filename_num = path + "importance_spherical_HO_3" + "D_num" +"_N_"+ str(n) + "_MC_"+ str(MC) +"_stringID_0" + ".txt";
            filename_ana = path + "importance_spherical_HO_3" + "D_ana" +"_N_"+ str(n) + "_MC_"+ str(MC) +"_stringID_0" + ".txt";

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

            plt.title("Number of particles = %i" % n)
            plt.errorbar(alphas, energies_a,yerr=std_a, fmt= "or",capsize=5, elinewidth=1, label = "Local energy (analytical)")
            plt.errorbar(alphas, energies_n,yerr=std_n, fmt= "ok",capsize=5, elinewidth=1, label = "Local energy (numerical)")
            plt.axhline(y=analytical_E, color='b', linestyle='-', label = "Local energy (exact)")
            plt.xlabel("alpha")
            plt.ylabel("Local energy")
            plt.legend()
            plt.show()

            alphas = []
            energies_a = []
            variances_a = []
            times_a = []
            energies_n = []
            variances_n = []
            times_n = []

    if task == "e":
        """
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

        GD_2 = np.loadtxt(GD[0])[1:]
        GD_16 = np.loadtxt(GD[1])[1:]
        GD_32 = np.loadtxt(GD[2])[1:]
        alphas= [GD_2,GD_16,GD_32]

        it_2 = np.linspace(1,len(GD_2),len(GD_2))
        it_16 = np.linspace(1,len(GD_16),len(GD_16))
        it_32 = np.linspace(1,len(GD_32),len(GD_32))
        its = [it_2,it_16,it_32]

        for i in range(3):
            plt.plot(its[i], alphas[i], "*",label = "N = %i" %numbers[i])
            plt.plot(its[i][-1], alphas[i][-1],"*", markersize = 12)

        plt.ylabel("Alpha value")
        plt.xlabel("Number of iterations")
        plt.title("Gradient descent for non-interacting bosons")
        plt.legend()
        plt.show()
        """

        print("Blocking results from GD with repulsion")

        path = "./Results/1e_implementing_gradient_descent_and_blocking/"
        os.chdir(path)
        runs = os.listdir()

        t2 = []
        t16 = []
        t32 = []

        N32 = [s for s in runs if "32_N_" in s]
        runs = [s for s in runs if "32_N_" not in s]
        N2 = [s for s in runs if "2_N_" in s]
        N16 = [s for s in runs if "16_N_" in s]

        thread0 = np.loadtxt(N2[0]);
        t2.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N2[1]);
        t2.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N2[2]);
        t2.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N2[3]);
        t2.append(thread3[0])
        thread3 = thread3[1:]

        joined_N2 = np.concatenate((thread0,thread1,thread2,thread3))

        thread0 = np.loadtxt(N16[0]);
        t2.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N16[1]);
        t2.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N16[2]);
        t2.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N16[3]);
        t2.append(thread3[0])
        thread3 = thread3[1:]

        joined_N16 = np.concatenate((thread0,thread1,thread2,thread3))

        thread0 = np.loadtxt(N32[0]);
        t2.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N32[1]);
        t2.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N32[2]);
        t2.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N32[3]);
        t2.append(thread3[0])
        thread3 = thread3[1:]

        joined_N32 = np.concatenate((thread0,thread1,thread2,thread3))

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


    if task == "g":

        print("Blocking results from GD with repulsion")

        path = "./Results/1g_implementing_repulsion/"
        os.chdir(path)
        runs = os.listdir()

        t2 = []
        t16 = []
        t32 = []

        N32 = [s for s in runs if "32_N_" in s]
        runs = [s for s in runs if "32_N_" not in s]
        N2 = [s for s in runs if "2_N_" in s]
        N16 = [s for s in runs if "16_N_" in s]

        thread0 = np.loadtxt(N2[0]);
        t2.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N2[1]);
        t2.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N2[2]);
        t2.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N2[3]);
        t2.append(thread3[0])
        thread3 = thread3[1:]

        joined_N2 = np.concatenate((thread0,thread1,thread2,thread3))

        thread0 = np.loadtxt(N16[0]);
        t2.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N16[1]);
        t2.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N16[2]);
        t2.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N16[3]);
        t2.append(thread3[0])
        thread3 = thread3[1:]

        joined_N16 = np.concatenate((thread0,thread1,thread2,thread3))

        thread0 = np.loadtxt(N32[0]);
        t2.append(thread0[0])
        thread0 = thread0[1:]
        thread1 = np.loadtxt(N32[1]);
        t2.append(thread1[0])
        thread1 = thread1[1:]
        thread2 = np.loadtxt(N32[2]);
        t2.append(thread2[0])
        thread2 = thread2[1:]
        thread3 = np.loadtxt(N32[3]);
        t2.append(thread3[0])
        thread3 = thread3[1:]

        joined_N32 = np.concatenate((thread0,thread1,thread2,thread3))

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

    if task == "h":

        path = "./Results/1h_one_body_densities/"
        os.chdir(path)
        OBD_list = os.listdir()

        OBD_list_I = [s for s in OBD_list if "Interaction" in s]
        OBD_list = [s for s in OBD_list if "Interaction" not in s]

        OBD_N2 = [s for s in OBD_list if "_N_2_" in s]
        OBD_N16 = [s for s in OBD_list if "_N_16_" in s]
        OBD_N32 = [s for s in OBD_list if "_N_32_" in s]

        OBD_N2_I = [s for s in OBD_list_I if "_N_2_" in s]
        OBD_N16_I = [s for s in OBD_list_I if "_N_16_" in s]
        OBD_N32_I = [s for s in OBD_list_I if "_N_32_" in s]

        thread0 = np.loadtxt(OBD_N2[0]);
        thread1 = np.loadtxt(OBD_N2[1]);
        thread2 = np.loadtxt(OBD_N2[2]);
        thread3 = np.loadtxt(OBD_N2[3]);

        joined_OBD_N2 = (thread0 + thread1 + thread2 + thread3)/4

        thread0 = np.loadtxt(OBD_N16[0]);
        thread1 = np.loadtxt(OBD_N16[1]);
        thread2 = np.loadtxt(OBD_N16[2]);
        thread3 = np.loadtxt(OBD_N16[3]);

        joined_OBD_N16 = (thread0 + thread1 + thread2 + thread3)/4

        thread0 = np.loadtxt(OBD_N32[0]);
        thread1 = np.loadtxt(OBD_N32[1]);
        thread2 = np.loadtxt(OBD_N32[2]);
        thread3 = np.loadtxt(OBD_N32[3]);

        joined_OBD_N32 = (thread0 + thread1 + thread2 + thread3)/4

        thread0 = np.loadtxt(OBD_N2_I[0]);
        thread1 = np.loadtxt(OBD_N2_I[1]);
        thread2 = np.loadtxt(OBD_N2_I[2]);
        thread3 = np.loadtxt(OBD_N2_I[3]);

        joined_OBD_N2_I = (thread0 + thread1 + thread2 + thread3)/4

        thread0 = np.loadtxt(OBD_N16_I[0]);
        thread1 = np.loadtxt(OBD_N16_I[1]);
        thread2 = np.loadtxt(OBD_N16_I[2]);
        thread3 = np.loadtxt(OBD_N16_I[3]);

        joined_OBD_N16_I = (thread0 + thread1 + thread2 + thread3)/4

        thread0 = np.loadtxt(OBD_N32_I[0]);
        thread1 = np.loadtxt(OBD_N32_I[1]);
        thread2 = np.loadtxt(OBD_N32_I[2]);
        thread3 = np.loadtxt(OBD_N32_I[3]);

        joined_OBD_N32_I = (thread0 + thread1 + thread2 + thread3)/4

        r = np.linspace(0,8,50)

        plt.title("N = 2")
        A = 2/np.trapz(joined_OBD_N2*r**2,r)
        plt.plot(r,A*joined_OBD_N2,label="No interaction")
        #plt.plot(r,joined_OBD_N2_I,label="Interaction")

        better_ana_rho = 2 * (np.sqrt(np.pi) ** (-3) * 4 * np.pi * np.exp(-r ** 2))
        plt.plot(r,better_ana_rho, label = "GOOD? rho")


        #plt.plot(r,ana_rho, label = "Analytical density")
        plt.xlabel("r")
        plt.ylabel("rho(r)")
        plt.legend()
        plt.show()




        plt.title("N = 16")
        A = 16/np.trapz(joined_OBD_N16*r**2,r)
        plt.plot(r,A*joined_OBD_N16,label="No interaction")
        #plt.plot(r,joined_OBD_N16_I,label="Interaction")

        better_ana_rho = 16 * (np.sqrt(np.pi) ** (-3) * 4 * np.pi * np.exp(-r ** 2))
        plt.plot(r,better_ana_rho, label = "GOOD? rho")

        #plt.plot(r,ana_rho, label = "Analytical density")
        plt.xlabel("r")
        plt.ylabel("rho(r)")
        plt.legend()
        plt.show()




        plt.title("N = 32")
        A = 32/np.trapz(joined_OBD_N32*r**2,r)
        plt.plot(r,A*joined_OBD_N32,label="No interaction")
        #plt.plot(r,joined_OBD_N32_I,label="Interaction")

        better_ana_rho = 32 * (np.sqrt(np.pi) ** (-3) * 4 * np.pi * np.exp(-r ** 2))
        plt.plot(r,better_ana_rho, label = "GOOD? rho")

        #plt.plot(r,ana_rho, label = "Analytical density")
        plt.xlabel("r")
        plt.ylabel("rho(r)")
        plt.legend()
        plt.show()

        """
        for i in range(len(OBD_list)):
            filename = OBD_list[i]
            N = N_list_new[i]
            alpha = alpha_list[i]
            OBD = np.loadtxt(filename)
            plt.plot(r,OBD,label="Non_interaction")

            OBD_integrate = np.trapz(OBD,r)
            A = OBD_integrate*np.sqrt(np.pi)/4

            filename_I = OBD_I_list[i]
            N_I = N_list_I_new[i]
            alpha_I = alpha_list_I[i]
            OBD_I = np.loadtxt(filename_I)
            plt.plot(r,OBD_I, label="Interacting")

            real_alpha = 0.5
            ana_rho = A*(r**2)*np.exp(-r**2)
            rho_integrate = np.trapz(ana_rho,r)
            print(rho_integrate)

            plt.plot(r,ana_rho, label = "Analytical density")


            plt.legend()
            plt.title("N = %i "% N_I)
            plt.xlabel("r")
            plt.ylabel("rho(r)")
            plt.show()

            """



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









if __name__ == "__main__":
    task = input("Which task to make plots for? \n"
                + "Choices are 'b' (simplest), 'c' (importance sampling), 'd' (gradient descent), 'f' (repulsion) or 'g' (one body densities): ")
    make_plots(task)
