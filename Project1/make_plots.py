import numpy as np
import matplotlib.pyplot as plt
import os


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

        #CONCAT
        """
        1. Hente opp GD utviklingen til alpha for string 0
        2. Plotte utviklingen for de forskjellige N verdiene


        1. Hente alle filer UTEN INTERACTION FORRAN
        2. Dele de opp etter N
        3. Merge alle filer med samme N
        4. Bootstrap/blocking
        5. Flotte estimater
        """
        lol = 1

    if task == "g":
        #CONCAT
        """
        1. Hente alle filer UTEN INTERACTION FORRAN
        2. Dele de opp etter N
        3. Merge alle filer med samme N
        4. Bootstrap/blocking
        5. Flotte estimater
        """
        lol=1

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

        r = np.linspace(0,5,100)

        plt.title("N = 2")
        plt.plot(r,joined_OBD_N2,label="No interaction")
        plt.plot(r,joined_OBD_N2_I,label="Interaction")

        OBD_integrate = np.trapz(joined_OBD_N2,r)
        A = OBD_integrate*4/np.sqrt(np.pi)
        ana_rho = A*(r**2)*np.exp(-r**2)

        plt.plot(r,ana_rho, label = "Analytical density")
        plt.xlabel("r")
        plt.ylabel("rho(r)")
        plt.legend()
        plt.show()

        plt.title("N = 16")
        plt.plot(r,joined_OBD_N16,label="No interaction")
        plt.plot(r,joined_OBD_N16_I,label="Interaction")
        OBD_integrate = np.trapz(joined_OBD_N16,r)
        A = OBD_integrate*4/np.sqrt(np.pi)
        ana_rho = A*(r**2)*np.exp(-r**2)

        plt.plot(r,ana_rho, label = "Analytical density")
        plt.xlabel("r")
        plt.ylabel("rho(r)")
        plt.legend()
        plt.show()

        plt.title("N = 32")
        plt.plot(r,joined_OBD_N32,label="No interaction")
        plt.plot(r,joined_OBD_N32_I,label="Interaction")
        OBD_integrate = np.trapz(joined_OBD_N32,r)
        A = OBD_integrate*4/np.sqrt(np.pi)
        ana_rho = A*(r**2)*np.exp(-r**2)

        plt.plot(r,ana_rho, label = "Analytical density")
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










if __name__ == "__main__":
    task = input("Which task to make plots for? \n"
                + "Choices are 'b' (simplest), 'c' (importance sampling), 'd' (gradient descent), 'f' (repulsion) or 'g' (one body densities): ")
    make_plots(task)
