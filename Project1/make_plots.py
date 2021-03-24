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

    if task == "h":
        path = "./Results/1h_one_body_densities/"
        os.chdir(path)
        OBD_list = os.listdir()
        OBD_I_list = []


        OBD_list = [s for s in OBD_list if "stringID_0" in s]
        OBD_I_list = [s for s in OBD_list if "Interaction" in s]
        OBD_list = [s for s in OBD_list if "Interaction" not in s]


        alpha_list = [OBD_list[i].split("_")[-1] for i in range(len(OBD_list))]
        alpha_list_I = [OBD_I_list[i].split("_")[-1] for i in range(len(OBD_I_list))]

        alpha_list = [float(alpha.split(".txt")[0]) for alpha in alpha_list]
        alpha_list_I = [float(alpha.split(".txt")[0]) for alpha in alpha_list_I]


        N_list = [float(OBD_list[i].split("_")[4]) for i in range(len(OBD_list))]
        N_list_I = [float(OBD_I_list[i].split("_")[4]) for i in range(len(OBD_I_list))]


        r = 0
        for i in range(len(OBD_list)):
            filename = OBD_list[i]
            N = N_list[i]
            alpha = alpha_list[i]
            OBD = np.loadtxt(filename)
            if i == 0:
                r = np.linspace(0,5,len(OBD))
            plt.plot(r,OBD)
            plt.show()
            """
            filename = OBD_I_list[i]
            N = N_list_I[i]
            alpha = alpha_list_I[i]
            OBD = np.loadtxt(filename)
            plt.plot(r,OBD)
            plt.show()


            """


        # Vil kunne plotte OBD for både interacting og ikke int med alpha verdi og N i plottet FOR EN STRING.






if __name__ == "__main__":
    task = input("Which task to make plots for? \n"
                + "Choices are 'b' (simplest), 'c' (importance sampling), 'd' (gradient descent), 'f' (repulsion) or 'g' (one body densities): ")
    make_plots(task)
