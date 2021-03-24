import os, sys
from make_plots import make_plots


def move_files_ask_plot(task, path, filenames):
    # Function for moving files and asking if plotting is wanted
    input("Press Enter to move data ")
    os.system("rm " + path + "*.txt")
    os.system("mv " + filenames + " " + path)

    # Take special care of one body density files
    if task == "e":
        os.system("rm ./Results/1h_one_body_densities/One_body_density*.txt")
        os.system("mv  One_body_density*.txt ./Results/1h_one_body_densities/")
    # if type is d or g: move OBD

    if task == "g":
        os.system("rm ./Results/1h_one_body_densities/Interaction_One_body_density*.txt")
        os.system("mv  Interaction_One_body_density*.txt ./Results/1h_one_body_densities/")

    if input("Press 'y' to make plots ") in ['y', 'Y']:
        make_plots(task)


print("Note: The data produced will overwrite previous results")

""" Compilation """

mac_linux_prompt = input("Run on macOS ('m') or Linux ('l')? ")

if mac_linux_prompt == 'l':
    os.system("g++  -o main.out Solver.cpp Psi.cpp main.cpp -std=c++11 -fopenmp")
elif mac_linux_prompt == 'm':
    os.system("g++-10 -o main.out Solver.cpp Psi.cpp main.cpp -std=c++11 -fopenmp")
else:
    print("Input not recognized. Aborting.")
    sys.exit(1)


""" Choose task """
task_prompt = input("Which task to run? Choices are 'b' (simplest), 'c' (importance sampling), 'e' (gradient descent + blocking + one body density) or 'g' (everything + repulsion): ")
if task_prompt not in ['b' , 'c', 'e', 'g']:
    print("Input not recognized. Aborting.")
    sys.exit(1)

print("Start producing data...")

""" Start producing and moving data """

if task_prompt == "b":

    num_particles = [1,10,100,500]
    dimentions = 3
    type_energy = [0, 1]
    mc_cycles = 10000
    type_sampling = 0
    num_threads = 1
    OBD_check = 0
    for c in type_energy:
        for n in num_particles:
            os.system("./main.out " + str(n) + " " + str(dimentions) + " " + str(mc_cycles) + " " + str(c) + " " \
                                    + str(type_sampling) + " " + str(num_threads) + " " + str(OBD_check))

    path = "./Results/1b_simple_noninteracting/"
    filenames = "spherical*.txt"
    move_files_ask_plot("b", path, filenames)


if task_prompt == "c":

    num_particles = [1,10,100,500]
    dimentions = 3
    type_energy = [0, 1]
    mc_cycles = 10000
    type_sampling = 1
    num_threads = 1
    OBD_check = 0

    for c in type_energy:
        for n in num_particles:
            os.system("./main.out " + str(n) + " " + str(dimentions) + " " + str(mc_cycles) + " " + str(c) + " " \
                                    + str(type_sampling) + " " + str(num_threads) + " " + str(OBD_check))

    path = "./Results/1c_implementing_importance_sampling/"
    filenames = "importance_spherical*.txt"
    move_files_ask_plot("c", path, filenames)

    """
    Husk å fikse delta t !!!!!!!!!!!!!!!
    """

if task_prompt == "e":

    num_particles = [2,16,64,128]
    num_etas = [0.08,0.01,0.001,0.0003]
    dimentions = 3
    mc_cycles = 1000
    type_energy = 0
    type_sampling = 2
    num_threads = 2
    OBD_check = 1
    for i in range(len(num_particles)):
        mc_cycles_optimal_run = 2**19/num_particles[i]
        os.system("./main.out " + str(num_particles[i]) + " " + str(dimentions) + " " + str(mc_cycles) + " " + str(type_energy) + " " \
                                    + str(type_sampling) + " " + str(num_threads) + " " + str(OBD_check) + " " + str(mc_cycles_optimal_run) + " " + str(num_etas[i]))

    path = "./Results/1e_implementing_gradient_descent_and_blocking/"
    filenames = "OPTIMAL_ALPHA*.txt"
    move_files_ask_plot("e", path, filenames)

if task_prompt == "g":
    num_particles = [2,16,32,128]
    num_etas = [0.08,0.0015,0.001,0.0003]
    dimentions = 3
    mc_cycles = 1000
    type_energy = 0
    type_sampling = 3
    num_threads = 2
    OBD_check = 1
    for i in range(len(num_particles)):
        mc_cycles_optimal_run = 2**19/num_particles[i]
        os.system("./main.out " + str(num_particles[i]) + " " + str(dimentions) + " " + str(mc_cycles) + " " + str(type_energy) + " " \
                                    + str(type_sampling) + " " + str(num_threads) + " " + str(OBD_check) + " " + str(mc_cycles_optimal_run) + " " + str(num_etas[i]))

    path = "./Results/1g_implementing_repulsion/"
    filenames = "INTERACTION_OPTIMAL_ALPHA*.txt"
    move_files_ask_plot("g", path, filenames)
