import os, sys
from make_plots import make_plots


def move_files_ask_plot(task, path, filenames):
    # Function for moving files and asking if plotting is wanted
    input("Press Enter to move data ")
    os.system("rm " + path + "*.txt")
    os.system("mv " + filenames + " " + path)

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
task_prompt = input("Which task to run? Choices are 'b' (simplest), 'c' (importance sampling), 'e' (gradient descent + blocking), 'g' (repulsion) or 'h' (one body densities): ")
if task_prompt not in ['b' , 'c', 'e', 'g', 'h']:
    print("Input not recognized. Aborting.")
    sys.exit(1)

print("Start producing data...")

""" Start producing and moving data """

if task_prompt == "b":

    num_particles = [1,10,50,100]
    dimentions = [1,2,3]
    type_energy = [0, 1]
    mc_cycles = 1000
    type_sampling = 0
    num_threads = 1
    OBD_check = 0
    for c in type_energy:
        for n in num_particles:
            for d in dimentions:
                os.system("./main.out " + str(n) + " " + str(d) + " " + str(mc_cycles) + " " + str(c) + " " \
                                        + str(type_sampling) + " " + str(num_threads) + " " + str(OBD_check))

    path = "./Results/1b_simple_noninteracting/"
    filenames = "spherical*.txt"
    move_files_ask_plot("b", path, filenames)


if task_prompt == "c":

    num_particles = 1
    dimentions = 1
    mc_cycles = 1000
    type_energy = 1
    type_sampling = 1
    num_threads = 1
    OBD_check = 0

    os.system("./main.out " + str(num_particles) + " " + str(dimentions) + " " + str(mc_cycles) + " " + str(type_energy) + " " \
                                    + str(type_sampling) + " " + str(num_threads) + " " + str(OBD_check))

    path = "./Results/1c_implementing_importance_sampling/"
    filenames = "importance_spherical*.txt"
    move_files_ask_plot("c", path, filenames)


if task_prompt == "e":

    num_particles = 5
    dimentions = 2
    mc_cycles = 1000
    type_energy = 1
    type_sampling = 2
    num_threads = 1
    OBD_check = 0
    mc_cycles_optimal_run = 1e17

    os.system("./main.out " + str(num_particles) + " " + str(dimentions) + " " + str(mc_cycles) + " " + str(type_energy) + " " \
                                    + str(type_sampling) + " " + str(num_threads) + " " + str(OBD_check) + " " + str(mc_cycles_optimal_run))

    path = "./Results/1e_implementing_gradient_descent_and_blocking/"
    filenames = "OPTIMAL_ALPHA*.txt"
    move_files_ask_plot("e", path, filenames)

if task_prompt == "g":
    pass

if task_prompt == "h":
    pass
