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
task_prompt = input("Which task to run? Choices are 'b' (simplest), 'c' (importance sampling), 'd' (gradient descent), 'f' (repulsion) or 'g' (one body densities): ")
if task_prompt not in ['b' , 'c', 'd', 'f', 'g']:
    print("Input not recognized. Aborting.")
    sys.exit(1)

print("Start producing data...")

""" Start producing and moving data """

if task_prompt == "b":

    num_particles = [1,4,6]
    dimentions = [1, 2, 3]

    mc_cycles = 100000
    type_energy = 1
    type_sampling = 0
    num_threads = 1
    OBD_check = 0

    for n in num_particles:
        for d in dimentions:
            os.system("./main.out " + str(n) + " " + str(d) + " " + str(mc_cycles) + " " + str(type_energy) + " " \
                                    + str(type_sampling) + " " + str(num_threads) + " " + str(OBD_check))

    path = "./Results/1b_simple_noninteracting/"
    filenames = "spherical*.txt"
    move_files_ask_plot("b", path, filenames)


if task_prompt == "c":

    # Fyll på her
    
    input("Press Enter to move data ")
    path = "./Results/1c_implementing_importance_sampling/"
    os.system("rm " + path + "*.txt")
    if input("Press 'y' to make plots ") in ['y', 'Y']:
        make_plots("c")

    # os.system("mv spherical*.txt " + path)

if task_prompt == "d":
    pass

if task_prompt == "e":
    pass

if task_prompt == "f":
    pass

if task_prompt == "g":
    pass 
