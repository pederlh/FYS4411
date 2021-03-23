import os, sys
import numpy as np
import matplotlib.pyplot as plt

""" Compilation """

mac_linux_prompt = input("Run on macOS ('m') or Linux ('l') ? ")
if mac_linux_prompt == 'l':
    os.system("g++  -o main.out Solver.cpp Psi.cpp main.cpp -std=c++11 -fopenmp")
elif mac_linux_prompt == 'm':
    os.system("g++-10 -o main.out Solver.cpp Psi.cpp main.cpp -std=c++11 -fopenmp")
else:
    print("Input not recognized. Aborting.")
    sys.exit(1)

os.system("echo Compiling finished. Executing...")

""" Execute """

num_particles = 10
dimentions = 3
mc_cycles = 1000
type_energy = 0
type_sampling = 3
num_threads = 2
OBD_check = 0
mc_cycles_optimal_run = 2**17

os.system("./main.out " + str(num_particles) + " " + str(dimentions) + " " + str(mc_cycles) + " " + str(type_energy) + " " \
                        + str(type_sampling) + " " + str(num_threads) + " " + str(OBD_check) + " " + str(mc_cycles_optimal_run))
