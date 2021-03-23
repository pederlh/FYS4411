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

