import numpy as np 
import matplotlib.pyplot as plt

def make_plots(task):
    if task == "c":
        print("2 1  2 3 2 1 2 3 2 11")
    else:
        print('lol')

if __name__ == "__main__":
    task = input("Which task to make plots for? \n" 
                + "Choices are 'b' (simplest), 'c' (importance sampling), 'd' (gradient descent), 'f' (repulsion) or 'g' (one body densities): ")
    make_plots(task)
