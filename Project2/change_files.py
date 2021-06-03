import os
import numpy as np

files  = os.listdir()

for file in files:
    if ".txt" in file:
        if not file == "a_joined.txt" and not file == "b_joined.txt" and not file == "w_joined.txt":
            vec = np.loadtxt(file)
            os.remove(file)
            file = file.replace(".txt",".npy")
            np.save(file,vec)
