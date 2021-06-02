import numpy as np
import os

eta_test = os.listdir()
eta_test.remove("join_files.py")
eta_test.remove("merged_files")
print(eta_test)
etas = [0.100000,0.150000,0.200000,0.250000,0.300000,0.350000,0.400000,0.450000,0.500000,0.550000,0.600000,0.650000,0.700000,0.750000,0.800000,0.850000,0.900000,0.950000,1.000000,1.050000,1.100000]

def is_float(n):
    try:
        float_n = float(n)
    except ValueError:
        return False
    else:
        return True

nums = []
its = []
for files in eta_test:
    for words in files.split("_"):
        if not words.isdigit():
            if is_float(words):
                nums.append(float(words))
        else:
            its.append(int(words))


nums, eta_test,its = zip(*sorted(zip(nums, eta_test,its)))

new_filenames = ["e_01","e_015","e_02","e_025","e_03","e_035","e_04","e_045","e_05","e_055","e_06","e_065", "e_07","e_075","e_08","e_085","e_09","e_095","e_10","e_105","e_11"]

j=0
for i in range(0,len(nums),2):
    data = data2 = ""
    # Reading data from file1
    with open(eta_test[i]) as fp:
        data = fp.read()

    # Reading data from file2
    with open(eta_test[i+1]) as fp:
        data2 = fp.read()

    mean_its = (its[i]+its[i+1])/2
    print(i)
    # Merging 2 files
    # To add the data of file2
    # from next line
    data += "\n"
    data += data2

    outfile = new_filenames[j] +  "_i_" + str(mean_its) + "_.txt"

    with open(outfile, 'w') as fp:
        fp.write(data)

    j+=1
