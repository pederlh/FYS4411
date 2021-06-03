import numpy as np
import os
"""
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

new_filenames = ["er_01","er_015","er_02","er_025","er_03","er_035","er_04","er_045","er_05","er_055","er_06","er_065", "er_07","er_075","er_08","er_085","er_09","er_095","er_10","er_105","er_11"]

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


eta_test = os.listdir()
eta_test.remove("join_files.py")
eta_test.remove("merged_files")

#etas = [0.100000,0.150000,0.200000,0.250000,0.300000,0.350000,0.400000,0.450000,0.500000,0.550000,0.600000,0.650000,0.700000,0.750000,0.800000,0.850000,0.900000,0.950000,1.000000,1.050000,1.100000]

sort_eta = []
nums = []
its = []

for i in range(len(eta_test)):
    w = eta_test[i].split("_")
    for j in range(i+1,len(eta_test)):
        w2 = eta_test[j].split("_")
        if w[1] == w2[1]:
            sort_eta.append(eta_test[i])
            sort_eta.append(eta_test[j])
            nums.append(w[1])
            nums.append(w2[1])
            its.append(float(w[3]))
            its.append(float(w2[3]))

sort_n = np.zeros(len(nums))
for i in range(len(nums)):
    sort_n[i] = float(nums[i])

sort_n, sort_eta,its = zip(*sorted(zip(sort_n, sort_eta,its)))

new_filenames = ["enn_01","enn_015","enn_02","enn_025","enn_03","enn_035","enn_04","enn_045","enn_05","enn_055","enn_06","enn_065", "enn_07","enn_075","enn_08","enn_085","enn_09","enn_095","enn_10","enn_105","enn_11"]


j=0
for i in range(0,len(sort_n),2):
    data =[]
    print(i)
    # Reading data from file1
    v1 = np.load(sort_eta[i])
    v2 = np.load(sort_eta[i+1])

    mean_its = (its[i]+its[i+1])/2
    outfile = new_filenames[j] +  "_i_" + str(mean_its) + "_.npy"

    vec = np.concatenate((v1,v2))
    np.save(outfile,vec)

    j+=1


eta_test = os.listdir()
eta_test.remove("join_files.py")
eta_test.remove("merged_files")

for i in range(len(eta_test)):
    w1 = eta_test[i].split("_")
    for j in range(i+1,len(eta_test)):
        w2 = eta_test[j].split("_")
        if w1[1] == w2[1]:
            i1 = float(w1[-2])
            i2 = float(w2[-2])
            im = (i1+i2)/2
            new_name = "e_"+w1[1]+"_i_" + str(im)+"_.npy"
            v1 = np.load(eta_test[i])
            v2 = np.load(eta_test[j])
            v = np.concatenate((v1,v2))
            np.save(new_name,v)

for file in eta_test:
    vec = np.loadtxt(file)
    file = file.replace(".txt",".npy")
    np.save(file,vec)
"""
