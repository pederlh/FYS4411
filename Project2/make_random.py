import numpy as np

np.random.seed(1)
nums = 1000
N = 2


q_factors       = [str(np.random.normal()) + "\n" for i in range(nums)]
rs              = [str(np.random.rand()) + "\n" for i in range(nums)]
idxs            = [str(np.random.randint(2)) + "\n" for i in range(nums)]
random_gausses  = [str(np.random.normal()) + "\n" for i in range(nums)]
acceptances     = [str(np.random.rand()) + "\n" for i in range(nums)]

with open("q_factors.txt", "w") as file1:
    file1.writelines(q_factors)

with open("rs.txt", "w") as file2:
    file2.writelines(rs)

with open("idxs.txt", "w") as file3:
    file3.writelines(idxs)

with open("random_gausses.txt", "w") as file4:
    file4.writelines(random_gausses)

with open("acceptances.txt", "w") as file5:
    file5.writelines(acceptances)
