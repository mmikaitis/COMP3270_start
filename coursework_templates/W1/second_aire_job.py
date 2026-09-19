import numpy as np

# Number of samples
N=10000

# Uniform random numbers in (0,1)
uniform_arr = np.random.random(N)

with open("uniform_numpy_aire.txt", "w") as file:
    for item in uniform_arr:
        file.write(str(item) + "\n")

# Normal distribution
normal_arr = np.random.normal(loc=0, scale=1, size=N)

with open("normal_numpy_aire.txt", "w") as file:
    for item in normal_arr:
        file.write(str(item) + "\n")

