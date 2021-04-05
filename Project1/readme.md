# VMC studies of Bose gas in harmonic oscillator trap
## Project 1 - Computational Physics II
The report containing background theory, explanation of the methods and results can be found [here](https://github.com/pederlh/FYS4411/blob/main/Project1/Article/Project_1___VMC_studies_of_QM_systems_NOTFIN.pdf) and the problem set for the project can be found [here](http://compphysics.github.io/ComputationalPhysics2/doc/Projects/2021/Project1/pdf/Project1.pdf).

### Summary
In this project we tackled the problem of both non-interacting and interacting boson gas particles in spherical and elliptical harmonic oscillator well traps using variational Monte Carlo simulations. In summary the code contained in this repository can be used to

- Run VMC simulations on N-particle systems using both brute force Metropolis and Metropolis-Hastings sampling.
- Estimate the expectation value of the local energy for N-particle systems
- Estimate the one body density of N-particle systems
- Perform gradient descent on the variational parameter alpha of the trial wave function

### Structure
The main structure of our code lies in the two classes `Solver.cpp` and `Psi.cpp`. The first contains methods for VMC calculations, gradient descent and Metropolis sampling. The second contains functions for the trial wave function, calculation of local energy and methods for evaluating the positions of the particle and the quantum force acting on them. The `main.cpp` file is there to create instances of the Solver class and write results to file. In the scripts folder there are Python files used to generate data and plots found in the report linked at the top.

### Running instructions
Instruction to produce new data and/or results for the tasks listed in the problem set are as follows. The `main.cpp` file is organized to take 7-9 input arguments to initiate the class `Solver.cpp`, so when running the makefile one has to specify at least the first seven. The (ordered) input arguments are
- `N` number of particles :  integer > 0
-  `D` spacial dimentions : either 1, 2 or 3
- `MC` number of Monte Carlo cycles : integer > 0
- `E` type of energy calculation : 0 for analytical expression, 1 for numerical differentiation
- `S` type of Metropolis sampling : 0 for brute force, 1 for importance sampling, 2 for importance sampling + gradient descent (non-interacting particles) and 3 for importance sampling + gradient descent (interacting particles)
- `T` number of threads : 1 or larger (depends on the number of cores your computer has)
- `O` one body density calculations : 1 is true, 0 is false

and if necessary for the task you wish to perform
 - `OptMC` number of Monte Carlo cycles for optimal run after gradient descent : integer > 0
 - `Eta` learning rate of gradient descent : usually 0.001, but can take any number.

The following example makefile compiles and executes the program

``` Ruby
linux: compile_par execute
mac: compile_par_mac execute

compile_par:
	g++  -o main.out Solver.cpp Psi.cpp main.cpp -std=c++11 -fopenmp

compile_par_mac:
	g++-10 -o main.out Solver.cpp Psi.cpp main.cpp -std=c++11 -fopenmp

execute:
	./main.out N D MC E S T O OptMC Eta

```

by running from your terminal

```console
make operative_system
```

where your choose either `mac` or `linux` depending on your computers operative system. If the number of threads is set to 1, everything is run on a single thread. If not, the OpenMP API is used to parallelize and run the calculations on the chosen number of threads.

### Reproducing data
We will now show examples of how to reproduce data for some of the tasks found in the problem set.

#### 1. Non-interacting bosons for various values of alpha using VMC without importance sampling.
For `N = 100` particles in `D = 3` dimentions with `MC = 10 000` Monte Carlo cycles using `E = 0` (the analytical expression for calculating the local energy) and `T = 4` threads the execute line in your makefile should look something like

``` Ruby
execute:
  ./main.out 100 3 10000 0 0 4 0
```

#### 2. Interacting bosons using gradient descent and calculating one body density  
For `N = 100` particles in `D = 3` dimentions with `MC = 1000` Monte Carlo cycles for the gradient descent using `E = 0` (the analytical expression for calculating the local energy) and `T = 4` threads. We also specify the learning rate `Eta = 0.001` and `OptMC = 2**18` cycles to be performed by each thread using the optimal variational parameter found by the gradient descent. The execute line in your makefile should look something like

``` Ruby
execute:
  ./main.out 100 3 1000 0 3 4 1 2**18 0.001
```
