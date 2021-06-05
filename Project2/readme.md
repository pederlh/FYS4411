# The restricted Boltzmannn machine applied to the quantum many body problem
## Project 2 - Computational Physics II
 The final report, containing background theory, explanation of the methods and results, can be found [here](https://github.com/pederlh/FYS4411/tree/main/Project2/Article).

#### Summary
In this project we apply a Gaussian-Binary RBM to the problem of two interacting bosons confined to a harmonic oscillator trap. We used Monte-Carlo simulations together with Metropolis-Hastings sampling to estimate the expectation value of the local energy. We implemented methods for performing gradient descent on the variational parameters of the NQS wave function. In summary the code contained in this repository can be used to

* Run RBM-VMC simulations on 2-particle systems using both brute force Metropolis and Metropolis- Hastings sampling.
* Estimate the expectation value of the local energy for 2-particle systems
* Perform gradient descent on the variational parameteres of the NQS-trial wave function


#### Try me
In order to test the code easily we have provided a python script which can be run to produce samples for the following systems

* 1 particle, 1 dimension
* 2 particles, 2 dimensions (with/without Coloumb interaction)

For both cases ADAM is used for performing gradient descent on the parameters and the HO-frequency is set to 1. For other specifications see the script ``TryMe.py`` where all parameters can easily be changed if other combinations are of interest.

For running the script simply enter in your terminal

```console
python3 TryMe.py
```  
You will then be prompted with questions concerning the number of particles and dimensions as well as interaction. The script will then run the RBM class and write samples of the local energy to a file. The blocking method is then applied and the mean and standard deviation of the dataset is printed to screen.


### Useful links for the eager reader

[Lecture notes and theory](http://compphysics.github.io/ComputationalPhysics2/doc/LectureNotes/_build/html/intro.html)

[Project description](http://compphysics.github.io/ComputationalPhysics2/doc/Projects/2021/Project2/Project2ML/pdf/Project2ML.pdf)
