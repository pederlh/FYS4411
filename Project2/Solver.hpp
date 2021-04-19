#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iomanip>
#include <cstdlib>
#include "omp.h"

#include "Psi.hpp"
using namespace std;

/*
        Class for executing Montecarlo simulations using metropolis sampling.
        Handles interacting and non-interacting particles.
        Includes functions for:
            - Monte Carlo simulation
            - Metropolis sampling
            - Metropolis-Hastings sampling
            - Gradient descent for variational parameter estimation
*/


class Solver {

private:

int N_, num_alphas_, MC_, MC_optimal_run_, D_, type_energy_, type_sampling_, thread_ID_, equi_cycles_;
double h_, step_, D_diff_, tol_GD_, eta_GD_, radi_, time_, start_time_, end_time_;
double *alphas_,*energies_, *variances_, *E_L_to_file_, *alpha_list_times_;
bool OBD_check_;
random_device rd_;

//Pointer to member function
void (Solver::*main_method)(double *shared_alphas);
void (Solver::*metropolis_sampling)(double alpha);
void (Solver::*Interaction_or_not_GD)(double *values, double alpha);
void (Solver::*Interaction_or_not_optimal)(double alpha, double *energies);

public:

//Constructor
Solver(int N, int MC, int MC_optimal_run, int D, int type_energy, int type_sampling, int num_threads, double learning_rate, double*shared_alphas);

//Makes use of functions from class Psi;
Psi wave;

// Sampling methods
void Metropolis(double alpha);
void Metropolis_importance(double alpha);

// Main methods: Find optimal alpha by looping over list or gradient descent
void Gradient_descent(double *shared_alphas);

// During GD process: Interacting or noninteracting case
void MonteCarlo_GD_noninteracting(double *values, double alpha);

// Final run after GD process: Interacting or noninteracting case
void MonteCarlo_optval_noninteracting(double alpha, double *energies);

// Other methods
void Equilibrate(double alpha);                     // Perform equilibration cycles
double Greens_function(int idx);                    // Evaluate Green's function
void Write_to_file(string outfilename);// Write results to file

};

#endif
