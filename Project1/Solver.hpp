#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iomanip>
#include <cstdlib>
#include "time.h"
#include <valarray>
#include "omp.h"

#include "Psi.hpp"
using namespace std;


class Solver {

private:

int N_, num_alphas_, MC_, MC_optimal_run_, D_, type_energy_, type_sampling_, thread_ID_, equi_cycles_, OBD_;
double h_, sum_, step_, D_diff_, tol_GD_, eta_GD_, radi_;
double *alphas_,*energies_, *variances_, *E_L_to_file_;
random_device rd_;

//Pointer to member function
void (Solver::*MC_method)();
void (Solver::*metropolis_sampling)(double alpha);
void (Solver::*Interaction_or_not_GD)(double *values, double alpha);
void (Solver::*Interaction_or_not_optimal)(double alpha, double *energies);


public:

//Constructor
Solver(int N, int num_alphas, int MC, int MC_optimal_run, int D, int type_energy, int type_sampling, int num_threads);

//Makes use of functions from class Psi;
Psi wave;

void Metropolis(double alpha);
void Metropolis_importance(double alpha);
void Metropolis_interaction(double alpha);

void MonteCarlo_alpha_list();
void MonteCarlo(double alpha, double *energies);
void MonteCarlo_interaction(double alpha, double *energies);

void MonteCarlo_GD(double *values, double alpha);
void MonteCarlo_GD_interaction(double *values, double alpha);

void Gradient_descent();
double Greens_function(int idx);
void One_body_density(double *bins);
// void ADAM();

void Write_to_file(string outfilename, double time);
// void Write_array_to_file(string outfilename, double *array, int len);
};

#endif
