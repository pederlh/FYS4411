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

int N_, num_alphas_, MC_, D_, type_energy_, type_sampling_, thread_ID_;
double h_, sum_, step_, D_diff_;
double *alphas_,*energies_, *variances_, *E_L_to_file_;
random_device rd_;

//Pointer to member function
void (Solver::*MC_method)();
void (Solver::*metropolis_sampling)(double alpha);


public:

//Constructor
Solver(int N, int num_alphas, int MC, int D, int type_energy, int type_sampling, int num_threads);

//Makes use of functions from class Psi;
Psi wave;

void Metropolis(double alpha);
void Metropolis_importance(double alpha);
void Metropolis_importance_interaction(double alpha);
//void MonteCarlo_burn();
void MonteCarlo();
void MonteCarlo2(double alpha, double *energies);
void MonteCarlo_GD(double *values, double alpha);
void MonteCarlo_GD_interaction(double *values, double alpha, string path);
void Gradient_descent();
void Gradient_descent_interaction();
double Greens_function(int idx);
// void ADAM();

void Write_to_file(string outfilename, double time);
// void Write_array_to_file(string outfilename, double *array, int len);
};

#endif
