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

using namespace std;


class Solver {
private:

int N_, num_alphas_, MC_, D_;
double h_, sum_;
double* alphas_, *energies_, *variances_;

double tf_forward_, tf_backward_, tf_middle_, laplace_tf_;
double *r_forward_, *r_backward_;

//Pointer to member function
double (Solver::*energy_calculation)(double alpha, double *r);

double trial_func_1D(double alpha, double* r); //Test wave function
double local_energy_1D_analytical(double alpha, double *r); //Local energy in one dimension
void MonteCarlo();

double local_energy_1D_brute_force(double alpha, double *r);

public:

//Constructor
Solver(int N, int num_alphas, int MC, int D, int type_energy);
void write_to_file(string outfilename);

};

#endif
