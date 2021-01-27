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

int N_, num_alphas_, MC_;
double* alphas_, *energies_, *variances_;

double trial_func_1D(double alpha, double* r); //Test wave function
double local_energy_1D_analytical(); //Local energy in one dimension
void MonteCarlo();



double local_energy_1D(double alpha, double r);

public:

//Constructor
Solver(int N);

void write_to_file(string outfilename);

};

#endif
