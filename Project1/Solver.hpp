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


#include "Psi.hpp"
using namespace std;


class Solver {
private:



int N_, num_alphas_, MC_, D_, type_energy_, type_sampling_;
double h_, sum_, step_, D_diff_;
double *alphas_,*energies_, *variances_;

random_device rd_;


//Pointer to member function
void (Solver::*MC_method)();
void (Solver::*metropolis_sampling)(double alpha);



public:

//Constructor
Solver(int N, int num_alphas, int MC, int D, int type_energy, int type_sampling);

//Makes use of functions from class Psi;
Psi wave;

void Metropolis(double alpha);
void Metropolis_importance(double alpha);
void Write_to_file(string outfilename, double time);
//void MonteCarlo_burn();
void MonteCarlo();
void MonteCarlo_SGD(double *values, double alpha);
void Gradient_descent();
double Greens_function(int idx);
void ADAM();

};

#endif
