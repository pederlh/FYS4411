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

using namespace std;


class Solver {
private:

int N_, num_alphas_, MC_, D_;
double h_, sum_, tmp_, step_;
double *alphas_,*energies_, *variances_;

double tf_middle_, laplace_tf_;
double **r_old_, *r_new_;
double **numerical_matrix_;


//Pointer to member function
double (Solver::*energy_calculation)(double alpha, double **r, double r_sum);

double trial_func(double alpha, double sum_r_squared);                 //Test wave function




double init_r_sum(double **r);                                  //Calculates the sum of the square of all posistions
double update_r_sum(double sum, double r_init, double r_move);  //Updates the sum of the square of all positions

double local_energy_brute_force(double alpha, double **r, double r_sum);
double local_energy_analytical(double alpha, double **r, double r_sum);     //Local energy analytical


public:

//Constructor
Solver(int N, int num_alphas, int MC, int D, int type_energy);
void write_to_file(string outfilename, double time);
void MonteCarlo();

};

#endif
