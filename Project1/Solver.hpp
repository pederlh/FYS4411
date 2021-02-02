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
double h_, sum_, tmp_;
double *alphas_,*energies_, *variances_;

double tf_middle_, laplace_tf_;
double **r_old_, *r_new_;
double **numerical_matrix_;


//Pointer to member function
double (Solver::*energy_calculation)(double alpha, double **r, double r_sum, int idx);

double trial_func(double alpha, double sum_r_squared);                 //Test wave function
double local_energy_1D_analytical(double alpha, double **r, double r_sum, int idx);     //Local energy in one dimension (analytical)
double local_energy_2D_analytical(double alpha, double **r, double r_sum, int idx);     //Local energy in two dimensions (analytical)
double local_energy_3D_analytical(double alpha, double **r, double r_sum, int idx);     //Local energy in three dimensions (analytical)

void MonteCarlo();

double init_r_sum(double **r);                                  //Calculates the sum of the square of all posistions
double update_r_sum(double sum, double r_init, double r_move);  //Updates the sum of the square of all positions

double local_energy_brute_force(double alpha, double **r, double r_sum, int idx);



public:

//Constructor
Solver(int N, int num_alphas, int MC, int D, int type_energy);
void write_to_file(string outfilename);

};

#endif
