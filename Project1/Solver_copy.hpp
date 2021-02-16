#ifndef SOLVER_COPY_HPP
#define SOLVER_COPY_HPP

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


class Solver_copy {
private:

friend class WaveFunction; //Class "Wavefunction" can access private members from Solver

int N_, num_alphas_, MC_, D_;
double h_, tmp_, step_;
double *alphas_,*energies_, *variances_;

double tf_middle_, laplace_tf_;
double **r_old_, *r_new_;
random_device rd_;


//Pointer to member function
double (Solver_copy::*energy_calculation)(double alpha, double r_sum);

double Trial_func(double alpha, double sum_r_squared);                 //Test wave function
double Update_r_sum(double sum, double r_init, double r_move);  //Updates the sum of the square of all positions

double Metropolis(double r2_new, double r2_old, double alpha);

double Local_energy_brute_force(double alpha, double r_sum);
double Local_energy_analytical(double alpha, double r_sum);     //Local energy analytical


public:

//Constructor
Solver_copy(int N, int num_alphas, int MC, int D, int type_energy);
void Write_to_file(string outfilename, double time);
void MonteCarlo();

};

#endif
