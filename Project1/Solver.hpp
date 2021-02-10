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

friend class WaveFunction; //Class "Wavefunction" can access private members from Solver

int N_, num_alphas_, MC_, D_;
double h_, sum_, tmp_, step_, D_diff_;
double *alphas_,*energies_, *variances_;

double tf_middle_, laplace_tf_;
double **r_old_, *r_new_, **quantum_force_old_, *quantum_force_new_;
double **numerical_matrix_;
random_device rd_;


//Pointer to member function
double (Solver::*energy_calculation)(double alpha, double r_sum);
void (Solver::*MC_method)();
double (Solver::*init_positions)(double r2_sum);
double(Solver::*metropolis_sampling)(double r2_new, double r2_old, double alpha, int move_idx, double move_P, double acc_P);

double Trial_func(double alpha, double sum_r_squared);                 //Test wave function


double Init_r_sum(double **r);                                  //Calculates the sum of the square of all posistions
double Update_r_sum(double sum, double r_init, double r_move);  //Updates the sum of the square of all positions

double Metropolis(double r2_new, double r2_old, double alpha, int move_idx, double move_P, double acc_P);
double Metropolis_importance(double r2_new, double r2_old, double alpha, int move_idx, double move_P, double acc_P);

double Local_energy_brute_force(double alpha, double r_sum);
double Local_energy_analytical(double alpha, double r_sum);     //Local energy analytical

double Initialize_positions(double r2_sum);

void Initialize_quantum_force(double alpha, double **positions, double **q_force);
void Update_quantum_force(double alpha);
double Greens_function(int idx);

public:

//Constructor
Solver(int N, int num_alphas, int MC, int D, int type_energy, int type_sampling);
void Write_to_file(string outfilename, double time);
void MonteCarlo();
void MonteCarlo_importance();
double MonteCarlo_SGD(double alpha);
void Gradient_descent();

};

#endif
