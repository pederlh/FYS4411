#ifndef PSI_HPP
#define PSI_HPP

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


class Psi {
private:
/*
double h_, sum_, tmp_, step_, D_diff_;

double tf_middle_, laplace_tf_;
double **r_old_, *r_new_, **quantum_force_old_, *quantum_force_new_;
random_device rd_;


//Pointer to member function
double (Psi::*energy_calculation)(double alpha, double r_sum);
void (Psi::*QF)(double alpha, double **positions, double **q_force);

double Trial_func(double alpha, double sum_r_squared);                 //Test wave function

double Update_r_sum(double sum, double r_init, double r_move);  //Updates the sum of the square of all positions

double Local_energy_brute_force(double alpha, double r_sum);
double Local_energy_analytical(double alpha, double r_sum);     //Local energy analytical

void Initialize_quantum_force(double alpha, double **positions, double **q_force);
void No_quantum_force(double alpha, double **positions, double **q_force);
void Update_quantum_force(double alpha);
double Greens_function(int idx);
*/
public:

    double h_, step_, D_diff_;
    int N_, D_;
    random_device rd_;

    double **r_old_, *r_new_;
    double **quantum_force_old_, *quantum_force_new_;
    double r2_sum_old_, r2_sum_new_;
    double tf_middle_, laplace_tf_;


//Constructor
    Psi();

    void Declare_positions(int N, int D);
    void Declare_quantum_force();
    void Initialize_positions();
    void Initialize_quantum_force(double alpha);

    double Greens_function(int idx);

    double Local_energy_analytical(double alpha);
    double Local_energy_brute_force(double alpha);
    double Update_r_sum(double sum, double r_init, double r_move);
    double Trial_func(double alpha, double sum_r_squared);
    void Update_quantum_force(double alpha);

};

#endif
