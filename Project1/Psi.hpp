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


//Pointer to member function
double (Psi::*energy_calculation)(double alpha, double r_sum);
void (Psi::*QF)(double alpha, double **positions, double **q_force);


void No_quantum_force(double alpha, double **positions, double **q_force);

*/
public:

    double h_, step_, D_diff_, beta_, a_;
    int N_, D_, case_;
    random_device rd_;

    double **r_old_, *r_new_;
    double *quantum_force_old_, *quantum_force_new_, *rkl_;
    double r2_sum_old_, r2_sum_new_;
    double tf_middle_, laplace_tf_;



    void Declare_position(int N, int D, double h, double step, int case_type);
    void Declare_position_interaction(int N, int D, double h, double step, int case_type);

    void Declare_quantum_force(double D_diff);
    double Initialize_positions();
    void Initialize_quantum_force(double alpha, int idx);
    void Initialize_quantum_force_interaction(double alpha, int idx);

    double Proposed_move(int idx);
    double Proposed_move_importance(int idx);
    double Proposed_move_interaction(int idx);

    double Local_energy_analytical(double alpha);
    double Local_energy_brute_force(double alpha);
    double Local_energy_interaction_analytical(double alpha);
    double Local_energy_interaction_brute_force(double alpha);


    double Update_r_sum(double sum, double r_init, double r_move);
    double Update_r_sum_interaction(double sum, double r_init, double r_move, double coord);


    double Trial_func(double alpha, double sum_r_squared);
    double Trial_func_interaction(double alpha, double sum_r_squared, string version, int idx);
    void Update_quantum_force(double alpha);
    void Update_quantum_force_interaction(double alpha, int idx);


    double Laplace_phi(int idx, double d2_phi, double alpha);

};

#endif
