#ifndef PSI_HPP
#define PSI_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iomanip>
#include <cstdlib>

using namespace std;

/*
        Class for trial wave function and hamiltonian.
        Contains functions for:
            - Evaluation of local energy
            - Trial wave function
            - Particle distribution
            - Quantum force

*/

class Psi {

public:

    double h_, step_, D_diff_, beta_, a_, move_step_;
    int N_, D_, case_;
    random_device rd_;

    double **r_old_, **r_new_;
    double *quantum_force_old_, *quantum_force_new_;
    double r2_sum_old_, r2_sum_new_;
    double tf_middle_, laplace_tf_;

    //Non-interacting case
    double Initialize_positions();
    void Declare_position(int N, int D, double h, double step, int case_type);
    double Proposed_move(int idx);
    double Local_energy_analytical(double alpha);
    double Local_energy_brute_force(double alpha);
    double Update_r_sum(double sum, double r_init, double r_move);
    double Trial_func(double alpha, double sum_r_squared);

    //Importance sampling
    void Declare_quantum_force(double D_diff);
    void Initialize_quantum_force(double alpha, int idx);
    double Proposed_move_importance(int idx);
    void Update_quantum_force(double alpha);

};

#endif
