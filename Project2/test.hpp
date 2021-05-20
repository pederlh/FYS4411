#ifndef TEST_HPP
#define TEST_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iomanip>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;


class test {

    /*
    Class for NQS wave function.
    */

private:
    void (test::*optimizer)();
    void (test::*MetropolisMethod)();


public:
    cube w_, dw_, E_dw_;
    mat a_, da_, E_da_;
    vec b_, db_, E_db_;
    vec Q_;
    int D_, N_, H_;
    double sigma_, sigma2_;
    bool interaction_;
    vec new_p,acc;
    ivec index;
    int lil_c, lil_n;

    mat r_old_, r_new_, quantum_force_, quantum_force_old_,quantum_force_new_;

    test(int num_particles,int dimentions, double eta, int MC, int type_sampling, bool interaction);
    double WaveFunction(mat r);
    void Q_factor(mat r);

    void Metropolis();
    void Metropolis_Hastings();
    double tf_old_, tf_new_, P_, D_diff_, t_step_;   //Parameters for Metropolis algorithm

    double MonteCarlo();
    double LocalEnergy();
    void Derivate_wavefunction();
    mat QuantumForce(mat r);
    double GreensFunction(int idx);

    void SGD();
    //For ADAM
    void ADAM();
    double eta_;
    cube mom_w_, second_mom_w_;
    vec mom_b_, second_mom_b_;
    mat mom_a_, second_mom_a_;
    int MC_;



};

#endif