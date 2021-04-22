#ifndef BOLTZMANNMACHINE_HPP
#define BOLTZMANNMACHINE_HPP

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


class BoltzmannMachine {

    /*
    Class for NQS wave function.
    */

private:
    void (BoltzmannMachine::*optimizer)();

    //For ADAM;
    void ADAM();
    double beta1_, beta2_, epsilon_, alpha_batch_, epsilon_batch_;
    double sigma_, eta_;
    cube mom_w_, second_mom_w_;
    vec mom_b_, second_mom_b_;
    mat mom_a_, second_mom_a_;
    int MC_;

public:
    cube w_, dw_, Energy_dw_;
    mat a_, da_, Energy_da_;
    vec b_, db_, Energy_db_;
    vec Q_;
    int D_, N_, H_;
    random_device rd_;

    mat r_old_, r_new_;

    BoltzmannMachine(int num_particles,int dimentions, double eta, int MC);
    void Initialize_SGD();
    void Initialize();
    double WaveFunction(mat r);
    void Q_factor(mat r);
    void Metropolis();
    double MonteCarlo();
    void SGD();
    double LocalEnergy();
    void Derivate_wavefunction();



};

#endif
