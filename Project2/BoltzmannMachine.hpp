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

public:
    cube w_, dw_;
    mat a_, da_;
    vec b_, db_;
    int D_, N_, hidden_nodes_;
    random_device rd_;

    void Initialize_Parameters(int dimentions, int num_particles, double eta);



};

#endif
