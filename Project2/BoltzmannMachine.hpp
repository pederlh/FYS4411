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

public:
    cube w_;
    mat a_;
    vec b_;
    int D_, N_, hidden_nodes_;
    random_device rd_;

void Initialize_Parameters(int dimentions, int num_particles);
};

#endif
