#ifndef BOLTZMANNMACHINE_HPP
#define BOLTZMANNMACHINE_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iomanip>
#include <cstdlib>

using namespace std;


class BoltzmannMachine {

public:
    double *a_, *b_, **w_;
    int D_, N_, hidden_layers_;
    random_device rd_;

void Initialize_Parameters();
};

#endif
