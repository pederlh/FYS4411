#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iomanip>
#include <cstdlib>
#include <armadillo>

#include "BoltzmannMachine.hpp"

using namespace std;
using namespace arma;

int main(int argc, char const *argv[]) {
    int num_part = 2;
    int dim = 2;
    int MC = 100000;
    double eta = 0.001;

    //typesampling: 0 for brute force, 1 for metropolis

    bool interaction = false;

    BoltzmannMachine solver(num_part, dim, eta, MC, 1, interaction);

    return 0;
}
