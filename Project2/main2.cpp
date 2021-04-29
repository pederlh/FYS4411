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
    int num_part = 6;
    int dim = 3;
    int MC = 10000;
    double eta = 0.001;

    //typesampling: 0 for brute force, 1 for metropolis

    BoltzmannMachine solver(num_part, dim, eta, MC, 1);

    return 0;
}
