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
    int num_part = 10;
    int dim = 3;
    int MC = 10000;
    double eta = 0.01;

    BoltzmannMachine solver(num_part, dim, eta, MC);

    return 0;
}
