#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iomanip>
#include <cstdlib>
#include <armadillo>

#include "BoltzmannMachine.hpp"
//#include "test.hpp"

using namespace std;
using namespace arma;


int main(int argc, char const *argv[]) {
    int num_part = 2;
    int dim = 2;
    int MC = 1000;
    double eta = 0.001;

    //typesampling: 0 for brute force, 1 for metropolis

    int interaction = 1;

    BoltzmannMachine solver(num_part, dim, eta, MC, 1, interaction);

    return 0;
}


/*
int main(int argc, char const *argv[]) {
    int num_part = 2;
    int dim = 2;
    int MC = 15;
    double eta = 0.001;

    //typesampling: 0 for brute force, 1 for metropolis

    bool interaction = true;

    test solver(num_part, dim, eta, MC, 1, interaction);

    return 0;
}
*/
