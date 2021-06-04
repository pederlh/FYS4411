#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include <list>
#include "omp.h"

#include "BoltzmannMachine.hpp"

using namespace std;
using namespace arma;


int main(int argc, char const *argv[]) {
    int num_part = 2;
    int dim = 2;
    int MC = pow(2,18);
    int num_hidden = 2;
    int num_threads = 1;
    double eta = 0.01;


    //typesampling: 0 for brute force, 1 for metropolis
    int typesampling = 1;
    //interaction: 0 for none, 1 for yes
    int interaction = 0;
    string optimizer = "ADAM";
    double omega = 1.0; //2D

    // Start parallelization
    if(num_threads > omp_get_max_threads()){
        cout << "Warning: requested number of threads (" << num_threads << ") is greater than omp_get_max_threads (" << omp_get_max_threads() << ")" << endl;
        cout << "Changing number of threads to omp_get_max_threads..." << endl;
        num_threads = omp_get_max_threads();
    }
    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
        int ID = omp_get_thread_num();

        // Initialize Solver object and perform calculations
        BoltzmannMachine solver(num_part, dim, eta, MC, typesampling, interaction, omega, num_hidden, optimizer,ID);

    }
    return 0;
}
