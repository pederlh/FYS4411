#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include "omp.h"

#include "BoltzPar.hpp"


using namespace std;
using namespace arma;


int main(int argc, char const *argv[]) {
    int num_part = 2;
    int dim = 2;
    int MC = pow(2,16);
    double eta = 0.001;

    //typesampling: 0 for brute force, 1 for metropolis

    int interaction = 1;
    double omega = 1.0; //2D
    //double omega = 1.0/4.0; //3D
    int num_threads = 2;

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
        BoltzPar solver(num_part, dim, eta, MC, 1, interaction, omega, ID);

    }


    return 0;
}
