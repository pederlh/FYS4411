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

    int num_part =  atoi(argv[1]);
    int dim =  atoi(argv[2]);
    int MC =  atoi(argv[3]);
    int num_hidden =  atoi(argv[4]);
    int num_threads =  atoi(argv[5]);
    double eta =  atof(argv[6]);


    //typesampling: 0 for brute force, 1 for metropolis hastings
    int typesampling =  atoi(argv[7]);
    //interaction: 0 for none, 1 for yes
    int interaction =  atoi(argv[8]);
    //opt: 0 for GD, 1 for ADAM
    int opt =  atoi(argv[9]);
    string optimizer;

    if (opt== 0){
        optimizer = "GD";
    }
    if (opt== 1){
        optimizer = "ADAM";
    }
    double omega = atoi(argv[10]);


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
        //list <double> omegas = {0.25, 0.08333333, 0.05, 0.0222417, 0.01826864, 0.01163857, 0.0086731, 0.00678601, 0.00595797, 0.00531003, 0.0047892};
        BoltzmannMachine solver(num_part, dim, eta, MC, typesampling, interaction, omega, num_hidden, optimizer,ID);


    }
    return 0;
}
