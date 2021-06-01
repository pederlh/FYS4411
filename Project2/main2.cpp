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
//#include "test.hpp"

using namespace std;
using namespace arma;


int main(int argc, char const *argv[]) {
    int num_part = 2;
    int dim = 2;
    int MC = pow(2,16);
    int num_hidden = 2;
    int num_threads = 1;
    double eta = 0.1;


    //typesampling: 0 for brute force, 1 for metropolis

    string optimizer = "GD";
    if (optimizer == "GD"){eta = 0.1;}
    if (optimizer == "ADAM"){eta = 0.1;}
    int interaction = 1;
    double omega = 1.0; //2D
    //double omega = 1.0/4.0; //3D

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
        list<int> layers = {1,2,3,4,5,6,7,8,9,10};
        for (int h : layers){
            BoltzmannMachine solver(num_part, dim, eta, MC, 1, interaction, omega, h, optimizer,ID);
        }

        //BoltzmannMachine solver(num_part, dim, eta, MC, 1, interaction, omega, num_hidden, optimizer,ID);

    }



    /*
    list<int> layers = {1,2,3,4,5,6,7,8,9, 10};
    for (int h : layers){
        BoltzmannMachine solver(num_part, dim, eta, MC, 1, interaction, omega, h);
    }

    list<double> omegas = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0};
    for (double w : omegas){
        BoltzmannMachine solver(num_part, dim, eta, MC, 1, interaction, w, num_hidden);
    }
    */


    return 0;
}

/*
int main(int argc, char const *argv[]) {
    int num_part = 2;
    int dim = 2;
    int MC = 30;
    double eta = 0.001;

    //typesampling: 0 for brute force, 1 for metropolis

    bool interaction = true;

    test solver(num_part, dim, eta, MC, 1, interaction);

    return 0;
}
*/
