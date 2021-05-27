#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include <list>

#include "BoltzmannMachine.hpp"
//#include "test.hpp"

using namespace std;
using namespace arma;


int main(int argc, char const *argv[]) {
    int num_part = 2;
    int dim = 2;
    int MC = pow(2,18);
    double eta = 0.1;
    int num_hidden = 2;


    //typesampling: 0 for brute force, 1 for metropolis

    string optimizer = "GD";
    int interaction = 0;
    double omega = 1.0; //2D
    //double omega = 1.0/4.0; //3D
    BoltzmannMachine solver(num_part, dim, eta, MC, 1, interaction, omega, num_hidden, optimizer);

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
