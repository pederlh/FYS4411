#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iomanip>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;



int main(int argc, char const *argv[]) {
    int N = 3;
    int D = 3;
    int hidden_nodes = 2;
    double std_norm_dist = 0.01;


    cube w = cube(N, D, hidden_nodes, fill::randn)*std_norm_dist;
    vec b = vec(hidden_nodes, fill::randn)*std_norm_dist;
    mat r = mat(N,D).fill(0.0);
    vec temp = vec(hidden_nodes).fill(0.0);
    vec Q = vec(hidden_nodes).fill(0.0);

    for (int h = 0; h < hidden_nodes; h++){
        mat W = w.slice(h);
        temp(h) = accu(r%W);
    }

    Q = b + temp;
    Q.print();




    return 0;
}
