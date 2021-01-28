#include "Solver.hpp"

int main(int argc, char const *argv[]) {

    int num_alphas = 10;
    int num_particles = 10;
    int mc_cycles = 1000;
    int dimentions = 1;
    Solver mysolver(num_particles, num_alphas, mc_cycles, dimentions);
    mysolver.MonteCarlo();
    mysolver.write_to_file("spherical_HO_1D.txt");

    return 0;
}
