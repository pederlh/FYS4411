#include "Solver.hpp"

int main(int argc, char const *argv[]) {

    int num_alphas = 10;
    int num_particles = 10;
    int mc_cycles = 1000000;
    int dimentions = 3;
    // type_energy = 0 for analytical, type_energy = 1 for brute force.
    int type_energy = 1;
    Solver mysolver(num_particles, num_alphas, mc_cycles, dimentions, type_energy);
    mysolver.write_to_file("spherical_HO_1D.txt");

    return 0;
}
