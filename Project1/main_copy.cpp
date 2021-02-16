#include "Solver_copy.hpp"

int main(int argc, char const *argv[]) {

    //int particles[] = {1, 10, 100, 500};
    //int dimentionss[] = {1,2,3};
    //int type_energies[] = {0,1};

    string outfilename, calc;
    clock_t start, stop;
    int num_alphas = 20;
    int num_particles = 10;
    int mc_cycles = 10000;
    int dimentions = 3;
    // type_energy = 0 for analytical, type_energy = 1 for brute force.
    int type_energy = 0;

    Solver_copy mysolver(num_particles, num_alphas, mc_cycles, dimentions, type_energy);

    start = clock();
    mysolver.MonteCarlo();
    stop = clock();
    double timeused = (double) (stop-start)/(CLOCKS_PER_SEC);

    if (type_energy==0){
        calc = "ana";
    }
    if (type_energy==1){
        calc = "num";
    }


    outfilename = "COPY" + to_string(dimentions)+ "D_" +calc+"_N_"+to_string(num_particles)+"_MC_"+ to_string(mc_cycles)+".txt";

    mysolver.Write_to_file(outfilename,timeused);


    return 0;
}
