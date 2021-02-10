#include "Solver.hpp"

int main(int argc, char const *argv[]) {

    //int particles[] = {1, 10, 100, 500};
    //int dimentionss[] = {1,2,3};
    //int type_energies[] = {0,1};

    string outfilename, calc;
    clock_t start, stop;
    int num_alphas = 10;
    int num_particles = 10;
    int mc_cycles = 100000;
    int dimentions = 3;
    // type_energy = 0 for analytical, type_energy = 1 for brute force.
    int type_energy = 0;

    //type_sampling = 0 for no importance sampling, type_sampling = 1 for importance sampling, =2 for SGD.
    int type_sampling = 2;


    Solver mysolver(num_particles, num_alphas, mc_cycles, dimentions, type_energy, type_sampling);
    /*
    start = clock();
    mysolver.MonteCarlo();
    stop = clock();
    double timeused = (double) (stop-start)/(CLOCKS_PER_SEC);
    */
    if (type_energy==0){
        calc = "ana";
    }
    if (type_energy==1){
        calc = "num";
    }



    double timeused = 10.10;

    if (type_sampling == 0){
        outfilename = "spherical_HO_" + to_string(dimentions)+ "D_" +calc+"_N_"+to_string(num_particles)+"_MC_"+ to_string(mc_cycles)+".txt";

        mysolver.Write_to_file(outfilename,timeused);
    }

    if (type_sampling == 1){
        outfilename = "importance_spherical_HO_" + to_string(dimentions)+ "D_" +calc+"_N_"+to_string(num_particles)+"_MC_"+ to_string(mc_cycles)+".txt";

        mysolver.Write_to_file(outfilename,timeused);
    }

    if (type_sampling == 2){
        outfilename = "SGD_spherical_HO_" + to_string(dimentions)+ "D_" +calc+"_N_"+to_string(num_particles)+"_MC_"+ to_string(mc_cycles)+".txt";

        //mysolver.Write_to_file_SGD(outfilename,timeused);
    }






    /*
    for (int p=0;p<4;p++){
        num_particles = particles[p];
        for(int d = 0; d<3;d++){
            dimentions = dimentionss[d];
            for (int e =0; e<2;e++){
                type_energy = type_energies[e];

                Solver mysolver(num_particles, num_alphas, mc_cycles, dimentions, type_energy);
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


                outfilename = "spherical_HO_" + to_string(dimentions)+ "D_" +calc+"_N_"+to_string(num_particles)+"_MC_"+ to_string(mc_cycles)+".txt";
                mysolver.write_to_file(outfilename,timeused);

            }
        }
    }
    */


    return 0;
}
