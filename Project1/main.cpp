#include "Solver.hpp"

int main(int argc, char const *argv[]) {

    //int particles[] = {1, 10, 100, 500};
    //int dimentionss[] = {1,2,3};
    //int type_energies[] = {0,1};

    string outfilename, calc;
    clock_t start, stop;
    int num_alphas = 15;
    int num_particles = 10;
    int mc_cycles = 1000;
    int dimentions = 3;

    // type_energy = 0 for analytical, type_energy = 1 for brute force.
    int type_energy = 0;

    //type_sampling = 0 for no importance sampling, type_sampling = 1 for importance sampling, = 2 for Gradient descent + importance sampling.
    int type_sampling = 2;

    // Number of threads
    int num_threads = 2;

    double start_time, end_time;

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
        int cycles_per_thread = mc_cycles/num_threads;         // Remainder will be added to master thread later
        #pragma omp master
        {
            if(num_threads > omp_get_num_threads()) cout << "Warning: Number of threads actually set to be" << omp_get_num_threads() << endl;
            cycles_per_thread += mc_cycles % num_threads;      // Add remaining cycles
        }
        

        // ------------------------------------------------
        // OBS OBS Når skal vi starte og slutte å måle tid?    
        start_time = omp_get_wtime();                          // Start recording time
        // ------------------------------------------------


        // Initialize Solver object and perform calculations
        Solver mysolver(num_particles, num_alphas, cycles_per_thread, dimentions, type_energy, type_sampling, ID);
    

        end_time = omp_get_wtime();
        double timeused = end_time - start_time;

        if (type_energy==0) calc = "ana";
        if (type_energy==1) calc = "num";

        if (type_sampling == 0){
            outfilename =   "spherical_HO_" + to_string(dimentions) + "D_" + calc + 
                            "_N_"+ to_string(num_particles) + "_MC_"+ to_string(mc_cycles) +
                            "_stringID_" + to_string(ID) + ".txt";

            mysolver.Write_to_file(outfilename,timeused);
        }

        if (type_sampling == 1){
            outfilename =   "importance_spherical_HO_" + to_string(dimentions) + "D_" + calc + 
                            "_N_" + to_string(num_particles) + "_MC_" + to_string(mc_cycles) + 
                            "_stringID_" + to_string(ID) + ".txt";

            mysolver.Write_to_file(outfilename,timeused);
        }
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