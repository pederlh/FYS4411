#include "Solver.hpp"

int main(int argc, char const *argv[]) {
    /*
    int num_alphas, num_particles, dimensions, mc_cycles, mc_cycles_optimal_run, type_energy, type_sampling, num_treads, OBD_check;
    if( argc != 10 ){
        cout << "-------------------------------------------------------" << endl;
        cout << "Bad Usage: " << argv[0] << " takes in nine arguments:" << endl;
        cout << "num_alphas, num_particles, dimensions, mc_cycles, mc_cycles_optimal_run, type_energy, type_sampling, num_treads, OBD_check" << endl;            
        cout << "-------------------------------------------------------" << endl;
        exit(1);
    }
    else{
        num_alphas              = atoi(argv[1]);    // Number of alpha values to be looped over (ignored when gradient descent used)
        num_particles           = atoi(argv[2]);    // Number of particles in system
        dimensions              = atoi(argv[3]);    // Number of spacial dimensions in system
        mc_cycles               = atoi(argv[4]);    // MC cycles used per run during gradient descent (GD) and in non-GD cases
        mc_cycles_optimal_run   = atoi(argv[5]);    // MC cycles used in big run after gradient descent
        type_energy             = atoi(argv[6]);    // type_energy = 0 for analytical, = 1 for brute force (Laplace operator in Hamiltonian).
        type_sampling           = atoi(argv[7]);    // Chooses how MC sampling and search for optimal alpha is done. See instructions below
        num_treads              = atoi(argv[8]);    // Number of threads (parallelization)
        OBD_check               = atoi(argv[9]);    // OBD_check = 0 skips calculation of one-body densites, = 1 does the calculation 

        // Choices for type_sampling parameter:
        // type_sampling = 0 for no importance sampling (loop through alpha-values)
        // type_sampling = 1 for importance sampling (loop through alpha-values)
        // type_sampling = 2 for importance sampling + gradient descent (non-interacting case)
        // type_sampling = 3 for importance sampling + gradient descent (interacting case)        
    }
    */     
    int num_alphas = 15;
    int num_particles = 10;
    int mc_cycles = 1000;
    int mc_cycles_optimal_run = pow(2,17);
    int dimentions = 3;
    int type_energy = 0;
    int type_sampling = 3;
    int OBD_check = true;

    // Certain combinations of input parameters are not implemented
    if(type_energy == 1 && type_sampling == 3){
        cout << "Brute force energy calculation not implemented for interactive case" << endl;
        return 1;
    }
    if(type_sampling <= 1 && OBD_check == true){
        cout << "Note: One body density calculation not implemented for non-gradient method." << endl;
    }
    // Number of threads
    int num_threads = 6;
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

        // For the cases without gradient descent, we time the calculations from initializing objects to before writing results to file.
        start_time = omp_get_wtime();                           // Start recording time

        // Initialize Solver object and perform calculations
        Solver mysolver(num_particles, num_alphas, mc_cycles, mc_cycles_optimal_run, dimentions, type_energy, type_sampling, ID);

        end_time = omp_get_wtime();                             // Stop recording time
        double timeused = end_time - start_time;


        // When type_sampling is 0 or 1, data is written to file with the functions below

        string calc, outfilename;
        if (type_energy==0) calc = "ana";   // Different strings appended to file name
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
    return 0;
}
