#include "Solver.hpp"

int main(int argc, char const *argv[]) {
    // /*
    int num_particles, dimentions, mc_cycles, type_energy, type_sampling, num_threads, OBD_check;
    double mc_cycles_optimal_run, learning_rate;
    if( argc < 8 || argc > 10 ){
        cout << "------------------------------------------------------" << endl;
        cout << "Bad Usage: " << argv[0] << " takes in eight or nine arguments:" << endl;
        cout << "num_particles, dimensions, mc_cycles, type_energy, type_sampling, num_treads, OBD_check, mc_cycles_optimal_run (when using gradient descent), learning rate (GD)" << endl;
        cout << "-------------------------------------------------------" << endl;

        cout << "got " << argc << endl;
        exit(1);
    }
    else{
        num_particles           = atoi(argv[1]);    // Number of particles in system
        dimentions              = atoi(argv[2]);    // Number of spacial dimensions in system
        mc_cycles               = atoi(argv[3]);    // MC cycles used per run during gradient descent (GD) and in non-GD cases
        type_energy             = atoi(argv[4]);    // type_energy = 0 for analytical, = 1 for brute force (Laplace operator in Hamiltonian).
        type_sampling           = atoi(argv[5]);    // Chooses how MC sampling and search for optimal alpha is done. See instructions below
        num_threads             = atoi(argv[6]);    // Number of threads (parallelization)
        OBD_check               = atoi(argv[7]);    // OBD_check = 0 skips calculation of one-body densites, = 1 does the calculation
        if (argc == 10){
        mc_cycles_optimal_run   = atoi(argv[8]);    // MC cycles used in big run after gradient descent
        learning_rate           = atof(argv[9]);    // Learning rate gradient descent
        }
        else{
        mc_cycles_optimal_run   = int(1e17);
        learning_rate = 0.01;
        }

        // Choices for type_sampling parameter:
        // type_sampling = 0 for no importance sampling (loop through alpha-values)
        // type_sampling = 1 for importance sampling (loop through alpha-values)
        // type_sampling = 2 for importance sampling + gradient descent (non-interacting case)
        // type_sampling = 3 for importance sampling + gradient descent (interacting case)
    }
    // */

    // Certain combinations of input parameters are not implemented
    if(type_energy == 1 && type_sampling == 3){
        cout << "Brute force energy calculation not implemented for interactive case" << endl;
        return 1;
    }
    if(type_sampling <= 1 && OBD_check == true){
        cout << "Note: One body density calculation not implemented for non-gradient method." << endl;
    }

    // Start parallelization
    if(num_threads > omp_get_max_threads()){
        cout << "Warning: requested number of threads (" << num_threads << ") is greater than omp_get_max_threads (" << omp_get_max_threads() << ")" << endl;
        cout << "Changing number of threads to omp_get_max_threads..." << endl;
        num_threads = omp_get_max_threads();
    }
    omp_set_num_threads(num_threads);

    double *shared_alphas = new double[num_threads];        // Array used in GD that all threads have access to.
                                                            // Used to find average alpha in GD.

    #pragma omp parallel
    {
        int ID = omp_get_thread_num();

        // Initialize Solver object and perform calculations
        Solver mysolver(num_particles, mc_cycles, mc_cycles_optimal_run, dimentions, type_energy, type_sampling, ID, learning_rate, shared_alphas);


        // When type_sampling is 0 or 1, data is written to file with the functions below

        string calc, outfilename;
        if (type_energy==0) calc = "ana";   // Different strings appended to file name
        if (type_energy==1) calc = "num";

        if (type_sampling == 0){
            outfilename =   "spherical_HO_" + to_string(dimentions) + "D_" + calc +
                            "_N_"+ to_string(num_particles) + "_MC_"+ to_string(mc_cycles) +
                            "_stringID_" + to_string(ID) + ".txt";

            mysolver.Write_to_file(outfilename);
        }

        if (type_sampling == 1){
            outfilename =   "importance_spherical_HO_" + to_string(dimentions) + "D_" + calc +
                            "_N_" + to_string(num_particles) + "_MC_" + to_string(mc_cycles) +
                            "_stringID_" + to_string(ID) + ".txt";

            mysolver.Write_to_file(outfilename);
        }
    }
    return 0;
}
