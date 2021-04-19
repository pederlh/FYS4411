#include "Solver.hpp"

int main(int argc, char const *argv[]) {
    int num_particles, dimentions, mc_cycles, type_energy, type_sampling, num_threads;
    double mc_cycles_optimal_run, learning_rate;


    num_particles = 10;
    dimentions = 3;
    mc_cycles = 10000;
    type_energy = 0;
    type_sampling = 2;
    num_threads = 2;
    mc_cycles_optimal_run = 10000;
    learning_rate = 0.01;

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

    }
    return 0;
}
