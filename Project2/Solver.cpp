#include "Solver.hpp"
#include "Psi.hpp"

Solver::Solver(int N, int MC, int MC_optimal_run, int D, int type_energy, int type_sampling, int thread_ID, double learning_rate, double*shared_alphas){
    D_ = D;                         // Dimentions
    N_ = N;                         // Number of particles
    h_ = 1.0;                       // Stepsize to determine the distributions space of particles
    step_ = h_*pow(10,-4);          // Stepsize used in numerical differentiation
    thread_ID_ = thread_ID;         // ID of thread (parallelization)
    equi_cycles_ = (int) 0.10*MC;   // Number of burn in cycles (used for equilibration)

    type_energy_ = type_energy;     // Analytical expression or numerical differentiation
    type_sampling_ = type_sampling; // Brute force = 0, importance sampling = 1, gradient descent = 2 or gradient descent with interaction = 3


    wave.Declare_position(N_, D_,h_, step_, 0);             // Non-interacting case


    MC_ = MC;                       // Number of monte carlo cycles
    OBD_check_ = true;              // "true" gives calculation of one body density


    if (type_sampling_ == 2){
        MC_optimal_run_ = MC_optimal_run;               // Number of MC cycles for run with optimal alpha value after GD
        D_diff_ = 0.5;
        wave.Declare_quantum_force(D_diff_);
        metropolis_sampling = &Solver::Metropolis_importance;

        main_method = &Solver::Gradient_descent;
        Interaction_or_not_GD = &Solver::MonteCarlo_GD_noninteracting;
        Interaction_or_not_optimal = &Solver::MonteCarlo_optval_noninteracting;
        tol_GD_ = 1e-9;//5*pow(10,-8);                                 // Acceptance tolerance for gradient descent
        eta_GD_ = learning_rate;                        // Learning rate for gradient descent
    }


    (this->*main_method)(shared_alphas);                               // Calls MC method and starts simulation
}



//Gradient Descent w/importance sampling for interacting and non-interacting bosons
void Solver::Gradient_descent(double * shared_alphas){
    string file;
    int iterations = 600;                                // Max number iterations of GD
    double *alpha_vals_GD = new double[iterations];     // Array to store alpha values from GD
    for (int z = 0; z < iterations; z++){
        alpha_vals_GD[z] = 0;
    }
    double *values = new double[3];                     // Array to contain energy, variance and energy derivative wrt alpha;
    double alpha_guess = 0.7;                          // Initial guess for alpha
    double alpha_guess_std;                             // Standard deviation of mean of alpha found by all threads
    int counter = 0;                                    // Counter to keep track of actual number of iterations
    double tol = tol_GD_;


     #pragma omp master
     {
         if (omp_get_num_threads() == 1) cout << "Start gradient descent" << endl;
         else cout << "Start gradient descent (showing progress master thread) N = "<< N_ << endl << endl;

         cout << setw(10) << "Alpha" << setw(12) << "Energy" << setw(16) << "Variance" << endl;
         cout << "--------------------------------------" << endl;
     }
    //Gradient descent
    for (int i = 0; i < iterations; i++){
        counter++;
        alpha_vals_GD[i] = alpha_guess;

        (this->*Interaction_or_not_GD)(values, alpha_guess);   //Points to correct Monte Carlo simulation (interacting/non-interacting)
        # pragma omp master
        {
            cout << setw(10) << setprecision(8) << alpha_guess << setw(12) << values[0] << setw(16) << values[1] << " ID: " << thread_ID_ << endl;
        }

        //Breaks GD if alpha provides acceptably low sample variance
        if (values[1] < tol && values[1]>0){
            break;
        }
        alpha_guess -= eta_GD_*values[2];
    }



    cout <<" Convergence after " << counter << " number of iterations  (thread " << to_string(thread_ID_) << ")" << endl;

    shared_alphas[thread_ID_] = alpha_guess;         // Write best alpha guess to shared array
    # pragma omp barrier
    # pragma omp master
    {
        cout << "Gradient descent finished, starting main MC calculations..." << endl;
        cout << "Best alphas found: " << endl;
        for (int i=0; i < omp_get_num_threads(); i++){
            cout << "Thread " << i << ": alpha = " << shared_alphas[i] << endl;
        }
    }
    # pragma omp barrier
    # pragma omp critical
    {
        double sum = 0.0;
        double sum_std = 0.0;
        for (int i =0; i < omp_get_num_threads();i++){
            sum += shared_alphas[i];
        }

        alpha_guess = sum/omp_get_num_threads();

        for (int i =0; i < omp_get_num_threads();i++){
            sum_std += pow(shared_alphas[i]-alpha_guess,2);
        }

        alpha_guess_std = sqrt((1./omp_get_num_threads())*sum_std);
    }

    # pragma omp master
    {
        cout <<" "<<endl;
        cout << "Mean alpha for N = " << N_ << " +- standard deviation:" <<endl;
        cout << alpha_guess << " +- " << alpha_guess_std << endl<<endl;
    }


    // Run large MC simulation for optimal value alpha obtained from GD
    double *optimal_energies = new double[N_*MC_optimal_run_];

    equi_cycles_ = (int) 0.10*MC_optimal_run_;                              // Adjust number of equilibration cycles

    (this->*Interaction_or_not_optimal)(alpha_guess, optimal_energies);     // Start big set of calculations

    # pragma omp barrier
    # pragma omp master
    {
        cout << "Finished calculations!" << endl;
    }

    if(type_sampling_==2){
        file = "OPTIMAL_ALPHA"+ to_string(N_) + "_N_stringID_" + to_string(thread_ID_) +
            "_alpha_" + to_string(alpha_guess) + "_E_L_samples.txt";
    }

    ofstream ofile;
    ofile.open(file);
    ofile << setprecision(15) << time_ << endl;
    for (int i = 0; i < N_*MC_optimal_run_; i++){
        ofile << setprecision(15) << optimal_energies[i] << endl;
    }
    ofile.close();
}

//Method for running Monte Carlo simulation with gradient descent (non-interacting bosons)
//  - Calculates the derivative of local energy w/regards to alpha to use in GD
void Solver::MonteCarlo_GD_noninteracting(double *values, double alpha){
    double energy, energy_squared, DeltaE, variance, DerivateE, Derivate_WF_E, sum_r, Derivate_WF;

    wave.Initialize_positions();         //Initialize random initial matrix of positions

    //Equilibration step: runs metropolis algorithm without sampling to equilibrate system
    Equilibrate(alpha);

    energy = 0;
    energy_squared = 0;
    Derivate_WF = 0;
    Derivate_WF_E = 0;

    //Monte Carlo simulation with metropolis sampling
    for (int cycle = 0; cycle < MC_; cycle++){
        for (int n = 0; n < N_; n++){

            (this->*metropolis_sampling)(alpha); //Metropolis test

            //Calculate local energy
            if (type_energy_ == 0){
                DeltaE = wave.Local_energy_analytical(alpha);
            }
            if (type_energy_==1){
                DeltaE = wave.Local_energy_brute_force(alpha);
            }

            energy += DeltaE;
            energy_squared += DeltaE*DeltaE;
            sum_r = -wave.r2_sum_old_;
            Derivate_WF += sum_r;
            Derivate_WF_E += sum_r*DeltaE;
            }
        }

    energy /= (MC_*N_);
    energy_squared /= (MC_*N_);
    Derivate_WF/=(MC_*N_);
    Derivate_WF_E/=(MC_*N_);
    variance = energy_squared - energy*energy;
    DerivateE = 2*(Derivate_WF_E - Derivate_WF*energy);  //The derivative of the local energy w/regards to alpha

    values[0] = energy;
    values[1] = variance;
    values[2] = DerivateE;

}

//Method for running Monte Carlo simulation for one (optimized) value of variational parameter alpha (non-interacting bosons)
void Solver::MonteCarlo_optval_noninteracting(double alpha, double *energies){
    double DeltaE;

    wave.Initialize_positions();     // Initialize random initial matrix of positions

    //Equilibration step: runs metropolis algorithm without sampling to equilibrate system
    Equilibrate(alpha);
    start_time_ = omp_get_wtime();

    //Monte Carlo simulation with metropolis sampling
    for (int cycle = 0; cycle < MC_optimal_run_; cycle++){
        for (int n = 0; n < N_; n++){

            (this->*metropolis_sampling)(alpha); //Metropolis test

            //Calculate local energy
            if (type_energy_ == 0){
                DeltaE = wave.Local_energy_analytical(alpha);
            }
            if (type_energy_ == 1){
                DeltaE = wave.Local_energy_brute_force(alpha);
            }
            energies[cycle*N_ + n] += DeltaE;
        }
    }

    end_time_ = omp_get_wtime();
    time_ = end_time_ - start_time_;

}



//Brute force metropolis test
void Solver::Metropolis(double alpha){

    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]
    uniform_int_distribution<int> RIG(0, N_-1);    //Random integer genererator [0,N_]
    double tf_old, tf_new, P;
    int idx = RIG(gen);

    wave.Proposed_move(idx);    //Moves particle with index "idx" and calculates new sum of r^2 (for trial wave function)

    tf_old = wave.Trial_func(alpha, wave.r2_sum_old_);             //Trial wave function of old position
    tf_new = wave.Trial_func(alpha, wave.r2_sum_new_);           //Trial wave function of new position
    P = (tf_new*tf_new)/(tf_old*tf_old);                         //Metropolis test
    if (RDG(gen) <= P){
        for (int k =0 ; k < D_;  k++){                   //Update initial position
            wave.r_old_(idx,k) = wave.r_new_(idx,k);
        }
        wave.r2_sum_old_ = wave.r2_sum_new_;
    }
    else{
        for (int k =0 ; k < D_;  k++){
            wave.r_new_(idx,k) = wave.r_old_(idx,k);
        }
        wave.r2_sum_new_ = wave.r2_sum_old_;
    }
}

//Metropolis with importance sampling (non-interacting bosons)
void Solver::Metropolis_importance(double alpha){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]
    uniform_int_distribution<int> RIG(0, N_-1);    //Random integer genererator [0,N_]

    double tf_old, tf_new, P, greensfunc;
    int idx = RIG(gen);

    wave.Initialize_quantum_force(alpha,idx);                   //Calculate quantum force for particles before proposed move
    wave.Proposed_move_importance(idx);     //Moves particle with index "idx" and calculates new sum of r^2 (for trial wave function)
    wave.Update_quantum_force(alpha,idx);                          //Calculate quantum force for particles after proposed move
    greensfunc = Greens_function(idx);                         //Calculate Greens function

    tf_old = wave.Trial_func(alpha, wave.r2_sum_old_);             //Trial wave function of old position
    tf_new = wave.Trial_func(alpha, wave.r2_sum_new_);           //Trial wave function of new position
    P = greensfunc*(tf_new*tf_new)/(tf_old*tf_old);            //Metropolis-Hastings test
    if (RDG(gen) <= P){
        for (int k =0 ; k < D_;  k++){                   //Update initial position
            wave.r_old_(idx,k) = wave.r_new_(idx,k);
        }
        wave.r2_sum_old_ = wave.r2_sum_new_;
    }
    else{
        for (int k =0 ; k < D_;  k++){
            wave.r_new_(idx,k)= wave.r_old_(idx,k);
        }
        wave.r2_sum_new_ = wave.r2_sum_old_;
    }
}


//Method for evaluating Greens function
double Solver::Greens_function(int idx){
    double old_2_new, new_2_old;
    for (int dd = 0; dd < D_; dd++){
        old_2_new += pow((wave.r_old_(idx,dd)-wave.r_new_(idx,dd)- D_diff_*step_*wave.quantum_force_new_(dd)),2)/(4*D_diff_*step_);
        new_2_old += pow((wave.r_new_(idx,dd)-wave.r_old_(idx,dd)- D_diff_*step_*wave.quantum_force_old_(dd)),2)/(4*D_diff_*step_);
    }
    old_2_new = exp(-old_2_new);
    new_2_old = exp(-new_2_old);

    return old_2_new/new_2_old;
}


// Method for performing equilibration cycles
void Solver::Equilibrate(double alpha){
    for (int equi_c= 0; equi_c < equi_cycles_; equi_c++){
        for (int n = 0; n< N_; n++){
            (this->*metropolis_sampling)(alpha); //Metropolis test
        }
    }
}


void Solver::Write_to_file(string outfilename){
    ofstream ofile;
    ofile << "Cycles this thread: " << MC_ << endl;
    ofile.open(outfilename);
    ofile << setw(5) << "alpha" << setw(15) << "energy" << setw(15) << "variance" << setw(15) << "time_used" << endl;
    for (int i = 0; i < num_alphas_; i++){
        ofile << setw(5) << setprecision(8) << alphas_[i];
        ofile << setw(15) << setprecision(8) << energies_[i];
        ofile << setw(15) << setprecision(8) << variances_[i];
        ofile << setw(15) << setprecision(8) << alpha_list_times_[i] << endl;
    }
    ofile.close();
}
