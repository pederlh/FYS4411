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


    // Initialize particle positions
    if (type_sampling_ == 3){
        wave.Declare_position_interaction(N_, D_,h_, step_, 1); // Interacting case
    }
    else{
        wave.Declare_position(N_, D_,h_, step_, 0);             // Non-interacting case
    }

    MC_ = MC;                       // Number of monte carlo cycles
    OBD_check_ = true;              // "true" gives calculation of one body density


    //Point to right member functions according to type of sampling
    if (type_sampling_ == 0){
        num_alphas_ = 8;                       // Number of alphas to test for
        metropolis_sampling = &Solver::Metropolis;

        main_method = &Solver::Alpha_list;
    }

    if (type_sampling_ == 1){
        num_alphas_ = 8;
        D_diff_ = 0.5;                                  // Diffusion constant in Greens function
        wave.Declare_quantum_force(D_diff_);
        metropolis_sampling = &Solver::Metropolis_importance;

        main_method = &Solver::Alpha_list;
    }

    if (type_sampling_ == 2){
        MC_optimal_run_ = MC_optimal_run;               // Number of MC cycles for run with optimal alpha value after GD
        D_diff_ = 0.5;
        wave.Declare_quantum_force(D_diff_);
        metropolis_sampling = &Solver::Metropolis_importance;

        main_method = &Solver::Gradient_descent;
        Interaction_or_not_GD = &Solver::MonteCarlo_GD_noninteracting;
        Interaction_or_not_optimal = &Solver::MonteCarlo_optval_noninteracting;
        tol_GD_ = 5*pow(10,-8);                                 // Acceptance tolerance for gradient descent
        eta_GD_ = learning_rate;                        // Learning rate for gradient descent
    }

    if (type_sampling_ == 3){
        MC_optimal_run_ = MC_optimal_run;
        D_diff_ = 0.5;
        wave.Declare_quantum_force(D_diff_);
        metropolis_sampling = &Solver::Metropolis_interacting;

        main_method = &Solver::Gradient_descent;
        Interaction_or_not_GD = &Solver::MonteCarlo_GD_interacting;
        Interaction_or_not_optimal = &Solver::MonteCarlo_optval_interacting;
        tol_GD_ = 1e-2;                                 // Acceptance tolerance for gradient descent
        eta_GD_ = learning_rate;                        // Learning rate for gradient descent
    }

    (this->*main_method)(shared_alphas);                               // Calls MC method and starts simulation
}

//Method for running Monte Carlo simulation for a list of values for variational parameter alpha
void Solver::Alpha_list(double * shared_alphas){
    double alpha, energy, energy_squared, DeltaE, variance;

    alphas_ = new double[num_alphas_];                 // Variational parameter
    for (int i = 0; i < num_alphas_; i++){
        alphas_[i] = 0.3 + 0.05*i; }                   // Fill array alphas
    energies_ = new double[num_alphas_];               // Array to hold energies for different values of alpha
    variances_ = new double[num_alphas_];              // Array to hold variances for different values of alpha
    alpha_list_times_ = new double[num_alphas_];

    for (int a=0; a < num_alphas_; a++){               //Loop over alpha values to perform MC simulations
        alpha = alphas_[a];
        energy = 0;
        energy_squared = 0;

        wave.r2_sum_old_ = wave.Initialize_positions();   // Initialize random initial matrix of positions

        //Equilibration step: runs metropolis algorithm without sampling to equilibrate system
        for (int equi_c= 0; equi_c < equi_cycles_; equi_c++){
            for (int n = 0; n< N_; n++){

                wave.r2_sum_new_ = wave.r2_sum_old_;

                (this->*metropolis_sampling)(alpha); //Metropolis test
            }
        }

        start_time_ = omp_get_wtime();
        //Monte Carlo simulation with metropolis sampling
        for (int cycle = 0; cycle < MC_; cycle++){
            for (int n = 0; n < N_; n++){

                wave.r2_sum_new_ = wave.r2_sum_old_;

                (this->*metropolis_sampling)(alpha);  //Metropolis test

                //Calculate local energy
                if (type_energy_ == 0){
                    DeltaE = wave.Local_energy_analytical(alpha);
                }
                if (type_energy_==1){
                    DeltaE = wave.Local_energy_brute_force(alpha);
                }

                energy += DeltaE;
                energy_squared += DeltaE*DeltaE;
                }

            }
        end_time_ = omp_get_wtime();
        time_ = end_time_ - start_time_;
        energy /= (MC_*N_);
        energy_squared /= (MC_*N_);
        variance = energy_squared - energy*energy;
        energies_[a] = energy;
        variances_[a] = variance;
        alpha_list_times_[a] = time_;
    }
}

//Gradient Descent w/importance sampling for interacting and non-interacting bosons
void Solver::Gradient_descent(double * shared_alphas){
    string file;
    int iterations = 500;                                // Max number iterations of GD
    double *alpha_vals_GD = new double[iterations];     // Array to store alpha values from GD
    for (int z = 0; z < iterations; z++){
        alpha_vals_GD[z] = 0;
    }
    double *values = new double[3];                     // Array to contain energy, variance and energy derivative wrt alpha;
    double alpha_guess = 0.55;                           // Initial guess for alpha
    int counter = 0;                                    // Counter to keep track of actual number of iterations
    double tol;

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
        if (type_sampling_ ==3){tol = tol_GD_*values[0];}
        else{tol = tol_GD_;}

        //Breaks GD if alpha provides acceptably low sample variance
        if (values[1] < tol ){
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
        for (int i =0; i < omp_get_num_threads();i++){
            sum += shared_alphas[i];
        }

        alpha_guess = sum/omp_get_num_threads();
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

    if(type_sampling_==3){
        file = "INTERACTION_OPTIMAL_ALPHA"+ to_string(N_) + "_N_stringID_" + to_string(thread_ID_) +
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

    wave.r2_sum_old_ = wave.Initialize_positions();         //Initialize random initial matrix of positions

    //Equilibration step: runs metropolis algorithm without sampling to equilibrate system
    Equilibrate(alpha);

    energy = 0;
    energy_squared = 0;
    Derivate_WF = 0;
    Derivate_WF_E = 0;

    //Monte Carlo simulation with metropolis sampling
    for (int cycle = 0; cycle < MC_; cycle++){
        for (int n = 0; n < N_; n++){

            wave.r2_sum_new_ = wave.r2_sum_old_;
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

//Method for running Monte Carlo simulation with gradient descent (interacting bosons)
//  - Caluculates the derivative of local energy w/regards to alpha to use in GD
void Solver::MonteCarlo_GD_interacting(double *values, double alpha){
    double energy, energy_squared, DeltaE, variance, DerivateE, Derivate_WF_E, sum_r, Derivate_WF;

    wave.r2_sum_old_ = wave.Initialize_positions(); //Initialize random initial matrix of positions

    //Equilibration step: runs metropolis algorithm without sampling to equilibrate system
    Equilibrate(alpha);

    energy = 0;
    energy_squared = 0;
    Derivate_WF = 0;
    Derivate_WF_E = 0;

    //Monte Carlo simulation with metropolis sampling
    for (int cycle = 0; cycle < MC_; cycle++){
        for (int n = 0; n < N_; n++){

            wave.r2_sum_new_ = wave.r2_sum_old_;

            (this->*metropolis_sampling)(alpha); //Metropolis test
            DeltaE = wave.Local_energy_interaction(alpha);  //Calculate local energy

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
    DerivateE = 2*(Derivate_WF_E - Derivate_WF*energy);     //The derivative of the local energy w/regards to alpha

    values[0] = energy;
    values[1] = variance;
    values[2] = DerivateE;

}

//Method for running Monte Carlo simulation for one (optimized) value of variational parameter alpha (non-interacting bosons)
void Solver::MonteCarlo_optval_noninteracting(double alpha, double *energies){
    double DeltaE;
    int num_bins = 100;                     // Number of bins OBD calculation
    double max_radi = 5.0;                  // Max radius OBD calculation
    double *bins = new double[num_bins];    // Array of bins where each bins[i] represents i*radius for one body density calculation

    //If one body density calculation is chosen: Initialize radius and bins array
    if (OBD_check_ == true){
        for (int i =0; i < num_bins; i++){
            bins[i] = 0.0;
        }
        radi_ = max_radi/num_bins;
    }

    wave.r2_sum_old_ = wave.Initialize_positions();     // Initialize random initial matrix of positions

    //Equilibration step: runs metropolis algorithm without sampling to equilibrate system
    Equilibrate(alpha);
    start_time_ = omp_get_wtime();

    //Monte Carlo simulation with metropolis sampling
    for (int cycle = 0; cycle < MC_optimal_run_; cycle++){
        for (int n = 0; n < N_; n++){

            if (OBD_check_ == true){One_body_density(bins);}       //Count particle positions for one body density calculation
            wave.r2_sum_new_ = wave.r2_sum_old_;
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

    //Write average particle distribution to file
    if (OBD_check_ == true){
        string OBD_file = "One_body_density_N_" + to_string(N_) + "_stringID_" + to_string(thread_ID_) + "_alpha_" + to_string(alpha) + ".txt";

        ofstream ofile2;
        double scale = (4./3)*M_PI*radi_*radi_;
        ofile2.open(OBD_file);
        for (int i = 0; i < num_bins; i ++){
            bins[i] /= (MC_optimal_run_*N_*scale*(pow(i+1,2)-pow(i,2)));
            //bins[i] /= (MC_optimal_run_*N_*pow(radi_,(D_-1)));
            //ELLER pow(r_old[i],D-1) ELLER volumet av kuleskallet V = 4/3*((i+1)**2-i**2)*radi**2;
            ofile2 << setprecision(15) <<bins[i]<<endl;
        }

        ofile2.close();
    }

}

//Method for running Monte Carlo simulation for one (optimized) value of variational parameter alpha (interacting bosons)
void Solver::MonteCarlo_optval_interacting(double alpha, double *energies){
    double DeltaE;
    int num_bins = 100;                     // Number of bins OBD calculation
    double max_radi = 5.0;                  // Max radius OBD calculation
    double *bins = new double[num_bins];    // Array of bins where each bins[i] represents i*radius for one body density calculation

    //If one body density calculation is chosen: Initialize radius and bins array
    if (OBD_check_ == true){
        for (int i =0; i < num_bins; i++){
            bins[i] = 0.0;
        }
        radi_ = max_radi/num_bins;
    }

    wave.r2_sum_old_ = wave.Initialize_positions();         //Initialize random initial matrix of positions

    //Equilibration step: runs metropolis algorithm without sampling to equilibrate system
    Equilibrate(alpha);

    //Monte Carlo simulation with metropolis sampling
    start_time_ = omp_get_wtime();
    for (int cycle = 0; cycle < MC_optimal_run_; cycle++){
        for (int n = 0; n < N_; n++){
            if (OBD_check_ == true){One_body_density(bins);}             //Count particle positions for one body density calculation
            wave.r2_sum_new_ = wave.r2_sum_old_;
            (this->*metropolis_sampling)(alpha);              //Metropolis test
            DeltaE = wave.Local_energy_interaction(alpha);    //Calculate local energy
            energies[cycle*N_ + n] += DeltaE;
        }
    }
    end_time_ = omp_get_wtime();
    time_ = end_time_ - start_time_;
    //Write average particle distribution to file
    if (OBD_check_ == true){
        string OBD_file = "Interaction_One_body_density_N_" + to_string(N_) + "_stringID_" + to_string(thread_ID_) + "_alpha_" + to_string(alpha) + ".txt";
        ofstream ofile2;
        double scale = (4./3)*M_PI*radi_*radi_;
        ofile2.open(OBD_file);
        for (int i = 0; i < num_bins; i ++){
            bins[i] /= (MC_optimal_run_*N_*scale*(pow(i+1,2)-pow(i,2)));
            //bins[i] /= (MC_optimal_run_*N_*pow(radi_,(D_-1))); //ELLER pow(r_old[i],D-1);
            ofile2 << setprecision(15) <<bins[i]<<endl;
        }
        ofile2.close();
    }
}


//Brute force metropolis test
void Solver::Metropolis(double alpha){

    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]
    uniform_int_distribution<int> RIG(0, N_-1);    //Random integer genererator [0,N_]
    double tf_old, tf_new, P;
    int idx = RIG(gen);

    wave.r2_sum_new_ = wave.Proposed_move(idx);    //Moves particle with index "idx" and calculates new sum of r^2 (for trial wave function)

    tf_old = wave.Trial_func(alpha, wave.r2_sum_old_);             //Trial wave function of old position
    tf_new = wave.Trial_func(alpha, wave.r2_sum_new_);           //Trial wave function of new position
    P = (tf_new*tf_new)/(tf_old*tf_old);                         //Metropolis test
    if (RDG(gen) <= P){
        for (int k =0 ; k < D_;  k++){                   //Update initial position
            wave.r_old_[idx][k] = wave.r_new_[k];
        }
        wave.r2_sum_old_ = wave.r2_sum_new_;
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
    wave.r2_sum_new_ = wave.Proposed_move_importance(idx);     //Moves particle with index "idx" and calculates new sum of r^2 (for trial wave function)
    wave.Update_quantum_force(alpha);                          //Calculate quantum force for particles after proposed move
    greensfunc = Greens_function(idx);                         //Calculate Greens function

    tf_old = wave.Trial_func(alpha, wave.r2_sum_old_);             //Trial wave function of old position
    tf_new = wave.Trial_func(alpha, wave.r2_sum_new_);           //Trial wave function of new position
    P = greensfunc*(tf_new*tf_new)/(tf_old*tf_old);            //Metropolis-Hastings test
    if (RDG(gen) <= P){
        for (int k =0 ; k < D_;  k++){                   //Update initial position
            wave.r_old_[idx][k] = wave.r_new_[k];
        }
        wave.r2_sum_old_ = wave.r2_sum_new_;
    }
}

//Metropolis with importance sampling (interacting bosons)
void Solver::Metropolis_interacting(double alpha){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]
    uniform_int_distribution<int> RIG(0, N_-1);    //Random integer genererator [0,N_]

    double tf_old, tf_new, P, greensfunc;
    int idx = RIG(gen);

    wave.Initialize_quantum_force_interaction(alpha, idx);          //Calculate quantum force for particles before proposed move
    wave.r2_sum_new_ = wave.Proposed_move_interaction(idx);         //Moves particle with index "idx" and calculates new sum of r^2 (for trial wave function)
    wave.Update_quantum_force_interaction(alpha, idx);              //Calculate quantum force for particles after proposed move
    greensfunc = Greens_function(idx);                              //Calculate Greens function

    tf_old = wave.Trial_func_interaction(alpha, wave.r2_sum_old_, "old",idx);             //Trial wave function of old position
    tf_new = wave.Trial_func_interaction(alpha, wave.r2_sum_new_, "new",idx);           //Trial wave function of new position
    P = greensfunc*(tf_new*tf_new)/(tf_old*tf_old);                                     //Metropolis-Hastings test
    if (RDG(gen) <= P){
        for (int k =0 ; k < D_;  k++){                   //Update initial posistion
            wave.r_old_[idx][k] = wave.r_new_[k];
        }
        wave.r2_sum_old_ = wave.r2_sum_new_;
    }
}


//Method for evaluating Greens function
double Solver::Greens_function(int idx){
    double old_2_new, new_2_old;
    for (int dd = 0; dd < D_; dd++){
        old_2_new += pow((wave.r_old_[idx][dd]-wave.r_new_[dd]- D_diff_*step_*wave.quantum_force_new_[dd]),2)/(4*D_diff_*step_);
        new_2_old += pow((wave.r_new_[dd]-wave.r_old_[idx][dd]- D_diff_*step_*wave.quantum_force_old_[dd]),2)/(4*D_diff_*step_);
    }
    old_2_new = exp(-old_2_new);
    new_2_old = exp(-new_2_old);

    return old_2_new/new_2_old;
}

// Method for counting numbers of particles in spherical shells centered at reference particle indexed 0. Shell radius is multiples of r0.
void Solver::One_body_density(double *bins){
    double r;                   // Distance from reference particles
    int bin_nr;                 // Corresponding bin in histogram
    for (int i = 1; i<N_;i++){
        r = 0.0;
        for (int j = 0; j < D_ ; j++){
            r += pow(wave.r_old_[0][j] - wave.r_old_[i][j],2);
        }
        r = sqrt(r);
        bin_nr = trunc(r/radi_);
        bins[bin_nr] +=1;
    }
}

// Method for performing equilibration cycles
void Solver::Equilibrate(double alpha){
    for (int equi_c= 0; equi_c < equi_cycles_; equi_c++){
        for (int n = 0; n< N_; n++){
            wave.r2_sum_new_ = wave.r2_sum_old_;

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
