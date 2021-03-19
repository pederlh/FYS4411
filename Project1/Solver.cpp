#include "Solver.hpp"
#include "Psi.hpp"

Solver::Solver(int N, int num_alphas, int MC, int MC_optimal_run, int D, int type_energy, int type_sampling, int thread_ID){
    D_ = D;
    N_ = N;
    h_ = 1.0;                       // Stepsize used for proposed move
    sum_ = 0;
    step_ = h_*pow(10,-4);         //Stepsize used in differentiation
    thread_ID_ = thread_ID;
    equi_cycles_ = 5000;

    type_energy_ = type_energy;
    type_sampling_ = type_sampling;

    if (type_sampling_ == 3){
        wave.Declare_position_interaction(N_, D_,h_, step_, 1);
    }
    else{
        wave.Declare_position(N_, D_,h_, step_, 0);
    }

    num_alphas_ = num_alphas;                          //
    MC_ = MC;
    MC_optimal_run_ = MC_optimal_run;
    OBD_ = 1;

    if (type_sampling_ == 0){
        MC_method = &Solver::MonteCarlo_alpha_list;
        metropolis_sampling = &Solver::Metropolis;
    }

    if (type_sampling_ == 1){
        D_diff_ = 0.5;                      //Diffusion constant in Greens function
        wave.Declare_quantum_force(D_diff_);
        metropolis_sampling = &Solver::Metropolis_importance;
        MC_method = &Solver::MonteCarlo_alpha_list;
    }

    if (type_sampling_ == 2){
        D_diff_ = 0.5;                      //Diffusion constant in Greens function
        wave.Declare_quantum_force(D_diff_);
        metropolis_sampling = &Solver::Metropolis_importance;
        MC_method = &Solver::Gradient_descent;
        Interaction_or_not_GD = &Solver::MonteCarlo_GD;
        Interaction_or_not_optimal = &Solver::MonteCarlo;
        tol_GD_ = 1e-9;
        eta_GD_ = 0.015;
    }

    if (type_sampling_ == 3){
        D_diff_ = 0.5;                      //Diffusion constant in Greens function
        wave.Declare_quantum_force(D_diff_);
        MC_method = &Solver::Gradient_descent;
        metropolis_sampling = &Solver::Metropolis_interaction;
        Interaction_or_not_GD = &Solver::MonteCarlo_GD_interaction;
        Interaction_or_not_optimal = &Solver::MonteCarlo_interaction;
        tol_GD_ = 1e-3;
        eta_GD_ = 0.01;
    }

    (this->*MC_method)();

}



void Solver::One_body_density(double *bins){
    double r;
    int bin_nr;
    for (int i =1; i<N_;i++){
        r = 0.0;
        for (int j = 0; j < D_ ; j++){
            r += pow(wave.r_old_[0][j] - wave.r_old_[i][j],2);
        }
        r = sqrt(r);
        bin_nr = trunc(r/radi_);
        bins[bin_nr] +=1;
    }
}
void Solver::MonteCarlo_alpha_list(){
    double alpha, energy, energy_squared, DeltaE, variance;

    alphas_ = new double[num_alphas_];                 //Variational parameter
    for (int i = 0; i < num_alphas_; i++) alphas_[i] = 0.1 + 0.05*i;
    energies_ = new double[num_alphas_];               //Array to hold energies for different values of alpha
    variances_ = new double[num_alphas_];              //Array to hold variances for different values of alpha

    for (int a=0; a < num_alphas_; a++){                //Loop over alpha values
        alpha = alphas_[a];
        energy = 0;
        energy_squared = 0;

        wave.r2_sum_old_ = wave.Initialize_positions();

        //Equilibration step: runs metropolis algorithm without sampling to equibrate system
        for (int equi_c= 0; equi_c < equi_cycles_; equi_c++){
            for (int n = 0; n< N_; n++){
                wave.r2_sum_new_ = wave.r2_sum_old_;

                (this->*metropolis_sampling)(alpha); //Metropolis test
            }
        }

        for (int cycle = 0; cycle < MC_; cycle++){
            for (int n = 0; n < N_; n++){

                wave.r2_sum_new_ = wave.r2_sum_old_;

                (this->*metropolis_sampling)(alpha); //Metropolis test

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
        energy /= (MC_*N_);
        energy_squared /= (MC_*N_);
        variance = energy_squared - energy*energy;
        energies_[a] = energy;
        variances_[a] = variance;
    }
}

void Solver::MonteCarlo(double alpha, double *energies){
    double DeltaE;
    int num_bins = 100;
    double max_radi = 5.0;
    double *bins = new double[num_bins];

    if (OBD_ == 1){
        for (int i =0; i < num_bins; i++){
            bins[i] = 0.0;
        }
        radi_ = max_radi/num_bins;
    }

    wave.r2_sum_old_ = wave.Initialize_positions();

    //Equilibration step: runs metropolis algorithm without sampling to equibrate system
    for (int equi_c= 0; equi_c < equi_cycles_; equi_c++){
        for (int n = 0; n< N_; n++){
            wave.r2_sum_new_ = wave.r2_sum_old_;

            (this->*metropolis_sampling)(alpha); //Metropolis test
        }
    }

    for (int cycle = 0; cycle < MC_optimal_run_; cycle++){
        for (int n = 0; n < N_; n++){

            // Trenger: Radius, antall bins
            // Finn avstand per partikkel
            // Heltallsdivisjon -> plasser i bins
            // Ta snitt til slutt
            if (OBD_ ==1){One_body_density(bins);}

            wave.r2_sum_new_ = wave.r2_sum_old_;

            (this->*metropolis_sampling)(alpha); //Metropolis test


            if (type_energy_ == 0){
                DeltaE = wave.Local_energy_analytical(alpha);
            }
            if (type_energy_ == 1){
                DeltaE = wave.Local_energy_brute_force(alpha);
            }


            energies[cycle*N_ + n] += DeltaE;
        }
    }

    if (OBD_ == 1){
        string OBD_file = "One_body_density_N_" + to_string(N_) + "_stringID_" + to_string(thread_ID_) + "_alpha_" + to_string(alpha) + ".txt";
        ofstream ofile2;
        ofile2.open(OBD_file);
        for (int i = 0; i < num_bins; i ++){
            bins[i] /= (MC_optimal_run_*N_*pow(i*radi_,(D_-1))); //ELLER pow(r_old[i],D-1);
            ofile2 << setprecision(15) <<bins[i]<<endl;
        }
        ofile2.close();
    }
}

void Solver::MonteCarlo_interaction(double alpha, double *energies){
    double DeltaE;
    int num_bins = 100;
    double max_radi = 5.0;
    double *bins = new double[num_bins];

    if (OBD_==1){
        for (int i =0; i < num_bins; i++){
            bins[i] = 0.0;
        }
        radi_ = max_radi/num_bins;
    }

    wave.r2_sum_old_ = wave.Initialize_positions();

    //Equilibration step: runs metropolis algorithm without sampling to equibrate system
    for (int equi_c= 0; equi_c < equi_cycles_; equi_c++){
        for (int n = 0; n< N_; n++){
            wave.r2_sum_new_ = wave.r2_sum_old_;

            (this->*metropolis_sampling)(alpha); //Metropolis test
        }
    }

    for (int cycle = 0; cycle < MC_optimal_run_; cycle++){
        for (int n = 0; n < N_; n++){

            if (OBD_==1){One_body_density(bins);}
            wave.r2_sum_new_ = wave.r2_sum_old_;
            (this->*metropolis_sampling)(alpha); //Metropolis test
            DeltaE = wave.Local_energy_interaction(alpha);
            energies[cycle*N_ + n] += DeltaE;
        }
    }

    if (OBD_==1){
        string OBD_file = "Interaction_One_body_density_N_" + to_string(N_) + "_stringID_" + to_string(thread_ID_) + "_alpha_" + to_string(alpha) + ".txt";
        ofstream ofile2;
        ofile2.open(OBD_file);
        for (int i = 0; i < num_bins; i ++){
            bins[i] /= (MC_optimal_run_*N_*pow(i*radi_,(D_-1))); //ELLER pow(r_old[i],D-1);
            ofile2 << setprecision(15) <<bins[i]<<endl;
        }
        ofile2.close();
    }
}

void Solver::MonteCarlo_GD(double *values, double alpha){
    double energy, energy_squared, DeltaE, variance, DerivateE, Derivate_WF_E, sum_r, Derivate_WF;

    energy = 0;
    energy_squared = 0;

    wave.r2_sum_old_ = wave.Initialize_positions();


    //Equilibration step: runs metropolis algorithm without sampling to equibrate system
    for (int equi_c= 0; equi_c < equi_cycles_; equi_c++){
        for (int n = 0; n< N_; n++){
            wave.r2_sum_new_ = wave.r2_sum_old_;

            (this->*metropolis_sampling)(alpha); //Metropolis test
        }
    }


    energy = 0;
    energy_squared = 0;

    for (int cycle = 0; cycle < MC_; cycle++){
        for (int n = 0; n < N_; n++){

            wave.r2_sum_new_ = wave.r2_sum_old_;

            (this->*metropolis_sampling)(alpha); //Metropolis test

            DeltaE = wave.Local_energy_brute_force(alpha);
            /*
            if (type_energy_ == 0){
                DeltaE = wave.Local_energy_analytical(alpha);
            }
            if (type_energy_==1){
                DeltaE = wave.Local_energy_brute_force(alpha);
            }
            */
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
    DerivateE = 2*(Derivate_WF_E - Derivate_WF*energy);

    values[0] = energy;
    values[1] = variance;
    values[2] = DerivateE;

}

void Solver::MonteCarlo_GD_interaction(double *values, double alpha){
    double energy, energy_squared, DeltaE, variance, DerivateE, Derivate_WF_E, sum_r, Derivate_WF;

    energy = 0;
    energy_squared = 0;

    wave.r2_sum_old_ = wave.Initialize_positions();

    //Equilibration step: runs metropolis algorithm without sampling to equibrate system
    for (int equi_c = 0; equi_c < equi_cycles_; equi_c++){
        for (int n = 0; n< N_; n++){
            wave.r2_sum_new_ = wave.r2_sum_old_;

            (this->*metropolis_sampling)(alpha); //Metropolis test
        }
    }

    for (int cycle = 0; cycle < MC_; cycle++){
        for (int n = 0; n < N_; n++){

            wave.r2_sum_new_ = wave.r2_sum_old_;

            (this->*metropolis_sampling)(alpha); //Metropolis test
            DeltaE = wave.Local_energy_interaction(alpha);

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
    DerivateE = 2*(Derivate_WF_E - Derivate_WF*energy);

    values[0] = energy;
    values[1] = variance;
    values[2] = DerivateE;

}


void Solver::Metropolis(double alpha){

    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]
    uniform_int_distribution<int> RIG(0, N_-1);    //Random integer genererator [0,N_]
    double tf_old, tf_new, P;
    int idx = RIG(gen);

    wave.r2_sum_new_ = wave.Proposed_move(idx);    //Moves particle with index "idx" and calculates new sum of r^2

    tf_old = wave.Trial_func(alpha, wave.r2_sum_old_);             //Trial wave function of old position
    tf_new = wave.Trial_func(alpha, wave.r2_sum_new_);           //Trial wave function of new position
    P = (tf_new*tf_new)/(tf_old*tf_old);                         //Metropolis test
    if (RDG(gen) <= P){
        for (int k =0 ; k < D_;  k++){                   //Update initial posistion
            wave.r_old_[idx][k] = wave.r_new_[k];
        }
        wave.r2_sum_old_ = wave.r2_sum_new_;
    }
}


void Solver::Metropolis_importance(double alpha){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]
    uniform_int_distribution<int> RIG(0, N_-1);    //Random integer genererator [0,N_]

    double tf_old, tf_new, P, greensfunc;
    int idx = RIG(gen);

    wave.Initialize_quantum_force(alpha,idx);
    wave.r2_sum_new_ = wave.Proposed_move_importance(idx);
    wave.Update_quantum_force(alpha);
    greensfunc = Greens_function(idx);

    tf_old = wave.Trial_func(alpha, wave.r2_sum_old_);             //Trial wave function of old position
    tf_new = wave.Trial_func(alpha, wave.r2_sum_new_);           //Trial wave function of new position
    P = greensfunc*(tf_new*tf_new)/(tf_old*tf_old);            //Metropolis test
    if (RDG(gen) <= P){
        for (int k =0 ; k < D_;  k++){                   //Update initial posistion
            wave.r_old_[idx][k] = wave.r_new_[k];
        }
        wave.r2_sum_old_ = wave.r2_sum_new_;
    }
}

void Solver::Metropolis_interaction(double alpha){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]
    uniform_int_distribution<int> RIG(0, N_-1);    //Random integer genererator [0,N_]

    double tf_old, tf_new, P, greensfunc;
    int idx = RIG(gen);

    wave.Initialize_quantum_force_interaction(alpha, idx);
    wave.r2_sum_new_ = wave.Proposed_move_interaction(idx);
    wave.Update_quantum_force_interaction(alpha, idx);
    greensfunc = Greens_function(idx);

    tf_old = wave.Trial_func_interaction(alpha, wave.r2_sum_old_, "old",idx);             //Trial wave function of old position
    tf_new = wave.Trial_func_interaction(alpha, wave.r2_sum_new_, "new",idx);           //Trial wave function of new position
    P = greensfunc*(tf_new*tf_new)/(tf_old*tf_old);            //Metropolis test
    if (RDG(gen) <= P){
        for (int k =0 ; k < D_;  k++){                   //Update initial posistion
            wave.r_old_[idx][k] = wave.r_new_[k];
        }
        wave.r2_sum_old_ = wave.r2_sum_new_;
    }
}


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


//Gradient Descent w/Importance sampling for interacting and non-interacting case
void Solver::Gradient_descent(){
    string file;                                        // path to where files should be stored
    int iterations = 50;                                // Max number iterations gradient descent
    double *alpha_vals_GD = new double[iterations];     // Array to store alpha values
    for (int z = 0; z < iterations; z++){
        alpha_vals_GD[z] = 0;
    }
    double *values = new double[3];
    double alpha_guess = 0.9;                           // Initial guess for alpha
    int counter = 0;                                    // Counter to keep track of actual number of iterations
    double E_old;


     #pragma omp master
     {
         if (omp_get_num_threads() == 1) cout << "Start gradient descent" << endl;
         else cout << "Start gradient descent (showing progress master thread)" << endl << endl;

         cout << setw(10) << "Alpha" << setw(12) << "Energy" << setw(16) << "Variance" << endl;
         cout << "--------------------------------------" << endl;
     }
    //cout << "Thread made it here: " << to_string(thread_ID_) << endl;
    for (int i = 0; i < iterations; i++){
        counter++;
        alpha_vals_GD[i] = alpha_guess;

        (this->*Interaction_or_not_GD)(values, alpha_guess);
        # pragma omp master
        {
            cout << setw(10) << setprecision(8) << alpha_guess << setw(12) << values[0] << setw(16) << values[1] << endl;
        }

        if (values[1] < tol_GD_){
            break;
        }
        alpha_guess -= eta_GD_*values[2];
    }

    cout <<"Number of iterations of gradient descent = " << counter << " for thread " << to_string(thread_ID_) << endl;

    # pragma omp barrier
    # pragma omp master
    {cout << "Gradient descent finished, starting main MC calculations..." << endl;}

    // Optimal run
    double *optimal_energies = new double[N_*MC_optimal_run_];

    if(type_sampling_==3){
    file = "INTERACTION_OPTIMAL_ALPHA"+ to_string(N_) + "_N_stringID_" + to_string(thread_ID_) +
            "_alpha_" + to_string(alpha_guess) + "_E_L_samples.txt";
    }

    if(type_sampling_==2){
    file = "OPTIMAL_ALPHA"+ to_string(N_) + "_N_stringID_" + to_string(thread_ID_) +
            "_alpha_" + to_string(alpha_guess) + "_E_L_samples.txt";
    }


    (this->*Interaction_or_not_optimal)(alpha_guess, optimal_energies);

    ofstream ofile;
    ofile.open(file);
    for (int i = 0; i < N_*MC_optimal_run_; i++){
        ofile << setprecision(15) << optimal_energies[i] << endl;
    }
    ofile.close();
}

void Solver::Write_to_file(string outfilename, double time){
    ofstream ofile;
    ofile.open(outfilename);
    ofile << setw(5) << "alpha" << setw(15) << "energy" << setw(15) << "variance" << endl;
    for (int i = 0; i < num_alphas_; i++){
        ofile << setw(5) << setprecision(8) << alphas_[i];
        ofile << setw(15) << setprecision(8) << energies_[i];
        ofile << setw(15) << setprecision(8) << variances_[i] << endl;
    }
    ofile << " " <<endl;
    ofile << "Cycles this thread: " << MC_ << endl;
    ofile << "Time used: " << time << " s" << endl;
    ofile.close();
}


/*
void Solver::ADAM(){
    double *values;
    values = new double[3];
    string path;
    double avg_1_mom = 0.0;         //Average first momentum
    double avg_2_mom = 0.0;         //Average second momentum
    double B_1 = 0.9;               //Decay rate for average first momentum
    double B_2 = 0.999;             //Decay rate for average second momentum
    double scaled_eta, update, eps;  //Scaled_eta = adaptive learning rate
    double Alphaa = 0.9;        //Initial guess for alpha
    double eta = 0.01;          //Initial learning rate

    int iterations = 50;
    cout <<"Alpha " << "Energy " << "Variance " << endl;
    for (int it = 0; it < iterations; it++){
        path = "./Results/StatisticalAnalysis/" + to_string(N_) + "_part/alpha_" + to_string(Alphaa) + "/E_L_samples.txt";
        MonteCarlo_GD(values, Alphaa, path);
        values[2] *= eta;
        avg_1_mom = B_1*avg_1_mom + (1-B_1)*values[2];
        avg_2_mom = B_2*avg_2_mom + (1-B_2)*values[2]*values[2];
        scaled_eta = eta*sqrt(1-pow(B_2,it+1))/(1-pow(B_1,it+1));
        eps = pow(10,-8)*sqrt(1-pow(B_2,it+1));

        update = scaled_eta*avg_1_mom/(sqrt(avg_2_mom)+eps);
        Alphaa -= update;
        cout <<Alphaa<<" " << values[0] << " " << values[1]<< " " << endl;
        if (values[1]< pow(10,-9)){
            break;
        }

    }
}
*/
