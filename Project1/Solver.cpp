#include "Solver.hpp"
#include "Psi.hpp"

Solver::Solver(int N, int num_alphas, int MC, int D,int type_energy, int type_sampling){
    D_ = D;
    N_ = N;
    h_ = 1.0;                                          //Stepsize
    sum_ = 0;
    step_ = h_*pow(10,-4);
    type_energy_ = type_energy;
    type_sampling_ = type_sampling;
    wave.Declare_positions(N_, D_,h_, step_);

    num_alphas_ = num_alphas;
    MC_ = MC;

    alphas_ = new double[num_alphas_];                  //Variational parameter
    for (int i = 0; i < num_alphas_; i++) alphas_[i] = 0.1 + 0.05*i;
    energies_ = new double[num_alphas_];               //Array to hold energies for different values of alpha
    variances_ = new double[num_alphas_];              //Array to hold variances for different values of alpha


    if (type_sampling_ == 0){
        MC_method = &Solver::MonteCarlo;
        metropolis_sampling = &Solver::Metropolis;
    }

    if (type_sampling_ == 1){
        MC_method = &Solver::MonteCarlo;
    }

    if (type_sampling_ == 2){
        MC_method = &Solver::Gradient_descent;
    }

    if (type_sampling_ == 3){
        MC_method = &Solver::ADAM;
    }

    if (type_sampling_ != 0){
        D_diff_ = 0.5;                      //Diffusion constant in Greens function
        wave.Declare_quantum_force(D_diff_);
        metropolis_sampling = &Solver::Metropolis_importance;
    }



    (this->*MC_method)();

}


void Solver::MonteCarlo(){
    double alpha, energy, energy_squared, DeltaE, variance;

    for (int a=0; a < num_alphas_; a++){                //Loop over alpha values
        alpha = alphas_[a];
        energy = 0;
        energy_squared = 0;

        wave.r2_sum_old_ = wave.Initialize_positions();

        if (type_sampling_ != 0){
            wave.Initialize_quantum_force(alpha);
        }

        //Equilibration step: runs metropolis algorithm without sampling to equibrate system
        int equi_cycles = 10000;
        for (int equi_c= 0; equi_c < equi_cycles; equi_c++){
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


void Solver::MonteCarlo_GD(double *values, double alpha, string path){
    double energy, energy_squared, DeltaE, variance, DerivateE, Derivate_WF_E, sum_r, Derivate_WF;

    double *E_L_to_file;
    E_L_to_file = new double[N_*MC_];
    energy = 0;
    energy_squared = 0;

    wave.r2_sum_old_ = wave.Initialize_positions();
    //wave.Initialize_quantum_force(alpha);


    //Equilibration step: runs metropolis algorithm without sampling to equibrate system
    int equi_cycles = 10000;
    for (int equi_c= 0; equi_c < equi_cycles; equi_c++){
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
            E_L_to_file[cycle*N_ + n] = DeltaE;
            energy += DeltaE;
            energy_squared += DeltaE*DeltaE;
            sum_r = -wave.r2_sum_old_;
            Derivate_WF += sum_r;
            Derivate_WF_E += sum_r*DeltaE;
            }
        }

    Write_array_to_file(path, E_L_to_file, N_*MC_);
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

    wave.r2_sum_new_ = wave.Proposed_move_importance(idx);
    wave.Update_quantum_force(alpha);
    greensfunc = Greens_function(idx);

    tf_old = wave.Trial_func(alpha, wave.r2_sum_old_);             //Trial wave function of old position
    tf_new = wave.Trial_func(alpha, wave.r2_sum_new_);           //Trial wave function of new position
    P = greensfunc*(tf_new*tf_new)/(tf_old*tf_old);            //Metropolis test
    if (RDG(gen) <= P){
        for (int k =0 ; k < D_;  k++){                   //Update initial posistion
            wave.r_old_[idx][k] = wave.r_new_[k];
            wave.quantum_force_old_[idx][k] = wave.quantum_force_new_[k];
        }
        wave.r2_sum_old_ = wave.r2_sum_new_;
    }
}


double Solver::Greens_function(int idx){
    double old_2_new, new_2_old;
    for (int dd = 0; dd < D_; dd++){
        old_2_new += pow((wave.r_old_[idx][dd]-wave.r_new_[dd]- D_diff_*step_*wave.quantum_force_new_[dd]),2)/(4*D_diff_*step_);

        new_2_old += pow((wave.r_new_[dd]-wave.r_old_[idx][dd]- D_diff_*step_*wave.quantum_force_old_[idx][dd]),2)/(4*D_diff_*step_);

    }
    old_2_new = exp(-old_2_new);
    new_2_old = exp(-new_2_old);

    return old_2_new/new_2_old;
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
    ofile<<" "<<endl;
    ofile << "Timeused: " << time <<endl;
    ofile.close();
}

void Solver::Write_array_to_file(string outfilename, double *array, int len){
    ofstream ofile1;
    ofile1.open(outfilename);
    for (int i = 0; i < len; i++){
        ofile1 <<array[i]<<endl;
    }
    ofile1.close();
}

//Gradient Descent w/Importance sampling
void Solver::Gradient_descent(){
    string file, file2;   //path to where files should be stored
    int iterations = 50;
    double *alp_vals =  new double[iterations]; //Array to store alpha values
    for (int z = 0; z < iterations; z++){
        alp_vals[z] = 0;
    }
    double *values = new double[3];
    double Alphaa = 0.9;        //Initial guess for alpha
    double eta = 0.015;
    int counter = 0;


    //Eller parallellsier dette også for forskjellige alphaer
    cout <<"Alpha " << "Energy " << "Variance " << endl;
    for (int it = 0; it < iterations; it++){
        counter += 1;
        alp_vals[it] = Alphaa;
        file = to_string(N_) + "_part_alpha_" + to_string(Alphaa) + "_E_L_samples.txt";
        MonteCarlo_GD(values, Alphaa, file);
        Alphaa -= eta*values[2];
        cout <<Alphaa<<" " << values[0] << " " << values[1]<< " " << endl;
        if (values[1]< pow(10,-9)){
            break;
        }
    }

    //Paralelliser dette!! Fra her
    MC_ = pow(2,19);
    file2 ="OPTIMAL_ALPHA"+ to_string(N_) + "_part_alpha_" + to_string(Alphaa) + "_E_L_samples.txt";
    MonteCarlo_GD(values, Alphaa, file2);
    //Til her!

    ofstream ofile2;
    ofile2.open("alpha_values_GD.txt");
    for (int p = 0; p < counter; p++){
        ofile2 <<alp_vals[p]<<endl;
    }
    ofile2.close();

}

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
