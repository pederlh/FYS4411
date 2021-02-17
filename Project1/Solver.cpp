#include "Solver.hpp"
//#include "Psi.hpp"

Solver::Solver(int N, int num_alphas, int MC, int D,int type_energy, int type_sampling){
    D_ = D;
    N_ = N;
    type_energy_ = type_energy;
    type_sampling_ = type_sampling;
    wave.Declare_positions(N_, D_);

    num_alphas_ = num_alphas;
    MC_ = MC;

    alphas_ = new double[num_alphas_];                  //Variational parameter
    for (int i = 0; i < num_alphas_; i++) alphas_[i] = 0.1 + 0.05*i;
    energies_ = new double[num_alphas_];               //Array to hold energies for different values of alpha
    variances_ = new double[num_alphas_];              //Array to hold variances for different values of alpha
    h_ = 1.0;                                          //Stepsize
    sum_ = 0;
    step_ = h_*pow(10,-4);

    if (type_sampling_ == 0){
        MC_method = &Solver::MonteCarlo;
        metropolis_sampling = &Solver::Metropolis;
    }

    if (type_sampling_ != 0){
        wave.Declare_quantum_force();
        MC_method = &Solver::MonteCarlo;
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

        wave.Initialize_positions();

        if (type_sampling_ != 0){
            wave.Initialize_quantum_force(alpha);
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


void Solver::Metropolis(double alpha){

    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]
    uniform_int_distribution<int> RIG(0, N_-1);    //Random integer genererator [0,N_]
    double tf_old, tf_new, P;
    int idx = RIG(gen);

    for (int k = 0; k < D_; k++){
        wave.r_new_[k] = wave.r_old_[idx][k] + h_ * (RDG(gen) - 0.5);
        wave.r2_sum_new_ = wave.Update_r_sum(wave.r2_sum_new_, wave.r_old_[idx][k], wave.r_new_[k]);
    }

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

    for (int k = 0; k < D_; k++){
        wave.r_new_[k] = wave.r_old_[idx][k] + h_ * (RDG(gen) - 0.5);
        wave.r2_sum_new_ = wave.Update_r_sum(wave.r2_sum_new_, wave.r_old_[idx][k], wave.r_new_[k]);
    }

    wave.Update_quantum_force(alpha);
    greensfunc = wave.Greens_function(idx);

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
