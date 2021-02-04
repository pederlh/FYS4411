#include "Solver.hpp"

Solver::Solver(int N, int num_alphas, int MC, int D,int type_energy){
    D_ = D;
    N_ = N;
    num_alphas_ = num_alphas;
    MC_ = MC;
    alphas_ = new double[num_alphas_];                  //Variational parameter
    for (int i = 0; i < num_alphas_; i++) alphas_[i] = 0.1 + 0.05*i;
    energies_ = new double[num_alphas_];               //Array to hold energies for different values of alpha
    variances_ = new double[num_alphas_];              //Array to hold variances for different values of alpha
    h_ = 1.0;                                          //Stepsize
    sum_ = 0;
    step_ = h_*pow(10,-2);

    if (type_energy == 0){
        energy_calculation = &Solver::local_energy_analytical;
    }
    if (type_energy == 1){
        energy_calculation = &Solver::local_energy_brute_force;
    }
}

void Solver::MonteCarlo(){

    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator (0,1)
    uniform_int_distribution<int> RIG(0,N_-1);    //Random integer genererator (0,1)

    double alpha, energy, energy_squared, P, DeltaE, variance;
    r_old_ = new double*[N_];
    r_new_ = new double[D_];
    for (int i= 0; i< N_ ; i++){
        r_old_[i] = new double[D_];
    }
    double tf_old, tf_new, r2_old, r2_new;            //New and old trial wave function (and r^2 sums)
    int move_idx;


    for (int a=0; a < num_alphas_; a++){
        alpha = alphas_[a];
        energy = 0;
        energy_squared = 0;

        for (int j = 0; j < N_; j++){                       //Initial posistion
            for (int k = 0; k < D_; k++){
                r_old_[j][k] = h_ * (RDG(gen) - 0.5);
            }
        }

        r2_old = init_r_sum(r_old_);               //Inital sum of all r_i^2
        tf_old = trial_func(alpha, r2_old);       //Initial trial wave function

        for (int cycle = 0; cycle < MC_; cycle++){
            for (int n = 0; n < N_; n++){
                move_idx = RIG(gen);                      //Index of proposed moved particle
                r2_new = r2_old;
                for (int k = 0; k < D_; k++){
                    r_new_[k] = r_old_[move_idx][k] + h_ * (RDG(gen) - 0.5);
                    r2_new = update_r_sum(r2_new, r_old_[move_idx][k], r_new_[k]);
                }

                tf_new = trial_func(alpha, r2_new);           //Trial wave function of new position
                P = (tf_new*tf_new)/(tf_old*tf_old);            //Metropolis test
                if (RDG(gen)<= P){
                    for (int k =0 ; k< D_;  k++){                   //Update initial posistion
                        r_old_[move_idx][k] = r_new_[k];
                    }
                    tf_old = tf_new;
                    r2_old = r2_new;
                }
                DeltaE = (this->*energy_calculation)(alpha,r_old_,r2_old); //Points to either analytical expression for local energy or numerical
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

double Solver::init_r_sum(double **r){
    sum_ = 0;
    for (int i=0; i<N_; i++){
        for (int j =0; j <D_; j++){
            sum_ += r[i][j]*r[i][j];
        }
    }
    return sum_;
}

double Solver::update_r_sum(double sum, double r_init, double r_move){
    sum -= r_init*r_init;
    sum += r_move*r_move;
    return sum;
}

double Solver::trial_func(double alpha, double sum_r_squared){
    return exp(-alpha*sum_r_squared);
}



double Solver::local_energy_analytical(double alpha, double **r, double r_sum){
    sum_ = 0;
    for (int i=0; i<N_; i++){
        tmp_= 0;
        for (int j = 0; j< D_; j++){
            tmp_ += r[i][j]*r[i][j];
        }
        sum_ += -2*alpha*alpha*tmp_ + D_*alpha + (1./2)*tmp_;
    }
    return sum_;

}


double Solver::local_energy_brute_force(double alpha, double **r, double r_sum){
    double dr_p, dr_m;
    laplace_tf_ = 0.0;

    for (int dd = 0; dd < D_; dd++){
        dr_p = r_sum;
        dr_m = r_sum;
        for (int nn = 0; nn < N_; nn++){
            dr_p = update_r_sum(dr_p, r[nn][dd], r[nn][dd] + step_);
            dr_m = update_r_sum(dr_m, r[nn][dd], r[nn][dd] - step_);
        }
        laplace_tf_  += trial_func(alpha,dr_p) + trial_func(alpha, dr_m);
    }

    tf_middle_ = trial_func(alpha,r_sum);
    laplace_tf_ -= 2*D_*tf_middle_;
    laplace_tf_ /= (step_*step_*tf_middle_);

    return (1./2)*(-laplace_tf_ + r_sum);
}

void Solver::write_to_file(string outfilename, double time){
    ofstream ofile;
    ofile.open(outfilename);
    ofile << "alpha" << " " << "energy" << " " << "variance" << endl;
    for (int i = 0; i < num_alphas_; i++){
        ofile << alphas_[i] << " " << energies_[i] << " " << variances_[i] << endl;
    }
    ofile<<" "<<endl;
    ofile << "Time used : " << time <<"sec" <<endl;
    ofile.close();
}
