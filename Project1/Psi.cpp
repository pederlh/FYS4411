#include "Psi.hpp"

Psi::Psi(){
    h_ = 1.0;
    step_ = h_*pow(10,-4);
}

void Psi::Declare_positions(int N, int D){

    D_ = D;
    N_ = N;

    r_new_ = new double[D_];
    r_old_ = new double*[N_];
    for (int i= 0; i< N_ ; i++){
        r_old_[i] = new double[D_];
    }

}

void Psi::Declare_quantum_force(){

    D_diff_ = 0.5; //Diffusion constant in Greens function

    quantum_force_new_ = new double[D_];
    quantum_force_old_ = new double*[N_];
    for (int i= 0; i< N_ ; i++){
        quantum_force_old_[i] = new double[D_];
    }

}

void Psi::Initialize_positions(){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]

    r2_sum_old_ = 0;

    for (int j = 0; j < N_; j++){                       //Initial posistion
        for (int k = 0; k < D_; k++){
            r_old_[j][k] = h_ * (RDG(gen) - 0.5);
        }
    }

    for (int i = 0; i < N_; i++){
        for (int j =0; j <D_; j++){
            r2_sum_old_ += r_old_[i][j]*r_old_[i][j];
        }
    }
}



void Psi::Initialize_quantum_force(double alpha){
    for (int j = 0; j < N_; j++){
        for (int k = 0; k < D_; k++){
            quantum_force_old_[j][k] = -4*alpha*r_old_[j][k];
        }
    }
}

void Psi::Update_quantum_force(double alpha){
    for (int dd =0; dd< D_; dd++){
        quantum_force_new_[dd] = -4*alpha*r_new_[dd];
    }
}

double Psi::Greens_function(int idx){
    double tmp1, tmp2;
    for (int dd = 0; dd < D_; dd++){
        tmp1 += pow((r_old_[idx][dd]-r_new_[dd]- D_diff_*step_*quantum_force_new_[dd]),2)*(1/(4*D_diff_*step_));
        tmp2 += pow((r_new_[dd]-r_old_[idx][dd]- D_diff_*step_*quantum_force_old_[idx][dd]),2)*(1/(4*D_diff_*step_));

    }
    tmp1 = exp(-tmp1);
    tmp2 = exp(-tmp2);
    return tmp1/tmp2;
    //return exp(tmp1/tmp2);
    //return 1.0;
}


double Psi::Update_r_sum(double sum, double r_init, double r_move){
    sum -= r_init*r_init;
    sum += r_move*r_move;
    return sum;
}

double Psi::Trial_func(double alpha, double sum_r_squared){
    return exp(-alpha*sum_r_squared);
}

double Psi::Local_energy_analytical(double alpha){
    return D_*N_*alpha + (1-4*alpha*alpha)*(1./2)*r2_sum_old_;
}

double Psi::Local_energy_brute_force(double alpha){
    double dr_p, dr_m;
    laplace_tf_ = 0.0;

    for (int nn = 0; nn < N_; nn++){
        for (int dd = 0; dd < D_; dd++){
            dr_p = Update_r_sum(r2_sum_old_, r_old_[nn][dd], r_old_[nn][dd] + step_);
            dr_m = Update_r_sum(r2_sum_old_, r_old_[nn][dd], r_old_[nn][dd] - step_);
            laplace_tf_  += Trial_func(alpha,dr_p) + Trial_func(alpha, dr_m);
        }
    }

    tf_middle_ = Trial_func(alpha,r2_sum_old_);
    laplace_tf_ -= 2*D_*N_*tf_middle_;
    laplace_tf_ /= (step_*step_*tf_middle_);

    return (1./2)*(-laplace_tf_ + r2_sum_old_);
}
