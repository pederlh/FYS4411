#include "Psi.hpp"


void Psi::Declare_position(int N, int D, double h, double step, int case_type){

    D_ = D;
    N_ = N;
    h_ = h;
    step_ = step;
    case_ = case_type;

    //Declaring position
    r_new_ = new double[D_];
    r_old_ = new double*[N_];
    for (int i= 0; i< N_ ; i++){
        r_old_[i] = new double[D_];
    }

}

void Psi::Declare_position_interaction(int N, int D, double h, double step, int case_type, double beta){

    D_ = D;
    N_ = N;
    h_ = h;
    step_ = step;
    beta_ = beta;
    case_ = case_type;
    a_ = 0.0043;

    //Declaring position
    r_new_ = new double[D_];
    r_old_ = new double*[N_];
    for (int i= 0; i< N_ ; i++){
        r_old_[i] = new double[D_];
    }

}

void Psi::Declare_quantum_force(double D_diff){

    D_diff_ = D_diff;
    quantum_force_new_ = new double[D_];
    quantum_force_old_ = new double[D_];

}

double Psi::Initialize_positions(){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]

    r2_sum_old_ = 0;

    for (int j = 0; j < N_; j++){                       //Initial posistion
        for (int k = 0; k < D_; k++){
            r_old_[j][k] = h_ * (RDG(gen) - 0.5);
        }
    }

    if (case_ == 1){
        for (int i = 0; i < N_; i++){
            for (int j =0; j <D_; j++){
                if(j==2){
                    r2_sum_old_ += r_old_[i][j]*r_old_[i][j]*beta_*beta_;
                }
                else{
                    r2_sum_old_ += r_old_[i][j]*r_old_[i][j];
                }
            }
        }
    }
    else{
        for (int i = 0; i < N_; i++){
            for (int j =0; j <D_; j++){
                r2_sum_old_ += r_old_[i][j]*r_old_[i][j];
            }
        }
    }

    return r2_sum_old_;
}


void Psi::Initialize_quantum_force(double alpha, int idx){

    for (int k = 0; k < D_; k++){
        quantum_force_old_[k] = -4*alpha*r_old_[idx][k];
    }

}


void Psi::Initialize_quantum_force_interaction(double alpha, int idx){
    double tmp =0.0;
    rkl_ = new double[N_];

    for(int n = 0; n < N_; n++){
        if (n != idx){
            for (int d = 0; d < D_; d++){
                tmp += pow(r_old_[idx][d] - r_old_[n][d],2);
            }
            rkl_[n] = sqrt(tmp);
            tmp = 0.0;
        }
        else{
            rkl_[n] = 0.0;
        }
    }


    for (int d = 0; d<D_;d++){
        if (d == 2){
            quantum_force_old_[d] = -4*alpha*beta_*r_old_[idx][d];
        }
        else{
            quantum_force_old_[d] = -4*alpha*r_old_[idx][d];
        }
        for (int n = 0; n<N_;n++){
            if (n != idx){
                quantum_force_old_[d] += (r_old_[idx][d]-r_old_[n][d])*(1/(rkl_[n]*rkl_[n])*(a_/(rkl_[n]-a_)));
            }
        }
    }

}


void Psi::Update_quantum_force(double alpha){

    for (int dim = 0; dim < D_; dim++){
        quantum_force_new_[dim] = -4*alpha*r_new_[dim];
    }
}


void Psi::Update_quantum_force_interaction(double alpha, int idx){
    double tmp =0.0;

    for(int n = 0; n < N_; n++){
        if (n != idx){
            for (int d = 0; d < D_; d++){
                tmp += pow(r_new_[d] - r_old_[n][d],2);
            }
            rkl_[n] = sqrt(tmp);
            tmp = 0.0;
        }
        else{
            rkl_[n] = 0.0;
        }
    }


    for (int d = 0; d<D_;d++){
        if (d == 2){
            quantum_force_new_[d] = -4*alpha*beta_*r_new_[d];
        }
        else{
            quantum_force_new_[d] = -4*alpha*r_new_[d];
        }
        for (int n = 0; n<N_;n++){
            if (n != idx){
                quantum_force_new_[d] += (r_new_[d]-r_old_[n][d])*(1/(rkl_[n]*rkl_[n])*(a_/(rkl_[n]-a_)));
            }
        }
    }
}

double Psi::Update_r_sum(double sum, double r_init, double r_move){
    sum -= r_init*r_init;
    sum += r_move*r_move;
    return sum;
}

double Psi::Update_r_sum_interaction(double sum, double r_init, double r_move, double coord){
    if (coord == 2){
        sum -= r_init*r_init*beta_*beta_;
        sum += r_move*r_move*beta_*beta_;
    }
    else{
        sum -= r_init*r_init;
        sum += r_move*r_move;
    }
    return sum;
}


double Psi::Proposed_move(int idx){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]

    for (int k = 0; k < D_; k++){
        r_new_[k] = r_old_[idx][k] + h_ * (RDG(gen) - 0.5);
        r2_sum_new_ = Update_r_sum(r2_sum_new_, r_old_[idx][k], r_new_[k]);
    }
    return r2_sum_new_;
}

double Psi::Proposed_move_importance(int idx){
    mt19937_64 gen(rd_());
    normal_distribution<double> NDG(0.0,1.0);   //Random number generated from gaussian distribution with mean = 0, std = 1;
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]

    double move_step = 0.005;
    for (int k = 0; k < D_; k++){
        r_new_[k] = r_old_[idx][k] + D_diff_*quantum_force_old_[k]*move_step + NDG(gen)*sqrt(move_step);
        r2_sum_new_ = Update_r_sum(r2_sum_new_, r_old_[idx][k], r_new_[k]);
    }
    return r2_sum_new_;
}

double Psi::Proposed_move_interaction(int idx){
    mt19937_64 gen(rd_());
    normal_distribution<double> NDG(0.0,1.0);   //Random number generated from gaussian distribution with mean = 0, std = 1;
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]

    double move_step = 0.005;
    for (int k = 0; k < D_; k++){
        r_new_[k] = r_old_[idx][k] + D_diff_*quantum_force_old_[k]*move_step + NDG(gen)*sqrt(move_step);
        r2_sum_new_ = Update_r_sum_interaction(r2_sum_new_, r_old_[idx][k], r_new_[k],k);
    }
    return r2_sum_new_;

}


double Psi::Trial_func(double alpha, double sum_r_squared){
    return exp(-alpha*sum_r_squared);
}

double Psi::Trial_func_interaction(double alpha, double sum_r_squared, string version, int idx){
    double exp_prod = -alpha*sum_r_squared;
    double u = 0;
    double diff;

    double **r_copy = new double*[N_];
    for (int i = 0; i< N_; i++){
        r_copy[i] = new double[D_];
        for (int j = 0; j<D_; j++){
            r_copy[i][j] = r_old_[i][j];

        }
    }

    if (version == "new"){
        for (int d = 0; d < D_;d++){
            r_copy[idx][d] = r_new_[d];
        }
    }

    for (int l = 0; l < N_; l++){
        for (int m = l+1; m < N_; m++){
            for (int dim = 0; dim <D_; dim ++){
                diff += pow(r_copy[l][dim]-r_copy[m][dim],2);
            }
            diff = sqrt(diff);
            if (diff <= a_){
                u += 0;
            }
            else{
                u += 1-(a_/diff);
            }
        }
    }

    return exp(exp_prod+u);
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


double Psi::Local_energy_interaction_analytical(double alpha){
    double d2_phi = 4*alpha*alpha;

    return 0.0;
}

double Psi::Local_energy_interaction_brute_force(double alpha){
    double dr_p, dr_m;
    laplace_tf_ = 0.0;

    for (int nn = 0; nn < N_; nn++){
        for (int dd = 0; dd < D_; dd++){
            dr_p = Update_r_sum_interaction(r2_sum_old_, r_old_[nn][dd], r_old_[nn][dd] + step_, dd);
            dr_m = Update_r_sum_interaction(r2_sum_old_, r_old_[nn][dd], r_old_[nn][dd] - step_, dd);
            laplace_tf_  += Trial_func_interaction(alpha,dr_p,"old",0) + Trial_func_interaction(alpha, dr_m,"old",0);
        }
    }

    tf_middle_ = Trial_func_interaction(alpha,r2_sum_old_,"old",0);
    laplace_tf_ -= 2*D_*N_*tf_middle_;
    laplace_tf_ /= (step_*step_*tf_middle_);

    return (1./2)*(-laplace_tf_ + r2_sum_old_);
}
