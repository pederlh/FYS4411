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
    quantum_force_old_ = new double*[N_];
    for (int i = 0; i < N_ ; i++){
        quantum_force_old_[i] = new double[D_];
    }

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




void Psi::Initialize_quantum_force(double alpha){
    for (int j = 0; j < N_; j++){
        for (int k = 0; k < D_; k++){
            quantum_force_old_[j][k] = -4*alpha*r_old_[j][k];
        }
    }
}

void Psi::Initialize_quantum_force_interaction(double alpha){
    double ** distances = new double*[N_];
    for (int j =0; j <N_; j++){
        distances[j] = new double[D_];
        for (int dim = 0; dim <D_; dim++){
            distances[j][dim] = 0.0;
        }
    }

    double **rkl = new double*[N_];  //MÅ FYLLES MED 0
    for (int z = 0; z < N_; z++){
        rkl[z] = new double[N_];
        for (int p =0; p < N_; p++){
            rkl[z][p] = 0.0;
        }
    }

    for(int j = 0; j < N_; j++){
        for(int t = 0; t< N_; t++){
            if(t != j){
                for (int d = 0; d < D_; d++){
                    rkl[j][t] += pow(r_old_[j][d] - r_old_[t][d],2);
                }
            rkl[j][t] = sqrt(rkl[j][t]);
            }
        }
    }

    for (int k = 0; k < N_;k++){
        for (int dim= 0; dim<D_;dim++){
            for (int l = 0; l<N_ ; l++){
                if (l!=k){
                    distances[k][dim] += ((r_old_[k][dim]-r_old_[l][dim])/rkl[k][l])*(a_/(rkl[k][l]*(rkl[k][l]-a_)));
                }
            }
            distances[j][dim] = distances[j][dim]*2;
        }
    }

    for (int j = 0; j < N_ ; j++){
        for (int k = 0; k <D_; k++){
            if ( k==2){
                quantum_force_old_[j][k] = -4*alpha*r_old_[j][k]*beta_ + distances[j][k];
            }
            else{
                quantum_force_old_[j][k] = -4*alpha*r_old_[j][k] + distances[j][k];
            }
        }
    }

}


void Psi::Update_quantum_force(double alpha){
    for (int dim = 0; dim < D_; dim++){
        quantum_force_new_[dim] = -4*alpha*r_new_[dim];
    }
}


void Psi::Update_quantum_force_interaction(double alpha){
    for (int dim = 0; dim < D_; dim++){
        quantum_force_new_[dim] = -4*alpha*r_new_[dim];
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
        r_new_[k] = r_old_[idx][k] + D_diff_*quantum_force_old_[idx][k]*move_step + NDG(gen)*sqrt(move_step);
        //r_new_[k] = r_old_[idx][k] + h_ * (RDG(gen) - 0.5);
        r2_sum_new_ = Update_r_sum(r2_sum_new_, r_old_[idx][k], r_new_[k]);
    }
    return r2_sum_new_;
}


double Psi::Trial_func(double alpha, double sum_r_squared){
    return exp(-alpha*sum_r_squared);
}

double Psi::Trial_func_interaction(double alpha, double sum_r_squared){
    double exp_prod = -alpha*sum_r_squared;
    double u = 0;
    double diff;

    for (int l = 0; l < N_; l++){
        for (int m = l+1; m < N_; m++){
            for (int dim = 0; dim <D_; dim ++){
                diff += pow(r_old_[l][dim]-r_old_[m][dim],2);
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

    return exp_prod*u;
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
            laplace_tf_  += Trial_func_interaction(alpha,dr_p) + Trial_func_interaction(alpha, dr_m);
        }
    }

    tf_middle_ = Trial_func_interaction(alpha,r2_sum_old_);
    laplace_tf_ -= 2*D_*N_*tf_middle_;
    laplace_tf_ /= (step_*step_*tf_middle_);

    return (1./2)*(-laplace_tf_ + r2_sum_old_);
}
