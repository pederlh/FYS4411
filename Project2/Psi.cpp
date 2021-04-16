#include "Psi.hpp"


//Function to initialize an instance of the wavefunction/hamiltonian class for non-interacting bosons.
void Psi::Declare_position(int N, int D, double h, double step, int case_type){

    D_ = D;                 //Dimentions
    N_ = N;                 //Number of particles
    h_ = h;                 //Stepsize to determine the distributions space of particles
    step_ = step;           //Stepsize used in numerical differentiation
    case_ = case_type;      // =0 for non-interacting bosons, =1 for interacting bosons
    move_step_ = 0.005;    //Delta t in solution of Langevin equation

    //Declaring position
    r_new  = mat(N_*D_).fill(0.);
    r_old_ = mat(N_*D_).fill(0.);

}


//Function to declare quantum force, called if importance sampling is chosen.
void Psi::Declare_quantum_force(double D_diff){

    D_diff_ = D_diff;                       //Diffusion constant in Green's function and solution of Langevin eq.

    //Declaring quantum force
    quantum_force_new_ = vec(D_).fill(0.);
    quantum_force_old_ = vec(D_).fill(0.);
}


//Function to intialize matrix of particles and sum of r**2 for the trial wave function
double Psi::Initialize_positions(){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]

    r2_sum_old_ = 0;    //Variable to hold the r**2 sum of the old positions for use in the trial wave function

     //Initialize posistions
    for (int j = 0; j < N_; j++){
        for (int k = 0; k < D_; k++){
            r_old_[j][k] = h_ * (RDG(gen) - 0.5);
            r_new_[j][k] = r_old_[j][k];
        }
    }

    //Initialize r**2 sum for non-interacting bosons
    for (int i = 0; i < N_; i++){
        for (int j =0; j <D_; j++){
            r2_sum_old_ += r_old_[i][j]*r_old_[i][j];
            r2_sum_new_ += r_old_[i][j]*r_old_[i][j];
        }
    }

    return r2_sum_old_;
}

//Function to initialize quantum force in the case of non-interacting bosons
//  -Called if importance sampling is choosen
void Psi::Initialize_quantum_force(double alpha, int idx){

    for (int d = 0; d < D_; d++){
        quantum_force_old_[d] = -4*alpha*r_old_[idx][d];  //Quantum force for particles before proposed move
    }

}


//Function to update quantum force in the case of non-interacting bosons
void Psi::Update_quantum_force(double alpha, int idx){

    for (int d = 0; d < D_; d++){
        quantum_force_new_[d] = -4*alpha*r_new_[idx][d];   //Quantum force for the moved particle
    }
}


//Function to update the r**2 sum in the trial wave function in the case of non-interacting bosons
double Psi::Update_r_sum(double sum, double r_init, double r_move){
    sum -= r_init*r_init;
    sum += r_move*r_move;
    return sum;
}


//Function to calculate a proposed move of one particle in the case of non-interacting bosons
// - Called if brute force metropolis sampling is chosen
double Psi::Proposed_move(int idx){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]

    for (int k = 0; k < D_; k++){
        r_new_[k] = r_old_[idx][k] + h_ * (RDG(gen) - 0.5);
        r2_sum_new_ = Update_r_sum(r2_sum_new_, r_old_[idx][k], r_new_[idx][k]);
    }
    return r2_sum_new_;
}


//Function to calculate a proposed move of one particle in the case of non-interacting bosons using solution of Langevin eq.
// - Called if importance sampling is chosen
double Psi::Proposed_move_importance(int idx){
    mt19937_64 gen(rd_());
    normal_distribution<double> NDG(0.0,1.0);   //Random number generated from gaussian distribution with mean = 0, std = 1;
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]

    for (int k = 0; k < D_; k++){
        r_new_[k] = r_old_[idx][k] + D_diff_*quantum_force_old_[k]*move_step_ + NDG(gen)*sqrt(move_step_);
        r2_sum_new_ = Update_r_sum(r2_sum_new_, r_old_[idx][k], r_new_[idx][k]);
    }
    return r2_sum_new_;
}



//Function to evalueate the trial wave function in the case of non-interacting bosons
double Psi::Trial_func(double alpha, double sum_r_squared){
    return exp(-alpha*sum_r_squared);
}


//Function to evaluate local energy using the analytical expression in the case of non-interacting bosons
double Psi::Local_energy_analytical(double alpha){
    return D_*N_*alpha + (1-4*alpha*alpha)*(1./2)*r2_sum_old_;
}

//Function to evaluate local energy using numerical differentiation in the case of non-interacting bosons
double Psi::Local_energy_brute_force(double alpha){
    double dr_p, dr_m;
    laplace_tf_ = 0.0;

    for (int nn = 0; nn < N_; nn++){
        for (int dd = 0; dd < D_; dd++){
            dr_p = Update_r_sum(r2_sum_old_, r_old_[nn][dd], r_old_[nn][dd] + step_);   //Position + delta r
            dr_m = Update_r_sum(r2_sum_old_, r_old_[nn][dd], r_old_[nn][dd] - step_);   //Position - delta r
            laplace_tf_  += Trial_func(alpha,dr_p) + Trial_func(alpha, dr_m);
        }
    }

    tf_middle_ = Trial_func(alpha,r2_sum_old_);
    laplace_tf_ -= 2*D_*N_*tf_middle_;
    laplace_tf_ /= (step_*step_*tf_middle_);

    return (1./2)*(-laplace_tf_ + r2_sum_old_);
}
