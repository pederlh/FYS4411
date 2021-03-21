#include "Psi.hpp"


//Function to initialize an instance of the wavefunction/hamiltonian class for non-interacting bosons.
void Psi::Declare_position(int N, int D, double h, double step, int case_type){

    D_ = D;             //Dimentions
    N_ = N;             //Number of particles
    h_ = h;             //Stepsize to determine the distributions space of particles
    step_ = step;       //Stepsize used in numerical differentiation
    case_ = case_type;  // =0 for non-interacting bosons, =1 for interacting bosons

    //Declaring position
    r_new_ = new double[D_];
    r_old_ = new double*[N_];
    for (int i= 0; i< N_ ; i++){
        r_old_[i] = new double[D_];
    }

}

//Function to initialize an instance of the wavefunction/hamiltonian class for interacting bosons.
void Psi::Declare_position_interaction(int N, int D, double h, double step, int case_type){

    D_ = D;                 //Dimentions
    N_ = N;                 //Number of particles
    h_ = h;                 //Stepsize to determine the distributions space of particles
    step_ = step;           //Stepsize used in numerical differentiation
    beta_ = 2.82843;        //Scaling for the z-contribution in the trial wave function
    case_ = case_type;      // =0 for non-interacting bosons, =1 for interacting bosons
    a_ = 0.0043;            // Hard-core diameter of bosons

    //Declaring position
    r_new_ = new double[D_];
    r_old_ = new double*[N_];
    for (int i= 0; i< N_ ; i++){
        r_old_[i] = new double[D_];
    }

}

//Function to declare quantum force, called if importance sampling is chosen.
void Psi::Declare_quantum_force(double D_diff){

    D_diff_ = D_diff;                       //Diffusion constant in Green's function and solution of Langevin eq.

    //Declaring quantum force
    quantum_force_new_ = new double[D_];
    quantum_force_old_ = new double[D_];
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
        }
    }

    //Initialize r**2 sum for interacting bosons
    if (case_ == 1){
        for (int i = 0; i < N_; i++){
            for (int j =0; j <D_; j++){
                if(j==2){
                    r2_sum_old_ += r_old_[i][j]*r_old_[i][j]*beta_;
                }
                else{
                    r2_sum_old_ += r_old_[i][j]*r_old_[i][j];
                }
            }
        }
    }
    //Initialize r**2 sum for non-interacting bosons
    else{
        for (int i = 0; i < N_; i++){
            for (int j =0; j <D_; j++){
                r2_sum_old_ += r_old_[i][j]*r_old_[i][j];
            }
        }
    }

    return r2_sum_old_;
}

//Function to initialize quantum force in the case of non-interacting bosons
//  -Called if importance sampling is choosen
void Psi::Initialize_quantum_force(double alpha, int idx){

    for (int k = 0; k < D_; k++){
        quantum_force_old_[k] = -4*alpha*r_old_[idx][k];  //Quantum force for particles before proposed move
    }

}

//Function to initialize quantum force in the case of interacting bosons
void Psi::Initialize_quantum_force_interaction(double alpha, int idx){
    double tmp =0.0;
    rkl_ = new double[N_];    //Array to hold distances between the particle proposed to move and all other particles

    for(int n = 0; n < N_; n++){
        for (int d = 0; d < D_; d++){
            tmp += pow(r_old_[idx][d] - r_old_[n][d],2);
        }
        rkl_[n] = sqrt(tmp);
        tmp = 0.0;
    }


    //Calculate quantum force in accordance to the analytical expression for QF for interacting bosons
    for (int d = 0; d<D_;d++){
        if (d == 2){
            quantum_force_old_[d] = -4*alpha*beta_*r_old_[idx][d];
        }
        else{
            quantum_force_old_[d] = -4*alpha*r_old_[idx][d];
        }
        for (int n = 0; n<N_;n++){
            if (n != idx){
                quantum_force_old_[d] += (r_old_[idx][d] - r_old_[n][d]) * (1/(rkl_[n]*rkl_[n])) * (a_/(rkl_[n] - a_));
            }
        }
        quantum_force_old_[d] = 2*quantum_force_old_[d];
    }

}

//Function to update quantum force in the case of non-interacting bosons
void Psi::Update_quantum_force(double alpha){

    for (int dim = 0; dim < D_; dim++){
        quantum_force_new_[dim] = -4*alpha*r_new_[dim];   //Quantum force for the moved particle
    }
}

//Function to update quantum force (after propsed move) in the case of interacting bosons
void Psi::Update_quantum_force_interaction(double alpha, int idx){
    double tmp =0.0;

    //Updates distances between proposed moved particle and all other particles
    for(int n = 0; n < N_; n++){
        for (int d = 0; d < D_; d++){
            tmp += pow(r_new_[d] - r_old_[n][d],2);
        }
        rkl_[n] = sqrt(tmp);
        tmp = 0.0;
    }


    //Calculate quantum force in accordance to the analytical expression for quantum force for interacting bosons
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
        quantum_force_new_[d] = 2*quantum_force_new_[d];
    }
}

//Function to update the r**2 sum in the trial wave function in the case of non-interacting bosons
double Psi::Update_r_sum(double sum, double r_init, double r_move){
    sum -= r_init*r_init;
    sum += r_move*r_move;
    return sum;
}


//Function to update the r**2 sum in the trial wave function in the case of interacting bosons
double Psi::Update_r_sum_interaction(double sum, double r_init, double r_move, double coord){
    if (coord == 2){
        sum -= r_init*r_init*beta_;
        sum += r_move*r_move*beta_;
    }
    else{
        sum -= r_init*r_init;
        sum += r_move*r_move;
    }
    return sum;
}

//Function to calculate a proposed move of one particle in the case of non-interacting bosons
// - Called if brute force metropolis sampling is chosen
double Psi::Proposed_move(int idx){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]

    for (int k = 0; k < D_; k++){
        r_new_[k] = r_old_[idx][k] + h_ * (RDG(gen) - 0.5);
        r2_sum_new_ = Update_r_sum(r2_sum_new_, r_old_[idx][k], r_new_[k]);
    }
    return r2_sum_new_;
}


//Function to calculate a proposed move of one particle in the case of non-interacting bosons using solution of Langevin eq.
// - Called if importance sampling is chosen
double Psi::Proposed_move_importance(int idx){
    mt19937_64 gen(rd_());
    normal_distribution<double> NDG(0.0,1.0);   //Random number generated from gaussian distribution with mean = 0, std = 1;
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]

    double move_step = 0.005;                //Delta t in solution of Langevin equation
    for (int k = 0; k < D_; k++){
        r_new_[k] = r_old_[idx][k] + D_diff_*quantum_force_old_[k]*move_step + NDG(gen)*sqrt(move_step);
        r2_sum_new_ = Update_r_sum(r2_sum_new_, r_old_[idx][k], r_new_[k]);
    }
    return r2_sum_new_;
}


//Function to calculate a proposed move of one particle in the case of interacting bosons using solution of Langevin eq.
double Psi::Proposed_move_interaction(int idx){
    mt19937_64 gen(rd_());
    normal_distribution<double> NDG(0.0,1.0);   //Random number generated from gaussian distribution with mean = 0, std = 1;
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]

    double move_step = 0.005;               //Delta t in solution of Langevin equation
    for (int k = 0; k < D_; k++){
        r_new_[k] = r_old_[idx][k] + D_diff_*quantum_force_old_[k]*move_step + NDG(gen)*sqrt(move_step);
        r2_sum_new_ = Update_r_sum_interaction(r2_sum_new_, r_old_[idx][k], r_new_[k], k);
    }
    return r2_sum_new_;

}


//Function to evalueate the trial wave function in the case of non-interacting bosons
double Psi::Trial_func(double alpha, double sum_r_squared){
    return exp(-alpha*sum_r_squared);
}

//Function to evaluate the trial wave function in the case of interacting bosons
double Psi::Trial_func_interaction(double alpha, double sum_r_squared, string version, int idx = 0){
    double exp_prod = -alpha*sum_r_squared;    //First term in trial WF
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

    //Evaluate second term i trial WF
    for (int l = 0; l < N_; l++){
        for (int m = l+1; m < N_; m++){
            diff= 0.0;
            for (int dim = 0; dim <D_; dim ++){
                diff += pow(r_copy[l][dim]-r_copy[m][dim],2);
            }
            diff = sqrt(diff);
            if (diff <= a_){                        //Checks if particles lie closer than hard-core diameter
                return 0.0;                         //If so -> WF collapses
            }
            else{
                u += log(1-(a_/diff));
            }
        }
    }

    return exp(exp_prod+u);
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


//Function to calculate the laplacian of the phi term of the trial wave function for interacting bosons
double Psi::Laplace_phi(int idx, double d2_phi, double alpha){
    d2_phi = 0.0;
    for (int d = 0; d< D_; d++){
        if (d ==2){
            d2_phi += beta_*beta_*r_old_[idx][d]*r_old_[idx][d];
        }
        else{
            d2_phi += r_old_[idx][d]*r_old_[idx][d];
        }
    }

    d2_phi *= 4*alpha*alpha;
    d2_phi -= 2*alpha*(2+beta_);

    return d2_phi;

}


//Function to evaluate local energy using the analytical expression in the case of interacting bosons
double Psi::Local_energy_interaction(double alpha){

    double* nabla_phi = new double[D_];
    double* distance_vec = new double[D_];


    double d2_phi, term2, term3, term4, u_der, d2_psi, V_ext, E_L, V_int, tmp;
    d2_psi = 0.0;
    V_ext = 0.0;


    for (int n = 0; n < N_ ;n++){
        d2_phi = Laplace_phi(n, d2_phi, alpha);
        term2 = 0.0;
        term3 = 0.0;
        term4 = 0.0;
        for (int d = 0; d < D_; d++){
            distance_vec[d] = 0.0;
            nabla_phi[d] = -2*alpha*r_old_[n][d];
        }
        for(int n2 = 0; n2 < N_; n2++){
            tmp = 0.0;
            for (int d2 = 0; d2 < D_; d2++){
                tmp += pow(r_old_[n][d2] - r_old_[n2][d2],2);
            }
            rkl_[n2] = sqrt(tmp);
        }
        for (int n2 = 0; n2<N_;n2++){
            if (n2 != n){
                u_der = (1/(rkl_[n2]*rkl_[n2])) * (a_/(rkl_[n2] - a_));
                term4 += (a_*a_ - 2*a_*rkl_[n2])/(rkl_[n2]*rkl_[n2]*(rkl_[n2]-a_)*(rkl_[n2]-a_)) + 2*u_der;
                for (int d = 0; d<D_;d++){
                    distance_vec[d] += (r_old_[n][d] - r_old_[n2][d]) * u_der;
                }
            }
        }

        for (int d =0; d<D_; d++){
            term2 += 2*nabla_phi[d]*distance_vec[d];
            term3 += distance_vec[d]*distance_vec[d];

            if(d==2) {V_ext += beta_*beta_*r_old_[n][d]*r_old_[n][d];}
            else {V_ext += r_old_[n][d]*r_old_[n][d];}

        }

        d2_psi += d2_phi + term2 + term3 + term4;
    }

    V_int = 0.0;

    //double diff;
    /*
    for (int i = 0; i < N_; i++){
        for (int j = i+1; j <N_; j++){
            diff = 0.0;
            for (int d = 0; d < D_; d++){
                diff += pow(r_old_[i][d]-r_old_[j][d],2);
            }
            diff = sqrt(diff);
            if (diff <= a_){
                V_int += 1e16;
                cout <<"Enormt bidrag" <<endl;
            }
        }
    }
    */

    E_L = (1./2)*(-d2_psi+V_ext) + V_int;

    return E_L;
}
