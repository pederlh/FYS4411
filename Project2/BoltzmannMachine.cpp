#include "BoltzmannMachine.hpp"

BoltzmannMachine::BoltzmannMachine(int num_particles,int dimentions, double eta, int MC)
{
    mt19937_64 gen(rd_());
    normal_distribution<double> NDG(0.0,0.001);   //Random number generated from gaussian distribution with mean = 0, std = 0.001;
    D_ = dimentions;
    N_ = num_particles;
    H_ = 4;
    MC_ = MC;
    double std_norm_dist = 0.1;  //standard deviation of normal distribution used to initialize weights/biases.

   //fill::randn = set each element to a random value from a normal/Gaussian distribution with zero mean and unit variance
   //Initializing weihgts/biases
   w_ = cube(N_, D_, H_, fill::randn)*std_norm_dist;
   a_ = mat(N_, D_, fill::randn)*std_norm_dist;
   b_ = vec(H_, fill::randn)*std_norm_dist;
   Q_ = vec(H_).fill(0.0); //HUsk å sette denne til 0 i funksjonen sin.
   sigma_ = 1;
   eta_ = eta; //Learning rate SGD.

   //Derivatives of the wave function w/respect to weights/biases for stochastic gradient descent
   dw_ = cube(N_, D_, H_).fill(0.0);
   da_ = mat(N_, D_).fill(0.0);
   db_ = vec(H_).fill(0.0);

   E_dw_ = cube(N_, D_, H_).fill(0.0);
   E_da_ = mat(N_, D_).fill(0.0);
   E_db_ = vec(H_).fill(0.0);

   ADAM();
}

void BoltzmannMachine::Initialize()
{
    //Declaring position
    r_old_ = mat(N_,D_, fill::randu) - 0.5;
    r_new_ = r_old_;
}

double BoltzmannMachine::WaveFunction(mat r)
{
    Q_factor(r);
    double term1 = 0.0;
    double term2 = 1.0;

    for (int n = 0; n < N_; n++){
        for (int d = 0 ; d < D_; d++){
            term1 += pow((r(n,d)-a_(n,d)),2);
        }
    }
    term1 = exp(-term1/(2*sigma_*sigma_));

    for (int h = 0; h < H_; h++){
        term2 *= (1.0 + exp(Q_(h)));
    }

    return term1*term2;

}

void BoltzmannMachine::Q_factor(mat r)
{
    vec temp = vec(H_).fill(0.0);
    Q_.fill(0.0);

    for (int h = 0; h < H_; h++){
        mat W = w_.slice(h);
        temp(h) = accu(r%W);    //sum of all elements in matrix defined by elementwise product of two other matrices
    }
    Q_ = b_ + temp;
}


double BoltzmannMachine::MonteCarlo()
{
    double energy,DeltaE;
    Initialize();

    cube delta_Psi_w = cube(N_, D_, H_).fill(0.0);
    mat delta_Psi_a = mat(N_, D_).fill(0.0);
    vec delta_Psi_b = vec(H_).fill(0.0);

    cube derivative_Psi_w = cube(N_, D_, H_).fill(0.0);
    mat derivative_Psi_a = mat(N_, D_).fill(0.0);
    vec derivative_Psi_b = vec(H_).fill(0.0);

    energy = 0;

    for (int cycle = 0; cycle < MC_; cycle++){
        for (int n = 0; n < N_; n++){

            Metropolis();
            DeltaE = LocalEnergy();
            Derivate_wavefunction();

            delta_Psi_w += dw_;
            delta_Psi_a += da_;
            delta_Psi_b += db_;

            energy += DeltaE;

            derivative_Psi_w += dw_*DeltaE;
            derivative_Psi_a += da_*DeltaE;
            derivative_Psi_b += db_*DeltaE;
        }
    }
    energy/=MC_*N_;
    derivative_Psi_w /= MC_*N_;
    derivative_Psi_a /= MC_*N_;
    derivative_Psi_b /= MC_*N_;
    delta_Psi_w /= MC_*N_;
    delta_Psi_a /= MC_*N_;
    delta_Psi_b /= MC_*N_;
    E_dw_ = 2*(derivative_Psi_w - delta_Psi_w*energy);
    E_da_ = 2*(derivative_Psi_a - delta_Psi_a*energy);
    E_db_ = 2*(derivative_Psi_b - delta_Psi_b*energy);

    return energy;
}

void BoltzmannMachine::Derivate_wavefunction()
{
    Q_factor(r_old_);
    da_ = (r_old_ - a_)/pow(sigma_,2);
    db_ = 1.0/(1.0 + exp(-Q_));
    for (int h = 0; h< H_; h++){
        dw_.slice(h) = dw_.slice(h)/(pow(sigma_,2)*(1.0 +  exp(-Q_(h))));
    }
}

void BoltzmannMachine::Metropolis(){

    double tf_old, tf_new, P;
    int idx = randi<int>(distr_param(0,N_-1));

    //Moves particle with index "idx"
    for (int d = 0; d < D_; d++){
        r_new_(idx,d) = r_old_(idx,d) + (randu() - 0.5);
    }

    tf_old = WaveFunction(r_old_);             //Trial wave function of old position
    tf_new = WaveFunction(r_new_);           //Trial wave function of new position
    P = (tf_new*tf_new)/(tf_old*tf_old);                         //Metropolis test
    if (randu() <= P){
        for (int d =0 ; d < D_;  d++){                   //Update initial position
            r_old_(idx,d) = r_new_(idx,d);
        }
    }
    else{
        for (int d =0 ; d < D_;  d++){
            r_new_(idx,d) = r_old_(idx,d);
        }
    }
}

double BoltzmannMachine::LocalEnergy()
{
    double energy = 0;
    Q_factor(r_old_);
    double der1_ln_psi, der2_ln_psi;

    for (int n = 0; n < N_; n++){
        for (int d =0; d < D_; d++){
            double sum1 =0;
            double sum2 =0;
            for (int h = 0; h < H_; h++){
                sum1 += w_(n,d,h)/(1.0+exp(-Q_(h)));
                sum2 += pow(w_(n,d,h),2)*exp(Q_(h))/pow((1.0 + exp(Q_(h))),2);
            }
            der1_ln_psi = -(r_old_(n,d) - a_(n,d))/pow(sigma_,2) + sum1/pow(sigma_,2);
            der2_ln_psi = -1.0/pow(sigma_,2) + sum2/(pow(sigma_,2));
            energy += 0.5*(-pow(der1_ln_psi,2) - der2_ln_psi + pow(r_old_(n,d),2));
        }
    }
    return energy;
}

/* Method for Adam optimization */
void BoltzmannMachine::ADAM()
{
    double epsilon = 1e-8;   //Value to avoid division by zero.
    double beta1 = 0.9;      //Decay rate
    double beta2 = 0.999;   //Decay rate
    double alpha_it, epsilon_it;

    //First momentum of weights/biases for stochastic gradient descent
    mom_w_ = cube(N_, D_, H_).fill(0.0);
    mom_a_ = mat(N_, D_).fill(0.0);
    mom_b_ = vec(H_).fill(0.0);

    //Second momentum of weights/biases for stochastic gradient descent
    second_mom_w_ = cube(N_, D_, H_).fill(0.0);
    second_mom_a_ = mat(N_, D_).fill(0.0);
    second_mom_b_ = vec(H_).fill(0.0);

    double Energy = 0;
    int its = 100;
    vec Energies = vec(its).fill(0.0);

    for (int i = 0; i < its;  i++){
        Energy = MonteCarlo();
        Energies(i) = Energy;
        cout << "Energy = " << Energy << endl;
        E_da_*=eta_;
        E_db_*=eta_;
        E_dw_*=eta_;

        mom_w_ = beta1*mom_w_ + (1-beta1)*E_dw_;
        mom_b_ = beta1*mom_b_ + (1-beta1)*E_db_;
        mom_a_ = beta1*mom_a_ + (1-beta1)*E_da_;

        second_mom_w_ = beta2*second_mom_w_ + (1-beta2)*(E_dw_%E_dw_);
        second_mom_b_ = beta2*second_mom_b_ + (1-beta2)*(E_db_%E_db_);
        second_mom_a_ = beta2*second_mom_a_ + (1-beta2)*(E_da_%E_da_);

        alpha_it = eta_*sqrt(1-pow(beta2, i+1))/(1-pow(beta1, i+1));
        epsilon_it = epsilon*sqrt(1-pow(beta2, i+1));

        w_ -= mom_w_*alpha_it/(sqrt(second_mom_w_) + epsilon_it);
        b_ -= mom_b_*alpha_it/(sqrt(second_mom_b_) + epsilon_it);
        a_ -= mom_a_*alpha_it/(sqrt(second_mom_a_) + epsilon_it);

    }
}


void BoltzmannMachine::SGD()
{
    double Energy = 0;
    int its = 100;
    vec Energies = vec(its).fill(0.0);

    for (int i = 0; i < its;  i++){
        Energy = MonteCarlo();
        Energies(i) = Energy;
        cout << "Energy = " << Energy << endl;
        a_ -= eta_*E_da_;
        b_ -= eta_*E_db_;
        w_ -= eta_*E_dw_;
    }
}
