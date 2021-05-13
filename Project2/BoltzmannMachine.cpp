#include "BoltzmannMachine.hpp"

BoltzmannMachine::BoltzmannMachine(int num_particles,int dimentions, double eta, int MC, int type_sampling, int interaction)
{
    arma_rng::set_seed_random();
    D_ = dimentions;
    N_ = num_particles;
    H_ = 2;
    MC_ = MC;
    double std_norm_dist = 0.001;  //standard deviation of normal distribution used to initialize weights/biases.

   //fill::randn = set each element to a random value from a normal/Gaussian distribution with zero mean and unit variance
   //Initializing weihgts/biases
   w_ = cube(N_, D_, H_, fill::randn)*std_norm_dist;
   a_ = mat(N_, D_, fill::randn)*std_norm_dist;
   b_ = vec(H_, fill::randn)*std_norm_dist;
   Q_ = vec(H_).fill(0.0);
   sigma_ = 1.0;
   sigma2_ = sigma_*sigma_;
   eta_ = eta; //Learning rate SGD.
   t_step_ = 0.05;
   D_diff_ = 0.5;
   interaction_ = interaction;

   //Initializing position
   r_old_ = mat(N_,D_, fill::randn)*sqrt(t_step_);
   r_new_ = mat(N_,D_).fill(0.0);
   for (int n = 0; n < N_ ; n++){
       for (int d = 0; d < D_; d++){
           r_new_(n,d) = r_old_(n,d);
       }
   }


   if (type_sampling == 0)
   {
       MetropolisMethod = &BoltzmannMachine::Metropolis;
   }

   if (type_sampling == 1)
   {
       MetropolisMethod = &BoltzmannMachine::Metropolis_Hastings;
       //Parameters for metropolis-hastings algo
       quantum_force_ = mat(N_, D_).fill(0.0);
       quantum_force_old_ = QuantumForce(r_old_);
       quantum_force_new_ = QuantumForce(r_old_);
   }

   SGD();
}



double BoltzmannMachine::MonteCarlo()
{
    double energy = 0.0;
    double DeltaE = 0.0;


    //Derivatives of the wave function w/respect to weights/biases for stochastic gradient descent
    dw_ = cube(N_, D_, H_).fill(0.0);
    da_ = mat(N_, D_).fill(0.0);
    db_ = vec(H_).fill(0.0);

    E_dw_ = cube(N_, D_, H_).fill(0.0);
    E_da_ = mat(N_, D_).fill(0.0);
    E_db_ = vec(H_).fill(0.0);

    cube delta_Psi_w = cube(N_, D_, H_).fill(0.0);
    mat delta_Psi_a = mat(N_, D_).fill(0.0);
    vec delta_Psi_b = vec(H_).fill(0.0);

    cube derivative_Psi_w = cube(N_, D_, H_).fill(0.0);
    mat derivative_Psi_a = mat(N_, D_).fill(0.0);
    vec derivative_Psi_b = vec(H_).fill(0.0);

    for (int cycle = 0; cycle < MC_; cycle++){
        for (int out_n = 0; out_n < N_; out_n++){
            (this->*MetropolisMethod)();
        }

        DeltaE = LocalEnergy();
        //cout << DeltaE << endl;
        Derivate_wavefunction();
        delta_Psi_w += dw_;
        delta_Psi_a += da_;
        delta_Psi_b += db_;

        energy += DeltaE;

        derivative_Psi_w += dw_*DeltaE;
        derivative_Psi_a += da_*DeltaE;
        derivative_Psi_b += db_*DeltaE;
    }

    energy/=MC_;
    derivative_Psi_w /= MC_;
    derivative_Psi_a /= MC_;
    derivative_Psi_b /= MC_;
    delta_Psi_w /= MC_;
    delta_Psi_a /= MC_;
    delta_Psi_b /= MC_;
    E_dw_ = 2*(derivative_Psi_w - energy*delta_Psi_w);
    E_da_ = 2*(derivative_Psi_a - energy*delta_Psi_a);
    E_db_ = 2*(derivative_Psi_b - energy*delta_Psi_b);

    return energy;
}

void BoltzmannMachine::Derivate_wavefunction()
{
    Q_factor(r_old_);
    da_ = (r_old_ - a_)/sigma2_;
    db_ = 1.0/(1.0 + exp(-Q_));
    for (int h = 0; h < H_; h++){
        dw_.slice(h) = w_.slice(h)/(sigma2_*(1.0 +  exp(-Q_(h))));
    }
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
    term1 = exp(-term1/(2*sigma2_));

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
        temp(h) = accu(r%w_.slice(h));    //sum of all elements in matrix defined by elementwise product of two other matrices
    }
    Q_ = b_ + temp;
}

void BoltzmannMachine::Metropolis()
{
    int idx = randi<int>(distr_param(0,N_-1));

    //Moves particle with index "idx"
    for (int d = 0; d < D_; d++){
        r_new_(idx,d) = r_old_(idx,d) + (randu() - 0.5);
    }

    tf_old_ = WaveFunction(r_old_);             //Trial wave function of old position
    tf_new_ = WaveFunction(r_new_);           //Trial wave function of new position
    P_ = (tf_new_*tf_new_)/(tf_old_*tf_old_);                         //Metropolis test
    if (randu() <= P_){
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

void BoltzmannMachine::Metropolis_Hastings()
{
    int idx = randi<int>(distr_param(0,N_-1));    //Random integer generated from uniform distribution [0,N-1];
    double rand_gauss = randn<double>();         //Random double generated from gaussian distribution with mean = 0, std = 1;
    double GreensFunc;

    //Proposed move of particle
    for (int d = 0; d < D_; d++){
        r_new_(idx,d) = r_old_(idx,d) + D_diff_*quantum_force_old_(idx,d)*t_step_ + rand_gauss*sqrt(t_step_);
    }
    quantum_force_new_ = QuantumForce(r_new_);

    GreensFunc = GreensFunction(idx);
    tf_old_ = WaveFunction(r_old_);             //Trial wave function of old position
    tf_new_ = WaveFunction(r_new_);           //Trial wave function of new position
    P_ = GreensFunc*(tf_new_*tf_new_)/(tf_old_*tf_old_);                         //Metropolis test
    if (randu() <= P_){
        for (int d = 0; d < D_;  d++){                   //Update initial position
            r_old_(idx,d) = r_new_(idx,d);
            quantum_force_old_(idx,d) = quantum_force_new_(idx,d);
        }
    }
    else{
        for (int d = 0; d < D_;  d++){
            r_new_(idx,d) = r_old_(idx,d);
            quantum_force_new_(idx,d) = quantum_force_old_(idx,d);
        }
    }
}

mat BoltzmannMachine::QuantumForce(mat r)
{
    mat temp = mat(N_,D_).fill(0.0);
    Q_factor(r);
     for(int h = 0; h < H_; h++){
         temp += w_.slice(h)/(1.0 + exp(-Q_(h)));
     }

     quantum_force_ = 2.0*((-(r-a_)/sigma2_) + (temp/sigma2_));
     return quantum_force_;
}


double BoltzmannMachine::GreensFunction(int idx)
{
    double G = 0.0;
    for (int d = 0; d < D_; d++){
        G += 0.5*(quantum_force_old_(idx,d) + quantum_force_new_(idx,d))*(D_diff_*t_step_*0.5*(quantum_force_old_(idx,d)-quantum_force_new_(idx,d)) - r_new_(idx,d) + r_old_(idx,d));
    }

    G = exp(G);
    return G;
}

double BoltzmannMachine::LocalEnergy()
{
    double delta_energy = 0.0;
    Q_factor(r_old_);
    double der1_ln_psi, der2_ln_psi;
    
    for (int n = 0; n < N_; n++){
        for (int d =0; d < D_; d++){
            double sum1 =0.0;
            double sum2 =0.0;
            for (int h = 0; h < H_; h++){
                sum1 += w_(n,d,h)/(1.0+exp(-Q_(h)));
                sum2 += pow(w_(n,d,h),2)*exp(Q_(h))/pow((1.0 + exp(Q_(h))),2);
            }
            der1_ln_psi = -(r_old_(n,d) - a_(n,d))/sigma2_ + sum1/sigma2_;
            der2_ln_psi = -1.0/sigma2_ + sum2/sigma2_;
            delta_energy += 0.5*(-pow(der1_ln_psi,2) - der2_ln_psi + pow(r_old_(n,d),2));
        }
    }


    if (interaction_==1)
    {
        double r_norm;
        for (int n1 = 0; n1 < N_; n1++){
            for (int n2 = 0; n2 < n1; n2++){
                r_norm = 0.0;
                for (int d = 0; d < D_; d++){
                    r_norm += pow((r_old_(n1,d) - r_old_(n2,d)), 2);
                }
                //cout <<r_norm << endl;
                delta_energy += 1.0/sqrt(r_norm);
            }
        }
    }
    //cout <<delta_energy<<endl;

    return delta_energy;
}

/* Method for Adam optimization */
void BoltzmannMachine::ADAM()
{
    double epsilon = 1e-12;   //Value to avoid division by zero.
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

    double Energy = 0.0;
    int its = 100;
    vec Energies = vec(its).fill(0.0);

    for (int i = 0; i < its;  i++){
        Energy = MonteCarlo();
        Energies(i) = Energy;
        cout << "Energy = " << setprecision(15) << Energy << endl;
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
    int its = 1000;
    vec Energies = vec(its).fill(0.0);

    for (int i = 0; i < its;  i++){
        Energy = MonteCarlo();
        Energies(i) = Energy;
        cout << "Energy = " << setprecision(15) << Energy << endl;
        a_ -= eta_*E_da_;
        b_ -= eta_*E_db_;
        w_ -= eta_*E_dw_;
    }
}
