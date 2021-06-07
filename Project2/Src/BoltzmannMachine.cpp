#include "BoltzmannMachine.hpp"

/* Constructor */
BoltzmannMachine::BoltzmannMachine(int num_particles,int dimentions, double eta, int MC, int type_sampling, int interaction, double omega, int num_hidden, string opt, int thread_ID)
{
    arma_rng::set_seed_random();
    thread_ID_ = thread_ID;
    D_ = dimentions;
    N_ = num_particles;
    H_ = num_hidden;
    MC_ = MC;
    double scale = 0.1;

   //Initializing weights/biases
   //fill::randn = set each element to a random value from a normal/Gaussian distribution with zero mean and unit variance
   w_ = cube(N_, D_, H_, fill::randn)*scale;
   a_ = mat(N_, D_, fill::randn)*scale;
   b_ = vec(H_, fill::randn)*scale;
   Q_ = vec(H_).fill(0.0);

   convergence_ = false;

   //Initializing various parameters
   omega_ = omega;
   omega2_ = omega*omega;
   sigma_ = 1.0/sqrt(omega_);
   sigma2_ = sigma_*sigma_;
   eta_ = eta; //Learning rate for gradient descent.
   t_step_ = 0.005;
   D_diff_ = 0.5;
   interaction_ = interaction;

   //Initializing positions
   r_old_ = mat(N_,D_, fill::randn)*sqrt(t_step_);
   r_new_ = mat(N_,D_).fill(0.0);
   for (int n = 0; n < N_ ; n++){
       for (int d = 0; d < D_; d++){
           r_new_(n,d) = r_old_(n,d);
       }
   }

   string ts;

   if (type_sampling == 0)
   {
       MetropolisMethod = &BoltzmannMachine::Metropolis;
       ts = "BF";
   }

   if (type_sampling == 1)
   {
       MetropolisMethod = &BoltzmannMachine::Metropolis_Hastings;
       //Parameters for metropolis-hastings algorithm
       quantum_force_ = mat(N_, D_).fill(0.0);
       quantum_force_old_ = QuantumForce(r_old_);
       quantum_force_new_ = QuantumForce(r_old_);
       ts = "Hast";
   }

   if (opt == "GD")
   {
       optimizer = &BoltzmannMachine::GD;
       filename_ = "EnergySamples_GD";
   }

   if (opt == "ADAM")
   {
       optimizer = &BoltzmannMachine::ADAM;
       filename_ = "EnergySamples_ADAM";
   }

   filename_ = filename_ +"_N_"+ to_string(N_) + "_D_"+ to_string(D_)+ "_H_" + to_string(H_) + "_eta_" + to_string(eta_) + "_MC_" + to_string(MC_) + "_sigma_" + to_string(sigma_) + "_ID_" + to_string(thread_ID_) + "_interaction_" + to_string(interaction_) +  "_omega_" + to_string(omega_) + "_ts_" + ts;

   (this->*optimizer)();
}


/* Method for performing Monte-Carlo simulations */
double BoltzmannMachine::MonteCarlo()
{
    double energy = 0.0;
    double DeltaE = 0.0;
    vec DeltaE_ = vec(MC_).fill(0.0);

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

    //Begin Monte-Carlo cycles
    for (int cycle = 0; cycle < MC_; cycle++){
        for (int out_n = 0; out_n < N_; out_n++){
            (this->*MetropolisMethod)();
        }

        DeltaE = LocalEnergy();                            //Calculate local energy
        DeltaE_(cycle) = DeltaE;
        Derivate_wavefunction();                        //Calculate derivative of wave function with respect to variational parameters
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
    if (convergence_){
        DeltaE_.save(filename_ + "_.txt", raw_ascii);                               //Save samples if convergence has been reached
    }
    return energy;
}


/* Method for calculating derivative of NQS trial wave function w/ reagards to variational parameters (PEDER HELP) */
void BoltzmannMachine::Derivate_wavefunction()
{
    Q_factor(r_old_);
    da_ = (r_old_ - a_)/sigma2_;
    db_ = 1.0/(1.0 + exp(-Q_));
    for (int h = 0; h < H_; h++){
        dw_.slice(h) = r_old_/(sigma2_*(1.0 +  exp(-Q_(h))));
    }
}


/* Method for calculating NQS wavefunction given set of particle positions */
double BoltzmannMachine::WaveFunction(mat r)
{
    Q_factor(r);
    double term1 = 0.0;
    double term2 = 1.0;

    for (int n = 0; n < N_; n++){
        for (int d = 0 ; d < D_; d++){
            term1 += pow((r(n,d) - a_(n,d)),2);
        }
    }

    term1 = exp(-term1/(2*sigma2_));

    for (int h = 0; h < H_; h++){
        term2 *= (1.0 + exp(Q_(h)));
    }

    return term1*term2;

}

/* Method for calculating term in exponentials */
void BoltzmannMachine::Q_factor(mat r)
{
    vec temp = vec(H_).fill(0.0);
    Q_.fill(0.0);

    for (int h = 0; h < H_; h++){
        temp(h) = accu(r%w_.slice(h));
    }

    Q_ = b_ + temp;
}


/* Method for performing metropolis sampling */
void BoltzmannMachine::Metropolis()
{
    int idx = randi<int>(distr_param(0,N_-1));  //Random integer generated from uniform distribution [0,N-1];

    //Proposed move of particle with randomly sampled index "idx"
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
}

/* Method for performing metropolis-hastings sampling */
void BoltzmannMachine::Metropolis_Hastings()
{
    int idx = randi<int>(distr_param(0,N_-1));    //Random integer generated from uniform distribution [0,N-1];
    double rand_gauss = randn<double>();         //Random double generated from gaussian distribution with zero mean and unit variance;
    double GreensFunc;

    //Proposed move of particle with randomly sampled index "idx"
    for (int d = 0; d < D_; d++){
        r_new_(idx,d) = r_old_(idx,d) + D_diff_*quantum_force_old_(idx,d)*t_step_ + rand_gauss*sqrt(t_step_);
    }

    quantum_force_new_ = QuantumForce(r_new_);

    GreensFunc = GreensFunction(idx);           //Calculate Greensfunction
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

/* Method for calculating quantum force (used in metropolis-hastings sampling) */
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

/* Method for calculating greensfunction (used in metropolis-hastings sampling) */
double BoltzmannMachine::GreensFunction(int idx)
{
    double G = 0.0;
    for (int d = 0; d < D_; d++){
        G += 0.5*(quantum_force_old_(idx,d) + quantum_force_new_(idx,d))*(D_diff_*t_step_*0.5*(quantum_force_old_(idx,d)-quantum_force_new_(idx,d)) - r_new_(idx,d) + r_old_(idx,d));
    }
    G = exp(G);
    return G;
}

/* Method for calculating the local energy (HELP PEDER) */
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
                sum2 += pow(w_(n,d,h),2)*exp(-Q_(h))/pow((1.0 + exp(-Q_(h))),2);
            }
            der1_ln_psi = -(r_old_(n,d) - a_(n,d))/sigma2_ + sum1/sigma2_;
            der2_ln_psi = -1.0/sigma2_ + sum2/sigma2_;
            delta_energy += 0.5*(-pow(der1_ln_psi,2) - der2_ln_psi + omega2_*pow(r_old_(n,d),2));
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
                if (r_norm > 6.5e-2){                                           //Avoid contributions from particles that are very close
                    delta_energy += 1.0/sqrt(r_norm);
                }
            }
        }
    }

    return delta_energy;
}

/* Method for Adam optimization */
void BoltzmannMachine::ADAM()
{
    cube w_new = w_;
    mat a_new = a_;
    vec b_new = b_;
    double tol = 1e-4;
    double Ana_e = D_*N_*omega_/2.0;
    if (N_ == 2  && interaction_ == 1 && D_ == 2){
        Ana_e = 3.0;
        tol = 1e-2;
    }

    cube w_joined = cube(N_, D_, H_).fill(0.0);
    mat a_joined = mat(N_, D_).fill(0.0);
    vec b_joined = vec(H_).fill(0.0);
    #pragma omp critical
    {
    w_joined.save("w_joined.txt");
    a_joined.save("a_joined.txt");
    b_joined.save("b_joined.txt");
    }

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


    double Energy = 0.0;
    its = 1000;
    vec Energies = vec(its).fill(0.0);
    double final_E = 0.0;

    #pragma omp master
     {
         if (omp_get_num_threads() == 1) cout << "Start ADAM" << endl;
         else cout << "Start ADAM (showing progress master thread) Dim = "<< D_ << endl << endl;
         cout << "--------------------------------------" << endl;
     }

    for (int i = 0; i < its;  i++){
        Energy = MonteCarlo();
        Energies(i) = Energy;

        #pragma omp master
        {
            cout << "Iteration: " << i << ", Energy = " << setprecision(15) << Energy << endl;
        }
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

        w_new -= mom_w_*alpha_it/(sqrt(second_mom_w_) + epsilon_it);
        b_new -= mom_b_*alpha_it/(sqrt(second_mom_b_) + epsilon_it);
        a_new -= mom_a_*alpha_it/(sqrt(second_mom_a_) + epsilon_it);


        if ((abs(accu(a_new - a_)) < tol && abs(accu(b_new-b_)) < tol && abs(accu(w_new-w_)) < tol) || abs(Ana_e-Energy) < tol){
            convergence_ = true;
            a_ = a_new;
            b_ = b_new;
            w_ = w_new;
            # pragma omp barrier
            # pragma omp critical
            {
                w_joined.load("w_joined.txt");
                a_joined.load("a_joined.txt");
                b_joined.load("b_joined.txt");

                a_joined += a_;
                b_joined += b_;
                w_joined += w_;

                w_joined.save("w_joined.txt");
                a_joined.save("a_joined.txt");
                b_joined.save("b_joined.txt");

            }
            # pragma omp barrier
            # pragma omp critical
            {
                w_joined.load("w_joined.txt");
                a_joined.load("a_joined.txt");
                b_joined.load("b_joined.txt");

                a_joined/=omp_get_num_threads();
                b_joined/=omp_get_num_threads();
                w_joined/=omp_get_num_threads();

                a_ = a_joined;
                b_ = b_joined;
                w_ = w_joined;
            }
            filename_ = filename_ + "_its_" + to_string(i);
            MC_ *=pow(2,2);
            final_E = MonteCarlo();
            break;
        }

        a_ = a_new;
        b_ = b_new;
        w_ = w_new;
    }
}


/* Method for vanilla Gradient Descent optimization */
void BoltzmannMachine::GD()
{
    double Energy = 0.0;
    its = 200;
    vec Energies = vec(its).fill(0.0);
    cube w_new = w_;
    mat a_new = a_;
    vec b_new = b_;
    double final_E = 0.0;

    cube w_joined = cube(N_, D_, H_).fill(0.0);
    mat a_joined = mat(N_, D_).fill(0.0);
    vec b_joined = vec(H_).fill(0.0);
    #pragma omp critical
    {
    w_joined.save("w_joined.txt");
    a_joined.save("a_joined.txt");
    b_joined.save("b_joined.txt");
    }
    double tol = 1e-3;
    double Ana_e = D_*N_*omega_/2.0;
    if (N_ == 2  && interaction_ == 1 && D_ == 2){
        Ana_e = 3.0;
        tol = 1e-2;
    }

    #pragma omp master
     {
         if (omp_get_num_threads() == 1) cout << "Start GD" << endl;
         else cout << "Start GD (showing progress master thread) Dim = "<< D_ << endl << endl;
         cout << "--------------------------------------" << endl;
     }

    for (int i = 0; i < its;  i++){
        Energy = MonteCarlo();
        Energies(i) = Energy;
        #pragma omp master
        {
        cout << "Iteration:" << i << ", Energy = " << setprecision(15) << Energy << endl;
        }

        a_new -= eta_*E_da_;
        b_new -= eta_*E_db_;
        w_new -= eta_*E_dw_;

        if ((abs(accu(a_new - a_)) < tol && abs(accu(b_new-b_)) < tol && abs(accu(w_new-w_)) < tol) || abs(Ana_e-Energy) < tol){
            convergence_ = true;
            a_ = a_new;
            b_ = b_new;
            w_ = w_new;
            # pragma omp barrier
            # pragma omp critical
            {
                w_joined.load("w_joined.txt");
                a_joined.load("a_joined.txt");
                b_joined.load("b_joined.txt");

                a_joined += a_;
                b_joined += b_;
                w_joined += w_;

                w_joined.save("w_joined.txt");
                a_joined.save("a_joined.txt");
                b_joined.save("b_joined.txt");

            }
            # pragma omp barrier
            # pragma omp critical
            {
                w_joined.load("w_joined.txt");
                a_joined.load("a_joined.txt");
                b_joined.load("b_joined.txt");

                a_joined/=omp_get_num_threads();
                b_joined/=omp_get_num_threads();
                w_joined/=omp_get_num_threads();

                a_ = a_joined;
                b_ = b_joined;
                w_ = w_joined;
            }
            filename_ = filename_ + "_its_" + to_string(i);
            MC_ *=pow(2,2);
            final_E = MonteCarlo();
            break;
        }

        a_ = a_new;
        b_ = b_new;
        w_ = w_new;
    }
}
