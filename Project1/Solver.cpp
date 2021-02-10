#include "Solver.hpp"

Solver::Solver(int N, int num_alphas, int MC, int D,int type_energy, int type_sampling){
    D_ = D;
    N_ = N;
    D_diff_ = 0.5; //Diffusion constant in Greens function
    num_alphas_ = num_alphas;
    MC_ = MC;
    alphas_ = new double[num_alphas_];                  //Variational parameter
    for (int i = 0; i < num_alphas_; i++) alphas_[i] = 0.1 + 0.05*i;
    energies_ = new double[num_alphas_];               //Array to hold energies for different values of alpha
    variances_ = new double[num_alphas_];              //Array to hold variances for different values of alpha
    h_ = 1.0;                                          //Stepsize
    sum_ = 0;
    step_ = h_*pow(10,-4);

    if (type_energy == 0){
        energy_calculation = &Solver::Local_energy_analytical;
    }
    if (type_energy == 1){
        energy_calculation = &Solver::Local_energy_brute_force;
    }

    if (type_sampling == 0){
        MC_method = &Solver::MonteCarlo;
        init_positions = &Solver::Initialize_positions;
        metropolis_sampling = &Solver::Metropolis;
    }

    if (type_sampling == 1){
        MC_method = &Solver::MonteCarlo_importance;
        init_positions = &Solver::Initialize_positions;
        metropolis_sampling = &Solver::Metropolis_importance;
    }

    if (type_sampling == 2){
        MC_method = &Solver::Gradient_descent;
        init_positions = &Solver::Initialize_positions;
        metropolis_sampling = &Solver::Metropolis_importance;

    }

    r_old_ = new double*[N_];
    quantum_force_old_ = new double*[N_];
    quantum_force_new_ = new double[D_];
    r_new_ = new double[D_];
    for (int i= 0; i< N_ ; i++){
        r_old_[i] = new double[D_];
        quantum_force_old_[i] = new double[D_];
    }

    (this->*MC_method)();


}

void Solver::Gradient_descent(){
    double E, V, dalpha;
    double Alphaa = 0.9;        //Initial guess for alpha
    double eta = 0.01;
    int iterations = 50;
    cout <<"Alpha " << "Energy " << "Variance " << endl;
    for (int it = 0; it < iterations; it++){
        E, V, dalpha = MonteCarlo_SGD(Alphaa);
        Alphaa -= eta*dalpha;
        cout <<Alphaa<<" " << E << " " << V<< " " << endl;

    }
}

double Solver::MonteCarlo_SGD(double alpha){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator (0,1)
    uniform_int_distribution<int> RIG(0,N_-1);    //Random integer genererator (0,1)
    double energy, energy_squared, DeltaE, variance, mp, ap, DerivateE, Derivate_WF_E, sum_r, Derivate_WF;
    double r2_sum_old, r2_sum_new;            //New and old r^2 sums (to go in the expression trial function)
    int idx;
    int burn_cycles = 10000;

    energy = 0;
    energy_squared = 0;

    r2_sum_old = Initialize_positions(r2_sum_old);   //Initialize random starting posistions of particle
    Initialize_quantum_force(alpha, r_old_, quantum_force_old_);

    for (int cycle = 0; cycle < MC_; cycle++){
        for (int n = 0; n < N_; n++){
            r2_sum_new = r2_sum_old;

            idx = RIG(gen);   //Index of proposed moved particle
            mp = RDG(gen);    //Magnitude of move
            ap = RDG(gen);    //Acceptance probability

            r2_sum_old = (this->*metropolis_sampling)(r2_sum_new, r2_sum_old, alpha, idx, mp, ap); //Metropolis test, updates position according to accepted/non accepted move

            DeltaE = (this->*energy_calculation)(alpha,r2_sum_old); //Points to either analytical expression for local energy or numerical
            energy += DeltaE;
            energy_squared += DeltaE*DeltaE;
            sum_r = -Init_r_sum(r_old_);
            Derivate_WF += sum_r;
            Derivate_WF_E += sum_r*DeltaE;
            }
        }
    energy /= (MC_*N_);
    energy_squared /= (MC_*N_);
    Derivate_WF/=(MC_*N_);
    Derivate_WF_E/=(MC_*N_);
    variance = energy_squared - energy*energy;
    DerivateE = 2*(Derivate_WF_E - Derivate_WF*energy);

    return energy, variance, DerivateE;
}

void Solver::MonteCarlo(){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator (0,1)
    uniform_int_distribution<int> RIG(0,N_-1);    //Random integer genererator (0,1)
    double alpha, energy, energy_squared, DeltaE, variance, mp, ap;
    double r2_sum_old, r2_sum_new;            //New and old r^2 sums (to go in the expression trial function)
    int idx;
    int burn_cycles = 10000;

    for (int a=0; a < num_alphas_; a++){                //Loop over alpha values
        alpha = alphas_[a];
        energy = 0;
        energy_squared = 0;

        r2_sum_old = Initialize_positions(r2_sum_old);   //Initialize random starting posistions of particle


        if (alpha == 0.5){
            double *b_energies, *b_var;
            int nums = burn_cycles/10;
            b_energies = new double[nums-1];
            b_var = new double[nums-1];
            int k = 1;
            for (int b = 0; b < burn_cycles; b++){
                for (int n = 0; n <N_; n++){
                    r2_sum_new = r2_sum_old;
                    idx = RIG(gen);
                    mp = RDG(gen);
                    ap = RDG(gen);

                    r2_sum_old = Metropolis(r2_sum_new, r2_sum_old, alpha, idx, mp, ap);

                    DeltaE = (this->*energy_calculation)(alpha,r2_sum_old); //Points to either analytical expression for local energy or numerical
                    energy += DeltaE;
                    energy_squared += DeltaE*DeltaE;
                }
                if (b==(10*k)){
                    energy /= (b*N_);
                    energy_squared /= (b*N_);
                    variance = energy_squared - energy*energy;

                    b_energies[k-1] = energy;
                    b_var[k-1] = variance;
                    k+=1;
                    energy = 0;
                    energy_squared = 0;
                }
            }
            ofstream ofile1;
            ofile1.open("./burn_in.txt");
            for (int i = 0; i < nums-1; i++){
                ofile1 << b_energies[i] <<" "<<b_var[i]<< endl;
            }
            ofile1.close();

        }



        energy = 0;
        energy_squared = 0;
        for (int cycle = 0; cycle < MC_; cycle++){
            for (int n = 0; n < N_; n++){
                r2_sum_new = r2_sum_old;

                idx = RIG(gen);   //Index of proposed moved particle
                mp = RDG(gen);    //Magnitude of move
                ap = RDG(gen);    //Acceptance probability

                r2_sum_old = Metropolis(r2_sum_new, r2_sum_old, alpha, idx, mp, ap); //Metropolis test, updates position according to accepted/non accepted move

                DeltaE = (this->*energy_calculation)(alpha,r2_sum_old); //Points to either analytical expression for local energy or numerical
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

void Solver::MonteCarlo_importance(){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator (0,1)
    uniform_int_distribution<int> RIG(0,N_-1);    //Random integer genererator (0,1)
    double alpha, energy, energy_squared, DeltaE, variance, mp, ap;
    double r2_sum_old, r2_sum_new;            //New and old r^2 sums (to go in the expression trial function)
    int idx;
    int burn_cycles = 10000;

    for (int a=0; a < num_alphas_; a++){                //Loop over alpha values
        alpha = alphas_[a];
        energy = 0;
        energy_squared = 0;

        r2_sum_old = Initialize_positions(r2_sum_old);   //Initialize random starting posistions of particle
        Initialize_quantum_force(alpha, r_old_, quantum_force_old_);

        for (int cycle = 0; cycle < MC_; cycle++){
            for (int n = 0; n < N_; n++){
                r2_sum_new = r2_sum_old;

                idx = RIG(gen);   //Index of proposed moved particle
                mp = RDG(gen);    //Magnitude of move
                ap = RDG(gen);    //Acceptance probability

                r2_sum_old = (this->*metropolis_sampling)(r2_sum_new, r2_sum_old, alpha, idx, mp, ap); //Metropolis test, updates position according to accepted/non accepted move

                DeltaE = (this->*energy_calculation)(alpha,r2_sum_old); //Points to either analytical expression for local energy or numerical
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



double Solver::Metropolis(double r2_new, double r2_old, double alpha, int move_idx, double move_P, double acc_P){

    double tf_old, tf_new, P;

    for (int k = 0; k < D_; k++){
        r_new_[k] = r_old_[move_idx][k] + h_ * (move_P - 0.5);
        r2_new = Update_r_sum(r2_new, r_old_[move_idx][k], r_new_[k]);
    }

    tf_old = Trial_func(alpha, r2_old);             //Trial wave function of old position
    tf_new = Trial_func(alpha, r2_new);           //Trial wave function of new position
    P = (tf_new*tf_new)/(tf_old*tf_old);            //Metropolis test
    if (acc_P <= P){
        for (int k =0 ; k< D_;  k++){                   //Update initial posistion
            r_old_[move_idx][k] = r_new_[k];
        }
        r2_old = r2_new;
    }
    return r2_old;
}

double Solver::Metropolis_importance(double r2_new, double r2_old, double alpha, int move_idx, double move_P, double acc_P){
    double tf_old, tf_new, P, greensfunc;

    for (int k = 0; k < D_; k++){
        r_new_[k] = r_old_[move_idx][k] + h_ * (move_P - 0.5);
        r2_new = Update_r_sum(r2_new, r_old_[move_idx][k], r_new_[k]);
    }
    Update_quantum_force(alpha);
    greensfunc = Greens_function(move_idx);

    tf_old = Trial_func(alpha, r2_old);             //Trial wave function of old position
    tf_new = Trial_func(alpha, r2_new);           //Trial wave function of new position
    P = greensfunc*(tf_new*tf_new)/(tf_old*tf_old);            //Metropolis test
    if (acc_P <= P){
        for (int k =0 ; k< D_;  k++){                   //Update initial posistion
            r_old_[move_idx][k] = r_new_[k];
            quantum_force_old_[move_idx][k] = quantum_force_new_[k];
        }
        r2_old = r2_new;
    }
    return r2_old;
}



void Solver::Initialize_quantum_force(double alpha, double **positions, double **q_force){
    for (int j = 0; j < N_; j++){                       //Initial posistion
        for (int k = 0; k < D_; k++){
            q_force[j][k] = -4*alpha*positions[j][k];
        }
    }
}

void Solver::Update_quantum_force(double alpha){
    for (int dd =0; dd< D_; dd++){
        quantum_force_new_[dd] = -4*alpha*r_new_[dd];
    }
}

double Solver::Greens_function(int idx){
    double tmp1, tmp2;
    for (int dd = 0; dd < D_; dd++){
        tmp1 += pow((r_new_[dd]-r_old_[idx][dd]- D_diff_*step_*quantum_force_old_[idx][dd]),2);
        tmp2 += pow((r_old_[idx][dd]-r_new_[dd]- D_diff_*step_*quantum_force_new_[dd]),2);
    }

    return exp(tmp1/tmp2);;
}

double Solver::Initialize_positions(double r2_sum){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator (0,1)

    for (int j = 0; j < N_; j++){                       //Initial posistion
        for (int k = 0; k < D_; k++){
            r_old_[j][k] = h_ * (RDG(gen) - 0.5);
        }
    }
    r2_sum = Init_r_sum(r_old_);               //Inital sum of all r_i^2

    return r2_sum;
}


double Solver::Init_r_sum(double **r){
    sum_ = 0;
    for (int i=0; i<N_; i++){
        for (int j =0; j <D_; j++){
            sum_ += r[i][j]*r[i][j];
        }
    }
    return sum_;
}

double Solver::Update_r_sum(double sum, double r_init, double r_move){
    sum -= r_init*r_init;
    sum += r_move*r_move;
    return sum;
}



double Solver::Trial_func(double alpha, double sum_r_squared){
    return exp(-alpha*sum_r_squared);
}

double Solver::Local_energy_analytical(double alpha, double r_sum){
    return D_*N_*alpha + (1-4*alpha*alpha)*(1./2)*r_sum;
}

double Solver::Local_energy_brute_force(double alpha, double r_sum){
    double dr_p, dr_m;
    laplace_tf_ = 0.0;

    for (int dd = 0; dd < D_; dd++){
        dr_p = r_sum;
        dr_m = r_sum;
        for (int nn = 0; nn < N_; nn++){
            dr_p = Update_r_sum(dr_p, r_old_[nn][dd], r_old_[nn][dd] + step_);
            dr_m = Update_r_sum(dr_m, r_old_[nn][dd], r_old_[nn][dd] - step_);
        }
        laplace_tf_  += Trial_func(alpha,dr_p) + Trial_func(alpha, dr_m);
    }

    tf_middle_ = Trial_func(alpha,r_sum);
    laplace_tf_ -= 2*D_*tf_middle_;
    laplace_tf_ /= (step_*step_*tf_middle_);

    return (1./2)*(-laplace_tf_ + r_sum);
}




void Solver::Write_to_file(string outfilename, double time){
    ofstream ofile;
    ofile.open(outfilename);
    ofile << "alpha" << " " << "energy" << " " << "variance" << endl;
    for (int i = 0; i < num_alphas_; i++){
        ofile << alphas_[i] << " " << energies_[i] << " " << variances_[i] << endl;
    }
    ofile<<" "<<endl;
    ofile << "Timeused: " << time <<endl;
    ofile.close();
}
