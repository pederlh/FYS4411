#include "Solver_copy.hpp"

Solver_copy::Solver_copy(int N, int num_alphas, int MC, int D,int type_energy){
    D_ = D;
    N_ = N;
    num_alphas_ = num_alphas;
    MC_ = MC;
    alphas_ = new double[num_alphas_];                  //Variational parameter
    for (int i = 0; i < num_alphas_; i++) alphas_[i] = 0.1 + 0.05*i;
    energies_ = new double[num_alphas_];               //Array to hold energies for different values of alpha
    variances_ = new double[num_alphas_];              //Array to hold variances for different values of alpha
    h_ = 1.0;                                          //Stepsize
    step_ = h_*pow(10,-4);

    if (type_energy == 0){
        energy_calculation = &Solver_copy::Local_energy_analytical;
    }
    if (type_energy == 1){
        energy_calculation = &Solver_copy::Local_energy_brute_force;
    }

    r_old_ = new double*[N_];
    r_new_ = new double[D_];
    for (int row = 0; row < N_ ; row++){
        r_old_[row] = new double[D_];
    }

}


void Solver_copy::MonteCarlo(){

    double alpha, energy, energy_squared, DeltaE, variance;
    double r2_sum_old, r2_sum_new;            //New and old r^2 sums (to go in the expression trial function)
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator (0,1)



    for (int a=0; a < num_alphas_; a++){                //Loop over alpha values
        alpha = alphas_[a];
        energy = 0;
        energy_squared = 0;
        r2_sum_old = 0;

        for (int j = 0; j < N_; j++){                       //Initial posistion
            for (int k = 0; k < D_; k++){
                r_old_[j][k] = h_ * (RDG(gen) - 0.5);
            }
        }

        for (int i = 0; i < N_; i++){
            for (int j = 0; j < D_; j++){
                r2_sum_old += r_old_[i][j]*r_old_[i][j];
            }
        }

        for (int cycle = 0; cycle < MC_; cycle++){
            for (int n = 0; n < N_; n++){
                r2_sum_new = r2_sum_old;

                r2_sum_old = Metropolis(r2_sum_new, r2_sum_old, alpha); //Metropolis test, updates position according to accepted/non accepted move

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


double Solver_copy::Metropolis(double r2_new, double r2_old, double alpha){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]
    uniform_int_distribution<int> RIG(0, N_-1);    //Random integer genererator [0,N_]

    double tf_old, tf_new, P;
    int idx = RIG(gen);

    for (int k = 0; k < D_; k++){
        r_new_[k] = r_old_[idx][k] + h_ * (RDG(gen) - 0.5);
        r2_new = Update_r_sum(r2_new, r_old_[idx][k], r_new_[k]);
    }

    tf_old = Trial_func(alpha, r2_old);             //Trial wave function of old position
    tf_new = Trial_func(alpha, r2_new);           //Trial wave function of new position

    P = (tf_new*tf_new)/(tf_old*tf_old);            //Metropolis test

    if (RDG(gen) <= P){
        for (int k =0 ; k < D_;  k++){                   //Update initial posistion
            r_old_[idx][k] = r_new_[k];
        }
        r2_old = r2_new;
    }

    return r2_old;
}


double Solver_copy::Update_r_sum(double sum, double r_init, double r_move){
    sum -= r_init*r_init;
    sum += r_move*r_move;
    return sum;
}

double Solver_copy::Trial_func(double alpha, double sum_r_squared){
    return exp(-alpha*sum_r_squared);
}

double Solver_copy::Local_energy_analytical(double alpha, double r_sum){
    return D_*N_*alpha + (1-4*alpha*alpha)*(1./2)*r_sum;
}

double Solver_copy::Local_energy_brute_force(double alpha, double r_sum){
    double dr_p, dr_m;
    laplace_tf_ = 0.0;

    for (int nn = 0; nn < N_; nn++){
        for (int dd = 0; dd < D_; dd++){
            dr_p = Update_r_sum(r_sum, r_old_[nn][dd], r_old_[nn][dd] + step_);
            //cout << dr_p << endl;
            dr_m = Update_r_sum(r_sum, r_old_[nn][dd], r_old_[nn][dd] - step_);
            laplace_tf_  += Trial_func(alpha,dr_p) + Trial_func(alpha, dr_m);
        }
    }

    tf_middle_ = Trial_func(alpha,r_sum);
    laplace_tf_ -= 2*D_*N_*tf_middle_;
    laplace_tf_ /= (step_*step_*tf_middle_);

    return (1./2)*(-laplace_tf_ + r_sum);
}


void Solver_copy::Write_to_file(string outfilename, double time){
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



//monte carlo som GD, kjør ytre alpha loop som håndterer både importance og ikke
void Solver::Gradient_descent(){
    double *values;
    values = new double[3];
    //double E, V, dalpha;
    double Alphaa = 0.9;        //Initial guess for alpha
    double eta = 0.01;
    int iterations = 50;
    cout <<"Alpha " << "Energy " << "Variance " << endl;
    for (int it = 0; it < iterations; it++){
        MonteCarlo_SGD(values, Alphaa);
        Alphaa -= eta*values[2];
        cout <<Alphaa<<" " << values[0] << " " << values[1]<< " " << endl;
        if (values[1]< pow(10,-9)){
            break;
        }

    }
}


void Solver::MonteCarlo_SGD(double *values, double alpha){
    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator (0,1)
    uniform_int_distribution<int> RIG(0,N_-1);    //Random integer genererator (0,1)
    double energy, energy_squared, DeltaE, variance, mp, ap, DerivateE, Derivate_WF_E, sum_r, Derivate_WF;
    double r2_sum_old, r2_sum_new;            //New and old r^2 sums (to go in the expression trial function)

    energy = 0;
    energy_squared = 0;
    r2_sum_old = 0;

    for (int j = 0; j < N_; j++){                       //Initial posistion
        for (int k = 0; k < D_; k++){
            wave.r_old_[j][k] = h_ * (RDG(gen) - 0.5);
        }
    }

    for (int i = 0; i < N_; i++){
        for (int j =0; j <D_; j++){
            r2_sum_old += wave.r_old_[i][j]*wave.r_old_[i][j];
        }
    }
    Initialize_quantum_force(alpha, wave.r_old_, quantum_force_old_);

    for (int cycle = 0; cycle < MC_; cycle++){
        for (int n = 0; n < N_; n++){
            r2_sum_new = r2_sum_old;

            r2_sum_old = (this->*metropolis_sampling)(r2_sum_new, r2_sum_old, alpha); //Metropolis test, updates position according to accepted/non accepted move

            if (type_energy_ == 0){
                DeltaE = wave.Local_energy_analytical(alpha);
            }
            if (type_energy_==1){
                DeltaE = wave.Local_energy_brute_force(alpha);
            }
            energy += DeltaE;
            energy_squared += DeltaE*DeltaE;
            sum_r = -wave.r2_sum_old_;
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

    values[0] = energy;
    values[1] = variance;
    values[2] = DerivateE;


}

//DENNE ER IKKE I BRUK PEDER :-) Skal forhåoentligvis brukes når jeg får svar på varianse problemet, så bare ignorer den for nå
void Solver::MonteCarlo_burn(){
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
        r2_sum_old = 0;

        for (int j = 0; j < N_; j++){                       //Initial posistion
            for (int k = 0; k < D_; k++){
                wave.r_old_[j][k] = h_ * (RDG(gen) - 0.5);
            }
        }

        for (int i = 0; i < N_; i++){
            for (int j =0; j <D_; j++){
                r2_sum_old += wave.r_old_[i][j]*wave.r_old_[i][j];
            }
        }


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

                    r2_sum_old = Metropolis(r2_sum_new, r2_sum_old, alpha);

                    if (type_energy_ == 0){
                        DeltaE = wave.Local_energy_analytical(alpha);
                    }
                    if (type_energy_==1){
                        DeltaE = wave.Local_energy_brute_force(alpha);
                    }
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

                r2_sum_old = (this->*metropolis_sampling)(r2_sum_new, r2_sum_old, alpha); //Metropolis test, updates position according to accepted/non accepted move

                if (type_energy_ == 0){
                    DeltaE = wave.Local_energy_analytical(alpha);
                }
                if (type_energy_== 1){
                    DeltaE = wave.Local_energy_brute_force(alpha);
                }
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


double Solver::Metropolis_importance(double r2_new, double r2_old, double alpha){

    mt19937_64 gen(rd_());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator [0,1]
    uniform_int_distribution<int> RIG(0, N_-1);    //Random integer genererator [0,N_]
    double tf_old, tf_new, P, greensfunc;
    int idx = RIG(gen);

    for (int k = 0; k < D_; k++){
        wave.r_new_[k] = wave.r_old_[idx][k] + h_ * (RDG(gen) - 0.5);
        r2_new = Update_r_sum(r2_new, wave.r_old_[idx][k], wave.r_new_[k]);
    }
    Update_quantum_force(alpha);
    greensfunc = Greens_function(idx);

    tf_old = wave.Trial_func(alpha, r2_old);             //Trial wave function of old position
    tf_new = wave.Trial_func(alpha, r2_new);           //Trial wave function of new position
    P = greensfunc*(tf_new*tf_new)/(tf_old*tf_old);            //Metropolis test
    if (RDG(gen) <= P){
        for (int k =0 ; k < D_;  k++){                   //Update initial posistion
            wave.r_old_[idx][k] = wave.r_new_[k];
            quantum_force_old_[idx][k] = quantum_force_new_[k];
        }
        r2_old = r2_new;
    }
    return r2_old;
}

void Solver::Update_quantum_force(double alpha){
    for (int dd =0; dd< D_; dd++){
        quantum_force_new_[dd] = -4*alpha*wave.r_new_[dd];
    }
}

double Solver::Greens_function(int idx){
    double tmp1, tmp2;
    for (int dd = 0; dd < D_; dd++){
        tmp1 += pow((wave.r_old_[idx][dd]-wave.r_new_[dd]- D_diff_*step_*quantum_force_new_[dd]),2)*(1/(4*D_diff_*step_));
        tmp2 += pow((wave.r_new_[dd]-wave.r_old_[idx][dd]- D_diff_*step_*quantum_force_old_[idx][dd]),2)*(1/(4*D_diff_*step_));

    }
    tmp1 = exp(-tmp1);
    tmp2 = exp(-tmp2);
    return tmp1/tmp2;
    //return exp(tmp1/tmp2);
}

//If brute force, this function (which does nothing) is called instead of Initialize quantum force
void Solver::No_quantum_force(double alpha, double **positions, double **q_force){}

void Solver::Initialize_quantum_force(double alpha, double **positions, double **q_force){
    for (int j = 0; j < N_; j++){                       //Initial posistion
        for (int k = 0; k < D_; k++){
            q_force[j][k] = -4*alpha*positions[j][k];
        }
    }
}
