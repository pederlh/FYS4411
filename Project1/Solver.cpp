#include "Solver.hpp"

Solver::Solver(int N, int num_alphas, int MC){
    N_ = N;
    num_alphas_ = num_alphas;
    MC_ = MC;
    alphas_ = new double[num_alphas_];                 //Variational parameter
    for (int i = 0; i < num_alphas_; i++) alphas_[i] = 0.1 + 0.05*i;
    energies_ = new double[num_alphas_];               //Array to hold energies for different values of alpha
    variances_ = new double[num_alphas_];              //Array to hold variances for different values of alpha
    h_ = 1.0;                                          //Stepsize

    MonteCarlo();
}


void Solver::MonteCarlo(){

    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator (0,1)

    double alpha, energy, energy_squared, P, DeltaE, variance;
    double *r_old_, *r_new_;               //New and old posistion
    r_old_ = new double[N_];
    r_new_ = new double[N_];

    double tf_old, tf_new;            //New and old trial wave function

    for (int i=0; i < num_alphas_; i++){
        alpha = alphas_[i];
        energy = 0;
        energy_squared = 0;
        for (int i = 0; i<N_;i++) r_old_[i] = h_ * (RDG(gen) - 0.5);                //Initial posistion
        tf_old = trial_func_1D(alpha, r_old_);       //Initial trial wave function
        for (int sycle = 0; sycle < MC_; sycle++){

            for (int i = 0; i<N_;i++) r_new_[i] = r_old_[i] + h_ * (RDG(gen) - 0.5);             //Proposed move to new posistion
            tf_new = trial_func_1D(alpha, r_new_);           //Trial wave function of new position
            P = (tf_new*tf_new)/(tf_old*tf_old);
            if (P > 1){
                for (i = 0; i<N_;i++) r_old_[i] = r_new_[i];
                tf_old = tf_new;
            }
            else if (RDG(gen) <= P){
                for (i = 0; i<N_;i++) r_old_[i] = r_new_[i];
                tf_old = tf_new;
            }

            DeltaE = local_energy_1D(alpha,r_old_);
            energy += DeltaE;
            energy_squared += DeltaE*DeltaE;
        }
        energy /= MC_;
        energy_squared /= MC_;
        variance = energy_squared - energy*energy;
        energies_[i] = energy;
        variances_[i] = variance;
    }
}

double Solver::trial_func_1D(double alpha, double*r){
    double sum;                                             //VIL IKKE MÅTTE DEKLARERE HVER GANG INNI HER, DETTE ER EN QUICK FIX
    for (int i=0; i<N_; i++){
        sum += r[i]*r[i];
    }
    return exp(-alpha*sum);
}

double Solver::local_energy_1D_analytical(){
    return (1/2)*N_;
}

double Solver::local_energy_1D(double alpha, double*r){
    double *r_forward = new double[N_];
    double *r_backward = new double[N_];
    for (int i = 0; i <N_ ; i++){
        r_forward[i] = r[i] + h_;
        r_backward[i] = r[i] - h_;
    }

    double tf_forward = trial_func_1D(alpha,r_forward);
    double tf_backward = trial_func_1D(alpha,r_backward);
    double tf_middle = trial_func_1D(alpha,r);

    double laplace_tf = (tf_forward - 2*tf_middle + tf_backward)/(h_*h_*tf_middle);

    double sum;                                             //VIL IKKE MÅTTE DEKLARERE HVER GANG INNI HER, DETTE ER EN QUICK FIX
    for (int i=0; i<N_; i++){
        sum += r[i]*r[i];
    }


    return (1/2)*(-laplace_tf + sum);
}

void Solver::write_to_file(string outfilename){
    ofstream ofile;
    ofile.open(outfilename);
    ofile << "alpha" << " " << "energy" << " " << "variance" << endl;
    for (int i = 0; i < num_alphas_; i++){
        ofile << alphas_[i] << " " << energies_[i] << " " << variances_[i] << endl;
    }
    ofile.close();
}
