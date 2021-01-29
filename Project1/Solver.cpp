#include "Solver.hpp"

Solver::Solver(int N, int num_alphas, int MC, int D,int type_energy){
    D_ = D;
    N_ = N;
    num_alphas_ = num_alphas;
    MC_ = MC;
    alphas_ = new double[num_alphas_];                 //Variational parameter
    for (int i = 0; i < num_alphas_; i++) alphas_[i] = 0.1 + 0.05*i;
    energies_ = new double[num_alphas_];               //Array to hold energies for different values of alpha
    variances_ = new double[num_alphas_];              //Array to hold variances for different values of alpha
    h_ = 1.0;                                          //Stepsize
    sum_ = 0;
    r_forward_ = new double[D_*N_];
    r_backward_ = new double[D_*N_];

    if (type_energy == 0){
        energy_calculation = &Solver::local_energy_1D_analytical;
    }
    if (type_energy == 1){
        energy_calculation = &Solver::local_energy_1D_brute_force;
    }
    MonteCarlo();

}

void Solver::MonteCarlo(){

    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RDG(0,1);    //Random double genererator (0,1)

    double alpha, energy, energy_squared, P, DeltaE, variance;
    double *r_old_, *r_new_;               //New and old posistion
    r_old_ = new double[D_*N_];
    r_new_ = new double[D_*N_];
    double tf_old, tf_new;            //New and old trial wave function

    for (int i=0; i < num_alphas_; i++){
        alpha = alphas_[i];
        energy = 0;
        energy_squared = 0;
        for (int j = 0; j<D_*N_;j++) r_old_[j] = h_ * (RDG(gen) - 0.5);                //Initial posistion
        tf_old = trial_func_1D(alpha, r_old_);       //Initial trial wave function
        for (int cycle = 0; cycle < MC_; cycle++){

            for (int k = 0; k<D_*N_;k++) r_new_[k] = r_old_[k] + h_ * (RDG(gen) - 0.5);             //Proposed move to new posistion
            tf_new = trial_func_1D(alpha, r_new_);           //Trial wave function of new position
            P = (tf_new*tf_new)/(tf_old*tf_old);
            if (RDG(gen)<= P){
                for (int l = 0; l<D_*N_;l++) r_old_[l] = r_new_[l];
                tf_old = tf_new;
            }
            DeltaE = (this->*energy_calculation)(alpha,r_old_);   //Points to either analytical expression for local energy or numerical
            energy += DeltaE;
            energy_squared += DeltaE*DeltaE;
        }
        energy /= MC_;

        energy_squared /= MC_;
        variance = energy_squared - energy*energy;
        energies_[i] = energy;
        variances_[i] = variance;
        cout << energies_[i] << endl;
    }

}

double Solver::trial_func_1D(double alpha, double*r){
    sum_ = 0;
    for (int i=0; i<D_*N_; i++){
        sum_ += r[i]*r[i];
    }
    return exp(-alpha*sum_);
}

double Solver::local_energy_1D_analytical(double alpha, double *r){
    sum_ = 0;
    for (int i=0; i<D_*N_; i++){
        sum_ += -(4*alpha*alpha*r[i]*r[i] - 2*alpha) + r[i]*r[i];
    }

    return (1./2)*sum_;
}

double Solver::local_energy_1D_brute_force(double alpha, double*r){
    for (int i = 0; i <D_*N_ ; i++){
        r_forward_[i] = r[i] + h_;
        r_backward_[i] = r[i] - h_;
    }

    tf_forward_ = trial_func_1D(alpha,r_forward_);
    tf_backward_ = trial_func_1D(alpha,r_backward_);
    tf_middle_ = trial_func_1D(alpha,r);


    laplace_tf_ = (tf_forward_ - 2*tf_middle_ + tf_backward_)/(h_*h_*tf_middle_);
    sum_ = 0;
    for (int i=0; i < D_*N_; i++){
        sum_ += r[i]*r[i];
    }
    return (1./2)*(-laplace_tf_ + sum_);
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
