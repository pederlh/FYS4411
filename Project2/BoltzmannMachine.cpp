#include "BoltzmannMachine.hpp"

void Initialize_Parameters(int dimentions, int num_paricles){
    mt19937_64 gen(rd_());
    normal_distribution<double> NDG(0.0,0.001);   //Random number generated from gaussian distribution with mean = 0, std = 0.001;
    D_ = dimentions;
    N_ = num_particles;
    hidden_layers_ = 4;

    w_ = new double**[N_];
    for (int n = 0; n < N_; n++){
       w_[n] = new double*[D_];
       for (int d = 0; d < D_; d++){
           w_[n][d] = new double[hidden_layers_];
           for (int h = 0; h < hidden_layers_; h++){
               w_[n][d][h] = NDG(gen);
           }
       }
   }

   a_ = new double *[N_];
   for (int n = 0; n < N_; n++){
       a_[n] = new double[D_];
       for (int d = 0; d < D_; d++){
           a_[n][d] = NDG(gen);
       }
   }

   b_ = new double*[hidden_layers_];
   for (int h = 0; h < hidden_layers_; h++){
       b_[h] = NDG(gen);
   }
}
