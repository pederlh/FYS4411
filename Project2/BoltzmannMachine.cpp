#include "BoltzmannMachine.hpp"

void BoltzmannMachine::Initialize_Parameters(int dimentions, int num_particles){
    mt19937_64 gen(rd_());
    normal_distribution<double> NDG(0.0,0.001);   //Random number generated from gaussian distribution with mean = 0, std = 0.001;
    D_ = dimentions;
    N_ = num_particles;
    hidden_nodes_ = 4;
    double std_norm_dist = 0.001;

   //fill::randn = set each element to a random value from a normal/Gaussian distribution with zero mean and unit variance
   w_ = cube(N_, D_, hidden_nodes_, fill::randn)*std_norm_dist;
   a_ = mat(N_, D_, fill::randn)*std_norm_dist;
   b_ = vec(hidden_nodes_, fill::randn)*std_norm_dist;


}
