#include "BoltzmannMachine.hpp"

void BoltzmannMachine::Initialize_Parameters(int dimentions, int num_particles, double eta)
{
    mt19937_64 gen(rd_());
    normal_distribution<double> NDG(0.0,0.001);   //Random number generated from gaussian distribution with mean = 0, std = 0.001;
    D_ = dimentions;
    N_ = num_particles;
    hidden_nodes_ = 4;
    double std_norm_dist = 0.001;  //standard deviation of normal distribution used to initialize weights/biases.

   //fill::randn = set each element to a random value from a normal/Gaussian distribution with zero mean and unit variance
   //Initializing weihgts/biases
   w_ = cube(N_, D_, hidden_nodes_, fill::randn)*std_norm_dist;
   a_ = mat(N_, D_, fill::randn)*std_norm_dist;
   b_ = vec(hidden_nodes_, fill::randn)*std_norm_dist;


}


void BoltzmannMachine::Initialize_SGD()
{
    eta_ = eta; //Learning rate SGD.
    epsilon_ = 1e-8;   //Value to avoid division by zero.

    //Derivatives of weights/biases for stochastic gradient descent
    dw_ = cube(N_, D_, hidden_nodes_).fill(0.0);
    da_ = mat(N_, D_).fill(0.0);
    db_ = vec(hidden_nodes_).fill(0.0);

    //First momentum of weights/biases for stochastic gradient descent
    mom_w_ = cube(N_, D_, hidden_nodes_).fill(0.0);
    mom_a_ = mat(N_, D_).fill(0.0);
    mom_b_ = vec(hidden_nodes_).fill(0.0);

    //Second momentum of weights/biases for stochastic gradient descent
    second_mom_w_ = cube(N_, D_, hidden_nodes_).fill(0.0);
    second_mom_a_ = mat(N_, D_).fill(0.0);
    second_mom_b_ = vec(hidden_nodes_).fill(0.0);

    sigma_ = 0.01;

}

/* Method for Adam optimization */
void BoltzmannMachine::ADAM()
{
    dw_ *= eta_;
    db_ *= eta_;
    da_ *= eta_;

    mom_w_ = beta1_*mom_w_ + (1-beta1_)*dw_;
    mom_b_ = beta1_*mom_b_ + (1-beta1_)*db_;
    mom_a_ = beta1_*mom_a_ + (1-beta1_)*da_;

    second_mom_w_ = beta2_*second_mom_w_ + (1-beta2_)*(dw_%dw_);
    second_mom_b_ = beta2_*second_mom_b_ + (1-beta2_)*(db_%db_);
    second_mom_a_ = beta2_*second_mom_a_ + (1-beta2_)*(da_%da_);

    //alpha_batch_ = eta_*sqrt(1-pow(beta2_, batch_+1))/(1-pow(beta1_, batch_+1));
    //epsilon_batch_ = epsilon_*sqrt(1-pow(beta2_, batch_+1));

    //weights_ -= mom_w_*alpha_batch_/(sqrt(second_mom_w_) + epsilon_batch_);
    //bias_ -= mom_b_*alpha_batch_/(sqrt(second_mom_b_) + epsilon_batch_);

    dw_.fill(0.);
    db_.fill(0.);
}
