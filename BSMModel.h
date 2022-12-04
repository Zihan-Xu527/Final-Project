//
// Created by Zihan Xu on 11/29/22.
//

#ifndef FINAL_PROJECT_BSMMODEL_H
#define FINAL_PROJECT_BSMMODEL_H
#include <random>
#include <vector>
#include <cmath>
#include "Payoff.h"

class BSMModel {
private:
    std::random_device rd; // for generating a random seed
    std::mt19937 rng;
    std::normal_distribution<double> normal;
    double interest_rate;  // interest rate
    double sigma;  // volatility coefficient
    int batch_size;
    double strike;


public:
//    BSMModel(){}
    BSMModel(double r, double sigma, int n) : interest_rate(r), sigma(sigma), batch_size(n), rng(rd()), normal(0.0, 1.0) { };
    void set_parameter(double mu_new, double sigma_new) { interest_rate = mu_new; sigma = sigma_new; };
    void set_seed(unsigned seed_new) { rng.seed(seed_new); };
    std::vector<std::vector<double>> sim_path(double maturity, double S_init, double strike, unsigned length_out);


};


#endif //FINAL_PROJECT_BSMMODEL_H
