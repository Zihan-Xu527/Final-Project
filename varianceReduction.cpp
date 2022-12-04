//
// Created by Zihan Xu on 12/3/22.
//

#include "varianceReduction.h"
#include "statistic_tool.h"
#include "BSMModel.h"
#include "Payoff.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "omp.h"

varianceReduction::varianceReduction(){

}

varianceReduction::varianceReduction(double r, double volatility, double s_init, double maturity, double strike){



}

varianceReduction::varianceReduction(double r, double volatility, double s_init, double maturity, double strike, int batch_size){
    interest_rate = r;
    sigma = volatility;
    S_init = s_init;
    T = maturity;
    K = strike;
    n = batch_size; //default batch size = 1000
    b = 0.;
    spot_maturity.assign(n, 0.);
    discounted_payoff.assign(n, 0.);
    sol.assign(n, 0.);
}

// simulate n trajectories based on BSM model
void varianceReduction::simulate(){
    BSMModel my_model(interest_rate, sigma, n);
    sim_traj = my_model.sim_path(T, S_init, K, 2);
}

// get terminal spot prices and discounted payoffs from n simulated trajectories
void varianceReduction::set_data(){
//    spot_maturity.resize(n);
//    discounted_payoff.resize(n);
    Call my_call(K);
//#pragma omp parallel for
    for (int i = 0; i < n ; i++){
        spot_maturity[i] = sim_traj[1][i];  // only terminal price is needed
        discounted_payoff[i] = my_call(spot_maturity[i]) * std::exp(- interest_rate * T) ;
    }
}

// get optimal coefficient b from terminal spot prices and discounted payoffs
double varianceReduction::getb(){
    double mean_spot = average(spot_maturity);
    double mean_payoff = average(discounted_payoff);
    double numerator = 0., denominator = 0.;
//#pragma omp parallel for
    for (int i = 0; i < n; i++){
        numerator += (spot_maturity[i] - mean_spot) * (discounted_payoff[i] - mean_payoff);
        denominator += pow(spot_maturity[i] - mean_spot, 2);
    }
    return numerator / denominator;
}
void varianceReduction::init(){
    simulate();
    set_data();
    b = getb();
//    sol.resize(n);
}

// get solution Y(b) = Y - b * (S(T) - E[S(T)])
void varianceReduction::exec(){
#pragma omp parallel for
    for(int i = 0; i < n; i++){
        sol[i] = discounted_payoff[i] - b * (spot_maturity[i] - std::exp(interest_rate * T) * S_init);
    }
}