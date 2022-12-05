//
// Created by Zihan Xu on 12/3/22.
//

#include "controlVariates.h"
#include "statistic_tool.h"
#include "BSMModel.h"
#include "Payoff.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "omp.h"

controlVariates::controlVariates(){

}


controlVariates::controlVariates(double r, double volatility, double s_init, double maturity, double strike, int batch_size){
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
void controlVariates::simulate(){
    BSMModel my_model(interest_rate, sigma, n);
    sim_traj = my_model.sim_path(T, S_init, K, 2);
}

// get terminal spot prices S(T) and discounted payoffs Y from n simulated trajectories
void controlVariates::set_data(){
    Call my_call(K);
    for (int i = 0; i < n ; i++){
        spot_maturity[i] = sim_traj[1][i];  // only terminal price is needed S(T)
        discounted_payoff[i] = my_call(spot_maturity[i]) * std::exp(- interest_rate * T) ; // e^{-rT}(S(T) - K)+
    }
}

// get optimal coefficient b from terminal spot prices and discounted payoffs
void controlVariates::set_b(){
    double mean_spot = average(spot_maturity);
    double mean_payoff = average(discounted_payoff);
    double numerator = 0., denominator = 0.;
// cov(S(T), Y(0)) / var(S(T))
    for (int i = 0; i < n; i++){
        numerator += (spot_maturity[i] - mean_spot) * (discounted_payoff[i] - mean_payoff);
        denominator += pow(spot_maturity[i] - mean_spot, 2);
    }
    b = numerator / denominator;
}

// simulate to get S(T), set_data to get Y, and use S(T), Y to get optimal coefficient b
void controlVariates::init(){
    simulate();
    set_data();
    set_b();
}

// control variates method: Y(b) = Y - b * (S(T) - E[S(T)])
void controlVariates::exec(){
#pragma omp parallel for
    for(int i = 0; i < n; i++){
        sol[i] = discounted_payoff[i] - b * (spot_maturity[i] - std::exp(interest_rate * T) * S_init);
    }
}