//
// Created by Zihan Xu on 12/3/22.
//

#ifndef FINAL_PROJECT_CONTROLVARIATES_H
#define FINAL_PROJECT_CONTROLVARIATES_H
#include "statistic_tool.h"
#include "BSMModel.h"
#include "Payoff.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "omp.h"

class controlVariates {
private:
    double interest_rate; // interest rate: r
    double sigma; // volatility coefficient: \sigma
    double S_init; // initial stock price: S(0)
    double T; // time to maturity: T
    double K; // strike: K
    double b; // optimal coefficient
    int n; // batch size
    std::vector<std::vector<double>> sim_traj; // simulated trajectories from BSM model
    std::vector<double> spot_maturity; // S(T)
    std::vector<double> discounted_payoff; // e^{-rT}*Y(T)
    std::vector<double> sol; // Y(b)
public:
    controlVariates();
    controlVariates(double r, double volatility, double s_init, double maturity, double strike);
    controlVariates(double r, double volatility, double s_init, double maturity, double strike, int newN);
    void simulate();
    void set_data();
    void set_b();
    void exec();
    void set_n(int newN){ n = newN; };
    void set_r(double newr){ interest_rate = newr; };
    void set_init(double newS0){ S_init = newS0; };
    void set_maturity(double newT){ T = newT; };
    void set_strike(double newK) { K = newK; };
    std::vector<double> get_payoff() { return discounted_payoff;};
    std::vector<double> get_sol(){ return sol;};
    void init();
};


#endif //FINAL_PROJECT_CONTROLVARIATES_H
