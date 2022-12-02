//
// Created by Zihan Xu on 11/29/22.
//
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "omp.h"
#include "statistic_tool.h"
#include "BSMModel.h"
#include "Payoff.h"
void check_reduction(double interest_rate, double sigma, double S_init, double maturity, Call call){
    BSMModel my_model(interest_rate, sigma);
    double exp_factor = std::exp(interest_rate * maturity);
    double expected_spot = S_init * exp_factor; //S(0)*exp(rT)


    std::vector<int> nn = {10, 100, 1000, 10000};
    std::vector<double> Yb_average(nn.size());


    for (int n = 0; n < nn.size(); n++){
        std::cout << "n = " << nn[n] <<", ";

        int m = 10000;
        std::vector<double> Y_bar(m);
        std::vector<double> Yb_bar(m);
#pragma omp parallel for
        for(int i = 0; i < m; i++) {

            std::vector<double> Yb(nn[n]);
            std::vector<double> spot_maturity(nn[n]);
            std::vector<double> discounted_payoff(nn[n]);

            for (int k=0; k<nn[n]; k++){
                spot_maturity[k] = (my_model.sim_path(maturity, S_init, 2))[1];  // only terminal price is needed
                discounted_payoff[k] = call(spot_maturity[k]) ;
            }
            Yb = getYb(spot_maturity, discounted_payoff, expected_spot);
            Yb_bar[i] = average(Yb);
            Y_bar[i] = average(discounted_payoff);
        }
        std::cout<<"variance of Y is "<<variance(Y_bar) <<", ";
        std::cout<<"variance of Y(b) is "<<variance(Yb_bar) <<", ";
        std::cout<<"rho^2 is "<<1-variance(Yb_bar)/variance(Y_bar)<<std::endl;

    }
}