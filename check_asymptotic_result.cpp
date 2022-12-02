//
// Created by Zihan Xu on 12/2/22.
//
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "omp.h"
#include "statistic_tool.h"
#include "BSMModel.h"
#include "Payoff.h"

void check_asymptotic_result(double interest_rate, double sigma, double S_init, double maturity, double strike){
    BSMModel my_model(interest_rate, sigma);
    Call call(strike);
    std::vector<int> nn = {10, 100, 1000, 10000};
    double EY = expected_Payoff(0, S_init, interest_rate, sigma, strike, maturity);
    std::cout<<"E[Y] = "<< EY << std::endl;
    double exp_factor = std::exp(interest_rate * maturity);
    double expected_spot = S_init * exp_factor; //S(0)*exp(rT)


    int m = 10000;



    for(int i = 0; i < nn.size(); i++){

        int n = nn[i];
        std::cout << "n = "<< n <<", ";
        int count = 0;

#pragma omp parallel for
        for(int j = 0; j < m; j++) {
            std::vector<double> spot_maturity(n);
            std::vector<double> discounted_payoff(n);

            for (int k = 0 ; k < n ; k++){
                spot_maturity[k] = (my_model.sim_path(maturity, S_init, 2))[1];  // only terminal price is needed
                discounted_payoff[k] = call(spot_maturity[k]) / exp_factor;
            }

            std::vector<double> Yb = getYb(spot_maturity, discounted_payoff, expected_spot);

            std::vector<double> interval = confidenceInterval(Yb);
            if( interval[0] < EY && EY < interval[1])
            {
                count += 1;
            }
        }

        std::cout << "the probability that E[Y] lies in the 95% confidence interval is " << count / (m + 0.)<< std::endl;

    }
}