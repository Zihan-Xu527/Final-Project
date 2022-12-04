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
#include "controlVariates.h"

void check_asymptotic_result(double r, double sigma, double S_0, double T, double K){

    std::vector<int> nn = {10, 100, 1000, 10000};
    double EY = expected_Payoff(0, S_0, r, sigma, K, T);
    std::cout<<"E[Y] = "<< EY << std::endl;
    int m = 10000;


    for(int i = 0; i < nn.size(); i++){

        int n = nn[i];
        std::cout << "n = "<< n <<", ";
        int count = 0;

#pragma omp parallel for
        for(int j = 0; j < m; j++) {
            std::vector<double> spot_maturity(n);
            std::vector<double> discounted_payoff(n);

            controlVariates VR(r, sigma, S_0, T, K, nn[i]);

            VR.init();
            VR.exec();

            std::vector<double> interval = confidenceInterval(VR.get_sol());
            if( interval[0] < EY && EY < interval[1])
            {
                count += 1;
            }
        }

        std::cout << "the probability that E[Y] lies in the 95% confidence interval is " << count / (m + 0.)<< std::endl;

    }
}