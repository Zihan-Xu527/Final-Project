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
#include "varianceReduction.h"

void exec(varianceReduction &VR, int m, std::vector<double> &Y_bar, std::vector<double> &Yb_bar);

void check_reduction(double r, double sigma, double S_0, double T, double K){
    varianceReduction VR(r, sigma, S_0, T, K);
    std::vector<int> nn = {10, 100, 1000, 10000};
//    omp_set_num_threads(20);

    for (int i = 0; i < nn.size(); i++){
        std::cout << "n = " << nn[i] <<", ";
        VR.set_n(nn[i]);

        int m = 10000;
        std::vector<double> Y_bar(m);
        std::vector<double> Yb_bar(m);

//#pragma omp parallel for
        exec(VR, m, Y_bar, Yb_bar);
//        double var_Y = variance(Y_bar), var_Yb = variance(Yb_bar);
//        double rho_squared = 1-var_Yb/var_Y;
//        std::cout << "variance of Y is " << var_Y << ", ";
//        std::cout <<"variance of Y(b) is "<< var_Yb << ", ";
//        std::cout <<"rho^2 is " << rho_squared << std::endl;

    }
}

void exec(varianceReduction &VR, int m, std::vector<double> &Y_bar, std::vector<double> &Yb_bar) {
//#pragma omp parallel for
    for(int j = 0; j < m; j++) {
        VR.init();
        VR.exec();
        Yb_bar[j] = average(VR.get_sol());
        Y_bar[j] = average(VR.get_payoff());
    }
}
