//
// Created by Zihan Xu on 12/6/22.
//
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "omp.h"
#include "statistic_tool.h"
#include "BSMModel.h"
#include "Payoff.h"
#include <chrono>
#include "controlVariates.h"
#include "statistic_tool.h"
#include <fstream>
void analysis(){
    double r = 0.05, sigma = 0.3;
    double S_0 = 50.0, T = 0.25;

    //Table 1: variance reduction for different strikes
    std::vector<double> strikes = {40., 45., 50., 55., 60., 65., 70.};
    int m = 10000;
    std::vector<double> Y_bar(m);
    std::vector<double> Yb_bar(m);
    for (int i=0; i<strikes.size(); i++){
        std::cout << "K = " << strikes[i] <<", ";
#pragma omp parallel for
        for(int j = 0; j < m; j++) {
            controlVariates cv(r, sigma, S_0, T, strikes[i]);
            cv.init();
            cv.exec();
            Yb_bar[j] = average(cv.get_sol());
            Y_bar[j] = average(cv.get_payoff());
        }
        double var_Y = variance(Y_bar), var_Yb = variance(Yb_bar);
        std::cout<<"variance of Y is " << var_Y <<", ";
        std::cout<<"variance of Y(b) is " << var_Yb <<", ";
        std::cout<<"rho^2 is "<< 1 - var_Yb/var_Y << std::endl;

    }

    // figure 2: linear regressions for strikes 40, 55, 70, n = 50
    std::vector<double> strike = {40., 55., 70.};
    for (int i=0; i<strike.size(); i++){

        controlVariates cv(r, sigma, S_0, T, strike[i], 50);
        cv.init();
        cv.exec();
        std::vector<double> sol = cv.get_sol();
        std::vector<double> spot_maturity = cv.get_spot();
        std::vector<double> discounted_payoff = cv.get_payoff();

        writeResultToFile(spot_maturity, discounted_payoff, sol, strike[i]);

    }
    
    // figure 3: K = 55, n = 10
    controlVariates cv(r, sigma, S_0, T, 55., 10);
    cv.init();
    cv.exec();
    std::vector<double> sol = cv.get_sol();
    std::vector<double> spot_maturity = cv.get_spot();
    std::vector<double> discounted_payoff = cv.get_payoff();

    writeResultToFile(spot_maturity, discounted_payoff, sol);


}