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


int main() {
    void check_reduction(double r, double sigma, double S_0, double T, double K);
    void check_asymptotic_result(double r, double sigma, double S_0, double T, double K);
    double r = 0.05, sigma = 0.3;
    double S_0 = 50.0, T = 0.25, K = 55.;

    auto startTime = std::chrono::high_resolution_clock::now();

    check_reduction(r, sigma, S_0, T, K);

    auto middleTime = std::chrono::high_resolution_clock::now();

    check_asymptotic_result(r, sigma, S_0, T, K);

    auto endTime = std::chrono::high_resolution_clock::now();


    std::chrono::duration<double, std::milli> fp_ms1 = middleTime - startTime;
    std::chrono::duration<double, std::milli> fp_ms2 = endTime - middleTime;


    std::cout << "Elapsed time in test 1: "<< fp_ms1.count()/1000 << "s." << std::endl;
    std::cout << "Elapsed time in test 2: "<< fp_ms2.count()/1000 << "s." << std::endl;
    std::cout << "Total time: "<< fp_ms1.count() + fp_ms2.count() << std::endl;

    //    //check reduction of variance
//    std::vector<double> KK = {40., 45., 50., 55., 60., 65., 70.};
//
//
//
//    for (int i=0; i<KK.size(); i++){
//        VR.set_strike(KK[i]);
//        std::cout << "K = " << KK[i] <<", ";
//        int m = 10000;
//        std::vector<double> Y_bar(m);
//        std::vector<double> Yb_bar(m);
////#pragma omp parallel for
//        for(int j = 0; j < m; j++) {
//            VR.exec();
//            Yb_bar[j] = average(VR.get_sol());
//            Y_bar[j] = average(VR.get_payoff());
//        }
//        double var_Y = variance(Y_bar), var_Yb = variance(Yb_bar);
//        std::cout<<"variance of Y is " << var_Y <<", ";
//        std::cout<<"variance of Y(b) is " << var_Yb <<", ";
//        std::cout<<"rho^2 is "<< 1 - var_Yb/var_Y << std::endl;
//
//    }

    return 0;
}