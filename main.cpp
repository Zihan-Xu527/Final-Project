#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "omp.h"
#include "statistic_tool.h"
#include "BSMModel.h"
#include "Payoff.h"
#include <chrono>
#include "varianceReduction.h"


int main() {
    void check_reduction(double r, double sigma, double S_0, double T, double K);
    double r = 0.05, sigma = 0.3;
    double S_0 = 50.0, T = 0.25, K = 55.;

//    BSMModel model(r, sigma, 20);
//    std::vector<std::vector<double>> sim=model.sim_path(T, S_0, 3);
//
//    for(int i=0; i<3; i++){
//        for(int j=0; j<20; j++){
//            std::cout<<sim[i][j]<<", ";
//
//    }
//        std::cout<<std::endl;
//
//    }
//    std::vector<double> result(20);
//    for

    auto startTime = std::chrono::high_resolution_clock::now();

//    varianceReduction VR(r, sigma, S_0, T, K);
//    VR.simulate();

    check_reduction(r, sigma, S_0, T, K);




//    auto startTime = std::chrono::high_resolution_clock::now();
//    void check_reduction(double interest_rate, double sigma, double S_init, double maturity, double strike);
//    void check_asymptotic_result(double interest_rate, double sigma, double S_init, double maturity, double strike);
//
//
//    double r = 0.05, sigma = 0.3;
//    double S_0 = 50.0, T = 0.25, K = 65.;
//
//    check_reduction(r, sigma, S_0, T, K);
//    auto middleTime = std::chrono::high_resolution_clock::now();
//    check_asymptotic_result(r, sigma, S_0, T, K);






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




//    int m = 10000;
//    std::vector<double> Y_bar(m);
//    std::vector<double> Yb_bar(m);
////#pragma omp parallel for
//    for(int j = 0; j < m; j++) {
//        VR.exec();
//        Yb_bar[j] = average(VR.get_sol());
//        Y_bar[j] = average(VR.get_payoff());
//    }
//    double var_Y = variance(Y_bar), var_Yb = variance(Yb_bar);
//    std::cout<<"variance of Y is " << var_Y <<", ";
//    std::cout<<"variance of Y(b) is " << var_Yb <<", ";
//    std::cout<<"rho^2 is "<< 1 - var_Yb/var_Y << std::endl;







    auto middleTime = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> fp_ms1 = middleTime - startTime;
//    std::chrono::duration<double, std::milli> fp_ms1 = middleTime - startTime;
//    std::chrono::duration<double, std::milli> fp_ms2 = endTime - middleTime;
//
//
    std::cout << "Part I: "<< fp_ms1.count() << std::endl;
//    std::cout << "Part II: "<< fp_ms2.count() << std::endl;
//    std::cout << "Total time: "<< fp_ms1.count() + fp_ms2.count() << std::endl;

    return 0;
}