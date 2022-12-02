#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "omp.h"
#include "statistic_tool.h"
#include "BSMModel.h"
#include "Payoff.h"
#include <chrono>



int main() {
    auto startTime = std::chrono::high_resolution_clock::now();
    void check_reduction(double interest_rate, double sigma, double S_init, double maturity, double strike);
    void check_asymptotic_result(double interest_rate, double sigma, double S_init, double maturity, double strike);


    double r = 0.05, sigma = 0.3;
    double S_0 = 50.0, T = 0.25, K = 55.;

    check_reduction(r, sigma, S_0, T, K);
    check_asymptotic_result(r, sigma, S_0, T, K);

    //check reduction of variance
//    std::vector<double> KK = {40., 45., 50., 55., 60., 65., 70.};


//    for (int i=0; i<KK.size(); i++){
//        std::cout<< "K = "<< KK[i]<< std::endl;
//        Call my_call(KK[i]);
//        check_reduction(r, sigma, S_0, T, my_call);
//    }











    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = endTime - startTime;


    std::cout << fp_ms.count() << std::endl;

    return 0;
}