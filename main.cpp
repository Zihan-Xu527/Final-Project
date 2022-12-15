#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "omp.h"
#include "src/statistic_tool.h"
#include "src/BSMModel.h"
#include "src/Payoff.h"
#include <chrono>
#include "src/controlVariates.h"
#include <fstream>

void check_reduction(double r, double sigma, double S_0, double T, double K){

    std::vector<int> nn = {10, 100, 1000, 10000};
    int m = 10000;
    std::vector<double> Y_mean(m);
    std::vector<double> Yb_mean(m);
    std::vector<double> X_mean(m);


    for (int i = 0; i < nn.size(); i++){
        std::cout << "n = " << nn[i] <<", ";
#pragma omp parallel for
        for(int j = 0; j < m; j++) {
            controlVariates cv(r, sigma, S_0, T, K, nn[i]);

            cv.init();
            cv.exec();

            Yb_mean[j] = average(cv.get_sol());
            Y_mean[j] = average(cv.get_payoff());
            X_mean[j] = average(cv.get_spot());
        }

        double var_Y = variance(Y_mean), var_Yb = variance(Yb_mean);
        double rho_squared = 1-var_Yb/var_Y;
        std::cout << "variance of Y is " << var_Y << ", ";
        std::cout <<"variance of Y(b) is "<< var_Yb << ", ";
        std::cout <<"rho^2 is " << rho_squared << std::endl;

    }
    writeResultToFile(X_mean, Y_mean, Yb_mean);
}

void check_asymptotic_result(double r, double sigma, double S_0, double T, double K, double EY){

    std::vector<int> nn = {10, 100, 1000, 10000};
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
int main() {

    double r = 0.05, sigma = 0.3;
    double S_0 = 50.0, T = 0.25, K = 55.; // modify K to 40 or 70 for figure 2, correspond to prob_distribution.m

    double S_T = S_0 * std::exp( r * T );
    std::cout << "E[S(T)] = "<< S_T ;
    double EY = expected_Payoff(0, S_0, r, sigma, K, T);
    std::cout <<", E[Y] = "<< EY << std::endl;

    auto startTime = std::chrono::high_resolution_clock::now();

    check_reduction(r, sigma, S_0, T, K);

    auto middleTime = std::chrono::high_resolution_clock::now();

    check_asymptotic_result(r, sigma, S_0, T, K, EY);

    auto endTime = std::chrono::high_resolution_clock::now();


    std::chrono::duration<double, std::milli> fp_ms1 = middleTime - startTime;
    std::chrono::duration<double, std::milli> fp_ms2 = endTime - middleTime;


    std::cout << "Elapsed time in test 1: "<< fp_ms1.count()/1000 << "s." << std::endl;
    std::cout << "Elapsed time in test 2: "<< fp_ms2.count()/1000 << "s." << std::endl;
    std::cout << "Total time: "<< (fp_ms1.count() + fp_ms2.count())/1000 << "s."<< std::endl;

    void analysis();
    analysis();

    return 0;
}