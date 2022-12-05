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
#include <fstream>

void writeResultToFile(std::vector<double> & vec1, std::vector<double> & vec2, std::vector<double> & vec3)
{
    std::ofstream output;
    int n = vec1.size();
    std::string fileName = "../sample_size";
    fileName.append("_").append(std::to_string(n)).append(".csv");
    output.open(fileName);
    for (int i = 0; i < n; i++)//save stock price and hedge payoff into a csv file
    {
        output << vec1[i] << "," << vec2[i] << "," << vec3[i] << "\n";
    }
    output.close();

}

int main() {
    void check_reduction(double r, double sigma, double S_0, double T, double K);
    void check_asymptotic_result(double r, double sigma, double S_0, double T, double K);
    double r = 0.05, sigma = 0.3;
    double S_0 = 50.0, T = 0.25, K = 55.;
    double S_T = S_0 * std::exp( r * T );

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

////Table 1: variance reduction for different strikes
//    std::vector<double> KK = {40., 45., 50., 55., 60., 65., 70.};
//    int m = 10000;
//    std::vector<double> Y_bar(m);
//    std::vector<double> Yb_bar(m);
//    for (int i=0; i<KK.size(); i++){
//        std::cout << "K = " << KK[i] <<", ";
//#pragma omp parallel for
//        for(int j = 0; j < m; j++) {
//            controlVariates cv(r, sigma, S_0, T, KK[i]);
//            cv.init();
//            cv.exec();
//            Yb_bar[j] = average(cv.get_sol());
//            Y_bar[j] = average(cv.get_payoff());
//        }
//        double var_Y = variance(Y_bar), var_Yb = variance(Yb_bar);
//        std::cout<<"variance of Y is " << var_Y <<", ";
//        std::cout<<"variance of Y(b) is " << var_Yb <<", ";
//        std::cout<<"rho^2 is "<< 1 - var_Yb/var_Y << std::endl;
//
//    }
//
//    // figure of linear regression
//    std::vector<int> nn = {10, 50, 100, 1000};
//    for (int i=0; i<nn.size(); i++){
////        std::cout << "n = " << nn[i] <<", ";
//        controlVariates cv(r, sigma, S_0, T, K, nn[i]);
//        cv.init();
//        cv.exec();
//        std::vector<double> sol = cv.get_sol();
//        std::vector<double> spot_maturity = cv.get_spot();
//        std::vector<double> discounted_payoff = cv.get_payoff();
//        double Yb_bar = average(sol);
//        double Y_bar = average(discounted_payoff);
//        double X_bar = average(spot_maturity);
//        writeResultToFile(spot_maturity, discounted_payoff, sol);
////        std::cout << " mean[Yb] = "<< Yb_bar << ", mean[Y] = "<<Y_bar<< ", mean[X] = "<< X_bar<< ", E[X] = "<< S_T<<std::endl;
//    }

    return 0;
}