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

void writeResultToFile(std::vector<double> & vec1, std::vector<double> & vec2, std::vector<double> & vec3)
{
    std::ofstream output;
    int n = vec1.size();
    std::string fileName = "../output/sample_size";
    fileName.append("_").append(std::to_string(n)).append(".csv");
    output.open(fileName);
    for (int i = 0; i < n; i++)
    {
        output << vec1[i] << "," << vec2[i] << "," << vec3[i] << "\n";
    }
    output.close();

}

void writeResultToFile(std::vector<double> & vec1, std::vector<double> & vec2, std::vector<double> & vec3, double strike)
{
    std::ofstream output;
    int n = vec1.size();
    std::string fileName = "../output/strike";
    fileName.append("_").append(std::to_string(int(strike))).append(".csv");
    output.open(fileName);
    for (int i = 0; i < n; i++)
    {
        output << vec1[i] << "," << vec2[i] << "," << vec3[i] << "\n";
    }
    output.close();

}

void check_reduction(double r, double sigma, double S_0, double T, double K){

    std::vector<int> nn = {10, 100, 1000, 10000};
//    omp_set_num_threads(8);
    int m = 10000;
    std::vector<double> Y_bar(m);
    std::vector<double> Yb_bar(m);


    for (int i = 0; i < nn.size(); i++){
        std::cout << "n = " << nn[i] <<", ";
#pragma omp parallel for
        for(int j = 0; j < m; j++) {
            controlVariates VR(r, sigma, S_0, T, K, nn[i]);

            VR.init();
            VR.exec();

            Yb_bar[j] = average(VR.get_sol());
            Y_bar[j] = average(VR.get_payoff());
        }

        double var_Y = variance(Y_bar), var_Yb = variance(Yb_bar);
        double rho_squared = 1-var_Yb/var_Y;
        std::cout << "variance of Y is " << var_Y << ", ";
        std::cout <<"variance of Y(b) is "<< var_Yb << ", ";
        std::cout <<"rho^2 is " << rho_squared << std::endl;

    }
}

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
int main() {

//    void check_asymptotic_result(double r, double sigma, double S_0, double T, double K);
    double r = 0.05, sigma = 0.3;
    double S_0 = 50.0, T = 0.25, K = 55.;
    double S_T = S_0 * std::exp( r * T );

    auto startTime = std::chrono::high_resolution_clock::now();

//    check_reduction(r, sigma, S_0, T, K);

    auto middleTime = std::chrono::high_resolution_clock::now();

//    check_asymptotic_result(r, sigma, S_0, T, K);

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
    // figure of linear regression
    std::vector<int> nn = {10, 50, 100, 1000};
    std::vector<double> strikes = {40., 45., 55., 70.};
    int sampleSize = 50;
    for (int i=0; i<strikes.size(); i++){

        controlVariates cv(r, sigma, S_0, T, strikes[i], sampleSize);
        cv.init();
        cv.exec();
        std::vector<double> sol = cv.get_sol();
        std::vector<double> spot_maturity = cv.get_spot();
        std::vector<double> discounted_payoff = cv.get_payoff();

        writeResultToFile(spot_maturity, discounted_payoff, sol, strikes[i]);

    }
    for (int i=0; i<nn.size(); i++){

        controlVariates cv(r, sigma, S_0, T, K, nn[i]);
        cv.init();
        cv.exec();
        std::vector<double> sol = cv.get_sol();
        std::vector<double> spot_maturity = cv.get_spot();
        std::vector<double> discounted_payoff = cv.get_payoff();

        writeResultToFile(spot_maturity, discounted_payoff, sol);

    }

    return 0;
}