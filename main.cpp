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
    void check_reduction(double interest_rate, double sigma, double S_init, double maturity, Call call);


    double r = 0.05, sigma = 0.3;
    double S_0 = 50.0, T = 0.25;
    Call my_call(55.);
    check_reduction(r, sigma, S_0, T, my_call);

    //check reduction of variance
//    std::vector<double> KK = {40., 45., 50., 55., 60., 65., 70.};


//    for (int i=0; i<KK.size(); i++){
//        std::cout<< "K = "<< KK[i]<< std::endl;
//        Call my_call(KK[i]);
//        check_reduction(r, sigma, S_0, T, my_call);
//    }




    double interest_rate = 0.05;
    double S_init = 50.0, maturity = 0.25;
    double strike = 55.0;
    Call call(strike);
    BSMModel my_model(interest_rate, sigma);
    double exp_factor = std::exp(interest_rate * maturity);
    double expected_spot = S_init * exp_factor; //S(0)*exp(rT)


    std::vector<int> nn = {10, 100, 1000, 10000};
    double EY = expected_Payoff(0, S_init, interest_rate, sigma, strike, maturity);
    std::cout<<"E[Y] = "<< EY << std::endl;
    int m = 10000;


    for(int i = 0; i < nn.size(); i++){
        int n = nn[i];
        std::cout << "n = "<< n <<", ";
        int count = 0;

#pragma omp parallel for
        for(int j = 0; j < m; j++) {
            std::vector<double> spot_maturity(n);
            std::vector<double> payoff(n);

            for (int k = 0 ; k < n ; k++){
                spot_maturity[k] = (my_model.sim_path(maturity, S_init, 2))[1];  // only terminal price is needed
                payoff[k] = call(spot_maturity[k]);
            }

            std::vector<double> Yb = getYb(spot_maturity, payoff, expected_spot);

            std::vector<double> interval = confidenceInterval(Yb);
            if( interval[0] < EY && EY < interval[1])
            {
                count += 1;
            }
        }

        std::cout << "the probability that E[Y] lies in the 95% confidence interval is " << count / (m + 0.)<< std::endl;

    }


    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = endTime - startTime;


    std::cout << fp_ms.count() << std::endl;

    return 0;
}