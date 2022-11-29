#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "omp.h"
#include "statistic_tool.h"
#include "BSMModel.h"
#include "Payoff.h"
int main() {
    void check_reduction(double interest_rate, double sigma, double S_init, double maturity, Call call);
   // set up model

    double r = 0.05, sigma = 0.3;
    double S_0 = 50.0, T = 0.25;
    std::vector<double> K = {40., 45., 50., 55., 60., 65., 70.};


    for (int i=0; i<K.size(); i++){
        std::cout<< "K = "<< K[i]<< std::endl;
        Call my_call(K[i]);
        check_reduction(r, sigma, S_0, T, my_call);
    }





    return 0;
}