//
// Created by Zihan Xu on 11/28/22.
//

#include "statistic_tool.h"
#include <vector>
#include <cmath>
#include <numeric>
#include "omp.h"
#include <iostream>
#include "math.h"

double average(std::vector<double> const & vec){
    return std::accumulate(vec.begin(), vec.end(), 0.0)/(vec.size());
//    int n = vec.size();
//    double sum = 0.;
//#pragma omp parallel for
//    for (int i = 0; i < n; i++){
//        sum += vec[i];
//    }
//
//    return sum/n;

}

double variance(std::vector<double> const & vec){
    int n = vec.size();
    double mean = average(vec);
    double variance = 0.;
//#pragma omp parallel for
    for (int i = 0; i < n; i++)
        variance += pow(vec[i] - mean, 2);

    return variance / (n - 1);

}
double std_deviation(std::vector<double> const & vec){
    return std::sqrt(variance(vec)); //sqrt(variance)
}

double std_error(std::vector<double> const & vec){
    return std_deviation(vec) / std::sqrt(vec.size()); //sd/sqrt(n)
}



double get_b(std::vector<double> const & spot, std::vector<double> const & payoff){
    if ( int(spot.size()) != int(payoff.size()) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");

    int n = spot.size();
    double mean_spot = average(spot);
    double mean_payoff = average(payoff);
    double numerator = 0., denominator = 0.;
//#pragma omp parallel for
    for (int i = 0; i < n; i++){
        numerator += (spot[i] - mean_spot) * (payoff[i] - mean_payoff);
        denominator += pow(spot[i] - mean_spot, 2);
    }
    return numerator / denominator;
}

std::vector<double> getYb(std::vector<double> const & spot, std::vector<double> const & payoff, double expected_spot){
    double b = get_b(spot, payoff);
    int n = spot.size();
    std::vector<double> Yb(n);
    for(int i = 0; i < n; i++){
        Yb[i] = payoff[i] - b * (spot[i] - expected_spot);
    }
    return Yb;
}


std::vector<double> confidenceInterval(std::vector<double> const & vec){
    std::vector<double> result(2);
    double mean = average(vec);
    double ste = std_error(vec);
    double z = 1.96; // default z = 1.96 for n > 1000
    if( vec.size() == 10 ){
        //n=10,df=9,t_9=2.262
        z = 2.262;
    }
    else if ( vec.size() == 100 ){
        // n=100,t=1.984
        z = 1.984;
    }
    else if ( vec.size() == 1000 ){
        // n=1000,z=1.962
        z = 1.962;
    }
    result[0] = mean - z * ste;
    result[1] = mean + z * ste;
    return result;
}

double normalCDF(double x)
{
    return std::erfc(-x / std::sqrt(2)) / 2;
}

double expected_Payoff(double t, double spot_price, double interest_rate, double volatility, double strike, double maturity)
{
    double d1 = (std::log(spot_price / strike) + (interest_rate + pow(volatility, 2) / 2) * (maturity - t)) / (volatility * std::sqrt(maturity - t));
    double d2 = d1 - volatility * std::sqrt(maturity - t);
    return normalCDF(d1) * spot_price * std::exp(interest_rate * maturity) - normalCDF(d2) * strike;
}