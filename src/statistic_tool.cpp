//
// Created by Zihan Xu on 11/28/22.
//

#include <fstream>
#include <chrono>
#include "controlVariates.h"
#include "Payoff.h"
#include "BSMModel.h"
#include "statistic_tool.h"
#include <vector>
#include <cmath>
#include <numeric>
#include "omp.h"
#include <iostream>
#include "math.h"



double average(std::vector<double> const & vec){
    return std::accumulate(vec.begin(), vec.end(), 0.0)/(vec.size());
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
    return normalCDF(d1) * spot_price  - normalCDF(d2) * strike / std::exp(interest_rate * maturity);
}

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