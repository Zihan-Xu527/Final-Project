//
// Created by Zihan Xu on 11/28/22.
//

#include "statistic_tool.h"
#include <vector>
#include <cmath>
#include <numeric>
#include "omp.h"

double average(std::vector<double> const & vec){
    return std::accumulate(vec.begin(), vec.end(), 0.0)/(vec.size());
}

double variance(std::vector<double> const & vec){
    int n = vec.size();
    double mean= average(vec);
    double variance = 0.0;
#pragma omp parallel for
    for (int i = 0; i < n; i++)
        variance += pow(vec[i] - mean, 2);

    return variance / (n - 1);

}

double covariance(std::vector<double> const & vec1, std::vector<double> const & vec2){
    if ( int(vec1.size()) != int(vec2.size()) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");

    int n = vec1.size();
    double mean1 = average(vec1);
    double mean2 = average(vec2);
    double covariance=0.;
#pragma omp parallel for
    for (int i = 0; i < n; i++)
        covariance += (vec1[i] - mean1) * (vec2[i] - mean2);

    return covariance / (n - 1);

}

double getb(std::vector<double> const & vec1, std::vector<double> const & vec2){
    if ( int(vec1.size()) != int(vec2.size()) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");

    int n = vec1.size();
    double mean1 = average(vec1);
    double mean2 = average(vec2);
    double numerator = 0., denominator = 0.;
#pragma omp parallel for
    for (int i = 0; i < n; i++){
        numerator += (vec1[i] - mean1) * (vec2[i] - mean2);
        denominator += pow(vec1[i] - mean1, 2);
    }
    return numerator / denominator;
}

std::vector<double> confidenceInterval(double mean, double sd, int n){
    std::vector<double> result(2);
    if (n < 100)
    { //n=10,df=9,t_9=2.262
        result[0] = mean - 2.262 * sd;
        result[1] = mean + 2.262 * sd;
    }
    else if (n < 1000)
    { //n=100,t=1.984
        result[0] = mean - 1.984 * sd;
        result[1] = mean + 1.984 * sd;
    }
    else
    { //z=1.96
        result[0] = mean - 1.96 * sd;
        result[1] = mean + 1.96 * sd;
    }
    return result;
}

double prob_interval(std::vector<double> const & interval, std::vector<double> const & samples){

    double count = 0.;
    int n = samples.size();

    for(int i = 0; i < n ; i++){
        if(samples[i] > interval[0] && samples[i] < interval[1])
            count += 1.;
    }
    return double( count / n );
}