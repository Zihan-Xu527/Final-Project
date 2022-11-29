//
// Created by Zihan Xu on 11/28/22.
//

#include "statistic_tool.h"
#include <vector>
#include <cmath>
#include <numeric>
#include "omp.h"

double average(std::vector<double> const & vec){
    int n = vec.size();
    double sum = 0.0;

#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        sum += vec[i];
    }
    return sum / n;
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

