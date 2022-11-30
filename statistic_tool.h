//
// Created by Zihan Xu on 11/28/22.
//

#ifndef FINAL_PROJECT_STATISTIC_TOOL_H
#define FINAL_PROJECT_STATISTIC_TOOL_H
#include <vector>
#include <numeric>


double average(std::vector<double> const & vec);
double variance(std::vector<double> const & vec);
double covariance(std::vector<double> const & vec1, std::vector<double> const & vec2);
double getb(std::vector<double> const & vec1, std::vector<double> const & vec2);
std::vector<double> confidenceInterval(double mean, double sd, int n);
double prob_interval(std::vector<double> const & interval, std::vector<double> const & samples);
#endif //FINAL_PROJECT_STATISTIC_TOOL_H
