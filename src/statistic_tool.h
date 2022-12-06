//
// Created by Zihan Xu on 11/28/22.
//

#ifndef FINAL_PROJECT_STATISTIC_TOOL_H
#define FINAL_PROJECT_STATISTIC_TOOL_H
#include <vector>
#include <numeric>


double average(std::vector<double> const & vec);
double variance(std::vector<double> const & vec);
double std_deviation(std::vector<double> const & vec);
double std_error(std::vector<double> const & vec);
std::vector<double> confidenceInterval(std::vector<double> const & vec);
double normalCDF(double x);
double expected_Payoff(double t, double spot_price, double interest_rate, double volatility, double strike, double maturity);
void writeResultToFile(std::vector<double> & vec1, std::vector<double> & vec2, std::vector<double> & vec3);
void writeResultToFile(std::vector<double> & vec1, std::vector<double> & vec2, std::vector<double> & vec3, double strike);
#endif //FINAL_PROJECT_STATISTIC_TOOL_H
